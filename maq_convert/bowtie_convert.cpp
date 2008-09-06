/*
 * This source file is derived from Maq v0.6.6.  It is distributed
 * under the GNU GENERAL PUBLIC LICENSE v2 with no warranty.  See file
 * "COPYING" in this directory for the terms.
 */

#include <string>
#include <iostream>
#include <set>
#include <map>
#include <stdio.h>
#include <seqan/sequence.h>
#include <algorithm>
#include "maqmap.h"
#include "algo.hh"
#include "bfa.h"
#include "tokenize.h"
#include "formats.h"
#include "pat.h"

enum {TEXT_MAP, BIN_MAP};

using namespace std;
static bool verbose	= false;

// A lookup table for fast integer logs
static int log_n[256];

static void print_usage()
{
	cout << "Usage: bowtie-convert <in.bwtmap> <out.map> <chr.bfa>" << endl;
	cout << "    -v     verbose output" << endl;
}

// Some of these limits come from Maq headers, beware...
static const int max_read_name = MAX_NAMELEN;
static const int max_read_bp = MAX_READLEN;

// Maq's default mapping quality
static int DEFAULT_QUAL = 25;

// Number of bases consider "reliable" on the five prime end of each read
static int MAQ_FIVE_PRIME = 28;

static inline int operator < (const maqmap1_t &a, const maqmap1_t &b)
{
	return (a.seqid < b.seqid) || (a.seqid == b.seqid && a.pos < b.pos);
}

// A simple way to compute mapping quality.
// Reads are ranked in three tiers, 0 seed mismatches, 1 seed mismatches,
// 2 seed mismatches.  Within each rank, the quality is reduced by the log of
// the number of alternative mappings at that level.
static inline int cal_map_qual(int default_qual,
							   unsigned int seed_mismatches,
							   unsigned int other_occs)
{
	if (seed_mismatches == 0)
		return 3 * default_qual - log_n[other_occs];
	if (seed_mismatches == 1)
		return 2 * default_qual - log_n[other_occs];
	else
		return default_qual - log_n[other_occs];
}

int convert_bwt_to_maq(const string& bwtmap_fname,
					   const string& maqmap_fname,
					   const map<string, unsigned int>& names_to_ids)
{
	FILE* bwtf = fopen(bwtmap_fname.c_str(), "r");

	if (!bwtf)
	{
		fprintf(stderr, "Error: could not open Bowtie mapfile %s for reading\n", bwtmap_fname.c_str());
		exit(1);
	}

	void* maqf = gzopen(maqmap_fname.c_str(), "w");
	if (!maqf)
	{
		fprintf(stderr, "Error: could not open Maq mapfile %s for writing\n", maqmap_fname.c_str());
		exit(1);
	}

	std::map<string, int> seqid_to_name;
	char bwt_buf[2048];
	static const int buf_size = 256;
	char name[buf_size];
	int bwtf_ret = 0;
	uint32_t seqid = 0;

	// Initialize a new Maq map table
	maqmap_t *mm = maq_new_maqmap();

	/**
	 Fields are:
		1) name (or for now, Id)
		2) orientations ('+'/'-')
		3) text name
		4) text offset
		5) sequence of hit (i.e. trimmed read)
	    6) quality values of sequence (trimmed)
		7) # of other hits in EBWT
		8) mismatch positions - this is a comma-delimited list of positions
			w.r.t. the 5 prime end of the read.
	 */
	const char* bwt_fmt_str = "%s %c %s %d %s %s %d %s";


	char orientation;
	char text_name[buf_size];
	unsigned int text_offset;
	char sequence[buf_size];

	char qualities[buf_size];
	unsigned int other_occs;
	char mismatches[buf_size];

	int max = 0;

	while (fgets(bwt_buf, 2048, bwtf))
	{
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;

		memset(mismatches, 0, sizeof(mismatches));

		bwtf_ret = sscanf(bwt_buf,
						  bwt_fmt_str,
						  name,
						  &orientation,
						  text_name,   // name of reference sequence
						  &text_offset,
						  sequence,
						  qualities,
						  &other_occs,
						  mismatches);

		if (bwtf_ret > 0 && bwtf_ret < 6)
		{
			fprintf(stderr, "Warning: found malformed record, skipping\n");
			continue;
		}

		if (mm->n_mapped_reads == (bit64_t)max)
		{
			max += 0x100000;
			mm->mapped_reads = (maqmap1_t*)realloc(mm->mapped_reads,
												   sizeof(maqmap1_t) * (max));
		}

		maqmap1_t *m1 = mm->mapped_reads + mm->n_mapped_reads;
		strncpy(m1->name, name, max_read_name-1);
		m1->name[max_read_name-1] = 0;

		text_name[buf_size-1] = '\0';

		// Convert sequence into Maq's bitpacked format
		memset(m1->seq, 0, max_read_bp);
		m1->size = strlen(sequence);

		map<string, unsigned int>::const_iterator i_text_id =
			names_to_ids.find(text_name);
		if (i_text_id == names_to_ids.end())
		{
			fprintf(stderr, "Warning: read maps to text not in BFA, skipping\n");
			continue;
		}

		m1->seqid = i_text_id->second;

//		if (seqid_to_name.find(text_name) == seqid_to_name.end()) {
//			// Map the alignment id to the name of the reference sequence
//			m1->seqid = seqid_to_name.size(); // 'seqid' is a unique id for alignments
//			seqid_to_name[text_name] = m1->seqid;
//		} else {
//			m1->seqid = seqid_to_name[text_name];
//		}

		int qual_len = strlen(qualities);

		for (int i = 0; i != m1->size; ++i)
		{
			int tmp = nst_nt4_table[(int)sequence[i]];
			m1->seq[i] = (tmp > 3)? 0 :
			(tmp <<6 | (i < qual_len ? (qualities[i] - 33) & 0x3f : 0));
		}

		vector<string> mismatch_tokens;
		tokenize(mismatches, ",", mismatch_tokens);

		vector<int> mis_positions;
		int five_prime_mismatches = 0;
		int three_prime_mismatches = 0;
		int seed_mismatch_quality_sum = 0;
		for (unsigned int i = 0; i < mismatch_tokens.size(); ++i)
		{
			mis_positions.push_back(atoi(mismatch_tokens[i].c_str()));
			int pos = mis_positions.back();
			if (pos < MAQ_FIVE_PRIME)
			{
				// Mismatch positions are always with respect to the 5' end of
				// the read, regardless of whether the read mapping is antisense
				if (orientation == '+')
					seed_mismatch_quality_sum += qualities[pos] - 33;
				else
					seed_mismatch_quality_sum += qualities[m1->size - pos - 1] - 33;

				++five_prime_mismatches;
			}
			else
			{
				// Maq doesn't care about qualities of mismatching bases beyond
				// the seed
				++three_prime_mismatches;
			}
		}

		// Only have 256 entries in the log_n table, so we need to pin other_occs
		other_occs = min(other_occs, 255u);

		m1->c[0] = m1->c[1] = 0;
		if (three_prime_mismatches + five_prime_mismatches)
			m1->c[1] = other_occs + 1; //need to include this mapping as well!
		else
			m1->c[0] = other_occs + 1; //need to include this mapping as well!

		// Unused paired-end data
		m1->flag = 0;
		m1->dist = 0;

		m1->pos = (text_offset << 1) | (orientation == '+'? 0 : 1);

		// info1's high 4 bits are the # of seed mismatches/
		// the low 4 store the total # of mismatches
		m1->info1 = (five_prime_mismatches << 4) |
			(three_prime_mismatches + five_prime_mismatches);

		// Sum of qualities in seed mismatches
		seed_mismatch_quality_sum = ((seed_mismatch_quality_sum <= 0xff) ?
									 seed_mismatch_quality_sum : 0xff);
		m1->info2 = seed_mismatch_quality_sum;

		m1->map_qual = cal_map_qual(DEFAULT_QUAL,
									five_prime_mismatches,
									other_occs);

		m1->seq[MAX_READLEN-1] = m1->alt_qual = m1->map_qual;

		++mm->n_mapped_reads;
		seqid++; // increment unique id for reads
	}

	mm->n_ref = names_to_ids.size();
	mm->ref_name = (char**)malloc(sizeof(char*) * mm->n_ref);

	for (map<string, unsigned int>::const_iterator i = names_to_ids.begin();
		 i != names_to_ids.end(); ++i)
	{
		char* name = strdup(i->first.c_str());
		if (name)
			mm->ref_name[i->second] = name;
		//cerr << mm->ref_name[i->first] << endl;
	}

	algo_sort(mm->n_mapped_reads, mm->mapped_reads);

	// Write out the header
	maqmap_write_header(maqf, mm);
	// Write out the alignments
	gzwrite(maqf, mm->mapped_reads, sizeof(maqmap1_t) * mm->n_mapped_reads);

	maq_delete_maqmap(mm);

	fclose(bwtf);
	gzclose(maqf);
	return 0;
}

void init_log_n()
{
	log_n[0] = -1;
	for (int i = 1; i != 256; ++i)
		log_n[i] = (int)(3.434 * log(i) + 0.5);
}

void get_names_from_bfa(const string& bfa_filename,
				   map<string, unsigned int>& names_to_ids)
{
	FILE* bfaf = fopen(bfa_filename.c_str(), "r");

	if (!bfaf)
	{
		fprintf(stderr, "Error: could not open Binary FASTA file %s for reading\n", bfa_filename.c_str());
		exit(1);
	}

	unsigned int next_id = 0;
	nst_bfa1_t *l;

	while ((l = nst_load_bfa1(bfaf)) != 0)
	{
		names_to_ids[l->name] = next_id++;
	}
}

int main(int argc, char **argv)
{
	string bwtmap_filename;
	string maqmap_filename;
	string bfa_filename;
	const char *short_options = "v";
	int next_option;
	do {
		next_option = getopt(argc, argv, short_options);
		switch (next_option) {
	   		case 'v': /* verbose */
				verbose = true;
				break;

			case -1: /* Done with options. */
				break;
			default:
				print_usage();
				return 1;
		}
	} while(next_option != -1);

	if(optind >= argc) {
		print_usage();
		return 1;
	}

	// The Bowtie output text file to be converted
	bwtmap_filename = argv[optind++];

	if(optind >= argc) {
		print_usage();
		return 1;
	}

	// The name of the binary Maq map to be written
	maqmap_filename = argv[optind++];

	// An optional argument:
	// a two-column text file of [Bowtie ref id, reference name string] pairs
	if(optind >= argc)
	{
		print_usage();
		return 1;
	}
	bfa_filename = string(argv[optind++]);

	init_log_n();

	map<string, unsigned int> names_to_ids;
	get_names_from_bfa(bfa_filename, names_to_ids);

	int ret = convert_bwt_to_maq(bwtmap_filename, maqmap_filename, names_to_ids);

	return ret;
}
