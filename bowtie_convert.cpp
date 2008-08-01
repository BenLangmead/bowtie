/*
 *  bowtie_convert.cpp
 *  Bowtie
 *
 *  Created by Cole Trapnell on 6/30/08.
 * 
 *
 */

#include <string>
#include <iostream>
#include <map>
#include <stdio.h>
#include <seqan/sequence.h>
#include "maq/maqmap.h"
#include "maq/algo.hh"
#include "tokenize.h"
#include "params.h"
#include "pat.h"

enum {TEXT_MAP, BIN_MAP};

using namespace std;
static bool verbose	= false;

// A lookup table for fast integer logs
static int log_n[256];

static void print_usage() 
{
	cout << "Usage: bowtie_convert <in.bwtmap> <out.map> [refnames]" << endl;
	cout << "    -v     verbose output" << endl;
}

// Some of these limits come from Maq headers, beware...
static const int max_read_name = MAX_NAMELEN;
static const int max_read_bp = MAX_READLEN;

// Maq's default mapping quality
static int DEFAULT_QUAL = 25;
static int FIVE_PRIME_PHRED_QUAL = 'Z' - 33;
//static int THREE_PRIME_PHRED_QUAL = ';' - 33;

// Number of bases consider "reliable" on the five prime end of each read
static int MAQ_FIVE_PRIME = 24;

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
					   const string* refnames_fname = NULL)
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

	map<unsigned int, string> seqid_to_name;
	char bwt_buf[2048];
	static const int buf_size = 256;
	char name[buf_size];
	int bwtf_ret = 0;
	
	FILE* refnamef = NULL;
	if (refnames_fname)
	{
		refnamef = fopen(refnames_fname->c_str(), "r");	
		if (!refnamef)
		{
			fprintf(stderr, 
					"Error: could not open reference names file for reading\n");
			exit(1);
		}
		while (fgets(bwt_buf, 2048, refnamef))
		{
			//char* nl = strrchr(bwt_buf, '\n');
			//if (nl) *nl = 0;
			int seqid = 0;
			bwtf_ret = sscanf(bwt_buf, "%d %s\n", &seqid, name);
			if (bwtf_ret != 2)
			{
				fprintf(stderr, 
						"Warning: bad sequence id, name pair, skipping\n");
				continue;
			}
			seqid_to_name[seqid] = name;
		}
	}
	
	// Initialize a new Maq map table
	maqmap_t *mm = maq_new_maqmap();

	/**
	 Fields are:
		1) name (or for now, Id)
		2) orientations ('+'/'-')
		3) text id
		4) text offset
		5) sequence of hit (i.e. trimmed read)
	    6) quality values of sequence (trimmed)
		7) # of other hits in EBWT 
		8) mismatch positions - this is a comma-delimited list of positions
			w.r.t. the 5 prime end of the read.
	 */
	char* bwt_fmt_str = "%s %c %d %d %s %s %d %s";

	
	char orientation;
	unsigned int text_id;
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
						  &text_id, 
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
		
		// Convert sequence into Maq's bitpacked format
		memset(m1->seq, 0, max_read_bp);
		m1->size = strlen(sequence);
		
		m1->seqid = text_id;
		
		if (!refnames_fname && 
			seqid_to_name.find(text_id) == seqid_to_name.end())
		{
			char name_buf[1024];
			sprintf(name_buf, "%d", text_id);
			seqid_to_name[text_id] = name_buf;
			
		}
		
		//const char* default_quals = "EDCCCBAAAA@@@@?>===<;;9:99987776666554444";
		//int qual_len = sizeof(default_quals);
		
		int qual_len = strlen(qualities);
		
		if (orientation == '+')
		{
			for (int i = 0; i != m1->size; ++i) 
			{
				int tmp = nst_nt4_table[(int)sequence[i]];
				m1->seq[i] = (tmp > 3)? 0 : 
				(tmp <<6 | (i < qual_len ? qualities[i] - 33 : 0));
			}
		}
		else
		{
			for (int i = m1->size - 1; i >= 0; --i)
			{
				int tmp = nst_nt4_table[(int)sequence[i]];
				tmp = (m1->seq[i] == 0) ? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
				m1->seq[m1->size-i-1] = m1->seq[i] = (tmp > 3)? 0 : 
				 (tmp <<6 |  (i < qual_len ? qualities[i] - 33 : 0));
			}	
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
		
		// FIXME: this is a bullshit mapping quality, we need to consider
		// mismatches, etc.
		// m1->map_qual = m1->seq[MAX_READLEN-1] = m1->alt_qual = DEFAULT_QUAL;
		
		m1->map_qual = cal_map_qual(DEFAULT_QUAL, 
									five_prime_mismatches, 
									other_occs);
		
		m1->seq[MAX_READLEN-1] = m1->alt_qual = m1->map_qual;
		
		++mm->n_mapped_reads;
	}
	
	mm->n_ref = seqid_to_name.size();
	mm->ref_name = (char**)malloc(sizeof(char*) * mm->n_ref);
	int j = 0;
	for (map<unsigned int,string>::iterator i = seqid_to_name.begin();
		 i != seqid_to_name.end(); ++i)
	{
		char* name = strdup(i->second.c_str());
		if (name)
			mm->ref_name[j++] = name;
		//cerr << mm->ref_name[i->first] << endl;
	}
	
	
	algo_sort(mm->n_mapped_reads, mm->mapped_reads);
	maqmap_write_header(maqf, mm);
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

int main(int argc, char **argv) 
{
	string bwtmap_filename;
	string maqmap_filename;
	string* refnames_filename = NULL;
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
	if(optind < argc)
		refnames_filename = new string(argv[optind++]);
	
	init_log_n();
	
	int ret = convert_bwt_to_maq(bwtmap_filename, maqmap_filename, refnames_filename);
	
	delete refnames_filename;
	return ret;
}
