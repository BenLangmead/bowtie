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
static bool verbose				= false;
//static int format				= TEXT_MAP;



static void print_usage() 
{
	cout << "Usage: bowtie_convert <in.bwtmap> <out.map>" << endl;
	cout << "    -v     verbose output" << endl;
}

// Some of these limits come from Maq headers, beware...
static const int max_read_name = MAX_NAMELEN;
static const int max_read_bp = MAX_READLEN;

// Maq's default mapping quality
static int DEFAULT_QUAL = 25;

// Number of bases consider "reliable" on the five prime end of each read
static int MAQ_FIVE_PRIME = 24;

static inline int operator < (const maqmap1_t &a, const maqmap1_t &b)
{
	return (a.seqid < b.seqid) || (a.seqid == b.seqid && a.pos < b.pos);
}

int convert_bwt_to_maq(const string& bwtmap_fname, 
					   const string& maqmap_fname)
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

	// Initialize a new Maq map table
	maqmap_t *mm = maq_new_maqmap();

	/**
	 Fields are:
		1) name (or for now, Id)
		2) orientations ('+'/'-')
		3) text id
		4) text offset
		5) sequence of hit (i.e. trimmed read)
		6) # of other hits in EBWT (ignored for now)
		7) mismatch positions - this is a comma-delimited list of positions
			w.r.t. the 5 prime end of the read.
	 */
	char* bwt_fmt_str = "%s %c %d %d %s %s %s";
	static const int buf_size = 256;
	char name[buf_size];
	char orientation;
	unsigned int text_id;
	unsigned int text_offset;
	char sequence[buf_size];
	char other_occs[buf_size];
	char mismatches[buf_size];
	int bwtf_ret = 0;
	int max = 0;
	char bwt_buf[2048];
	
	map<unsigned int, string> seqid_to_name;
	
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
						  other_occs, 
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
		
		// FIXME: we need to get real text names, not just ascii conversions
		// of ids, and put them in this mapping.
		m1->seqid = text_id;
		if (seqid_to_name.find(text_id) == seqid_to_name.end())
		{
			char name_buf[1024];
			sprintf(name_buf, "%d", text_id);
			seqid_to_name[text_id] = name_buf;
			
		}
		
		if (orientation == '+')
		{
			for (int i = 0; i != m1->size; ++i) 
			{
				int tmp = nst_nt4_table[(int)sequence[i]];
				m1->seq[i] = (tmp > 3)? 0 : (tmp <<6 | DEFAULT_QUAL);
			}
		}
		else
		{
			for (int i = m1->size - 1; i >= 0; --i)
			{
				int tmp = nst_nt4_table[(int)sequence[i]];
				tmp = (m1->seq[i] == 0) ? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
				m1->seq[m1->size-i-1] = m1->seq[i] = (tmp > 3)? 0 : (tmp <<6 | DEFAULT_QUAL);
			}	
		}
		
		vector<string> mismatch_tokens;
		tokenize(mismatches, ",", mismatch_tokens);
		
		vector<int> mis_positions;
		int five_prime_mismatches = 0;
		int three_prime_mismatches = 0;
		for (unsigned int i = 0; i < mismatch_tokens.size(); ++i)
		{
			mis_positions.push_back(atoi(mismatch_tokens[i].c_str()));
			if (mis_positions.back() > MAQ_FIVE_PRIME)
				++five_prime_mismatches;
			else
				++three_prime_mismatches;
		}
		
		m1->c[0] = 0;
		m1->c[1] = 0;
		m1->flag = 0;
		m1->dist = 0;
		
		m1->pos = (text_offset)<<1 | (orientation == '+'? 0 : 1);
		m1->info1 = five_prime_mismatches << 4 | 
			(three_prime_mismatches + five_prime_mismatches);
		
		m1->info2 = (three_prime_mismatches + five_prime_mismatches) * DEFAULT_QUAL;
		
		// FIXME: this is a bullshit mapping quality, we need to consider
		// mismatches, etc.
		m1->map_qual = m1->seq[MAX_READLEN-1] = m1->alt_qual = DEFAULT_QUAL;
		++mm->n_mapped_reads;
	}
	
	mm->n_ref = seqid_to_name.size();
	mm->ref_name = (char**)malloc(sizeof(char*) * mm->n_ref);
	for (map<unsigned int,string>::iterator i = seqid_to_name.begin();
		 i != seqid_to_name.end(); ++i)
	{
		mm->ref_name[i->first] = strdup(i->second.c_str());
	}
	
	
	algo_sort(mm->n_mapped_reads, mm->mapped_reads);
	maqmap_write_header(maqf, mm);
	gzwrite(maqf, mm->mapped_reads, sizeof(maqmap1_t) * mm->n_mapped_reads);

	maq_delete_maqmap(mm);
	
	fclose(bwtf); 
	gzclose(maqf);
	return 0;
}

int main(int argc, char **argv) 
{
	string bwtmap_filename;
	string maqmap_filename;
	
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
	bwtmap_filename = argv[optind++];
	
	if(optind >= argc) {
		print_usage();
		return 1;
	}
	maqmap_filename = argv[optind++];
	
	return convert_bwt_to_maq(bwtmap_filename, maqmap_filename);
}