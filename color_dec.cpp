/*
 * color_dec.cpp
 *
 *  Created on: October 15, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include <string>
#include <stdlib.h>
#include "alphabet.h"
#include "color_dec.h"
#include "color.h"
#include "qual.h"

using namespace std;

// 4-bit pop count
static int alts[] = {
	-1, 1, 1, 2, 1, 2, 2, 3,
	 1, 2, 2, 3, 2, 3, 3, 4
};
static int firsts[] = {
	-1, 0, 1, 0, 2, 0, 1, 0,
	 3, 0, 1, 0, 2, 0, 1, 0
};

/**
 * Given a nucleotide mask, pick one matching nucleotide at random.
 */
static int randFromMask(int mask) {
	assert_gt(mask, 0);
	if(alts[mask] == 1) return firsts[mask];
	assert_gt(mask, 0);
	assert_lt(mask, 16);
	int r = rand() % alts[mask];
	assert_geq(r, 0);
	assert_lt(r, alts[mask]);
	for(int i = 0; i < 4; i++) {
		if((mask & (1 << i)) != 0) {
			if(r == 0) return i;
			r--;
		}
	}
	cerr << "Shouldn't get here" << endl;
	throw 1;
	return -1;
}

/**
 * Does a 2-bit-encoded base match any bit in a mask?
 */
static inline bool matches(int i, int j) {
	return ((1 << i) & j) != 0;
}

/**
 * Given the dynamic programming table, trace backwards from the last
 * column and populate the 's' and 'cmm' strings accordingly.  Whenever
 * there are multiple equally good ways of backtracking, choose one at
 * random.
 */
static void backtrack(int table[4][6][1025], // filled-in table
                      const char *read, size_t readi, size_t readf,
                      const char *ref, size_t refi, size_t reff,
                      char *s,    // final nucleotide string
                      char *cmm,  // color mismatches
                      char *nmm,  // nucleotide mismatches
                      int& cmms,  // # color mismatches
                      int& nmms)  // # nucleotide mismatches
{
	const size_t len = reff-refi;
	cmms = nmms = 0;
	int min = INT_MAX;
	int bests = 0;
	// Determine best base in final column of table
	for(int i = 0; i < 4; i++) {
		// Install minimum and backtrack info
		int m = table[i][4][len-1];
		if(m < min) {
			min = m;
			bests = (1 << i);
		} else if(m == min) {
			bests |= (1 << i);
		}
	}
	// i <- position of rightmost nucleotide
	int i = (int)len-1;
	// to <- rightmost nucleotide
	int to = randFromMask(bests);
	while(true) {
		bests = table[to][5][i]; // get next best mask
		s[i--] = to; // install best nucleotide
		if(i < 0) break; // done
		assert_gt(bests, 0);
		assert_lt(bests, 16);
		to = randFromMask(bests); // select
	}
	// Determine what reference nucleotides were matched against
	for(size_t i = 0; i < len; i++) {
		if(matches(s[i], ref[refi+i])) {
			assert_eq(1, alts[(int)ref[refi+i]]);
			// Just plain matched
			nmm[i] = 'M';
		} else {
			// If ref is ambiguous here, does it matter which one we
			// choose?  I don't think so.
			assert_eq(1, alts[(int)ref[refi+i]]);
			// SNP here
			nmm[i] = 'S';
			nmms++;
		}
	}
	for(size_t i = 0; i < len-1; i++) {
		int c1 = (int)read[readi+i]; // actual
		int c2 = dinuc2color[(int)s[i]][(int)s[i+1]]; // decoded
		assert_leq(c1, 4); assert_geq(c1, 0);
		if(c1 != c2 || c1 == 4) {
			// Actual != decoded
			assert_lt(c2, 4); assert_geq(c2, 0);
			cmm[i] = "0123."[c2];
			cmms++;
		} else {
			cmm[i] = 'M';
		}
	}
	// done
}

/**
 * Decode the colorspace read 'read' as aligned against the reference
 * string 'ref', assuming that it's a hit.
 */
void decodeHit(
		const char *read, // ASCII colors, '0', '1', '2', '3', '.'
		const char *qual, // ASCII quals, Phred+33 encoded
		size_t readi, // offset of first character within 'read' to consider
		size_t readf, // offset of last char (exclusive) in 'read' to consider
		const char *ref, // reference sequence, as masks
		size_t refi, // offset of first character within 'ref' to consider
		size_t reff, // offset of last char (exclusive) in 'ref' to consider
		int snpPhred, // penalty incurred by a SNP
		char *ns,  // decoded nucleotides are appended here
		char *cmm, // where the color mismatches are in the string
		char *nmm, // where nucleotide mismatches are in the string
		int& cmms, // number of color mismatches
		int& nmms) // number of nucleotide mismatches
{
	assert_lt(refi, reff);
	assert_lt(readi, readf);
	assert_eq(reff-refi-1, readf-readi);

	//
	// Dynamic programming table; good for colorspace reads up to 1024
	// colors in length.
	//
	int table[4][6][1025];
	// 0 -> A, 1 -> C, 2 -> G, 3 -> T, 4 -> min(A, C, G, T),
	// 5 -> backtrack mask, 6 -> min mismatches

	// The first column of the table just considers the first
	// nucleotide and whether it matches the ref nucleotide.
	for(int to = 0; to < 4; to++) {
		if(matches(to, ref[refi])) {
			// The assigned subject nucleotide matches the reference;
			// no penalty
			table[to][0][0] = 0;
			table[to][1][0] = 0;
			table[to][2][0] = 0;
			table[to][3][0] = 0;
			table[to][4][0] = 0;
			table[to][5][0] = 15;
		} else {
			// The assigned subject nucleotide does not match the
			// reference nucleotide, so we add a SNP penalty
			table[to][0][0] = snpPhred;
			table[to][1][0] = snpPhred;
			table[to][2][0] = snpPhred;
			table[to][3][0] = snpPhred;
			table[to][4][0] = snpPhred;
			table[to][5][0] = 15;
		}
	}

	// Successive columns examine successive alignment positions
	int omin = INT_MAX, t = 0;
	for(size_t c = readi; c < readf; c++) {
		const int readc = (int)read[c];
		assert_leq(readc, 4);
		assert_geq(readc, 0);
		omin = INT_MAX;
		// t <- index of column in dynamic programming table
		t = (int)(c - readi + 1);
		const int refc = ref[refi + t];
		int from[] = { table[0][4][t-1], table[1][4][t-1],
		               table[2][4][t-1], table[3][4][t-1] };
		// For each downstream nucleotide
		for(int to = 0; to < 4; to++) {
			// For each upstream nucleotide
			int min = INT_MAX;
			const int goodfrom = nuccol2nuc[to][readc];
			int q = qual[c];
			// Reward the preceding position
			if(goodfrom < 4) from[goodfrom] -= q;
			min = from[0];
			table[to][5][t] = 1;
			if(from[1] < min) {
				min = from[1];
				table[to][5][t] = 2;
			} else if(from[1] == min) {
				table[to][5][t] |= 2;
			}
			if(from[2] < min) {
				min = from[2];
				table[to][5][t] = 4;
			} else if(from[2] == min) {
				table[to][5][t] |= 4;
			}
			if(from[3] < min) {
				min = from[3];
				table[to][5][t] = 8;
			} else if(from[3] == min) {
				table[to][5][t] |= 8;
			}
			min += q;
			if(!matches(to, refc)) {
				min += snpPhred;
			}
			table[to][4][t] = min;
			if(min < omin) omin = min;
			if(goodfrom < 4) from[goodfrom] += q;
		}
	}

	t++;
	assert_eq(t, (int)(reff - refi));
	// Install the best backward path into ns, cmm, nmm
	backtrack(table,
	          read, readi, readi + t - 1,
	          ref, refi, refi + t,
	          ns, cmm, nmm, cmms, nmms);
}

#ifdef MAIN_COLOR_DEC

#include <sstream>
#include <getopt.h>

static const char *short_opts = "s:m:r:e:";
static struct option long_opts[] = {
	{(char*)"snppen",  required_argument, 0, 's'},
	{(char*)"misspen", required_argument, 0, 'm'},
	{(char*)"seed",    required_argument, 0, 'r'},
	{(char*)"maxpen",  required_argument, 0, 'e'}
};

template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

int main(int argc, char **argv) {
	int option_index = 0;
	int next_option;
	int snppen = 30;
	int misspen = 20;
	int maxPenalty = 70;
	unsigned seed = 0;
	do {
		next_option = getopt_long(argc, argv, short_opts, long_opts, &option_index);
		switch (next_option) {
			case 's': snppen = parse<int>(optarg); break;
			case 'm': misspen = parse<int>(optarg); break;
			case 'r': seed = parse<unsigned>(optarg); break;
			case 'e': maxPenalty = parse<int>(optarg); break;
			case -1: break;
			default: {
				cerr << "Unknown option: " << (char)next_option << endl;
				exit(1);
			}
		}
	} while(next_option != -1);
	srand(seed);
	if(argc - optind < 2) {
		cerr << "Not enough options" << endl;
		exit(1);
	}
	string read, ref;
	read = argv[optind];
	for(size_t i = 0; i < read.length(); i++) {
		read[i] = asc2col[(int)read[i]];
		assert_leq(read[i], 4);
		assert_geq(read[i], 0);
	}
	ref = argv[optind+1];
	for(size_t i = 0; i < ref.length(); i++) {
		int num = 0;
		int alts[] = {4, 4, 4, 4};
		decodeNuc(toupper(ref[i]), num, alts);
		assert_leq(num, 4);
		assert_gt(num, 0);
		ref[i] = 0;
		for(int j = 0; j < num; j++) {
			ref[i] |= (1 << alts[j]);
		}
	}
	string ns;
	string quals;
	quals.resize(read.length(), misspen);
	string cmm, nmm;
	int score = decode(read, quals, 0, read.length(),
	                   ref, 0, ref.length(), maxPenalty,
	                   snppen, ns, cmm, nmm);
	cout << " Score: " << score << " (max: " << maxPenalty << ")" << endl;
	cout << "   MMs:  ";
	for(size_t i = 0; i < cmm.length(); i++) {
		cout << cmm[i] << " ";
	}
	cout << endl;
	cout << "Colors:  ";
	for(size_t i = 0; i < read.length(); i++) {
		printColor((int)read[i]);
		cout << " ";
	}
	cout << endl;
	cout << " Bases: ";
	for(size_t i = 0; i < ns.length(); i++) {
		cout << "ACGTN"[(int)ns[i]] << " ";
	}
	cout << endl;
	cout << "   Ref: ";
	for(size_t i = 0; i < ref.length(); i++) {
		cout << mask2iupac[(int)ref[i]] << " ";
	}
	cout << endl;
	cout << "   MMs: ";
	for(size_t i = 0; i < ref.length(); i++) {
		cout << nmm[i] << " ";
	}
	cout << endl;
}
#endif
