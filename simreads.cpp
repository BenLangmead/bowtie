#ifdef SIMREADS_MAIN

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <getopt.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "params.h"
#include "sequence_io.h"
#include "tokenize.h"

using namespace std;
using namespace seqan;

// TODO: add option to scale error rate
// TODO: allow user to specify error profile
// TODO: use better model for qual values

static uint32_t numReads = 10;         // number of reads to generate
static int      readLen  = 35;         // length of reads to generate
static int64_t  cutoff   = 0xffffffff; // max # of reference bases
static int      seed     = 0;          // srandom seed
static int      errScale = 2;          // srandom seed
static bool     verbose  = false;      // be talkative
static bool     addzs    = true;       // add Zs from 3' end
static uint32_t snpRate  = 0;

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: simreads [options]* <reference_in> [<fa_outfile>] [<fq_outfile>]" << endl
	    << "    reference_in         comma-separated list of (multi-)FASTA reference files" << endl
	    << "    fa_outfile           FASTA output file (default: stdout)" << endl
	    << "    fq_outfile           FASTQ output file (default: stderr)" << endl
	    << "Options:" << endl
	    << "    -r/--numreads <int>  # reads to generate (default: 10)" << endl
	    << "    -l/--length <int>    length of reads to generate (up to 35)" << endl
	    << "    -c/--cutoff <int>    generate reads only up to <int>-length prefix of ref" << endl
	    << "    -e/--errscale <int>  increase errors by <int> times (0=disable errs)" << endl
	    << "    -p/--polyRate <int>  artificial polymorphism rate (0=disable)" << endl
	    << "    -s/--seed <int>      seed for random number generator" << endl
	    << "    -n/--nozs            don't add stretches of Ns from 3' end" << endl
	    << "    -v/--verbose         verbose output (for debugging)" << endl
	    ;
}

static const char *short_options = "r:l:c:s:e:p:nv";

static struct option long_options[] = {
	/* These options set a flag. */
	{"numreads", required_argument, 0, 'r'},
	{"length",   required_argument, 0, 'l'},
	{"cutoff",   required_argument, 0, 'c'},
	{"seed",     required_argument, 0, 's'},
	{"errscale", required_argument, 0, 'e'},
	{"polyRate", required_argument, 0, 'p'},
	{"nozs",     no_argument, 0, 'n'},
	{"verbose",  no_argument, 0, 'v'},
	{0, 0, 0, 0}
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', then output the given error message and
 * exit with an error and a usage message.
 */
template<typename T>
static int parseNumber(T lower, const char *errmsg) {
	char *endPtr= NULL;
	T t = (T)strtoll(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (t < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			exit(1);
		}
		return t;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	exit(1);
	return -1;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, char **argv) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
		switch (next_option) {
	   		case 'r':
	   			numReads = parseNumber<uint32_t>(1, "-r/--numreads arg must be at least 1");
	   			break;
	   		case 'l':
	   			readLen = parseNumber<int>(5, "-l/--length arg must be at least 5");
	   			break;
	   		case 'p':
	   			snpRate = parseNumber<uint32_t>(0, "-p/--polyRate arg must be at least 0");
	   			break;
	   		case 'c':
	   			cutoff = parseNumber<int64_t>(1, "-c/--cutoff arg must be at least 1");
	   			break;
	   		case 's':
	   			seed = parseNumber<int>(0, "-s/--seed arg must be at least 0");
	   			break;
	   		case 'e':
	   			errScale = parseNumber<int>(0, "-e/--errscale arg must be at least 0");
	   			break;
	   		case 'n': addzs = false; break;
	   		case 'v': verbose = true; break;
			case -1: /* Done with options. */ break;
			case 0: if (long_options[option_index].flag != 0) break;	
			default: 
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
		}
	} while(next_option != -1);
}

// Simple per-cycle error model from eyeballing error rate graphs from
// various sources.
static float solexaErrs[] = {
	0.0001, // cycle 0
	0.0001, // cycle 1
	0.0001, // cycle 2
	0.0001, // cycle 3
	0.0001, // cycle 4
	0.0001, // cycle 5
	0.0002, // cycle 6
	0.0003, // cycle 7
	0.0003, // cycle 8
	0.0004, // cycle 9
	0.0007, // cycle 10
	0.0012, // cycle 11
	0.0021, // cycle 12
	0.0031, // cycle 13
	0.0040, // cycle 14
	0.0050, // cycle 15
	0.0053, // cycle 16
	0.0056, // cycle 17
	0.0060, // cycle 18
	0.0064, // cycle 19
	0.0068, // cycle 20
	0.0073, // cycle 21
	0.0078, // cycle 22
	0.0084, // cycle 23
	0.0092, // cycle 24
	0.0100, // cycle 25 (1%)
	0.0109, // cycle 26
	0.0118, // cycle 27
	0.0125, // cycle 28
	0.0160, // cycle 30
	0.0201, // cycle 31 (2%)
	0.0250, // cycle 32
	0.0300, // cycle 33 (3%)
	0.0350, // cycle 34
	0.0400, // cycle 35 (4%)
};

// nProbs[i] is the probability that ~i Ns should be written to the
// read starting at the 3' end.  These numbers are based on the results
// of analyzing the SRR001115 1000-Genomes reads with
// summarize_solexa.pl.
/*
static float nProbs[] = {
	0.66000, // 0 Ns
	0.10700, // 1 Ns
	0.06900, // 2 Ns
	0.04300, // 3 Ns
	0.03200, // 4 Ns
	0.02000, // 5 Ns
	0.01100, // 6 Ns
	0.00400, // 7 Ns
	0.00100, // 8 Ns
	0.00001, // 9 Ns
	0.00001, // 10 Ns
	0.00001, // 11 Ns
	0.00001, // 12 Ns
	0.00001, // 13 Ns
	0.00001, // 14 Ns
	0.00001, // 15 Ns
	0.00001, // 16 Ns
	0.00001, // 17 Ns
	0.00001, // 18 Ns
	0.00001, // 19 Ns
	0.00001, // 20 Ns
	0.00001, // 21 Ns
	0.00001, // 22 Ns
	0.00001, // 23 Ns
	0.00001, // 24 Ns
	0.00001, // 25 Ns
	0.00001, // 26 Ns
	0.00001, // 27 Ns
	0.00001, // 28 Ns
	0.00001, // 29 Ns
	0.00080, // 30 Ns
	0.00100, // 31 Ns
	0.00300, // 32 Ns
	0.00800, // 33 Ns
	0.03100, // 34 Ns
	0.01000, // 35 Ns
};
*/

// nProbs[i] is the probability that ~i Ns should be written to the
// read starting at the 3' end.  These numbers are based on the results
// of analyzing the SRR001115 1000-Genomes reads with
// summarize_solexa.pl.  Out of 8839010 reads.  To make it fit into a
// range of 35 bases, I removed the following entries and subtracted
// those counts from the total.
// 962,     // 14 Ns
// 974,     // 15 Ns
// 940,     // 16 Ns
// 788,     // 17 Ns
// 756,     // 18 Ns
// 799,     // 19 Ns
// 780,     // 20 Ns
// 923,     // 21 Ns
// 923,     // 22 Ns
// 1049,    // 23 Ns
// 1144,    // 24 Ns
// 1568,    // 25 Ns
static uint32_t nProbs2Total = 8839010
	- 962
	- 974
	- 940
	- 788
	- 756
	- 799
	- 780
	- 923
	- 923
	- 1049
	- 1144
	- 1568;
static uint32_t nProbs2[] = {
	8386952, // 0 Ns
	8486,    // 1 Ns
	4439,    // 2 Ns
	3622,    // 3 Ns
	4764,    // 4 Ns
	6309,    // 5 Ns
	9069,    // 6 Ns
	8804,    // 7 Ns
	12100,   // 8 Ns
	4738,    // 9 Ns
	3028,    // 10 Ns
	1754,    // 11 Ns
	1362,    // 12 Ns
	1236,    // 13 Ns
	1041,    // 14 Ns
	1004,    // 15 Ns
	1260,    // 16 Ns
	1556,    // 17 Ns
	2878,    // 18 Ns
	13293,   // 19 Ns
	21642,   // 20 Ns
	3765,    // 21 Ns
	2388,    // 22 Ns
	2403,    // 23 Ns
	4039,    // 24 Ns
	5348,    // 25 Ns
	5246,    // 26 Ns
	6892,    // 27 Ns
	8914,    // 28 Ns
	30024,   // 29 Ns
	37498,   // 30 Ns
	8343,    // 31 Ns
	20949,   // 32 Ns
	58408,   // 33 Ns
	119359,  // 34 Ns
	14491,   // 35 Ns
};

// "Average" FASTQ quality characters as calculated over a FASTQ file
// corresponding to 1000-Genomes run SRR001113 using the
// fq_qual_freqs.pl script
static string fastqQuals = "EDCCCBAAAA@@@@?>===<;;9:9998777655443322";

template<typename TStr>
void driver(vector<string>& infiles, ostream& faout, ostream& fqout) {
	typedef typename Value<TStr>::Type TVal;
	typedef typename Value<TStr>::Type TVal;
	vector<TStr> ss;
	int64_t savedCutoff = cutoff;
	readSequenceFiles<TStr, Fasta>(infiles, ss, cutoff, -1, false);
	uint32_t bases = (uint32_t)(savedCutoff - cutoff);
	
	// Polymorphize the reference
	if(snpRate != 0) {
		uint32_t numSnps = bases / snpRate;
		for(uint32_t i = 0; i < numSnps; i++) {
			uint32_t r = rand() % bases;
			for(size_t j = 0; j < ss.size(); j++) {
				if(r < length(ss[j])) {
					ss[j][r] = (Dna)(((int)ss[j][r] + (rand()%3)) & 3);
					break;
				}
				r -= length(ss[j]);
			}
		}
	}
	
	// Chop characters off the end of the qual string as appropriate
	if((size_t)readLen < fastqQuals.length()) {
		fastqQuals = fastqQuals.substr(0, readLen);
	}
	// Scale error thresholds as appropriate
	if(errScale != 1) {
		for(size_t i = 0; i < 35; i++) {
			solexaErrs[i] *= errScale;
		}
	}
	// For each read
	for(uint32_t ri = 0; ri < numReads; ri++) {
		uint32_t boff = ((uint32_t)random()) % bases;
		TStr read;
		// Extract read from a text
		for(size_t ti = 0; ti < ss.size(); ti++) {
			uint32_t tlen = length(ss[ti]);
			if(boff < tlen) {
				boff = min(boff, tlen-readLen);
				if(verbose) {
					cout << "Grabbing from position " << boff << " in text " << ti << endl;
				}
				for(int i = 0; i < readLen; i++) {
					append(read, ss[ti][boff+i]);
				}
				if((random() % 2) == 0) {
					read = reverseComplement(read);
				}
				break;
			}
			boff -= tlen;
		}
		assert_gt(length(read), 0);
		// Mutate the read
		if(errScale > 0) {
			for(int i = 0; i < readLen; i++) {
				float f = (random() % 1000000) / 1000000.0;
				if(f < solexaErrs[i]) {
					// Mutate not just this base, but all subsequent
					// bases
					TVal old = read[i];
					read[i] = (TVal)((((int)read[i])+(random()&3))%ValueSize<TVal>::VALUE);
					// I tuned this to 4 in order to make the mapping
					// rates look like SRR001115
					if((random() % 5) != 0) {
						for(int j = i+1; j < readLen; j++) {
							read[j] = (TVal)((((int)read[j])+(random()&3))%ValueSize<TVal>::VALUE);
						}
						break;
					}
					if(verbose) {
						cout << "Mutated position " << i << " from " << old << " to " << read[i] << endl;
					}
				}
			}
		}
		// Possibly create a streak of Ns
		if(addzs) {
			uint32_t f = random() % nProbs2Total;
			for(int i = 0; i < 35; i++) {
				if(f < nProbs2[i]) {
					if(i == 0) break;
					for(int j = (int)readLen-1; j >= (int)readLen-i; j--) {
						if((random() % 10) != 0) {
							read[j] = 'N';
						}
					}
					break;
				}
				f -= nProbs2[i];
			}
		}
		
		// FASTA output
		faout << ">r" << ri << endl;
		faout << read << endl;
		// FASTQ output
		fqout << "@r" << ri << endl;
		fqout << read << endl;
		fqout << "+" << endl;
		fqout << fastqQuals << endl;
	}
}

int main(int argc, char **argv) {
	string infile;
	string faOutfile = "", fqOutfile = "";
	ostream *faout, *fqout;
	vector<string> infiles;
	parseOptions(argc, argv);
	srandom(seed);
	// Get input filename
	if(optind >= argc) {
		cerr << "No input sequence or sequence file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	infile = argv[optind++];
	tokenize(infile, ",", infiles);
	if(infiles.size() < 1) {
		cerr << "Tokenized input file list was empty!" << endl;
		printUsage(cerr);
		return 1;
	}
	// Get optional FASTA output filename
	if(optind < argc) {
		faOutfile = argv[optind++];
	}
	if(faOutfile.length() > 0) {
		faout = new ofstream(faOutfile.c_str());
	} else {
		faout = &cout;
	}
	// Get optional FASTQ output filename
	if(optind < argc) {
		fqOutfile = argv[optind++];
	}
	if(fqOutfile.length() > 0) {
		fqout = new ofstream(fqOutfile.c_str());
	} else {
		fqout = &cerr;
	}
	// Optionally summarize
	if(verbose) {
		cout << "Settings:" << endl
		     << "  Input files:" << endl;
		for(size_t i = 0; i < infiles.size(); i++) {
			cout << "    " << infiles[i] << endl;
		}
		cout << "  FASTA output file: " << faOutfile << endl;
		cout << "  FASTQ output file: " << fqOutfile << endl;
		cout << "  Num reads: " << numReads << endl;
		cout << "  Read length: " << readLen << endl;
		cout << "  Cutoff: " << cutoff << endl;
		cout << "  Seed: " << seed << endl;
	}
	driver<String<Dna5> >(infiles, *faout, *fqout); // allow Ns
}

#endif
