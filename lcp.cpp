#ifdef LCP_MAIN

#include <iostream>
#include <string>
#include <getopt.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include "sequence_io.h"
#include "params.h"
#include "tokenize.h"
#include "ebwt.h"

using namespace std;
using namespace seqan;

static int verbose           = 0;
static int sanityCheck       = 0;
static int format            = FASTA; // input sequence format
static int seqType           = DNA;   // input sequence alphabet

// Ebwt parameters
static int32_t chunkRate     = 11;

static const char *short_options = "drpfegcvsh:";

static struct option long_options[] = {
	/* These options set a flag. */
	{"verbose",      no_argument, 0, 'v'},
	{"sanity",       no_argument, 0, 's'},
	{"dna",          no_argument, &seqType, DNA},
	{"rna",          no_argument, &seqType, RNA},
	{"amino",        no_argument, &seqType, AMINO_ACID},
	{"protein",      no_argument, &seqType, AMINO_ACID},
	{"pro",          no_argument, &seqType, AMINO_ACID},
	{"chunkRate",    required_argument, 0, 'h'},
	{0, 0, 0, 0}
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: lcp [options]* <seq_in_list>" << endl
	    << "<seq_in_list> = Input sequences/sequence filenames separated by commas (no spaces)" << endl
	    << "    -d --dna                input uses 4-char DNA character set" << endl
	    << "    -r --rna                input uses 4-char RNA character set" << endl
	    << "    -p --protein            input uses 20+-char protein character set" << endl
	    //<< "    -k                      input file is packed (output of PackFasta)" << endl
	    << "    -f                      input file is Fasta (default)" << endl
	    << "    -e                      input file is Embl" << endl
	    << "    -g                      input file is Genbank" << endl
	    << "    -c                      input sequence is provided on command line (<seq_in>)" << endl
	    << "    -h --chunkRate <int>    # of characters in a text chunk" << endl
	    << "    -v --verbose            verbose output (for debugging)" << endl
	    << "    -s --sanity             enable sanity checks (much slower/increased memory usage)" << endl;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			exit(1);
		}
		return (int32_t)l;
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
			case 'd': seqType = DNA; break;
			case 'r': seqType = RNA; break;
	   		case 'p': seqType = AMINO_ACID; break;

	   		//case 'k': format = PACKED; break;
	   		case 'f': format = FASTA; break;
	   		case 'g': format = GENBANK; break;
	   		case 'e': format = EMBL; break;
	   		case 'a': format = RAW; break;
	   		case 'c': format = CMDLINE; break;

	   		case 'h':
	   			chunkRate = parseInt(1, "-h/--chunkRate arg must be at least 1");
	   			break;

	   		case 'v': verbose = true; break;
	   		case 's': sanityCheck = true; break;

			case -1: /* Done with options. */
				break;
			case 0:
				if (long_options[option_index].flag != 0)
					break;	
			default: 
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
		}
	} while(next_option != -1);
}

/**
 * 
 */
int main(int argc, char **argv) {
	typedef String<Dna> TStr;
	string infile;
	vector<string> infiles;
	
	parseOptions(argc, argv);

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

	// Optionally summarize
	if(verbose) {
		cout << "Input files (" << sequence_type_names[seqType] << ", "
		     << file_format_names[format] << "):" << endl;
		for(size_t i = 0; i < infiles.size(); i++) {
			cout << "  " << infiles[i] << endl;
		}
	}
	
	switch(seqType) {
		case DNA: {
			vector<String<Dna> > ss;
			if(verbose) cout << "About to read text files" << endl;
			switch(format) {
				case FASTA:   readSequenceFiles<TStr, Fasta>(infiles, ss);   break;
				case EMBL:    readSequenceFiles<TStr, Embl>(infiles, ss);    break;
			    case GENBANK: readSequenceFiles<TStr, Genbank>(infiles, ss); break;
				case RAW:     readSequenceFiles<TStr, Raw>(infiles, ss);     break;
				case CMDLINE: readSequenceString<TStr>(infile, ss);          break;
				//case PACKED:  readPackedSequenceFiles<TStr>(infiles, ss);    break;
				default: assert(false);
			}
			// Check that input is non-empty
			if(ss.size() == 0) {
				cerr << "Error: Empty input!  Check that file format is correct." << endl;
				return 1;
			}
			// Echo input strings
			if(verbose) {
				cout << "DNA input strings:" << endl;
				for(unsigned int i = 0; i < ss.size(); i++) {
					uint32_t ssz = length(ss[i]);
					if(ssz < 1000) {
						cout << "  " << i << ": " << ss[i] << endl;
					} else {
						cout << "  " << i << ": (" << ssz << " chars)" << endl;
					}
				}
			}
			String<Dna> joined = Ebwt<TStr>::join(ss, chunkRate, true);
			
			// Create suffix array
			String<uint32_t> sa;
			resize(sa, length(joined));
			createSuffixArray(sa, joined, Skew7());
			
			// Create LCP info
			String<uint32_t> lcp;
			resize(lcp, length(joined));
			createLCPTable(lcp, joined, sa, Kasai());
			
			// Analyze LCP info
			size_t min = 1 << 31;
			size_t max = 0;
			size_t ords[32] = {
				0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0,
			};
			for(size_t j = 0; j < length(joined); j++) {
				if(lcp[j] > max) max = lcp[j];
				if(lcp[j] < min) min = lcp[j];
				for(size_t k = 0; k < 32; k++) {
					if((lcp[j] >> k) > 0) {
						ords[k]++;
					} else {
						break;
					}
				}
			}
			cout << "min=" << min << ", max=" << max << endl;
			for(size_t k = 0; k < 32; k++) {
				if(ords[k] == 0) break;
				cout << "  >= 2^" << k << ": " << ords[k] << endl;
			}
			break;
		}
		case RNA: {
			cerr << "Error: No RNA support yet; only DNA support" << endl;
			break;
		}
		case AMINO_ACID: {
			cerr << "Error: No protein support yet; only DNA support" << endl;
			break;
		}
		default: {
			assert(false);
		}
	}
	
	return 0;
}
#endif
