#ifdef BLOCKWISE_SA_MAIN

#include <getopt.h>
#include "blockwise_sa.h"

static bool verbose     = false;
static bool quiet       = false;
static bool sanityCheck = false;
static int  bmax        = 10;
static int  dcv         = 7;

/// Accept -q, -v , -s, -b, -d
static const char *short_options = "qvsb:d:";

static struct option long_options[] = {
	{"verbose",  no_argument, 0, 'v'}, // treat --verbose like -v
	{"quiet",    no_argument, 0, 'q'}, // treat --quiet like -q
	{"sanity",   no_argument, 0, 's'}, // treat --sanity like -s
	{"bmax",     required_argument, 0, 'b'}, // treat --bmax like -b
	{"dcv",      required_argument, 0, 'd'}, // treat --dcv like -d
	{0, 0, 0, 0}
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: blockwise_sa [options]* <input_sequence>" << endl
	    << "Description: Builds a suffix array for a sequence blockwise and outputs it." << endl
        << "Options:" << endl
	    << "    -b --bmax <int>  set maximum bucket size" << endl
	    << "    -d --dcv <int>   set difference-cover sample periodicity" << endl
	    << "    -v --verbose     verbose output (for debugging)" << endl
	    << "    -s --sanity      enable sanity checking (increased runtime/memore usage)" << endl
	    << "    -q --quiet       output nothing unless there is a problem (for testing)" << endl;
}

/**
 * Read command-line args
 */
static void parseOptions(int argc, char **argv) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
	int next_option;
	do {
		char *endPtr = NULL;
		long l;
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
		switch (next_option) {
			case 'v':
				verbose = 1;
				break;
			case 'q':
				quiet = 1;
				break;
	   		case 's':
	   			sanityCheck = true;
				break;
	   		case 'b':
	   			l = strtol(optarg, &endPtr, 10);
	   			if(endPtr != NULL) {
	   				if(l < 1) {
	   					cerr << "-b/--bmax <int> argument must be at least 1" << endl;
	   					printUsage(cerr);
	   					exit(1);
	   				}
	   				bmax = (int)l;
	   			}
	   			break;
	   		case 'd':
	   			l = strtol(optarg, &endPtr, 10);
	   			if(endPtr != NULL) {
	   				if(l < 3) {
	   					cerr << "-d/--dcv <int> argument must be at least 3" << endl;
	   					printUsage(cerr);
	   					exit(1);
	   				}
	   				dcv = (int)l;
	   			}
	   			break;
			case -1: /* Done with options. */
				break;
			default: 
				cerr << "Unknown option: " << (char)next_option << endl;
				//printUsage(cerr);
				exit(1);
		}
	} while(next_option != -1);
}

/**
 * Driver for experimenting with blockwise SA construction.
 */
int main(int argc, char **argv) {
	
	parseOptions(argc, argv);
	if(verbose) {
		cout << "Difference-cover sample periodicity: " << dcv << endl
		     << "Maximum bucket size: " << bmax << endl;
	}
	
	// Get input
	if(optind >= argc) {
		cerr << "No input sequence specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	String<Dna> text = argv[optind++];
	if(length(text) < 1) {
		cout << "Text is empty!" << endl;
		return 1;
	}
	
	// Initialize blockwise builder
	KarkkainenBlockwiseSA<String<Dna> > bsa(text, bmax, dcv, 0, sanityCheck, verbose, cout);
	assert(bsa.hasMoreBlocks());
	assert_eq(bsa.size(), length(text)+1);
	
	// Keep appending blocks onto asmSA until it contains the entire SA
	String<uint32_t> asmSA;
	reserve(asmSA, length(text)+1);
	while(bsa.hasMoreSuffixes()) {
		appendNE(asmSA, bsa.nextSuffix());
	}
	assert_eq(length(asmSA), length(text)+1);

	if(!quiet) {
		// Print it
		for(size_t i = 0; i < length(asmSA); i++) {
			cout << asmSA[i];
			if(i < length(asmSA)-1) {
				cout << ",";
			}
		}
		cout << endl;
	}

	if(sanityCheck) {
		// Create entire suffix array the old fashioned way
		String<uint32_t> sa;
		String<Dna5> text5 = text;
		append(text5, 'N'); // Force end-of-string to be higher than other chars
		resize(sa, length(text5), Exact());
		createSuffixArray(sa, text5, Skew7()); // Use difference-cover method
		assert_eq(length(sa), length(text5));
		
		// Check whether it matches blockwise SA
		if(verbose) cout << "Checking against Old-fashioned SA" << endl;
		bool matches = length(sa) == length(asmSA);
		if(matches) {
			for(size_t i = 0; i < length(sa); i++) {
				if(sa[i] != asmSA[i]) matches = false;
			}
		}

		// Print it if it didn't match (even if --quiet - we want our
		// tester to know something's wrong)
		if(!matches) {
			cout << "Didn't match!  Old-fashioned SA builder gave:" << endl;
			for(size_t i = 0; i < length(sa); i++) {
				cout << sa[i];
				if(i < length(sa)-1) {
					cout << ",";
				}
			}
			cout << endl;
			return 1;
		}

		// Now generate the suffix array using a silly blockwise
		// builder and iterate through both it and the true blockwise
		// SA in tandem.
		if(verbose) cout << "Checking against SillyBlockwiseDnaSA" << endl;
		bsa.resetSuffixItr();
		SillyBlockwiseDnaSA<String<Dna> > sillyBsa(text, bmax, sanityCheck, verbose, cout);
		assert_eq(sillyBsa.size(), length(text)+1);
		assert_eq(sillyBsa.size(), length(sa));
		assert_eq(bsa.size(), sillyBsa.size());
		for(size_t i = 0; i < length(sa); i++) {
			uint32_t sillyBsaSuf = sillyBsa.nextSuffix();
			uint32_t bsaSuf = bsa.nextSuffix();
			assert_eq(bsaSuf, sillyBsaSuf);
		}
		assert(!sillyBsa.hasMoreBlocks());
		assert(!bsa.hasMoreBlocks());
	}
	
	return 0; // success
}
#endif
