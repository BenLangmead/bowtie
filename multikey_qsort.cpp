#include "multikey_qsort.h"

using namespace std;
using namespace seqan;

#ifdef MULTIKEY_QSORT_MAIN

#include <getopt.h>
#include "timer.h"
#include "params.h"
#include "tokenize.h"
#include "sequence_io.h"

static bool verbose     = false;
static bool quiet       = false;
static bool sanityCheck = false;
static bool suffixes    = false;
static int format       = FASTA; // input sequence format
static int seqType      = DNA;   // input sequence alphabet
static size_t upto      = 0xffffffff;
static uint32_t dcv     = 0;
static int seed         = 0;

static const int ARG_DCV       = 259;
static const int ARG_SEED      = 260;

static const char *short_options = "d:rpfegcxqvsu:";

static struct option long_options[] = {
	/* These options set a flag. */
	{"verbose",  no_argument, 0, 'v'},
	{"quiet",    no_argument, 0, 'q'},
	{"sanity",   no_argument, 0, 's'},
	{"suffixes", no_argument, 0, 'x'},
	{"dna",      no_argument, &seqType, DNA},
	{"rna",      no_argument, &seqType, RNA},
	{"amino",    no_argument, &seqType, AMINO_ACID},
	{"protein",  no_argument, &seqType, AMINO_ACID},
	{"pro",      no_argument, &seqType, AMINO_ACID},
	{"upto",     required_argument, 0, 'u'},
	{"dcv",      required_argument, 0, ARG_DCV},
	{"seed",     required_argument, 0, ARG_SEED},
	{0, 0, 0, 0}
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: multikey_qsort [options]* <input_sequences>" << endl
	    << "Description: Sorts a list of input sequences, or suffixes of a single sequence." << endl
        << "Options:" << endl
	    << "    -d --dna                input uses 4-char DNA character set" << endl
	    << "    -r --rna                input uses 4-char RNA character set" << endl
	    << "    -p --protein            input uses 20+-char protein character set" << endl
	    //<< "    -k                      input file is packed (output of PackFasta)" << endl
	    << "    -f                      input file is Fasta (default)" << endl
	    << "    -e                      input file is Embl" << endl
	    << "    -g                      input file is Genbank" << endl
	    << "    -c                      input sequence is provided on command line (<seq_in>)" << endl
	    << "    -x --suffixes           sort a random sample of the suffixes of the first input" << endl
	    << "    -u --upto <int>         sort up to <int> characters into each string" << endl
	    << "    -v --verbose            verbose output (for debugging)" << endl
	    << "    -s --sanity             enable sanity checking (increased runtime/memore usage)" << endl
	    << "    -q --quiet              output nothing unless there is a problem (for testing)" << endl
	    << "    -d --dcv                set difference-cover periodicity (and enable difference cover)" << endl;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static size_t parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			exit(1);
		}
		return (size_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	exit(1);
	return 0;
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
			case 'd': seqType = DNA;        break;
			case 'r': seqType = RNA;        break;
	   		case 'p': seqType = AMINO_ACID; break;
	
	   		//case 'k': format = PACKED;    break;
	   		case 'f': format = FASTA;     break;
	   		case 'g': format = GENBANK;   break;
	   		case 'e': format = EMBL;      break;
	   		case 'a': format = RAW;       break;
	   		case 'c': format = CMDLINE;   break;
	   		
			case 'v': verbose = 1;        break;
			case 'q': quiet = 1;          break;
	   		case 's': sanityCheck = true; break;
	   		case 'x': suffixes = true;    break;
	   		
	   		case 'u':
	   			upto = parseInt(1, "-u/--upto arg must be at least 1");
	   			break;
	   		case ARG_DCV:
	   			dcv  = parseInt(3, "--dcv arg must be at least 3");
	   			break;
	   		case ARG_SEED:
	   			seed = parseInt(0, "--seed arg must be at least 0");
	   			break;
			case -1: /* Done with options. */ break;
			default: 
				printUsage(cerr);
				exit(1);
		}
	} while(next_option != -1);
}

/**
 * 
 */
template<typename TStr>
static void driver(const char * type,
                   const string& infile,
                   const vector<string>& infiles)
{
	typedef String<uint32_t> TU32Str;
	typedef DifferenceCoverSample<TStr> TDC;
	typedef Value<TStr> TChar;
	vector<TStr> ss;
	if(verbose) {
		cout << "About to read text files (" << file_format_names[format] << ")" << endl;
	}
	switch(format) {
		case FASTA:   readSequenceFiles<TStr, Fasta>(infiles, ss);   break;
		case EMBL:    readSequenceFiles<TStr, Embl>(infiles, ss);    break;
	    case GENBANK: readSequenceFiles<TStr, Genbank>(infiles, ss); break;
		#ifdef PACKED_STRINGS
		case RAW:
			cerr << "RAW format not supported in packed-string mode" << endl;
			exit(1); break;
		#else
		case RAW:     readSequenceFiles<TStr, Raw>(infiles, ss);     break;
		#endif
		case CMDLINE: readSequenceString<TStr>(infile, ss);          break;
		//case PACKED:  readPackedSequenceFiles<TStr>(infiles, ss);    break;
		default: assert(false);
	}
	// Check that input is non-empty
	if(ss.size() == 0) {
		cerr << "Error: Empty input!  Check that file format is correct." << endl;
		exit(1);
	}
	if(verbose && !quiet) {
		cout << "Input strings:" << endl;
		for(size_t i = 0; i < ss.size(); i++) {
			cout << "  " << ss[i] << endl;
		}
	}
	if(suffixes) {
		// Select a random subset of the suffixes of ss[0] to sort
		uint32_t* sufs = new uint32_t[length(ss[0])];
		size_t sufslen = 0;
		uint32_t* idxs = new uint32_t[length(ss[0])];
		size_t idxslen = 0;
		for(uint32_t i = 0; i < length(ss[0]); i++) {
			if((random() % 20) >= 3) {
				sufs[sufslen++] = i;
			}
		}
		if(sufslen == 0) {
			sufs[sufslen++] = 0;
		}
		for(uint32_t i = 0; i < sufslen; i++) {
			idxs[idxslen++] = i;
		}
		if(dcv > 0) {
			TDC dc(ss[0], dcv, verbose, sanityCheck, cout);
			dc.build();
			mkeyQSortSufDc(ss[0],
			               sufs,
			               sufslen,
			               dc,
			               ValueSize<Dna>::VALUE,
			               verbose,
			               sanityCheck);
			if(!quiet) {
				cout << "Sorted strings:" << endl;
				printSuffixList(ss[0], sufs, sufslen, idxs, cout);
			}
		} else {
			mkeyQSortSuf2(ss[0],
			              sufs,
			              sufslen,
			              idxs,
			              ValueSize<Dna>::VALUE,
			              verbose,
			              sanityCheck,
			              upto);
			if(!quiet) {
				cout << "Sorted strings:" << endl;
				printSuffixList(ss[0], sufs, sufslen, idxs, cout);
			}
		}
	} else {
		// Convert ss to a String
		String<TStr> sss;
		for(size_t i = 0; i < ss.size(); i++) {
			append(sss, ss[i]);
		}
		// Sort all of the strings in ss
		mkeyQSort(sss, ValueSize<TChar>::VALUE, verbose, sanityCheck);
		if(!quiet) {
			cout << "Sorted strings:" << endl;
			printStringList(sss, cout);
		}
	}
}

int main(int argc, char** argv) {
	typedef String<Dna> TStr;
	string infile;
	vector<string> infiles;
	parseOptions(argc, argv);
	// Get input
	if(optind >= argc) {
		cerr << "No input sequence specified!" << endl;
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
	if(verbose) {
		cout << "Sanity checking: " << (sanityCheck? "enabled" : "disabled") << endl;
		cout << "Difference-cover period: " << dcv << endl;
		cout << "Random seed: " << seed << endl;
		cout << "Upto: " << upto << endl;
	}
	// Seed random number generator
	srandom(seed);
	switch(seqType) {
		case DNA: {
			Timer timer(cout, "Total time for call to driver(): ", verbose);
			#ifdef PACKED_STRINGS
			driver<String<Dna, Packed<Alloc<> > > >("DNA (packed)", infile, infiles);
			#else
			driver<String<Dna, Alloc<> > >("DNA", infile, infiles);
			#endif
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
