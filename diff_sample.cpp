#include <iostream>
#include <math.h>
#include <strings.h>
#include <seqan/sequence.h>
#include "diff_sample.h"
#include "diff_covers.h"

using namespace std;
using namespace seqan;

#ifdef DIFF_SAMPLE_MAIN

#include <getopt.h>
#include "params.h"
#include "sequence_io.h"
#include "tokenize.h"

static bool verbose          = false;
static bool quiet            = false;
static bool sanityCheck      = false;
static bool showColbournLing = false;
static bool exhaustive       = false;
static bool anchorMaps       = false;
static int format            = FASTA; // input sequence format
static int seqType           = DNA;   // input sequence alphabet

static const char *short_options = "drpfegvsceaq";

#define ARG_CL 256
#define ARG_EX 257
#define ARG_AN 258

static struct option long_options[] = {
	/* These options set a flag. */
	{"dna",        no_argument, &seqType, DNA},
	{"rna",        no_argument, &seqType, RNA},
	{"amino",      no_argument, &seqType, AMINO_ACID},
	{"protein",    no_argument, &seqType, AMINO_ACID},
	{"pro",        no_argument, &seqType, AMINO_ACID},
	{"verbose",    no_argument, 0, 'v'},
	{"quiet",      no_argument, 0, 'q'},
	{"sanity",     no_argument, 0, 's'},
	{"cl",         no_argument, 0, ARG_CL},
	{"exhaustive", no_argument, 0, ARG_EX},
	{"anchorMaps", no_argument, 0, ARG_AN},
	{0, 0, 0, 0}
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: diff_sample [options]* <int> [<int>] [<seq_in>]" << endl
	    << "  Calculates a difference cover for given periodicity." << endl
	    << "  If two ints are specified, will calcualte all difference covers in that range." << endl
	    << "    -d --dna                input uses 4-char DNA character set" << endl
	    << "    -r --rna                input uses 4-char RNA character set" << endl
	    << "    -p --protein            input uses 20+-char protein character set" << endl
	    //<< "    -k                      input file is packed (output of PackFasta)" << endl
	    << "    -f                      input file is Fasta (default)" << endl
	    << "    -e                      input file is Embl" << endl
	    << "    -a                      input file is raw text (single sequence)" << endl
	    << "    -g                      input file is Genbank" << endl
	    << "    -c                      input sequence is provided on command line (<seq_in>)" << endl
        << "    --cl                    show results of calculating Colbourn and Ling" << endl
        << "    --exhaustive            show results of calculating exhaustively" << endl
        << "    --anchorMaps            calculate anchor maps for difference-cover samples" << endl
	    << "    -v --verbose            verbose output (for debugging)" << endl
	    << "    -q --quiet              quiet output (only print warnings and errors)" << endl
	    << "    -s --sanity             enable sanity checks (much slower/increased memory usage)" << endl;
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
			case 'd':
				seqType = DNA;
				break;
			case 'r':
				seqType = RNA;
				break;
	   		case 'p':
	   			seqType = AMINO_ACID;
				break;
	
	   		//case 'k':
			//	format = PACKED;
			//	break;
	   		case 'f':
				format = FASTA;
				break;
	   		case 'g':
				format = GENBANK;
				break;
	   		case 'e':
				format = EMBL;
				break;
	   		case 'a':
				format = RAW;
				break;
	   		case 'c':
				format = CMDLINE;
				break;

	   		case 'v':
				verbose = true;
				break;
	   		case 'q':
				quiet = true;
				break;
	   		case 's':
				sanityCheck = true;
				break;
	   		case ARG_CL:
	   			showColbournLing = true;
				break;
	   		case ARG_EX:
	   			exhaustive = true;
				break;
	   		case ARG_AN:
	   			anchorMaps = true;
				break;
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

int main(int argc, char** argv) {
	typedef String<Dna> TStr;
	
	parseOptions(argc, argv);
	// Get first v
	if(optind >= argc) {
		cerr << "No periodicity specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	uint32_t v = atoi(argv[optind++]);
	uint32_t v2 = v;
	// Get optional second v
	if(optind < argc) {
		v2 = atoi(argv[optind++]);
	}
	if(verbose) { cout << "v: " << v << endl; }
	if(verbose && v != v2) { cout << "v2: " << v2 << endl; }

	TStr s;
	if(optind < argc) {
		vector<TStr> ss;
		string infile = argv[optind++];
		vector<string> infiles;
		tokenize(infile, ",", infiles);
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
			exit(1);
		}
		s = ss[0];
	}
	
	calcColbournAndLingDCs<uint32_t>(verbose, sanityCheck);
	if(exhaustive) {
		if(verbose) { cout << "Exhaustively calculating vs from " << v << " to " << v2 << endl; }
		for(size_t vv = v; vv <= v2; vv++) {
			calcExhaustiveDC(vv);
		}
	} else {
		for(size_t vv = v; vv <= v2; vv++) {
			String<uint32_t> dc = getDiffCover(vv, verbose, sanityCheck);
			if(empty(dc)) {
				cout << "Could not calculate difference cover for v=" << vv << endl;
				return 1;
			}
			if(!quiet) {
				cout << "v=" << vv << ", |D|=" << length(dc) << endl;
				// Print out difference cover (and optionally calculate
				// anchor map)
				for(size_t i = 0; i < length(dc); i++) {
					cout << dc[i];
					if(i < length(dc)-1) cout << ",";
				}
				cout << endl;
			}
			// Print out anchor map w/r/t left-hand offset
			if(anchorMaps && !quiet) {
				String<uint32_t> amap = getDeltaMap<uint32_t>(vv, dc);
				for(size_t i = 0; i < vv; i++) {
					assert_neq(0xffffffff, amap[i]);
					cout << amap[i];
					if(i < vv-1) cout << ",";
				}
				cout << endl;
			}
			
			if(!empty(s)) {
				DifferenceCoverSample<TStr> dcs(s, vv, verbose, sanityCheck, cout);
				dcs.build();
				if(!quiet) {
					dcs.print(cout);
				}
			}
		}
	}
	if(showColbournLing) {
		cout << "Colbourn and Ling DCs:" << endl;
		for(int i = 0; i < 16; i++) {
			cout << i << ": maxV=" << clDCs[i].maxV << ", numSamples=" << clDCs[i].numSamples << endl;
			cout << "  ";
			for(int j = 0; j < 128; j++) {
				cout << clDCs[i].samples[j];
				if(j < 127 && clDCs[i].samples[j+1] > 0) cout << ",";
				if(clDCs[i].samples[j+1] == 0) break;
			}
			cout << endl;
		}
	}
}
#endif
