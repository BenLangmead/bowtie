#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include "assert_helpers.h"
#include "tokenize.h"
#include "hit.h"

enum {
	FORMAT_DEFAULT = 1,
	FORMAT_BIN
};

// Build parameters
static bool verbose     = true;  // be talkative (default)
static bool sanityCheck = false; // do slow sanity checks
static bool showVersion = false; //
static int informat  = FORMAT_BIN;
static int outformat = FORMAT_DEFAULT;

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: bowtie-maptool [options]* <align_in> [<align_out>]" << endl
	    << "    align_in         alignments output by bowtie" << endl
	    << "    align_out        write output alignments to this file (default: stdout)" << endl
	    << "Options:" << endl
	    << "    -d/--indef       <align_in> is default bowtie output" << endl
	    << "    -b/--inbin       <align_in> is bowtie -b/--binout output (default)" << endl
	    << "    -D/--outdef      <align_out> will be default bowtie output  (default)" << endl
	    << "    -B/--outbin      <align_out> will be -b/--binout bowtie output" << endl
	    << "    -v/--verbose     verbose output (for debugging)" << endl
	  //<< "    -s/--sanity      enable sanity checks (much slower/increased memory usage)" << endl
	    << "    -h/--help        print detailed description of tool and its options" << endl
	    << "    --version        print version information and quit" << endl
	    ;
}

static const char *short_options = "hvsdbDB";

enum {
	ARG_VERSION = 256
};

static struct option long_options[] = {
	{"indef",   no_argument, 0, 'd'},
	{"inbin",   no_argument, 0, 'b'},
	{"outdef",  no_argument, 0, 'D'},
	{"outbin",  no_argument, 0, 'B'},
	{"verbose", no_argument, 0, 'v'},
	{"sanity",  no_argument, 0, 's'},
	{"help",    no_argument, 0, 'h'},
	{"version", no_argument, 0, ARG_VERSION},
	{0, 0, 0, 0} // terminator
};

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, char **argv) {
    int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
		switch (next_option) {
	   		case 'h':
				//printLongUsage(cout);
				printUsage(cout);
				exit(0);
				break;
	   		case '?':
				printUsage(cerr);
				exit(0);
				break;
	   		case 'v': verbose = false; break;
	   		case 's': sanityCheck = true; break;
	   		case ARG_VERSION: showVersion = true; break;
	   		case 'd': informat = FORMAT_DEFAULT; break;
	   		case 'b': informat = FORMAT_BIN; break;
	   		case 'D': outformat = FORMAT_DEFAULT; break;
	   		case 'B': outformat = FORMAT_BIN; break;
			case -1: /* Done with options. */ break;
			case 0: if (long_options[option_index].flag != 0) break;
			default:
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
		}
	} while(next_option != -1);
}

static char *argv0 = NULL;

/**
 * main function.  Parses command-line arguments.
 */
int main(int argc, char **argv) {
	string infile;
	vector<string> infiles;
	string outfile;
	parseOptions(argc, argv);
	ostream *out = &cout;
	argv0 = argv[0];
	if(showVersion) {
		cout << argv0 << " version " << BOWTIE_VERSION << endl;
		cout << "Built on " << BUILD_HOST << endl;
		cout << BUILD_TIME << endl;
		cout << "Compiler: " << COMPILER_VERSION << endl;
		cout << "Options: " << COMPILER_OPTIONS << endl;
		cout << "Sizeof {int, long, long long, void*, size_t}: {" << sizeof(int)
		     << ", " << sizeof(long) << ", " << sizeof(long long)
		     << ", " << sizeof(void *)
		     << ", " << sizeof(size_t) << "}" << endl;
		cout << "Source hash: " << EBWT_SEARCH_HASH << endl;
		return 0;
	}

	// Get input filename
	if(optind >= argc) {
		cerr << "No input alignments file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	infile = argv[optind++];
	tokenize(infile, ",", infiles);

	// Get output filename
	if(optind < argc) {
		out = new ofstream(argv[optind]);
	}

	for(size_t i = 0; i < infiles.size(); i++) {
		Hit h;
		cout << "Processing " << infiles[i] << endl;
		ifstream in(infiles[i].c_str(), ios_base::out | ios_base::binary);
		BinaryHitSink::readHit(h, in, verbose);
		in.close();
	}
	(*out) << "Output" << endl;

	// Close and delete output
	if(optind < argc) {
		((ofstream*)out)->close();
		delete out;
	}

	return 0;
}
