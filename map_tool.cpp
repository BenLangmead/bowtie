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
	FORMAT_BIN,
	FORMAT_CONCISE,
	FORMAT_FASTA,
	FORMAT_FASTQ
};

static bool verbose     = false;       // be talkative
static bool showVersion = false;       // show version info and exit
static bool refIdx      = false;       // print reference idxs instead of names
static int informat  = FORMAT_BIN;     // format of input alignments
static int outformat = FORMAT_DEFAULT; // format of output alignments

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: bowtie-maptool [options]* <align_in> [<align_out>]" << endl
	    << "    align_in         alignment file output by bowtie, or \"-\" for stdin" << endl
	    << "    align_out        write output alignments to this file (default: stdout)" << endl
	    << "Options:" << endl
	    << "    -d/--indef       <align_in> is default bowtie output" << endl
	    << "    -b/--inbin       <align_in> is bowtie -b/--binout output (default)" << endl
	    << "    -D/--outdef      <align_out> is default bowtie output  (default)" << endl
	    << "    -B/--outbin      <align_out> is -b/--binout bowtie output" << endl
	    << "    -C/--outconcise  <align_out> is --concise bowtie output" << endl
	    << "    -Q/--outfastq    <align_out> is aligned reads in FASTQ format" << endl
	    << "    -F/--outfasta    <align_out> is aligned reads in FASTA format" << endl
	    << "    -v/--verbose     verbose output (for debugging)" << endl
	    << "    -h/--help        print detailed description of tool and its options" << endl
	    << "    --version        print version information and quit" << endl
	    ;
}

static const char *short_options = "hvsdbDBCQF";

enum {
	ARG_VERSION = 256,
	ARG_REFIDX
};

static struct option long_options[] = {
	{"indef",      no_argument, 0, 'd'},
	{"inbin",      no_argument, 0, 'b'},
	{"outdef",     no_argument, 0, 'D'},
	{"outbin",     no_argument, 0, 'B'},
	{"outconcise", no_argument, 0, 'C'},
	{"outfastq",   no_argument, 0, 'Q'},
	{"outfasta",   no_argument, 0, 'F'},
	{"verbose",    no_argument, 0, 'v'},
	{"help",       no_argument, 0, 'h'},
	{"refidx",     no_argument, 0, ARG_REFIDX},
	{"version",    no_argument, 0, ARG_VERSION},
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
	   		case 'v': verbose = true; break;
	   		case ARG_VERSION: showVersion = true; break;
	   		case ARG_REFIDX: refIdx = true; break;
	   		case 'd': informat = FORMAT_DEFAULT; break;
	   		case 'b': informat = FORMAT_BIN; break;
	   		case 'D': outformat = FORMAT_DEFAULT; break;
	   		case 'B': outformat = FORMAT_BIN; break;
	   		case 'C': outformat = FORMAT_CONCISE; break;
	   		case 'Q': outformat = FORMAT_FASTQ; break;
	   		case 'F': outformat = FORMAT_FASTA; break;
			case -1: /* Done with options. */ break;
			case 0: if (long_options[option_index].flag != 0) break;
			default:
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
		}
	} while(next_option != -1);
}

/**
 * Print the read involved in an alignment as a FASTQ record.
 */
static void fastqAppend(ostream& out, Hit& h) {
	if(!h.fw) {
		reverseComplementInPlace(h.patSeq);
		reverseInPlace(h.quals);
	}
	out << "@" << h.patName << endl
	    << h.patSeq << endl
	    << "+" << endl
	    << h.quals << endl;
}

/**
 * Print the read involved in an alignment as a FASTA record.
 */
static void fastaAppend(ostream& out, Hit& h) {
	if(!h.fw) reverseComplementInPlace(h.patSeq);
	out << ">" << h.patName << endl << h.patSeq << endl;
}

/**
 * Parse command-line options, iterate through input alignments and
 * output appropriate converted alignment.
 */
int main(int argc, char **argv) {
	string infile;
	vector<string> infiles;
	string outfile;
	parseOptions(argc, argv);
	ostream *out = &cout;
	if(showVersion) {
		cout << argv[0] << " version " << BOWTIE_VERSION << endl;
		cout << "Built on " << BUILD_HOST << endl;
		cout << BUILD_TIME << endl;
		cout << "Compiler: " << COMPILER_VERSION << endl;
		cout << "Options: " << COMPILER_OPTIONS << endl;
		cout << "Sizeof {int, long, long long, void*, size_t}: {" << sizeof(int)
		     << ", " << sizeof(long) << ", " << sizeof(long long)
		     << ", " << sizeof(void *)
		     << ", " << sizeof(size_t) << "}" << endl;
		cout << "Source hash: " << EBWT_MAPTOOL_HASH << endl;
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
	} else if(outformat == FORMAT_BIN) {
		cerr << "If -B/--outbin is specified, <align_out> must also be specified" << endl;
		exit(1);
	}

	// Process all input files
	vector<string> refnames;
	for(size_t i = 0; i < infiles.size(); i++) {
		istream *inp;
		if(infiles[i] == "-") {
			inp = &cin;
		} else {
			inp = new ifstream(infiles[i].c_str(), ios_base::out | ios_base::binary);
		}
		istream& in = *inp;
		while(in.good() && !in.eof()) {
			Hit h;
			bool good;
			vector<string>* inrefnames = &refnames;
			vector<string>* outrefnames = NULL;
			if(!refIdx) {
				outrefnames = &refnames;
			}
			if(informat == FORMAT_BIN) {
				good = BinaryHitSink::readHit(h, in, inrefnames, verbose);
			} else {
				good = VerboseHitSink::readHit(h, in, inrefnames, verbose);
			}
			if(!good) continue; // bad alignment; skip it

			if(outformat == FORMAT_BIN) {
				BinaryHitSink::append(*out, h,outrefnames);
			} else if(outformat == FORMAT_DEFAULT) {
				VerboseHitSink::append(*out, h, outrefnames, 0 /* partition */);
			} else if(outformat == FORMAT_FASTQ) {
				fastqAppend(*out, h);
			} else if(outformat == FORMAT_FASTA) {
				fastaAppend(*out, h);
			} else {
				ConciseHitSink::append(*out, h, false /* reportOpps */);
			}
		}
		if(infiles[i] != "-") {
			((ifstream*)inp)->close();
			delete inp;
		}
	}

	// Close and delete output
	if(optind < argc) {
		((ofstream*)out)->close();
		delete out;
	}

	return 0;
}
