#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <getopt.h>
#include "assert_helpers.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "timer.h"
#include "ref_read.h"
#include "hit.h"
#include "rot_buf.h"

/**
 * \file Driver for the bowtie-asm assembly (consensus) tool.
 */

// Build parameters
static bool verbose     = true;  // be talkative (default)
static int sanityCheck  = 0;     // do slow sanity checks
static bool showVersion = false;
static bool partitioned = true;

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: bowtie-asm [options]* <alignments_in> [<reference_in>]" << endl
	    << "    alignments_in        output from bowtie" << endl
	    << "    reference_in         write Ebwt data to files with this dir/basename" << endl
	    << "Options:" << endl
	    << "    -v/--verbose         verbose output (for debugging)" << endl
	    //<< "    -s/--sanity          enable sanity checks (much slower/increased memory usage)" << endl
	    << "    -h/--help            print detailed description of tool and its options" << endl
	    << "    --version            print version information and quit" << endl
	    ;
}

static const char *short_options = "hvs";

enum {
	ARG_VERSION = 256
};

static struct option long_options[] = {
	{"verbose", no_argument, 0, 'v'},
	{"sanity",  no_argument, 0, 's'},
	{"help",    no_argument, 0, 'h'},
	{"version", no_argument, 0, ARG_VERSION},
	{0, 0, 0, 0} // terminator
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

static char *argv0 = NULL;

/**
 *
 */
template<typename T>
static void processAlignment(const char *line,
                             AlignmentSink<T, 1024>& asink,
                             ColumnAnalyzer<T>& analyzer)
{
	string s(line);
	if(!s.empty() && s.find_first_of("\t", 0) != string::npos) {
		Hit h = VerboseHitSink::parseHit(s);
		asink.addAlignment(h, &analyzer); // Send the alignment to the sink
	}
}

/**
 * Sort alignments that are partitioned, where partitions are sorted
 * but alignments within partitions are not.
 */
template<typename T>
static void processPartitionedAlignments(
		istream& alfile,
        AlignmentSink<T, 1024>& asink,
        ColumnAnalyzer<T>& analyzer)
{
	char buf[4096];
	char partbuf[4096];
	vector<Hit> bucket;
	while(true) {
		alfile.getline(buf, 4096);
		if(alfile.eof()) break;
		if(alfile.bad()) {
			cerr << "Alignment file set \"bad\" bit" << endl;
			exit(1);
		}
		if(alfile.fail()) {
			cerr << "A line from the alignment file was longer than 4K" << endl;
			exit(1);
		}
		bool samePart = true;
		size_t pos = 0;
		// Assumes all partition labels are the same length
		while(buf[pos] != '\t') {
			if(samePart && partbuf[pos] != buf[pos]) {
				samePart = false;
			}
			if(!samePart) {
				partbuf[pos] = buf[pos];
			}
			pos++;
		}
		if(!samePart) {
			sort(bucket.begin(), bucket.end());
			// Sort the bucket and release it
			for(size_t i = 0; i < bucket.size(); i++) {
				asink.addAlignment(bucket[i], &analyzer);
			}
			bucket.clear();
		}
		string line(buf);
		if(!line.empty() && line.find_first_of("\t", 0) != string::npos) {
			Hit h = VerboseHitSink::parseHit(string(buf));
			bucket.push_back(h);
		}
	}
	// Sort and flush the bucket if necessary
	if(bucket.size() > 0) {
		sort(bucket.begin(), bucket.end());
		// Sort the bucket and release it
		for(size_t i = 0; i < bucket.size(); i++) {
			asink.addAlignment(bucket[i], &analyzer);
		}
	}
}

/**
 *
 */
template<typename T>
static void processSortedAlignments(istream& alfile,
                                    AlignmentSink<T, 1024>& asink,
                                    ColumnAnalyzer<T>& analyzer)
{
	char buf[4096];
	while(true) {
		alfile.getline(buf, 4096);
		if(alfile.eof()) break;
		if(alfile.bad()) {
			cerr << "Alignment file set \"bad\" bit" << endl;
			exit(1);
		}
		if(alfile.fail()) {
			cerr << "A line from the alignment file was longer than 4K" << endl;
			exit(1);
		}
		processAlignment(buf, asink, analyzer);
	}
}

/**
 * main function.  Parses command-line arguments.
 */
int main(int argc, char **argv) {

	string infile;
	vector<string> infiles;
	parseOptions(argc, argv);
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
		cout << "Source hash: " << EBWT_ASM_HASH << endl;
		return 0;
	}

	// Get input filename
	if(optind >= argc) {
		cerr << "No alignments file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	infile = argv[optind++];
	tokenize(infile, ",", infiles);

	for(size_t i = 0; i < infiles.size(); i++) {
		SNPColumnCharPairAnalyzer ca(std::cout);
		RotatingCharPairAlignmentBuf<1024> cpBuf;
		cpBuf.reset(0, &ca);
		if(partitioned) {
			if(infiles[i] == "-") {
				processPartitionedAlignments(cin, cpBuf, ca);
			} else {
				ifstream inin(infiles[i].c_str());
				if(!inin.good()) {
					cerr << "Error opening alignment file \""
					     << infiles[i] << "\" for reading" << endl;
				}
				processPartitionedAlignments(inin, cpBuf, ca);
			}
		} else {
			if(infiles[i] == "-") {
				processSortedAlignments(cin, cpBuf, ca);
			} else {
				ifstream inin(infiles[i].c_str());
				if(!inin.good()) {
					cerr << "Error opening alignment file \""
					     << infiles[i] << "\" for reading" << endl;
				}
				processSortedAlignments(inin, cpBuf, ca);
			}
		}
		cpBuf.finalize(&ca);
		if(verbose) {
			cout << "Finished processing alignment file " << infiles[i] << endl;
		}
	}

	return 0;
}
