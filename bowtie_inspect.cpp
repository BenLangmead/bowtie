#include <string>
#include <vector>
#include <iostream>
#include <getopt.h>
#include <seqan/find.h>

#include "endian_swap.h"
#include "ebwt.h"

using namespace std;
using namespace seqan;

static bool showVersion = false; // just print version and quit?
static int verbose     = 0;  // be talkative
static int names_only  = 0;  // just print the sequence names in the index
static int across      = 60; // number of characters across in FASTA output

static const char *short_options = "vh?na:";

static const int ARG_VERSION = 256;

static struct option long_options[] = {
	{(char*)"verbose", no_argument,       0, 'v'},
	{(char*)"version", no_argument,       0, ARG_VERSION},
	{(char*)"names",   no_argument,       0, 'n'},
	{(char*)"help",    no_argument,       0, 'h'},
	{(char*)"across",  required_argument, 0, 'a'},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out
	<< "Usage: bowtie-inspect [options]* <ebwt_base>" << endl
	<< "  <ebwt_base>        ebwt filename minus trailing .1.ebwt/.2.ebwt" << endl
	<< "Options:" << endl
	<< "  -a/--across <int>  Number of characters across in FASTA output (default: 60)" << endl
	<< "  -n/--names         Print reference sequence names only" << endl
	<< "  -v/--verbose       Verbose output (for debugging)" << endl
	<< "  -h/--help          print detailed description of tool and its options" << endl
	;
}

/**
 * Print a detailed usage message to the provided output stream.
 *
 * Manual text converted to C++ string with something like:
 * cat MANUAL  | tail -100 | sed -e 's/\"/\\\"/g' | \
 *   sed -e 's/^/"/' | sed -e 's/$/\\n"/'
 */
static void printLongUsage(ostream& out) {
	out <<
	"\n"
	" Using the 'bowtie-inspect' Index Inspector\n"
	" ------------------------------------------\n"
	"\n"
	" 'bowtie-inspect' extracts information from a Bowtie index about the\n"
	" original reference sequences that were used to build it.  By default,\n"
	" the tool will output a FASTA file containing the sequences of the\n"
	" original references (with all non-A/C/G/T characters converted to Ns).\n"
	" It can also be used to extract just the reference sequence names using\n"
	" the -n option.\n"
	"\n"
	"  Command Line\n"
	"  ------------\n"
	" \n"
	" Usage: bowtie-inspect [options]* <ebwt_base>\n"
	"\n"
	"  <ebwt_base>        The basename of the index to be inspected.  The\n"
	"                     basename is the name of any of the four index\n"
	"                     files up to but not including the first period.\n"
	"                     bowtie first looks in the current directory for\n"
	"                     the index files, then looks in the 'indexes'\n"
	"                     subdirectory under the directory where the\n"
	"                     currently-running 'bowtie' executable is located,\n"
	"                     then looks in the directory specified in the\n"
	"                     BOWTIE_INDEXES environment variable.\n"
	"\n"
	" Options:\n"
	"\n"
	"  -a/--across <int>  When printing FASTA output, output a newline\n"
	"                     character every <int> bases (default: 60).\n"
	"\n"
	"  -n/--names         Print reference sequence names only; ignore\n"
	"                     sequence.\n"
	"\n"
	"  -v/--verbose       Print verbose output (for debugging).\n"
	"\n"
	"  -h/--help          Print detailed description of tool and its options\n"
	"                     (from MANUAL).\n"
	"\n"
	;
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
    int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
		switch (next_option) {
	   		case 'h':
	   		case '?':
	   			printLongUsage(cout); exit(0); break;
	   		case 'v': verbose = true; break;
	   		case ARG_VERSION: showVersion = true; break;
			case 'n': names_only = true; break;
			case 'a': across = parseInt(-1, "-a/--across arg must be at least 1"); break;
			case -1: break; /* Done with options. */
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

void print_fasta_record(ostream& fout,
						const string& defline,
						const string& seq)
{
	fout << ">";
	fout << defline << endl;

	if(across > 0) {
		size_t i = 0;
		while (i + across < seq.length())
		{
			fout << seq.substr(i, across) << endl;
			i += across;
		}
		if (i < seq.length())
			fout << seq.substr(i) << endl;
	} else {
		fout << seq << endl;
	}
}

template<typename TStr>
void print_index_sequences(ostream& fout, Ebwt<TStr>& ebwt)
{
	vector<string>* refnames = &(ebwt.refnames());

	TStr cat_ref;
	ebwt.restore(cat_ref);

	uint32_t curr_ref = 0xffffffff;
	string curr_ref_seq = "";
	uint32_t curr_ref_len = 0xffffffff;
	uint32_t last_text_off = 0;
	size_t orig_len = seqan::length(cat_ref);
	uint32_t tlen = 0xffffffff;
	bool first = true;
	for(size_t i = 0; i < orig_len; i++) {
		uint32_t tidx = 0xffffffff;
		uint32_t textoff = 0xffffffff;
		tlen = 0xffffffff;

		ebwt.joinedToTextOff(1 /* qlen */, i, tidx, textoff, tlen);

		if (tidx != 0xffffffff && textoff < tlen)
		{
			if (curr_ref != tidx)
			{
				if (curr_ref != 0xffffffff)
				{
					// Add trailing gaps, if any exist
					if(curr_ref_seq.length() < curr_ref_len) {
						curr_ref_seq += string(curr_ref_len - curr_ref_seq.length(), 'N');
					}
					print_fasta_record(fout, (*refnames)[curr_ref], curr_ref_seq);
				}
				curr_ref = tidx;
				curr_ref_seq = "";
				curr_ref_len = tlen;
				last_text_off = 0;
				first = true;
			}

			uint32_t textoff_adj = textoff;
			if(first && textoff > 0) textoff_adj++;
			if (textoff_adj - last_text_off > 1)
				curr_ref_seq += string(textoff_adj - last_text_off - 1, 'N');

			curr_ref_seq.push_back(getValue(cat_ref,i));
			last_text_off = textoff;
			first = false;
		}
	}
	if (curr_ref < refnames->size())
	{
		// Add trailing gaps, if any exist
		if(curr_ref_seq.length() < curr_ref_len) {
			curr_ref_seq += string(curr_ref_len - curr_ref_seq.length(), 'N');
		}
		print_fasta_record(fout, (*refnames)[curr_ref], curr_ref_seq);
	}

}

static char *argv0 = NULL;

void print_index_sequence_names(const string& fname, ostream& fout)
{
	vector<string> p_refnames;
	string adjust = adjustEbwtBase(argv0, fname, verbose);
	readEbwtRefnames(adjust, p_refnames);
	for(size_t i = 0; i < p_refnames.size(); i++) {
		cout << p_refnames[i] << endl;
	}
}

typedef Ebwt<String<Dna, Packed<Alloc<> > > > TPackedEbwt;

static void driver(const string& ebwtFileBase, const string& query) {
	// Adjust
	string adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, verbose);

	if (names_only) {
		print_index_sequence_names(adjustedEbwtFileBase, cout);
	} else {
		// Initialize Ebwt object
		TPackedEbwt ebwt(adjustedEbwtFileBase, true, -1, -1, false, verbose);
	    // Load whole index into memory
		ebwt.loadIntoMemory();
		print_index_sequences(cout, ebwt);
		// Evict any loaded indexes from memory
		if(ebwt.isInMemory()) {
			ebwt.evictFromMemory();
		}
	}
}

/**
 * main function.  Parses command-line arguments.
 */
int main(int argc, char **argv) {
	string ebwtFile;  // read serialized Ebwt from this file
	string query;   // read query string(s) from this file
	vector<string> queries;
	string outfile; // write query results to this file
	argv0 = argv[0];
	parseOptions(argc, argv);
	if(showVersion) {
		cout << argv0 << " version " << BOWTIE_VERSION << endl;
		if(sizeof(void*) == 4) {
			cout << "32-bit" << endl;
		} else if(sizeof(void*) == 8) {
			cout << "64-bit" << endl;
		} else {
			cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
		}
		cout << "Built on " << BUILD_HOST << endl;
		cout << BUILD_TIME << endl;
		cout << "Compiler: " << COMPILER_VERSION << endl;
		cout << "Options: " << COMPILER_OPTIONS << endl;
		cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
		     << sizeof(int)
		     << ", " << sizeof(long) << ", " << sizeof(long long)
		     << ", " << sizeof(void *) << ", " << sizeof(size_t)
		     << ", " << sizeof(off_t) << "}" << endl;
		return 0;
	}

	// Get input filename
	if(optind >= argc) {
		cerr << "No index name given!" << endl;
		printUsage(cerr);
		return 1;
	}
	ebwtFile = argv[optind++];

	// Optionally summarize
	if(verbose) {
		cout << "Input ebwt file: \"" << ebwtFile << "\"" << endl;
		cout << "Output file: \"" << outfile << "\"" << endl;
		cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
#ifdef NDEBUG
		cout << "Assertions: disabled" << endl;
#else
		cout << "Assertions: enabled" << endl;
#endif
	}
	driver(ebwtFile, query);
	return 0;
}
