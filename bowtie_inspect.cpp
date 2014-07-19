#include <string>
#include <vector>
#include <iostream>
#include <getopt.h>
#include <stdexcept>
#include <seqan/find.h>

#include "assert_helpers.h"
#include "endian_swap.h"
#include "ebwt.h"
#include "reference.h"

using namespace std;
using namespace seqan;

static bool showVersion = false; // just print version and quit?
int verbose             = 0;  // be talkative
static int names_only   = 0;  // just print the sequence names in the index
static int summarize_only = 0; // just print summary of index and quit
static int across       = 60; // number of characters across in FASTA output
static bool extra       = false; // print extra summary info
static bool exclAllGaps = false; // print extra summary info
static bool refFromEbwt = false; // true -> when printing reference, decode it from Ebwt instead of reading it from BitPairReference
static string wrapper;
static const char *short_options = "vhnsea:";

enum {
	ARG_VERSION = 256,
	ARG_USAGE,
	ARG_EXTRA,
	ARG_EXCL_AMBIG,
	ARG_WRAPPER
};

static struct option long_options[] = {
	{(char*)"verbose",  no_argument,        0, 'v'},
	{(char*)"version",  no_argument,        0, ARG_VERSION},
	{(char*)"usage",    no_argument,        0, ARG_USAGE},
	{(char*)"extra",    no_argument,        0, ARG_EXTRA},
	{(char*)"excl-ambig",no_argument,       0, ARG_EXCL_AMBIG},
	{(char*)"names",    no_argument,        0, 'n'},
	{(char*)"summary",  no_argument,        0, 's'},
	{(char*)"help",     no_argument,        0, 'h'},
	{(char*)"across",   required_argument,  0, 'a'},
	{(char*)"ebwt-ref", no_argument,        0, 'e'},
	{(char*)"wrapper",  required_argument,  0, ARG_WRAPPER},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out
	<< "Usage: bowtie-inspect [options]* <ebwt_base>" << endl
	<< "  <ebwt_base>        ebwt filename minus trailing .1." + gEbwt_ext + "/.2." + gEbwt_ext << endl
	<< endl
	<< "  By default, prints FASTA records of the indexed nucleotide sequences to" << endl
	<< "  standard out.  With -n, just prints names.  With -s, just prints a summary of" << endl
	<< "  the index parameters and sequences.  With -e, preserves colors if applicable." << endl
	<< endl
	<< "Options:" << endl;
	if(wrapper == "basic-0") {
		out << "  --large-index      force inspection of the 'large' index, even if a" << endl
			<< "                     'small' one is present." << endl;
	}
	out << "  -a/--across <int>  Number of characters across in FASTA output (default: 60)" << endl
	<< "  -n/--names         Print reference sequence names only" << endl
	<< "  -s/--summary       Print summary incl. ref names, lengths, index properties" << endl
	<< "  -e/--ebwt-ref      Reconstruct reference from ebwt (slow, preserves colors)" << endl
	<< "  -v/--verbose       Verbose output (for debugging)" << endl
	<< "  -h/--help          print detailed description of tool and its options" << endl
	<< "  --help             print this usage message" << endl
	;
	if(wrapper.empty()) {
		cerr << endl
		     << "*** Warning ***" << endl
			 << "'boowtie-inspect' was run directly.  It is recommended "
			 << "to use the wrapper script instead."
			 << endl << endl;
	}
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
			throw 1;
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
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
			case ARG_WRAPPER:
				wrapper = optarg;
				break;
			case ARG_USAGE:
			case 'h':
				printUsage(cout);
				throw 0;
				break;
			case 'v': verbose = true; break;
			case ARG_VERSION: showVersion = true; break;
			case ARG_EXCL_AMBIG: exclAllGaps = true; break;
			case ARG_EXTRA: extra = true; break;
			case 'e': refFromEbwt = true; break;
			case 'n': names_only = true; break;
			case 's': summarize_only = true; break;
			case 'a': across = parseInt(-1, "-a/--across arg must be at least 1"); break;
			case -1: break; /* Done with options. */
			case 0:
				if (long_options[option_index].flag != 0)
					break;
			default:
				printUsage(cerr);
				throw 1;
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

/**
 * Given output stream, name and length, print a string of Ns with the
 * appropriate number of columns.
 */
void print_alln_ref_sequence(
	ostream& fout,
	const string& name,
	size_t len)
{
	fout << ">" << name << "\n";
	size_t j = 0;
	for(size_t i = 0; i < len; i += across) {
		while(j < len && j < i+across) {
			fout << 'N';
			j++;
		}
		fout << "\n";
	}
}

/**
 * Given output stream, BitPairReference, reference index, name and
 * length, print the whole nucleotide reference with the appropriate
 * number of columns.
 */
void print_ref_sequence(
	ostream& fout,
	BitPairReference& ref,
	const string& name,
	size_t refi,
	size_t len)
{
	bool newlines = across > 0;
	int myacross = across > 0 ? across : 60;
	size_t incr = myacross * 1000;
	uint32_t *buf = new uint32_t[(incr + 128)/4];
	fout << ">" << name << "\n";
	for(size_t i = 0; i < len; i += incr) {
		size_t amt = min(incr, len-i);
		assert_leq(amt, incr);
		int off = ref.getStretch(buf, refi, i, amt);
		uint8_t *cb = ((uint8_t*)buf) + off;
		for(size_t j = 0; j < amt; j++) {
			if(newlines && j > 0 && (j % myacross) == 0) fout << "\n";
			assert_range(0, 4, (int)cb[j]);
			fout << "ACGTN"[(int)cb[j]];
		}
		fout << "\n";
	}
	delete buf;
}

/**
 * Create a BitPairReference encapsulating the reference portion of the
 * index at the given basename.  Iterate through the reference
 * sequences, sending each one to print_ref_sequence to print.
 */
void print_ref_sequences(
	ostream& fout,
	bool color,
	const vector<string>& refnames,
	const TIndexOffU* plen,
	const string& adjustedEbwtFileBase)
{
	BitPairReference ref(
		adjustedEbwtFileBase, // input basename
		color,                // true -> expect colorspace reference
		false,                // sanity-check reference
		NULL,                 // infiles
		NULL,                 // originals
		false,                // infiles are sequences
		true,                 // load sequence
		false,                // memory-map
		false,                // use shared memory
		false,                // sweep mm-mapped ref
		verbose,              // be talkative
		verbose);             // be talkative at startup
#ifdef ACCOUNT_FOR_ALL_GAP_REFS
	for(size_t i = 0; i < ref.numRefs(); i++) {
		if(ref.isAllGaps(i) && !exclAllGaps) {
			print_alln_ref_sequence(
				fout,
				refnames[i],
				ref.len(i));
		} else {
			print_ref_sequence(
				fout,
				ref,
				refnames[i],
				ref.shrinkIdx(i),
				ref.len(i));
		}
	}
#else
	assert_eq(refnames.size(), ref.numNonGapRefs());
	for(size_t i = 0; i < ref.numNonGapRefs(); i++) {
		print_ref_sequence(
			fout,
			ref,
			refnames[i],
			i,
			plen[i] + (color ? 1 : 0));
	}
#endif
}

/**
 * Given an index, reconstruct the reference by LF mapping through the
 * entire thing.
 */
template<typename TStr>
void print_index_sequences(
	ostream& fout,
	Ebwt<TStr>& ebwt,
	const BitPairReference& refs)
{
	vector<string>* refnames = &(ebwt.refnames());

	TStr cat_ref;
	ebwt.restore(cat_ref);

	TIndexOffU curr_ref = OFF_MASK;
	string curr_ref_seq = "";
	TIndexOffU curr_ref_len = OFF_MASK;
	uint32_t last_text_off = 0;
	size_t orig_len = seqan::length(cat_ref);
	TIndexOffU tlen = OFF_MASK;
	bool first = true;
	for(size_t i = 0; i < orig_len; i++) {
		TIndexOffU tidx = OFF_MASK;
		TIndexOffU textoff = OFF_MASK;
		tlen = OFF_MASK;

		ebwt.joinedToTextOff(1 /* qlen */, (TIndexOffU)i, tidx, textoff, tlen);

		if (tidx != OFF_MASK && textoff < tlen)
		{
			if (curr_ref != tidx)
			{
				if (curr_ref != OFF_MASK)
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

			TIndexOffU textoff_adj = textoff;
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
	readEbwtRefnames(fname, p_refnames);
	for(size_t i = 0; i < p_refnames.size(); i++) {
		cout << p_refnames[i] << endl;
	}
}

typedef Ebwt<String<Dna, Packed<Alloc<> > > > TPackedEbwt;

/**
 * Print a short summary of what's in the index and its flags.
 */
void print_index_summary(
	const string& fname,
	ostream& fout,
	const BitPairReference& refs)
{
	int32_t flags = readFlags(fname);
	int32_t flagsr = readFlags(fname + ".rev");
	bool color = readEbwtColor(fname);
	bool entireReverse = readEntireReverse(fname + ".rev");
	TPackedEbwt ebwt(
		fname,
		color,                // index is colorspace
		-1,                   // don't require entire reverse
		true,                 // index is for the forward direction
		-1,                   // offrate (-1 = index default)
		-1,
		false,                // use memory-mapped IO
		false,                // use shared memory
		false,                // sweep memory-mapped memory
		true,                 // load names?
		//false,                // load SA sample?
		NULL,                 // no reference map
		verbose,              // be talkative?
		verbose,              // be talkative at startup?
		false,                // pass up memory exceptions?
		false);               // sanity check?
	vector<string> p_refnames;
	readEbwtRefnames(fname, p_refnames);
	if(extra) {
		cout << "Flags" << '\t' << (-flags) << endl;
		cout << "Reverse flags" << '\t' << (-flagsr) << endl;
	}
	cout << "Colorspace" << '\t' << (color ? "1" : "0") << endl;
	if(extra) {
		cout << "Concat then reverse" << '\t' << (entireReverse ? "1" : "0") << endl;
		cout << "Reverse then concat" << '\t' << (entireReverse ? "0" : "1") << endl;
		cout << "nPat" << '\t' << ebwt.nPat() << endl;
		cout << "refnames.size()" << '\t' << p_refnames.size() << endl;
		cout << "refs.numRefs()" << '\t' << refs.numRefs() << endl;
		cout << "refs.numNonGapRefs()" << '\t' << refs.numNonGapRefs() << endl;
	}
	cout << "SA-Sample" << "\t1 in " << (1 << ebwt.eh().offRate()) << endl;
	cout << "FTab-Chars" << '\t' << ebwt.eh().ftabChars() << endl;
	for(size_t i = 0; i < ebwt.nPat(); i++) {
		cout << "Sequence-" << (i+1)
		     << '\t' << p_refnames[refs.expandIdx((uint32_t)i)]
		     << '\t' << (ebwt.plen()[i] + (color ? 1 : 0))
		     << endl;
	}
	if(extra) {
		cout << "RefRecords:\n";
		for(size_t i = 0; i < refs.refRecords().size(); i++) {
			RefRecord r = refs.refRecords()[i];
			cout << r.first << "\t(" << r.off << ", " << r.len << ")" << endl;
		}
	}
}

static void driver(
	const string& ebwtFileBase,
	const string& query)
{
	// Adjust
	string adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, verbose);
	if (names_only) {
		print_index_sequence_names(adjustedEbwtFileBase, cout);
		return;
	}
	bool color = readEbwtColor(adjustedEbwtFileBase);
	BitPairReference refs(
		adjustedEbwtFileBase,
		color,
		false,
		NULL,
		NULL,
		false,
		false, // don't load sequence (yet)
		false,
		false,
		false, // mmSweep
		verbose,
		verbose);
	if(summarize_only) {
		print_index_summary(adjustedEbwtFileBase, cout, refs);
	} else {
		// Initialize Ebwt object
		TPackedEbwt ebwt(
			adjustedEbwtFileBase,
			color,                // index is colorspace
			-1,                   // don't care about entire-reverse
			true,                 // index is for the forward direction
			-1,                   // offrate (-1 = index default)
			-1,
			false,                // use memory-mapped IO
			false,                // use shared memory
			false,                // sweep memory-mapped memory
			true,                 // load names?
			//true,                 // load SA sample?
			NULL,                 // no reference map
			verbose,              // be talkative?
			verbose,              // be talkative at startup?
			false,                // pass up memory exceptions?
			false);               // sanity check?
		// Load whole index into memory
		if(refFromEbwt) {
			ebwt.loadIntoMemory(-1, -1, true, false);
			print_index_sequences(cout, ebwt, refs);
		} else {
			vector<string> refnames;
			readEbwtRefnames(adjustedEbwtFileBase, refnames);
			print_ref_sequences(
				cout,
				readEbwtColor(ebwtFileBase),
				refnames,
				ebwt.plen(),
				adjustedEbwtFileBase);
		}
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
	try {
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
	} catch(std::exception& e) {
		cerr << "Command: ";
		for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
		cerr << endl;
		return 1;
	} catch(int e) {
		if(e != 0) {
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
		}
		return e;
	}
}
