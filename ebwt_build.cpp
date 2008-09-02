#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <getopt.h>
#include "assert_helpers.h"
#include "endian_swap.h"
#include "ebwt.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "timer.h"
#include "ref_read.h"

/**
 * \file Driver for the bowtie-build indexing tool.
 */

// Build parameters
static int verbose           = 0;     // be talkative
static int sanityCheck       = 0;     // do slow sanity checks
static int format            = FASTA; // input sequence format
static uint32_t bmax         = 0xffffffff; // max blockwise SA bucket size
static uint32_t bmaxMultSqrt = 0xffffffff; // same, as multplier of sqrt(n)
static uint32_t bmaxDivN     = 8;     // same, as divisor of n
static int dcv               = 1024;  // bwise SA difference-cover sample sz
static int noDc              = 0;     // disable difference-cover sample
static int entireSA          = 0;     // 1 = disable blockwise SA
static int seed              = 0;     // srandom seed
static int showVersion       = 0;     // just print version and quit?
static bool doubleEbwt       = true;  // build forward and reverse Ebwts
static int64_t cutoff        = 0xffffffff; // max # of reference bases
//   Ebwt parameters
static int32_t lineRate      = 6;  // a "line" is 64 bytes
static int32_t linesPerSide  = 1;  // 1 64-byte line on a side
static int32_t offRate       = 5;  // sample 1 out of 32 SA elts
static int32_t ftabChars     = 10; // 10 chars in initial lookup table
static int32_t chunkRate     = 11; // each chunk is 2K
static int32_t bigEndian     = 0;  // little endian

// Argument constants for getopts
static const int ARG_BMAX      = 256;
static const int ARG_BMAX_MULT = 257;
static const int ARG_BMAX_DIV  = 258;
static const int ARG_DCV       = 259;
static const int ARG_SEED      = 260;
static const int ARG_CUTOFF    = 261;

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>" << endl
	    << "    reference_in            comma-separated list of files with ref sequences" << endl
	    << "    ebwt_outfile_base       write Ebwt data to files with this dir/basename" << endl
	    << "Options:" << endl
	    << "    -f                      reference files are Fasta (default)" << endl
	    << "    -c                      reference sequences given on cmd line (as <seq_in>)" << endl
	    //<< "    -d/--double             build forward and reverse Ebwts for fast 1-mismatch" << endl
	    //<< "    --entiresa              build whole suffix array at once; huge mem footprint" << endl
	    << "    --bmax <int>            max SA bucket sz for blockwise suffix-array builder" << endl
	    << "    --bmaxmultsqrt <int>    max SA bucket sz as multiple of sqrt(ref len)" << endl
	    << "    --bmaxdivn <int>        max SA bucket sz as divisor of ref len (default: 8)" << endl
	    << "    --dcv <int>             diff-cover period for blockwise (default: 1024)" << endl
	    << "    --nodc                  disable difference cover (blockwise is quadratic)" << endl
	    //<< "    -l/--linerate <int>     line rate (single line is 2^rate bytes)" << endl
	    //<< "    -i/--linesperside <int> # lines in a side" << endl
	    << "    -o/--offrate <int>      SA index is kept every 2^offRate BWT chars" << endl
	    << "    -t/--ftabchars <int>    # of characters in initial lookup table key" << endl
	    << "    -h/--chunkrate <int>    # of characters in a text chunk" << endl
	    << "    --big --little          endianness (default: little, this host: "
	    << (currentlyBigEndian()? "big":"little") << ")" << endl
	    << "    --seed <int>            seed for random number generator" << endl
	    << "    --cutoff <int>          truncate reference at prefix of <int> bases" << endl
	    << "    -v/--verbose            verbose output (for debugging)" << endl
	    //<< "    -s/--sanity             enable sanity checks (much slower/increased memory usage)" << endl
	    << "    --version               print version information and quit" << endl
	    ;
}

static const char *short_options = "vrpscfl:i:o:t:h:";

static struct option long_options[] = {
	{"verbose",      no_argument,       0,            'v'},
	{"sanity",       no_argument,       0,            's'},
	//{"double",       no_argument,       0,            'd'},
	{"little",       no_argument,       &bigEndian,   0},
	{"big",          no_argument,       &bigEndian,   1},
	{"bmax",         required_argument, 0,            ARG_BMAX},
	{"bmaxmultsqrt", required_argument, 0,            ARG_BMAX_MULT},
	{"bmaxdivn",     required_argument, 0,            ARG_BMAX_DIV},
	{"dcv",          required_argument, 0,            ARG_DCV},
	{"nodc",         no_argument,       &noDc,        1},
	{"seed",         required_argument, 0,            ARG_SEED},
	{"entiresa",     no_argument,       &entireSA,    1},
	{"version",      no_argument,       &showVersion, 1},
	{"linerate",     required_argument, 0,            'l'},
	{"linesperside", required_argument, 0,            'i'},
	{"offrate",      required_argument, 0,            'o'},
	{"ftabchars",    required_argument, 0,            't'},
	{"chunkrate",    required_argument, 0,            'h'},
	{"cutoff",       required_argument, 0,            ARG_CUTOFF},
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
	   		case 'f': format = FASTA; break;
	   		case 'c': format = CMDLINE; break;
	   		//case 'd': doubleEbwt = true; break;
	   		case 'l':
   				lineRate = parseNumber<int>(3, "-l/--lineRate arg must be at least 3");
	   			break;
	   		case 'i':
	   			linesPerSide = parseNumber<int>(1, "-i/--linesPerSide arg must be at least 1");
	   			break;
	   		case 'o':
	   			offRate = parseNumber<int>(0, "-o/--offRate arg must be at least 0");
	   			break;
	   		case 't':
	   			ftabChars = parseNumber<int>(1, "-t/--ftabChars arg must be at least 1");
	   			break;
	   		case 'h':
	   			chunkRate = parseNumber<int>(1, "-h/--chunkRate arg must be at least 1");
	   			break;
	   		case ARG_BMAX:
	   			bmax = parseNumber<uint32_t>(1, "--bmax arg must be at least 1");
	   			break;
	   		case ARG_BMAX_MULT:
	   			bmaxMultSqrt = parseNumber<uint32_t>(1, "--bmaxMultSqrt arg must be at least 1");
	   			break;
	   		case ARG_BMAX_DIV:
	   			bmaxDivN = parseNumber<uint32_t>(1, "--bmaxDivN arg must be at least 1");
	   			break;
	   		case ARG_DCV:
	   			dcv = parseNumber<int>(3, "--dcv arg must be at least 3");
	   			break;
	   		case ARG_SEED:
	   			seed = parseNumber<int>(0, "--seed arg must be at least 0");
	   			break;
	   		case ARG_CUTOFF:
	   			cutoff = parseNumber<int64_t>(1, "--cutoff arg must be at least 1");
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
 * Drive the Ebwt construction process and optionally sanity-check the
 * result.
 */
template<typename TStr>
static void driver(const char * type,
                   const string& infile,
                   const vector<string>& infiles,
                   const string& outfile,
                   bool reverse = false)
{
	vector<istream*> is;
	RefReadInParams refparams(cutoff, -1, reverse);
	if(format == CMDLINE) {
		// Adapt sequence strings to stringstreams open for input
		for(size_t i = 0; i < infiles.size(); i++) {
			stringstream ss;
			ss << ">" << i << endl << infiles[i] << endl;
			is.push_back(new stringstream(ss.str(), ios::in));
		}
	} else {
		assert_eq(FASTA, format);
		// Adapt sequence files to ifstreams
		for(size_t i = 0; i < infiles.size(); i++) {
			is.push_back(new ifstream(infiles[i].c_str()));
		}
	}
	vector<uint32_t> szs;
	uint32_t sztot;
	{
		if(verbose) cout << "Reading reference sizes" << endl;
		Timer _t(cout, "  Time reading reference sizes: ", verbose);
		sztot = fastaRefReadSizes(is, szs, refparams);
	}
	assert_gt(sztot, 0);
	assert_gt(szs.size(), 0);
	// Construct Ebwt from input strings and parameters
	Ebwt<TStr> ebwt(lineRate,
	                linesPerSide,
	                offRate,      // suffix-array sampling rate
	                ftabChars,    // number of chars in initial arrow-pair calc
	                chunkRate,    // alignment
	                outfile,      // basename for .?.ebwt files
	                !entireSA,    // useBlockwise
	                bmax,         // block size for blockwise SA builder
                    bmaxMultSqrt, // block size as multiplier of sqrt(len)
                    bmaxDivN,     // block size as divisor of len
	                noDc? 0 : dcv,// difference-cover period
	                is,           // list of input streams
	                szs,          // list of reference sizes
	                sztot,        // total size of all references
	                refparams,    // reference read-in parameters
	                seed,         // pseudo-random number generator seed
	                -1,           // override offRate
	                verbose,      // be talkative
	                sanityCheck); // verify results and internal consistency
	// Note that the Ebwt is *not* resident in memory at this time.  To
	// load it into memory, call ebwt.loadIntoMemory()
	if(verbose) {
		// Print Ebwt's vital stats
		ebwt.eh().print(cout);
	}
	if(sanityCheck) {
		// Try restoring the original string (if there were
		// multiple texts, what we'll get back is the joined,
		// padded string, not a list)
		ebwt.loadIntoMemory();
		TStr s2; ebwt.restore(s2);
		ebwt.evictFromMemory();
		{
			TStr joinedss = Ebwt<TStr>::join(is, szs, sztot, refparams, ebwt.eh().chunkRate(), seed);
			assert_eq(length(joinedss), length(s2));
			assert_eq(joinedss, s2);
		}
		if(verbose) {
			if(length(s2) < 1000) {
				cout << "Passed restore check: " << s2 << endl;
			} else {
				cout << "Passed restore check: (" << length(s2) << " chars)" << endl;
			}
		}
	}
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
	argv0 = argv[0];
	if(showVersion) {
		cout << argv0 << " version " << BOWTIE_VERSION << endl;
		cout << "Built on " << BUILD_HOST << endl;
		cout << BUILD_TIME << endl;
		cout << "Compiler: " << COMPILER_VERSION << endl;
		cout << "Options: " << COMPILER_OPTIONS << endl;
		cout << "Sizeof {int, long, long long, void*}: {" << sizeof(int)
		     << ", " << sizeof(long) << ", " << sizeof(long long)
		     << ", " << sizeof(void *) << "}" << endl;
		cout << "Source hash: " << EBWT_BUILD_HASH << endl;
		return 0;
	}

	// Get input filename
	if(optind >= argc) {
		cerr << "No input sequence or sequence file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	infile = argv[optind++];

	// Get output filename
	if(optind >= argc) {
		cerr << "No output file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	outfile = argv[optind++];

	tokenize(infile, ",", infiles);
	if(infiles.size() < 1) {
		cerr << "Tokenized input file list was empty!" << endl;
		printUsage(cerr);
		return 1;
	}

	// Optionally summarize
	if(verbose) {
		cout << "Settings:" << endl
		     << "  Output file: \"" << outfile << "\"" << endl
		     << "  Line rate: " << lineRate << " (line is " << (1<<lineRate) << " bytes)" << endl
		     << "  Lines per side: " << linesPerSide << " (side is " << ((1<<lineRate)*linesPerSide) << " bytes)" << endl
		     << "  Offset rate: " << offRate << " (one in " << (1<<offRate) << ")" << endl
		     << "  FTable chars: " << ftabChars << endl
		     << "  Chunk rate: " << chunkRate << " (chunk is " << (1 << chunkRate) << " bytes)" << endl;
		if(bmax == 0xffffffff) {
			cout << "  Max bucket size: default" << endl;
		} else {
			cout << "  Max bucket size: " << bmax << endl;
		}
		if(bmaxMultSqrt == 0xffffffff) {
			cout << "  Max bucket size, sqrt multiplier: default" << endl;
		} else {
			cout << "  Max bucket size, sqrt multiplier: " << bmaxMultSqrt << endl;
		}
		if(bmaxDivN == 0xffffffff) {
			cout << "  Max bucket size, len divisor: default" << endl;
		} else {
			cout << "  Max bucket size, len divisor: " << bmaxDivN << endl;
		}
		cout << "  Difference-cover sample period: " << dcv << endl;
		cout << "  Endianness: " << (bigEndian? "big":"little") << endl
	         << "  Actual local endianness: " << (currentlyBigEndian()? "big":"little") << endl
		     << "  Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
#ifdef NDEBUG
		cout << "  Assertions: disabled" << endl;
#else
		cout << "  Assertions: enabled" << endl;
#endif
		cout << "  Random seed: " << seed << endl;
		cout << "  Sizeofs: void*:" << sizeof(void*) << ", int:" << sizeof(int) << ", long:" << sizeof(long) << endl;
		cout << "Input files (Dna), " << file_format_names[format] << "):" << endl;
		for(size_t i = 0; i < infiles.size(); i++) {
			cout << "  " << infiles[i] << endl;
		}
	}
	// Seed random number generator
	srand(seed);
	int64_t origCutoff = cutoff; // save cutoff since it gets modified
	{
		Timer timer(cout, "Total time for call to driver(): ", verbose);
		#ifdef PACKED_STRINGS
		driver<String<Dna, Packed<Alloc<> > > >("DNA (packed)", infile, infiles, outfile);
		#else
		driver<String<Dna, Alloc<> > >("DNA", infile, infiles, outfile);
		#endif
	}
	cutoff = origCutoff; // reset cutoff for backward Ebwt
	if(doubleEbwt) {
		srand(seed);
		Timer timer(cout, "Total time for backward call to driver(): ", verbose);
		#ifdef PACKED_STRINGS
		driver<String<Dna, Packed<Alloc<> > > >("DNA (packed)", infile, infiles, outfile + ".rev", true);
		#else
		driver<String<Dna, Alloc<> > >("DNA", infile, infiles, outfile + ".rev", true);
		#endif
	}
	return 0;
}
