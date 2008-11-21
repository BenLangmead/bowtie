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
static bool verbose          = true;  // be talkative (default)
static int sanityCheck       = 0;     // do slow sanity checks
static int format            = FASTA; // input sequence format
static uint32_t bmax         = 0xffffffff; // max blockwise SA bucket size
static uint32_t bmaxMultSqrt = 0xffffffff; // same, as multplier of sqrt(n)
static uint32_t bmaxDivN     = 4;          // same, as divisor of n
static int dcv               = 1024;  // bwise SA difference-cover sample sz
static int noDc              = 0;     // disable difference-cover sample
static int entireSA          = 0;     // 1 = disable blockwise SA
static int seed              = 0;     // srandom seed
static int showVersion       = 0;     // just print version and quit?
static bool doubleEbwt       = true;  // build forward and reverse Ebwts
static int64_t cutoff        = -1;    // max # of reference bases
//   Ebwt parameters
static int32_t lineRate      = 6;  // a "line" is 64 bytes
static int32_t linesPerSide  = 1;  // 1 64-byte line on a side
static int32_t offRate       = 5;  // sample 1 out of 32 SA elts
static int32_t isaRate       = -1; // sample rate for ISA; default: don't sample
static int32_t ftabChars     = 10; // 10 chars in initial lookup table
//static int32_t chunkRate     = 12; // Now set automatically
static int     bigEndian     = 0;  // little endian
static bool    useBsearch    = true;
static bool    nsToAs        = false;
static bool    autoMem       = false;

// Argument constants for getopts
static const int ARG_BMAX      = 256;
static const int ARG_BMAX_MULT = 257;
static const int ARG_BMAX_DIV  = 258;
static const int ARG_DCV       = 259;
static const int ARG_SEED      = 260;
static const int ARG_CUTOFF    = 261;
static const int ARG_PMAP      = 262;
static const int ARG_ISARATE   = 263;
static const int ARG_NTOA      = 264;

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
	    //<< "    -n/--noblocks           disable blockwise - faster, uses more memory (default)" << endl
	    << "    -a/--auto               automatically pick bmax/dcv fitting available memory" << endl
	    << "    --bmax <int>            max bucket sz for blockwise suffix-array builder" << endl
	    //<< "    --bmaxmultsqrt <int>    max bucket sz as multiple of sqrt(ref len)" << endl
	    << "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)" << endl
	    << "    --dcv <int>             diff-cover period for blockwise (default: 1024)" << endl
	    << "    --nodc                  disable diff-cover (algorithm becomes quadratic)" << endl
	    //<< "    -l/--linerate <int>     line rate (single line is 2^rate bytes)" << endl
	    //<< "    -i/--linesperside <int> # lines in a side" << endl
	    << "    -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)" << endl
	    << "    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)" << endl
	    //<< "    --chunkrate <int>       # of characters in a text chunk" << endl
	    << "    --ntoa                  convert Ns in reference to As" << endl
	    << "    --big --little          endianness (default: little, this host: "
	    << (currentlyBigEndian()? "big":"little") << ")" << endl
	    << "    --seed <int>            seed for random number generator" << endl
	    << "    --cutoff <int>          truncate reference at prefix of <int> bases" << endl
	    << "    --oldpmap               use old reference mapping; yields larger index" << endl
	    << "    -q/--quiet              verbose output (for debugging)" << endl
	    //<< "    -s/--sanity             enable sanity checks (much slower/increased memory usage)" << endl
	    << "    -h/--help               print detailed description of tool and its options" << endl
	    << "    --version               print version information and quit" << endl
	    ;
}

/**
 * Print a detailed usage message to the provided output stream.
 *
 * Manual text converted to C++ string with something like:
 * cat MANUAL  | head -304 | tail -231 | sed -e 's/\"/\\\"/g' | \
 *   sed -e 's/^/"/' | sed -e 's/$/\\n"/'
 */
static void printLongUsage(ostream& out) {
	out <<
	" \n"
	" Using the 'bowtie-build' Indexer\n"
	" --------------------------------\n"
	" \n"
	" Use 'bowtie-build' to build a Bowtie index from a set of DNA\n"
	" sequences.  'bowtie-build' outputs a set of 4 files named\n"
	" NAME.1.ebwt, $NAME.2.ebwt, $NAME.rev.1.ebwt, and $NAME.rev.2.ebwt,\n"
	" where $NAME is <index_basename> supplied by the user on the command\n"
	" line.  These files together constitute the index: they are all that is\n"
	" needed to align reads to the reference sequences.  The sequence files\n"
	" themselves are not needed by Bowtie once the index is built.  \n"
	" \n"
	" The Bowtie index is based on the FM Index of Ferragina and Manzini,\n"
	" and the algorithm used to build it is the blockwise algorithm of\n"
	" Karkkainen.  For more information on these techniques, see these\n"
	" references:\n"
	" \n"
	"  * Ferragina, P. and Manzini, G. 2000. Opportunistic data structures\n"
	"    with applications. In Proceedings of the 41st Annual Symposium on\n"
	"    Foundations of Computer Science (November 12 - 14, 2000). FOCS\n"
	"  * Ferragina, P. and Manzini, G. 2001. An experimental study of an\n"
	"    opportunistic index. In Proceedings of the Twelfth Annual ACM-SIAM\n"
	"    Symposium on Discrete Algorithms (Washington, D.C., United States,\n"
	"    January 07 - 09, 2001). 269-278.\n"
	"  * Karkkainen, J. 2007. Fast BWT in small space by blockwise suffix\n"
	"    sorting. Theor. Comput. Sci. 387, 3 (Nov. 2007), 249-257\n"
	" \n"
	" bowtie-build allows the user to trade off between running time and\n"
	" memory footprint.  By default, bowtie-build will use an amount of\n"
	" memory equal to about 30% the size of the suffix array for the\n"
	" sequence: about 4 gigabytes for the human genome.\n"
	" \n"
	" The indexer provides many options pertaining to the \"shape\" of the\n"
	" index, e.g. --offrate governs the number of marked rows.  All of\n"
	" these options are potentially profitable trade-offs.  They have been\n"
	" set to defaults that are reasonable for most cases according to our\n"
	" experiments.\n"
	"\n"
	" If bowtie-build aborts with an out-of-memory error, try re-running it\n"
	" with a larger --bmaxdivn (the default is 4).\n"
	"\n"
	" Because bowtie-build uses 32-bit indexes internally, it can handle up\n"
	" to a maximum of 2^32-1 (somewhat more than 4 billion) characters in an\n"
	" index.  Because bowtie-build adds padding to separate distinct\n"
	" reference sequences in the index, bowtie-build actually supports fewer\n"
	" than 2^32-1 characters when there are multuple reference sequences.\n"
	" If your reference exceeds 2^32-1 characters either before or after\n"
	" padding, bowtie-build will print an error message and abort.  To\n"
	" resolve this, please divide your reference sequences into smaller\n"
	" batches and/or chunks and build a separate index for each.\n"
	"\n"
	"  Command Line\n"
	"  ------------\n"
	"\n"
	" Usage: bowtie-build [options]* <reference_in> <index_basename>\n"
	" \n"
	"    <reference_in>          A comma-separated list of files containing\n"
	"                            the reference sequences to be aligned to,\n"
	"                            or, if -c is specified, the sequences\n"
	"                            themselves. E.g., this might be\n"
	"                            \"chr1.fa,chr2.fa,chrX.fa,chrY.fa\", or, if\n"
	"                            -c is specified, this might be\n"
	"                            \"GGTCATCCT,ACGGGTCGT,CCGTTCTATGCGGCTTA\".\n"
	"\n"
	"    <index_basename>        The basename of the index files to write.\n"
	"\n"
	" Options:\n"
	"\n"
	"    -f                      The reference input files (specified as\n"
	"                            <reference_in>) are FASTA files (usually\n"
	"                            having extension .fa, .mfa, .fna or\n"
	"                            similar).\n"
	"\n"
	"    -c                      The reference sequences are given on the\n"
	"                            command line.  I.e. <reference_in> is a\n"
	"                            comma-separated list of sequences rather\n"
	"                            than a list of FASTA files.\n"
	"\n"
	"    -a                      Instead of taking a user-specified --bmax/\n"
	"                            --bmaxdivn and --dcv, bowtie-build\n"
	"                            automatically finds appropriate values for\n"
	"                            those parameters given the memory\n"
	"                            available.  bowtie-build will start at the\n"
	"                            specified --bmax/--bmaxdivn and --dcv (or\n"
	"                            the defaults, if not specified) and try\n"
	"                            more memory-economical settings until\n"
	"                            indexing completes with no out-of-memory\n"
	"                            exceptions.\n"
	"\n"
	"    --bmax <int>            The maximum number of suffixes allowed in a\n"
	"                            block.  Allowing more suffixes per block\n"
	"                            makes indexing faster, but increases memory\n"
	"                            overhead.  Overrides any previous\n"
	"                            specification of --bmax, or --bmaxdivn.\n"
	"                            Default: --bmaxdivn 4.\n"
	"\n"
	"    --bmaxdivn <int>        The maximum number of suffixes allowed in a\n"
	"                            block, expressed as a fraction of the\n"
	"                            length of the reference.  Overrides any\n"
	"                            previous specification of --bmax or\n"
	"                            --bmaxdivn. Default: --bmaxdivn 4.\n"
	"\n"
	"    --dcv <int>             Use <int> as the period for the difference-\n"
	"                            cover sample.  A larger period yields less\n"
	"                            memory overhead, but may make suffix\n"
	"                            sorting slower, especially if repeats are\n"
	"                            present.  Must be a power of 2 no greater\n"
	"                            than 4096.  Default: 1024.\n"
	"\n"
	"    --nodc                  Disable use of the difference-cover sample.\n"
	"                            Suffix sorting becomes quadratic worst-case\n"
	"                            time.\n"
	"\n"
	"    -o/--offrate <int>      To map alignments back to positions on the\n"
	"                            reference sequences, it's necessary to\n"
	"                            annotate (\"mark\") some or all of the\n"
	"                            Burrows-Wheeler rows with their\n"
	"                            corresponding location on the genome.  The\n"
	"                            offrate governs how many rows get marked:\n"
	"                            the indexer will mark every 2^<int> rows.\n"
	"                            Marking more rows makes reference-position\n"
	"                            lookups faster, but requires more memory to\n"
	"                            hold the annotations at runtime.  The\n"
	"                            default is 5 (every 32nd row is marked; for \n"
	"                            human genome, annotations occupy about 340\n"
	"                            megabytes).  \n"
	"\n"
	"    -t/--ftabchars <int>    The ftab is the lookup table used to\n"
	"                            calculate an initial Burrows-Wheeler range\n"
	"                            with respect to the first <int> characters\n"
	"                            of the query.  A larger <int> yields a\n"
	"                            larger lookup table but faster query times.\n"
	"                            The ftab has size 4^(<int>+1) bytes.  The\n"
	"                            default is 10 (ftab is 4MB).\n"
	"\n"
	"    --ntoa                  Convert Ns in the reference sequence to As\n"
	"                            before building the index.  By default, Ns\n"
	"                            are simply excluded from the index and\n"
	"                            'bowtie' will not find alignments that.\n"
	"                            overlap them.\n"
	"\n"
	"    --big --little          Endianness to use when serializing integers\n"
	"                            to the index file.  Default: little-endian\n"
	"                            (recommended for Intel- and AMD-based\n"
	"                            architectures).\n"
	"    \n"
	"    --seed <int>            Use <int> as the seed for pseudo-random\n"
	"                            number generator.\n"
	"    \n"
	"    --cutoff <int>          Index only the first <int> bases of the\n"
	"                            reference sequences (cumulative across\n"
	"                            sequences) and ignore the rest.\n"
	"    \n"
	"    -q/--quiet              ebwt-build is verbose by default.  With\n"
	"                            this option ebwt-build will print only\n"
	"                            error messages.\n"
	"    \n"
	"    -h/--help               Print detailed description of tool and its\n"
	"                            options (from MANUAL).\n"
	"    \n"
	"    --version               Print version information and quit.\n"
	" \n"
	;
}

static const char *short_options = "qraph?nscfl:i:o:t:h:";

static struct option long_options[] = {
	{"quiet",        no_argument,       0,            'q'},
	{"sanity",       no_argument,       0,            's'},
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
	{"auto",         no_argument,       0,            'a'},
	{"noblocks",     required_argument, 0,            'n'},
	{"linerate",     required_argument, 0,            'l'},
	{"linesperside", required_argument, 0,            'i'},
	{"offrate",      required_argument, 0,            'o'},
	{"isarate",      required_argument, 0,            ARG_ISARATE},
	{"ftabchars",    required_argument, 0,            't'},
	//{"chunkrate",    required_argument, 0,            'h'},
	{"help",         no_argument,       0,            'h'},
	{"cutoff",       required_argument, 0,            ARG_CUTOFF},
	{"ntoa",         no_argument,       0,            ARG_NTOA},
	{"oldpmap",      no_argument,       0,            ARG_PMAP},
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
	   		case ARG_ISARATE:
	   			isaRate = parseNumber<int>(0, "--isaRate arg must be at least 0");
	   			break;
	   		case 't':
	   			ftabChars = parseNumber<int>(1, "-t/--ftabChars arg must be at least 1");
	   			break;
//	   		case 'h':
//	   			chunkRate = parseNumber<int>(1, "-h/--chunkRate arg must be at least 1");
//	   			break;
	   		case 'n':
	   			// all f-s is used to mean "not set", so put 'e' on end
	   			bmax = 0xfffffffe;
	   			break;
	   		case 'h':
				printLongUsage(cout);
				exit(0);
				break;
	   		case '?':
				printUsage(cerr);
				exit(0);
				break;
	   		case ARG_BMAX:
	   			bmax = parseNumber<uint32_t>(1, "--bmax arg must be at least 1");
	   			bmaxMultSqrt = 0xffffffff; // don't use multSqrt
	   			bmaxDivN = 0xffffffff;     // don't use multSqrt
	   			break;
	   		case ARG_BMAX_MULT:
	   			bmaxMultSqrt = parseNumber<uint32_t>(1, "--bmaxmultsqrt arg must be at least 1");
	   			bmax = 0xffffffff;     // don't use bmax
	   			bmaxDivN = 0xffffffff; // don't use multSqrt
	   			break;
	   		case ARG_BMAX_DIV:
	   			bmaxDivN = parseNumber<uint32_t>(1, "--bmaxdivn arg must be at least 1");
	   			bmax = 0xffffffff;         // don't use bmax
	   			bmaxMultSqrt = 0xffffffff; // don't use multSqrt
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
	   		case ARG_NTOA: nsToAs = true; break;
	   		case ARG_PMAP: useBsearch = false; break;
	   		case 'a': autoMem = true; break;
	   		case 'q': verbose = false; break;
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
	if(bmax < 40) {
		cerr << "Warning: specified bmax is very small (" << bmax << ").  This can lead to" << endl
		      << "extremely slow performance and memory exhaustion.  Perhaps you meant to specify" << endl
		      << "a small --bmaxdivn?" << endl;
	}
}

/**
 * Drive the Ebwt construction process and optionally sanity-check the
 * result.
 */
template<typename TStr>
static void driver(const char * type,
                   const string& infile,
                   vector<string>& infiles,
                   const string& outfile,
                   bool reverse = false)
{
	vector<FileBuf*> is;
	RefReadInParams refparams(cutoff, -1, reverse, nsToAs);
	assert_gt(infiles.size(), 0);
	if(format == CMDLINE) {
		// Adapt sequence strings to stringstreams open for input
		stringstream *ss = new stringstream();
		for(size_t i = 0; i < infiles.size(); i++) {
			(*ss) << ">" << i << endl << infiles[i] << endl;
		}
		FileBuf *fb = new FileBuf(ss);
		assert(fb != NULL);
		assert(!fb->eof());
		assert(fb->get() == '>');
		ASSERT_ONLY(fb->reset());
		assert(!fb->eof());
		is.push_back(fb);
	} else {
		// Adapt sequence files to ifstreams
		for(size_t i = 0; i < infiles.size(); i++) {
			FILE *f = fopen(infiles[i].c_str(), "r");
			FileBuf *fb = new FileBuf(f);
			assert(fb != NULL);
			assert(!fb->eof());
			assert(fb->get() == '>');
			ASSERT_ONLY(fb->reset());
			assert(!fb->eof());
			is.push_back(fb);
		}
	}
	vector<RefRecord> szs;
	uint32_t sztot;
	int32_t chunkRate = 0;
	{
		if(verbose) cout << "Reading reference sizes" << endl;
		Timer _t(cout, "  Time reading reference sizes: ", verbose);
		sztot = fastaRefReadSizes(is, szs, refparams);
		if(!useBsearch) {
			chunkRate = EbwtParams::calcBestChunkRate(szs, offRate, lineRate, linesPerSide);
			if(verbose) cout << "  Choose best chunkRate: " << chunkRate << endl;
		}
	}
	assert_gt(sztot, 0);
	assert_gt(szs.size(), 0);
	// Construct Ebwt from input strings and parameters
	Ebwt<TStr> ebwt(lineRate,
	                linesPerSide,
	                offRate,      // suffix-array sampling rate
	                isaRate,      // ISA sampling rate
	                ftabChars,    // number of chars in initial arrow-pair calc
	                useBsearch? -1 : chunkRate, // alignment
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
	                -1,           // override isaRate
	                verbose,      // be talkative
	                autoMem,      // pass exceptions up to the toplevel so that we can adjust memory settings automatically
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
		cout << "Sizeof {int, long, long long, void*, size_t}: {" << sizeof(int)
		     << ", " << sizeof(long) << ", " << sizeof(long long)
		     << ", " << sizeof(void *)
		     << ", " << sizeof(size_t) << "}" << endl;
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
		     << "  Output files: \"" << outfile << ".*.ebwt\"" << endl
		     << "  Line rate: " << lineRate << " (line is " << (1<<lineRate) << " bytes)" << endl
		     << "  Lines per side: " << linesPerSide << " (side is " << ((1<<lineRate)*linesPerSide) << " bytes)" << endl
		     << "  Offset rate: " << offRate << " (one in " << (1<<offRate) << ")" << endl
		     << "  FTable chars: " << ftabChars << endl
		     ;
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
		cout << "  Reference base cutoff: ";
		if(cutoff == -1) cout << "none"; else cout << cutoff << " bases";
		cout << endl;
		cout << "  Endianness: " << (bigEndian? "big":"little") << endl
	         << "  Actual local endianness: " << (currentlyBigEndian()? "big":"little") << endl
		     << "  Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
#ifdef NDEBUG
		cout << "  Assertions: disabled" << endl;
#else
		cout << "  Assertions: enabled" << endl;
#endif
		cout << "  Random seed: " << seed << endl;
		cout << "  Sizeofs: void*:" << sizeof(void*) << ", int:" << sizeof(int) << ", long:" << sizeof(long) << ", size_t:" << sizeof(size_t) << endl;
		cout << "Input files DNA, " << file_format_names[format] << ":" << endl;
		for(size_t i = 0; i < infiles.size(); i++) {
			cout << "  " << infiles[i] << endl;
		}
	}
	// Seed random number generator
	srand(seed);
	int64_t origCutoff = cutoff; // save cutoff since it gets modified
	{
		Timer timer(cout, "Total time for call to driver() for forward index: ", verbose);
		#ifdef PACKED_STRINGS
		driver<String<Dna, Packed<Alloc<> > > >("DNA (packed)", infile, infiles, outfile);
		#else
		driver<String<Dna, Alloc<> > >("DNA", infile, infiles, outfile);
		#endif
	}
	cutoff = origCutoff; // reset cutoff for backward Ebwt
	if(doubleEbwt) {
		srand(seed);
		Timer timer(cout, "Total time for backward call to driver() for mirror index: ", verbose);
		#ifdef PACKED_STRINGS
		driver<String<Dna, Packed<Alloc<> > > >("DNA (packed)", infile, infiles, outfile + ".rev", true);
		#else
		driver<String<Dna, Alloc<> > >("DNA", infile, infiles, outfile + ".rev", true);
		#endif
	}
	return 0;
}
