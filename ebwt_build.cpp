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
#include "filebuf.h"
#include "reference.h"

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
static bool    autoMem       = true;
static bool    packed        = false;
static bool    writeRef      = true;  // write compact reference to .3.ebwt/.4.ebwt
static bool    justRef       = false; // *just* write compact reference, don't index

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
	    << "    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting" << endl
	    << "    -p/--packed             use packed strings internally; slower, uses less mem" << endl
	    << "    --bmax <int>            max bucket sz for blockwise suffix-array builder" << endl
	    //<< "    --bmaxmultsqrt <int>    max bucket sz as multiple of sqrt(ref len)" << endl
	    << "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)" << endl
	    << "    --dcv <int>             diff-cover period for blockwise (default: 1024)" << endl
	    << "    --nodc                  disable diff-cover (algorithm becomes quadratic)" << endl
	    << "    -r/--noref              don't build .3/.4.ebwt (packed reference) portion" << endl
	    << "    -3/--justref            just build .3/.4.ebwt (packed reference) portion" << endl
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
	"\n"
	" Using the 'bowtie-build' Indexer\n"
	" --------------------------------\n"
	"\n"
	" Use 'bowtie-build' to build a Bowtie index from a set of DNA\n"
	" sequences.  bowtie-build outputs a set of 6 files with suffixes\n"
	" .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt,\n"
	" where the prefix is the <ebwt_outfile_base> parameter supplied by the\n"
	" user on the command line.  These files together constitute the index:\n"
	" they are all that is needed to align reads to the reference sequences.\n"
	" The original sequence files are no longer used by Bowtie once the\n"
	" index is built.  \n"
	"\n"
	" Use of Karkkainen's blockwise algorithm (see reference #4 below)\n"
	" allows bowtie-build to trade off between running time and memory\n"
	" usage. bowtie-build has three options governing how it makes this\n"
	" trade: -p/--packed, --bmax/--bmaxdivn, and --dcv.  By default, bowtie-\n"
	" build will automatically search for the settings that yield the best\n"
	" running time without exhausting memory.  This behavior can be disabled\n"
	" using the -a/--noauto option.\n"
	"\n"
	" The indexer provides options pertaining to the \"shape\" of the index,\n"
	" e.g. --offrate governs the fraction of Burrows-Wheeler rows that are\n"
	" \"marked\" (i.e., the \"density\" of the suffix-array sample; see\n"
	" reference #2).  All of these options are potentially profitable trade-\n"
	" offs depending on the application.  They have been set to defaults\n"
	" that are reasonable for most cases according to our experiments.  See\n"
	" \"High Performance Tips\" in the \"Using the 'bowtie' Aligner\" section\n"
	" for additional details.\n"
	"\n"
	" Because bowtie-build uses 32-bit pointers internally, it can handle up\n"
	" to a maximum of 2^32-1 (somewhat more than 4 billion) characters in an\n"
	" index.  If your reference exceeds 2^32-1 characters, bowtie-build will\n"
	" print an error message and abort.  To resolve this, divide your\n"
	" reference sequences into smaller batches and/or chunks and build a\n"
	" separate index for each.\n"
	" \n"
	" If your computer has more than 3-4 GB of memory and you would like to\n"
	" exploit that fact to make index building faster, you must use a 64-bit\n"
	" version of the bowtie-build binary.  The 32-bit version of the binary\n"
	" is restricted to using less than 4 GB of memory.  If a 64-bit pre-\n"
	" built binary does not yet exist for your platform on the sourceforge\n"
	" download site, you will need to build one from source.\n"
	"\n"
	" The Bowtie index is based on the FM Index of Ferragina and Manzini,\n"
	" which in turn is based on the Burrows-Wheeler transform.  The\n"
	" algorithm used to build the index is based on the blockwise algorithm\n"
	" of Karkkainen.  For more information on these techniques, see these\n"
	" references:\n"
	"\n"
	" 1. Burrows M, Wheeler DJ: A block sorting lossless data compression\n"
	"    algorithm. Digital Equipment Corporation, Palo Alto, CA 1994,\n"
	"    Technical Report 124.\n"
	" 2. Ferragina, P. and Manzini, G. 2000. Opportunistic data structures\n"
	"    with applications. In Proceedings of the 41st Annual Symposium on\n"
	"    Foundations of Computer Science (November 12 - 14, 2000). FOCS\n"
	" 3. Ferragina, P. and Manzini, G. 2001. An experimental study of an\n"
	"    opportunistic index. In Proceedings of the Twelfth Annual ACM-SIAM\n"
	"    Symposium on Discrete Algorithms (Washington, D.C., United States,\n"
	"    January 07 - 09, 2001). 269-278.\n"
	" 4. Karkkainen, J. 2007. Fast BWT in small space by blockwise suffix\n"
	"    sorting. Theor. Comput. Sci. 387, 3 (Nov. 2007), 249-257\n"
	"\n"
	"  Command Line\n"
	"  ------------\n"
	"\n"
	" Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>\n"
	"\n"
	"    <reference_in>          A comma-separated list of FASTA files\n"
	"                            containing the reference sequences to be\n"
	"                            aligned to, or, if -c is specified, the\n"
	"                            sequences themselves. E.g., this might be\n"
	"                            \"chr1.fa,chr2.fa,chrX.fa,chrY.fa\", or, if\n"
	"                            -c is specified, this might be\n"
	"                            \"GGTCATCCT,ACGGGTCGT,CCGTTCTATGCGGCTTA\".\n"
	"\n"
	"    <ebwt_outfile_base>     The basename of the index files to write.\n"
	"                            By default, bowtie-build writes files named\n"
	"                            NAME.1.ebwt, NAME.2.ebwt, NAME.3.ebwt,\n"
	"                            NAME.4.ebwt, NAME.rev.1.ebwt, and\n"
	"                            NAME.rev.2.ebwt, where NAME is the\n"
	"                            basename.\n"
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
	"    -a/--noauto             Disable the default behavior whereby\n"
	"                            bowtie-build automatically selects values\n"
	"                            for --bmax/--dcv/--packed parameters\n"
	"                            according to the memory available.  User\n"
	"                            may specify values for those parameters.\n"
	"                            If memory is exhausted during indexing, an\n"
	"                            error message will be printed; it is up to\n"
	"                            the user to try new parameters.\n"
	"\n"
	"    -p/--packed             Use a packed (2-bits-per-nucleotide)\n"
	"                            representation for DNA strings.  This saves\n"
	"                            memory but makes indexing 2-3 times slower.\n"
	"                            Default: off.  This is configured\n"
	"                            automatically by default; use -a/--noauto\n"
	"                            to configure manually.\n"
	"\n"
	"    --bmax <int>            The maximum number of suffixes allowed in a\n"
	"                            block.  Allowing more suffixes per block\n"
	"                            makes indexing faster, but increases memory\n"
	"                            overhead.  Overrides any previous\n"
	"                            specification of --bmax, --bmaxmultsqrt or\n"
	"                            --bmaxdivn.  Default: --bmaxdivn 4.  This\n"
	"                            is configured automatically by default; use\n"
	"                            -a/--noauto to configure manually.\n"
	"\n"
	"    --bmaxdivn <int>        The maximum number of suffixes allowed in a\n"
	"                            block, expressed as a fraction of the\n"
	"                            length of the reference.  Overrides any\n"
	"                            previous specification of --bmax,\n"
	"                            --bmaxmultsqrt or --bmaxdivn. Default:\n"
	"                            --bmaxdivn 4.  This is configured\n"
	"                            automatically by default; use -a/--noauto\n"
	"                            to configure manually.\n"
	"\n"
	"    --dcv <int>             Use <int> as the period for the difference-\n"
	"                            cover sample.  A larger period yields less\n"
	"                            memory overhead, but may make suffix\n"
	"                            sorting slower, especially if repeats are\n"
	"                            present.  Must be a power of 2 no greater\n"
	"                            than 4096.  Default: 1024.  This is\n"
	"                            configured automatically by default; use\n"
	"                            -a/--noauto to configure manually.\n"
	"\n"
	"    --nodc                  Disable use of the difference-cover sample.\n"
	"                            Suffix sorting becomes quadratic-time in\n"
	"                            the worst case (where the worst case is an\n"
	"                            extremely repetitive reference).  Default:\n"
	"                            off.\n"
	"\n"
	"    -r/--noref              Do not build the NAME.3.ebwt and\n"
	"                            NAME.4.ebwt portions of the index, which\n"
	"                            contain a bitpacked version of the\n"
	"                            reference sequences and are (currently)\n"
	"                            only used for paired-end alignment.\n"
	"\n"
	"    -3/--justref            Build *only* the NAME.3.ebwt and\n"
	"                            NAME.4.ebwt portions of the index, which\n"
	"                            contain a bitpacked version of the\n"
	"                            reference sequences and are (currently)\n"
	"                            only used for paired-end alignment.\n"
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
	"                            'bowtie' will not find alignments that\n"
	"                            overlap them.\n"
	"\n"
	"    --big --little          Endianness to use when serializing integers\n"
	"                            to the index file.  Default: little-endian\n"
	"                            (recommended for Intel- and AMD-based\n"
	"                            architectures).\n"
	"\n"
	"    --seed <int>            Use <int> as the seed for pseudo-random\n"
	"                            number generator.\n"
	"\n"
	"    --cutoff <int>          Index only the first <int> bases of the\n"
	"                            reference sequences (cumulative across\n"
	"                            sequences) and ignore the rest.\n"
	"\n"
	"    --oldpmap               bowtie-build switched schemes for mapping\n"
	"                            \"joined\" reference locations to original\n"
	"                            reference locations in version 0.9.8.  The\n"
	"                            new scheme has the advantage that it does\n"
	"                            not use padding.  This option activates the\n"
	"                            old padding-based scheme used in versions\n"
	"                            prior to 0.9.8.  Versions of 'bowtie' prior\n"
	"                            to 0.9.8 can query only indexes that use\n"
	"                            the old scheme.  Version of 'bowtie'\n"
	"                            starting with 0.9.8 can query indexes using\n"
	"                            either scheme.  This option will be\n"
	"                            deprecated in version 1.0. \n"
	"\n"
	"    -q/--quiet              bowtie-build is verbose by default.  With\n"
	"                            this option ebwt-build will print only\n"
	"                            error messages.\n"
	"\n"
	"    -h/--help               Print detailed description of tool and its\n"
	"                            options (from MANUAL).\n"
	"\n"
	"    --version               Print version information and quit.\n"
	"\n"
	;
}

static const char *short_options = "qraph?nscfl:i:o:t:h:3";

static struct option long_options[] = {
	{"quiet",        no_argument,       0,            'q'},
	{"sanity",       no_argument,       0,            's'},
	{"packed",       no_argument,       0,            'p'},
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
	{"noauto",       no_argument,       0,            'a'},
	{"noblocks",     required_argument, 0,            'n'},
	{"linerate",     required_argument, 0,            'l'},
	{"linesperside", required_argument, 0,            'i'},
	{"offrate",      required_argument, 0,            'o'},
	{"isarate",      required_argument, 0,            ARG_ISARATE},
	{"ftabchars",    required_argument, 0,            't'},
	{"help",         no_argument,       0,            'h'},
	{"cutoff",       required_argument, 0,            ARG_CUTOFF},
	{"ntoa",         no_argument,       0,            ARG_NTOA},
	{"oldpmap",      no_argument,       0,            ARG_PMAP},
	{"justref",      no_argument,       0,            '3'},
	{"noref",        no_argument,       0,            'r'},
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
	   		case 'p': packed = true; break;
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
	   		case '3':
	   			justRef = true;
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
	   		case 'a': autoMem = false; break;
	   		case 'q': verbose = false; break;
	   		case 's': sanityCheck = true; break;
	   		case 'r': writeRef = false; break;

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
static void driver(const string& infile,
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
			if (f == NULL)
			{
				cerr << "Error: could not open "<< infiles[i] << endl;
				exit(1);
			}
			FileBuf *fb = new FileBuf(f);
			assert(fb != NULL);
			assert(!fb->eof());
			assert(fb->get() == '>');
			ASSERT_ONLY(fb->reset());
			assert(!fb->eof());
			is.push_back(fb);
		}
	}
	// Vector for the ordered list of "records" comprising the input
	// sequences.  A record represents a stretch of unambiguous
	// characters in one of the input sequences.
	vector<RefRecord> szs;
	uint32_t sztot;
	int32_t chunkRate = 0;
	{
		if(verbose) cout << "Reading reference sizes" << endl;
		Timer _t(cout, "  Time reading reference sizes: ", verbose);
		if(!reverse && (writeRef || justRef)) {
			string file3 = outfile + ".3.ebwt";
			string file4 = outfile + ".4.ebwt";
			BitpairOutFileBuf bpout(file4.c_str());
			// Read in the sizes of all the unambiguous stretches of the
			// genome into a vector of RefRecords
			sztot = fastaRefReadSizes(is, szs, refparams, &bpout);
			// Write the records back out to a '.3.ebwt' file.
			ofstream fout3(file3.c_str(), ios::binary);
			if(!fout3.good()) {
				cerr << "Could not open index file for writing: \"" << file3 << "\"" << endl
					 << "Please make sure the directory exists and that permissions allow writing by" << endl
					 << "Bowtie." << endl;
				exit(1);
			}
			bool be = currentlyBigEndian();
			writeU32(fout3, 1, be); // endianness sentinel
			writeU32(fout3, szs.size(), be); // write # records
			// Write the records themselves
			for(size_t i = 0; i < szs.size(); i++) {
				szs[i].write(fout3, be);
			}
			bpout.close();
			fout3.close();
#ifndef NDEBUG
			if(sanityCheck) {
				BitPairReference bpr(outfile, true, &infiles, NULL, format == CMDLINE);
			}
#endif
		} else {
			// Read in the sizes of all the unambiguous stretches of the
			// genome into a vector of RefRecords
			sztot = fastaRefReadSizes(is, szs, refparams);
		}
	}
	if(justRef) return;
	// If we have to insert padding between reference sequences...
	if(!useBsearch) {
		// Then calculate the space-optimal size for the padding
		chunkRate = EbwtParams::calcBestChunkRate(szs, offRate, lineRate, linesPerSide);
		if(verbose) cout << "  Choose best chunkRate: " << chunkRate << endl;
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
	                !reverse,     // fw
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
		     << "  Strings: " << (packed? "packed" : "unpacked") << endl
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
		if(!packed) {
			try {
				driver<String<Dna, Alloc<> > >(infile, infiles, outfile);
			} catch(bad_alloc& e) {
				if(autoMem) {
					cerr << "Switching to a packed string representation." << endl;
					packed = true;
				} else {
					throw e;
				}
			}
		}
		if(packed) {
			driver<String<Dna, Packed<Alloc<> > > >(infile, infiles, outfile);
		}
	}
	cutoff = origCutoff; // reset cutoff for backward Ebwt
	if(doubleEbwt) {
		srand(seed);
		Timer timer(cout, "Total time for backward call to driver() for mirror index: ", verbose);
		if(!packed) {
			try {
				driver<String<Dna, Alloc<> > >(infile, infiles, outfile + ".rev", true);
			} catch(bad_alloc& e) {
				if(autoMem) {
					cerr << "Switching to a packed string representation." << endl;
					packed = true;
				} else {
					throw e;
				}
			}
		}
		if(packed) {
			driver<String<Dna, Packed<Alloc<> > > >(infile, infiles, outfile + ".rev", true);
		}
	}
	return 0;
}
