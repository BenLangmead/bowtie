#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <seqan/find.h>
#include <getopt.h>
#include <vector>
#include <pthread.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "endian_swap.h"
#include "ebwt.h"
#include "formats.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "hit.h"
#include "pat.h"
#include "bitset.h"
#include "threading.h"
#include "range_cache.h"
#include "refmap.h"
#include "annot.h"
#include "aligner.h"
#include "aligner_0mm.h"
#include "aligner_1mm.h"
#include "aligner_23mm.h"
#include "aligner_seed_mm.h"
#include "aligner_metrics.h"
#include "sam.h"
#ifdef CHUD_PROFILING
#include <CHUD/CHUD.h>
#endif

using namespace std;
using namespace seqan;

static vector<string> mates1;  // mated reads (first mate)
static vector<string> mates2;  // mated reads (second mate)
static vector<string> mates12; // mated reads (1st/2nd interleaved in 1 file)
static string adjustedEbwtFileBase;
static bool verbose;      // be talkative
static bool startVerbose; // be talkative at startup
static bool quiet;        // print nothing but the alignments
static int sanityCheck;   // enable expensive sanity checks
static int format;        // default read format is FASTQ
static string origString; // reference text, or filename(s)
static int revcomp;       // search for reverse complements?
static int seed;          // srandom() seed
static int timing;        // whether to report basic timing data
static bool allHits;      // for multihits, report just one
static bool rangeMode;    // report BWT ranges instead of ref locs
static int showVersion;   // just print version and quit?
static int ipause;        // pause before maching?
static uint32_t qUpto;    // max # of queries to read
static int trim5;         // amount to trim from 5' end
static int trim3;         // amount to trim from 3' end
static int reportOpps;    // whether to report # of other mappings
static int offRate;       // keep default offRate
static int isaRate;       // keep default isaRate
static int mismatches;    // allow 0 mismatches by default
static char *patDumpfile; // filename to dump patterns to
static bool solexaQuals;  // quality strings are solexa quals, not phred, and subtract 64 (not 33)
static bool phred64Quals; // quality chars are phred, but must subtract 64 (not 33)
static bool integerQuals; // quality strings are space-separated strings of integers, not ASCII
static int maqLike;       // do maq-like searching
static int seedLen;       // seed length (changed in Maq 0.6.4 from 24)
static int seedMms;       // # mismatches allowed in seed (maq's -n)
static int qualThresh;    // max qual-weighted hamming dist (maq's -e)
static int maxBtsBetter;  // max # backtracks allowed in half-and-half mode
static int maxBts;        // max # backtracks allowed in half-and-half mode
static int nthreads;      // number of pthreads operating concurrently
static output_types outType;  // style of output
static bool randReadsNoSync;  // true -> generate reads from per-thread random source
static int numRandomReads;    // # random reads (see Random*PatternSource in pat.h)
static int lenRandomReads;    // len of random reads (see Random*PatternSource in pat.h)
static bool fullIndex;        // load halves one at a time and proceed in phases
static bool noRefNames;       // true -> print reference indexes; not names
static string dumpAlFaBase;   // basename of FASTA files to dump aligned reads to
static string dumpAlFqBase;   // basename of FASTQ files to dump aligned reads to
static string dumpAlBase;     // basename of same-format files to dump aligned reads to
static string dumpUnalFaBase; // basename of FASTA files to dump unaligned reads to
static string dumpUnalFqBase; // basename of FASTQ files to dump unaligned reads to
static string dumpUnalBase;   // basename of same-format files to dump unaligned reads to
static string dumpMaxFaBase;  // basename of FASTA files to dump reads with more than -m valid alignments to
static string dumpMaxFqBase;  // basename of FASTQ files to dump reads with more than -m valid alignments to
static string dumpMaxBase;    // basename of same-format files to dump reads with more than -m valid alignments to
static uint32_t khits;  // number of hits per read; >1 is much slower
static uint32_t mhits;  // don't report any hits if there are > mhits
static bool better;     // true -> guarantee alignments from best possible stratum
static bool oldBest;    // true -> guarantee alignments from best possible stratum (the old way)
static bool strata;     // true -> don't stop at stratum boundaries
static bool refOut;     // if true, alignments go to per-ref files
static int partitionSz; // output a partitioning key in first field
static bool noMaqRound; // true -> don't round quals to nearest 10 like maq
static bool forgiveInput; // let read input be a little wrong w/o complaining or dying
static bool useSpinlock;  // false -> don't use of spinlocks even if they're #defines
static bool fileParallel; // separate threads read separate input files in parallel
static bool useShmem;     // use shared memory to hold the index
static bool useMm;        // use memory-mapped files to hold the index
static bool mmSweep;      // sweep through memory-mapped files immediately after mapping
static bool stateful;     // use stateful aligners
static uint32_t prefetchWidth; // number of reads to process in parallel w/ --stateful
static uint32_t minInsert;     // minimum insert size (Maq = 0, SOAP = 400)
static uint32_t maxInsert;     // maximum insert size (Maq = 250, SOAP = 600)
static bool mate1fw;           // -1 mate aligns in fw orientation on fw strand
static bool mate2fw;           // -2 mate aligns in rc orientation on fw strand
static uint32_t mixedThresh;   // threshold for when to switch to paired-end mixed mode (see aligner.h)
static uint32_t mixedAttemptLim; // number of attempts to make in "mixed mode" before giving up on orientation
static bool dontReconcileMates;  // suppress pairwise all-versus-all way of resolving mates
static uint32_t cacheLimit;      // ranges w/ size > limit will be cached
static uint32_t cacheSize;       // # words per range cache
static int offBase;              // offsets are 0-based by default, but configurable
static bool tryHard;             // set very high maxBts, mixedAttemptLim
static uint32_t skipReads;       // # reads/read pairs to skip
static bool nofw; // don't align fw orientation of read
static bool norc; // don't align rc orientation of read
static bool strandFix;  // attempt to fix strand bias
static bool randomizeQuals; // randomize quality values
static bool stats; // print performance stats
static int chunkPoolMegabytes;    // max MB to dedicate to best-first search frames per thread
static int chunkSz;    // size of single chunk disbursed by ChunkPool
static bool chunkVerbose; // have chunk allocator output status messages?
static bool recal;
static int recalMaxCycle;
static int recalMaxQual;
static int recalQualShift;
static bool useV1;
static bool reportSe;
static const char * refMapFile;  // file containing a map from index coordinates to another coordinate system
static const char * annotMapFile;  // file containing a map from reference coordinates to annotations
static size_t fastaContLen;
static size_t fastaContFreq;
static bool hadoopOut; // print Hadoop status and summary messages
static bool fuzzy;
static bool fullRef;
static bool samNoHead; // don't print any header lines in SAM output
static bool samNoSQ;   // don't print @SQ header lines

static void resetOptions() {
	mates1.clear();
	mates2.clear();
	mates12.clear();
	adjustedEbwtFileBase	= "";
	verbose					= 0;
	startVerbose			= 0;
	quiet					= false;
	sanityCheck				= 0;  // enable expensive sanity checks
	format					= FASTQ; // default read format is FASTQ
	origString				= ""; // reference text, or filename(s)
	revcomp					= 1; // search for reverse complements?
	seed					= 0; // srandom() seed
	timing					= 0; // whether to report basic timing data
	allHits					= false; // for multihits, report just one
	rangeMode				= false; // report BWT ranges instead of ref locs
	showVersion				= 0; // just print version and quit?
	ipause					= 0; // pause before maching?
	qUpto					= 0xffffffff; // max # of queries to read
	trim5					= 0; // amount to trim from 5' end
	trim3					= 0; // amount to trim from 3' end
	reportOpps				= 0; // whether to report # of other mappings
	offRate					= -1; // keep default offRate
	isaRate					= -1; // keep default isaRate
	mismatches				= 0; // allow 0 mismatches by default
	patDumpfile				= NULL; // filename to dump patterns to
	solexaQuals				= false; // quality strings are solexa quals, not phred, and subtract 64 (not 33)
	phred64Quals			= false; // quality chars are phred, but must subtract 64 (not 33)
	integerQuals			= false; // quality strings are space-separated strings of integers, not ASCII
	maqLike					= 1;   // do maq-like searching
	seedLen					= 28;  // seed length (changed in Maq 0.6.4 from 24)
	seedMms					= 2;   // # mismatches allowed in seed (maq's -n)
	qualThresh				= 70;  // max qual-weighted hamming dist (maq's -e)
	maxBtsBetter			= 125; // max # backtracks allowed in half-and-half mode
	maxBts					= 800; // max # backtracks allowed in half-and-half mode
	nthreads				= 1;     // number of pthreads operating concurrently
	outType					= OUTPUT_FULL;  // style of output
	randReadsNoSync			= false; // true -> generate reads from per-thread random source
	numRandomReads			= 50000000; // # random reads (see Random*PatternSource in pat.h)
	lenRandomReads			= 35;    // len of random reads (see Random*PatternSource in pat.h)
	fullIndex				= true;  // load halves one at a time and proceed in phases
	noRefNames				= false; // true -> print reference indexes; not names
	dumpAlFaBase			= "";    // basename of FASTA files to dump aligned reads to
	dumpAlFqBase			= "";    // basename of FASTQ files to dump aligned reads to
	dumpAlBase				= "";    // basename of same-format files to dump aligned reads to
	dumpUnalFaBase			= "";    // basename of FASTA files to dump unaligned reads to
	dumpUnalFqBase			= "";    // basename of FASTQ files to dump unaligned reads to
	dumpUnalBase			= "";    // basename of same-format files to dump unaligned reads to
	dumpMaxFaBase			= "";    // basename of FASTA files to dump reads with more than -m valid alignments to
	dumpMaxFqBase			= "";    // basename of FASTQ files to dump reads with more than -m valid alignments to
	dumpMaxBase				= "";    // basename of same-format files to dump reads with more than -m valid alignments to
	khits					= 1;     // number of hits per read; >1 is much slower
	mhits					= 0xffffffff; // don't report any hits if there are > mhits
	better					= false; // true -> guarantee alignments from best possible stratum
	oldBest					= false; // true -> guarantee alignments from best possible stratum (the old way)
	strata					= false; // true -> don't stop at stratum boundaries
	refOut					= false; // if true, alignments go to per-ref files
	partitionSz				= 0;     // output a partitioning key in first field
	noMaqRound				= false; // true -> don't round quals to nearest 10 like maq
	forgiveInput			= false; // let read input be a little wrong w/o complaining or dying
	useSpinlock				= true;  // false -> don't use of spinlocks even if they're #defines
	fileParallel			= false; // separate threads read separate input files in parallel
	useShmem				= false; // use shared memory to hold the index
	useMm					= false; // use memory-mapped files to hold the index
	mmSweep					= false; // sweep through memory-mapped files immediately after mapping
	stateful				= false; // use stateful aligners
	prefetchWidth			= 1;     // number of reads to process in parallel w/ --stateful
	minInsert				= 0;     // minimum insert size (Maq = 0, SOAP = 400)
	maxInsert				= 250;   // maximum insert size (Maq = 250, SOAP = 600)
	mate1fw					= true;  // -1 mate aligns in fw orientation on fw strand
	mate2fw					= false; // -2 mate aligns in rc orientation on fw strand
	mixedThresh				= 4;     // threshold for when to switch to paired-end mixed mode (see aligner.h)
	mixedAttemptLim			= 100;   // number of attempts to make in "mixed mode" before giving up on orientation
	dontReconcileMates		= true;  // suppress pairwise all-versus-all way of resolving mates
	cacheLimit				= 5;     // ranges w/ size > limit will be cached
	cacheSize				= 0;     // # words per range cache
	offBase					= 0;     // offsets are 0-based by default, but configurable
	tryHard					= false; // set very high maxBts, mixedAttemptLim
	skipReads				= 0;     // # reads/read pairs to skip
	nofw					= false; // don't align fw orientation of read
	norc					= false; // don't align rc orientation of read
	strandFix				= true;  // attempt to fix strand bias
	randomizeQuals			= false; // randomize quality values
	stats					= false; // print performance stats
	chunkPoolMegabytes		= 32;    // max MB to dedicate to best-first search frames per thread
	chunkSz					= 16;    // size of single chunk disbursed by ChunkPool
	chunkVerbose			= false; // have chunk allocator output status messages?
	recal					= false;
	recalMaxCycle			= 64;
	recalMaxQual			= 40;
	recalQualShift			= 2;
	useV1					= true;
	reportSe				= false;
	refMapFile				= NULL;  // file containing a map from index coordinates to another coordinate system
	annotMapFile			= NULL;  // file containing a map from reference coordinates to annotations
	fastaContLen			= 0;
	fastaContFreq			= 0;
	hadoopOut				= false; // print Hadoop status and summary messages
	fuzzy					= false; // reads will have alternate basecalls w/ qualities
	fullRef					= false; // print entire reference name instead of just up to 1st space
	samNoHead				= false; // don't print any header lines in SAM output
	samNoSQ					= false; // don't print @SQ header lines
}

// mating constraints

static const char *short_options = "fF:qbzhcu:rv:s:at3:5:o:e:n:l:w:p:k:m:1:2:I:X:x:B:yS";

enum {
	ARG_ORIG = 256,
	ARG_SEED,
	ARG_DUMP_PATS,
	ARG_RANGE,
	ARG_CONCISE,
	ARG_SOLEXA_QUALS,
	ARG_MAXBTS,
	ARG_VERBOSE,
	ARG_STARTVERBOSE,
	ARG_QUIET,
	ARG_RANDOM_READS,
	ARG_RANDOM_READS_NOSYNC,
	ARG_NOOUT,
	ARG_FAST,
	ARG_AL,
	ARG_UN,
	ARG_MAXDUMP,
	ARG_REFIDX,
	ARG_DUMP_NOHIT,
	ARG_DUMP_HHHIT,
	ARG_SANITY,
	ARG_OLDBEST,
	ARG_BETTER,
	ARG_BEST,
	ARG_REFOUT,
	ARG_ISARATE,
	ARG_SEED_EXTEND,
	ARG_PARTITION,
	ARG_integerQuals,
	ARG_FORGIVE_INPUT,
	ARG_NOMAQROUND,
	ARG_USE_SPINLOCK,
	ARG_FILEPAR,
	ARG_SHMEM,
	ARG_MM,
	ARG_MMSWEEP,
	ARG_STATEFUL,
	ARG_PREFETCH_WIDTH,
	ARG_FF,
	ARG_FR,
	ARG_RF,
	ARG_MIXED_ATTEMPTS,
	ARG_NO_RECONCILE,
	ARG_CACHE_LIM,
	ARG_CACHE_SZ,
	ARG_NO_FW,
	ARG_NO_RC,
	ARG_SKIP,
	ARG_STRAND_FIX,
	ARG_RANDOMIZE_QUALS,
	ARG_STATS,
	ARG_ONETWO,
	ARG_PHRED64,
	ARG_PHRED33,
	ARG_CHUNKMBS,
	ARG_CHUNKSZ,
	ARG_CHUNKVERBOSE,
	ARG_RECAL,
	ARG_STRATA,
	ARG_PEV2,
	ARG_CHAINOUT,
	ARG_CHAININ,
	ARG_REFMAP,
	ARG_ANNOTMAP,
	ARG_REPORTSE,
	ARG_HADOOPOUT,
	ARG_FUZZY,
	ARG_FULLREF,
	ARG_USAGE,
	ARG_SAM_NOHEAD,
	ARG_SAM_NOSQ
};

static struct option long_options[] = {
	{(char*)"verbose",      no_argument,       0,            ARG_VERBOSE},
	{(char*)"startverbose", no_argument,       0,            ARG_STARTVERBOSE},
	{(char*)"quiet",        no_argument,       0,            ARG_QUIET},
	{(char*)"sanity",       no_argument,       0,            ARG_SANITY},
	{(char*)"exact",        no_argument,       0,            '0'},
	{(char*)"1mm",          no_argument,       0,            '1'},
	{(char*)"2mm",          no_argument,       0,            '2'},
	{(char*)"pause",        no_argument,       &ipause,      1},
	{(char*)"orig",         required_argument, 0,            ARG_ORIG},
	{(char*)"all",          no_argument,       0,            'a'},
	{(char*)"concise",      no_argument,       0,            ARG_CONCISE},
	{(char*)"binout",       no_argument,       0,            'b'},
	{(char*)"noout",        no_argument,       0,            ARG_NOOUT},
	{(char*)"solexa-quals", no_argument,       0,            ARG_SOLEXA_QUALS},
	{(char*)"integer-quals",no_argument,       0,            ARG_integerQuals},
	{(char*)"time",         no_argument,       0,            't'},
	{(char*)"trim3",        required_argument, 0,            '3'},
	{(char*)"trim5",        required_argument, 0,            '5'},
	{(char*)"seed",         required_argument, 0,            ARG_SEED},
	{(char*)"qupto",        required_argument, 0,            'u'},
	{(char*)"al",           required_argument, 0,            ARG_AL},
	{(char*)"un",           required_argument, 0,            ARG_UN},
	{(char*)"max",          required_argument, 0,            ARG_MAXDUMP},
	{(char*)"offrate",      required_argument, 0,            'o'},
	{(char*)"isarate",      required_argument, 0,            ARG_ISARATE},
	{(char*)"reportopps",   no_argument,       &reportOpps,  1},
	{(char*)"version",      no_argument,       &showVersion, 1},
	{(char*)"dumppats",     required_argument, 0,            ARG_DUMP_PATS},
	{(char*)"revcomp",      no_argument,       0,            'r'},
	{(char*)"maqerr",       required_argument, 0,            'e'},
	{(char*)"seedlen",      required_argument, 0,            'l'},
	{(char*)"seedmms",      required_argument, 0,            'n'},
	{(char*)"filepar",      no_argument,       0,            ARG_FILEPAR},
	{(char*)"help",         no_argument,       0,            'h'},
	{(char*)"threads",      required_argument, 0,            'p'},
	{(char*)"khits",        required_argument, 0,            'k'},
	{(char*)"mhits",        required_argument, 0,            'm'},
	{(char*)"minins",       required_argument, 0,            'I'},
	{(char*)"maxins",       required_argument, 0,            'X'},
	{(char*)"best",         no_argument,       0,            ARG_BEST},
	{(char*)"better",       no_argument,       0,            ARG_BETTER},
	{(char*)"oldbest",      no_argument,       0,            ARG_OLDBEST},
	{(char*)"strata",       no_argument,       0,            ARG_STRATA},
	{(char*)"nomaqround",   no_argument,       0,            ARG_NOMAQROUND},
	{(char*)"refidx",       no_argument,       0,            ARG_REFIDX},
	{(char*)"range",        no_argument,       0,            ARG_RANGE},
	{(char*)"maxbts",       required_argument, 0,            ARG_MAXBTS},
	{(char*)"randread",     no_argument,       0,            ARG_RANDOM_READS},
	{(char*)"randreadnosync", no_argument,     0,            ARG_RANDOM_READS_NOSYNC},
	{(char*)"phased",       no_argument,       0,            'z'},
	{(char*)"dumpnohit",    no_argument,       0,            ARG_DUMP_NOHIT},
	{(char*)"dumphhhit",    no_argument,       0,            ARG_DUMP_HHHIT},
	{(char*)"refout",       no_argument,       0,            ARG_REFOUT},
	{(char*)"seedextend",   no_argument,       0,            ARG_SEED_EXTEND},
	{(char*)"partition",    required_argument, 0,            ARG_PARTITION},
	{(char*)"forgive",      no_argument,       0,            ARG_FORGIVE_INPUT},
	{(char*)"nospin",       no_argument,       0,            ARG_USE_SPINLOCK},
	{(char*)"stateful",     no_argument,       0,            ARG_STATEFUL},
	{(char*)"prewidth",     required_argument, 0,            ARG_PREFETCH_WIDTH},
	{(char*)"ff",           no_argument,       0,            ARG_FF},
	{(char*)"fr",           no_argument,       0,            ARG_FR},
	{(char*)"rf",           no_argument,       0,            ARG_RF},
	{(char*)"mixthresh",    required_argument, 0,            'x'},
	{(char*)"pairtries",    required_argument, 0,            ARG_MIXED_ATTEMPTS},
	{(char*)"noreconcile",  no_argument,       0,            ARG_NO_RECONCILE},
	{(char*)"cachelim",     required_argument, 0,            ARG_CACHE_LIM},
	{(char*)"cachesz",      required_argument, 0,            ARG_CACHE_SZ},
	{(char*)"nofw",         no_argument,       0,            ARG_NO_FW},
	{(char*)"norc",         no_argument,       0,            ARG_NO_RC},
	{(char*)"offbase",      required_argument, 0,            'B'},
	{(char*)"tryhard",      no_argument,       0,            'y'},
	{(char*)"skip",         required_argument, 0,            's'},
	{(char*)"strandfix",    no_argument,       0,            ARG_STRAND_FIX},
	{(char*)"randquals",    no_argument,       0,            ARG_RANDOMIZE_QUALS},
	{(char*)"stats",        no_argument,       0,            ARG_STATS},
	{(char*)"12",           required_argument, 0,            ARG_ONETWO},
	{(char*)"phred33-quals", no_argument,      0,            ARG_PHRED33},
	{(char*)"phred64-quals", no_argument,      0,            ARG_PHRED64},
	{(char*)"solexa1.3-quals", no_argument,    0,            ARG_PHRED64},
	{(char*)"chunkmbs",     required_argument, 0,            ARG_CHUNKMBS},
	{(char*)"chunksz",      required_argument, 0,            ARG_CHUNKSZ},
	{(char*)"chunkverbose", no_argument,       0,            ARG_CHUNKVERBOSE},
	{(char*)"mm",           no_argument,       0,            ARG_MM},
	{(char*)"shmem",        no_argument,       0,            ARG_SHMEM},
	{(char*)"mmsweep",      no_argument,       0,            ARG_MMSWEEP},
	{(char*)"recal",        no_argument,       0,            ARG_RECAL},
	{(char*)"pev2",         no_argument,       0,            ARG_PEV2},
	{(char*)"chainout",     no_argument,       0,            ARG_CHAINOUT},
	{(char*)"chainin",      no_argument,       0,            ARG_CHAININ},
	{(char*)"refmap",       required_argument, 0,            ARG_REFMAP},
	{(char*)"annotmap",     required_argument, 0,            ARG_ANNOTMAP},
	{(char*)"reportse",     no_argument,       0,            ARG_REPORTSE},
	{(char*)"hadoopout",    no_argument,       0,            ARG_HADOOPOUT},
	{(char*)"fuzzy",        no_argument,       0,            ARG_FUZZY},
	{(char*)"fullref",      no_argument,       0,            ARG_FULLREF},
	{(char*)"usage",        no_argument,       0,            ARG_USAGE},
	{(char*)"sam",          no_argument,       0,            'S'},
	{(char*)"sam-nohead",   no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"sam-nosq",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"sam-noSQ",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: " << endl
        << "  bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]" << endl
        << endl
	    << "  <m1>    Comma-separated list of files containing upstream mates (or the" << endl
	    << "          sequences themselves, if -c is set) paired with mates in <m2>" << endl
	    << "  <m2>    Comma-separated list of files containing downstream mates (or the" << endl
	    << "          sequences themselves if -c is set) paired with mates in <m1>" << endl
	    << "  <r>     Comma-separated list of files containing Crossbow-style reads.  Can be" << endl
	    << "          a mixture of paired and unpaired.  Specify \"-\" for stdin." << endl
	    << "  <s>     Comma-separated list of files containing unpaired reads, or the" << endl
	    << "          sequences themselves, if -c is set.  Specify \"-\" for stdin." << endl
	    << "  <hit>   File to write hits to (default: stdout)" << endl
	    << "Input:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    << "  -r                 query input files are raw one-sequence-per-line" << endl
	    << "  -c                 query sequences given on cmd line (as <mates>, <singles>)" << endl
	    << "  -s/--skip <int>    skip the first <int> reads/pairs in the input" << endl
	    << "  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads" << endl
		<< "  --phred33-quals    input quals are Phred+33 (default)" << endl
		<< "  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)" << endl
		<< "  --solexa-quals     input quals are from GA Pipeline ver. < 1.3" << endl
		<< "  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3" << endl
		<< "  --integer-quals    qualities are given as space-separated integers (not ASCII)" << endl
	    << "Alignment:" << endl
	    << "  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)" << endl
	    << "  -e/--maqerr <int>  max sum of mismatch quals (rounds like maq; default: 70)" << endl
	    << "  -l/--seedlen <int> seed length (default: 28)" << endl
		<< "  --nomaqround       disable Maq-like quality rounding (to nearest 10 <= 30)" << endl
	    << "  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities" << endl
	    << "  -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)" << endl
	    << "  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)" << endl
	    << "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)" << endl
	    << "  --nofw/--norc      do not align to forward/reverse-complement reference strand" << endl
	    << "  --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)" << endl
	    << "  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)" << endl
	    << "  -y/--tryhard       try hard to find valid alignments, at the expense of speed" << endl
	    << "  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 32)" << endl
	    << "Reporting:" << endl
	    << "  -k <int>           report up to <int> good alignments per read (default: 1)" << endl
	    << "  -a/--all           report all alignments per read (much slower than low -k)" << endl
	    << "  -m <int>           suppress all alignments if > <int> exist (def: no limit)" << endl
	    //<< "  --better           alignments guaranteed best possible stratum (old --best)" << endl
	    << "  --best             hits guaranteed best stratum; ties broken by quality" << endl
	    << "  --strata           hits in sub-optimal strata aren't reported (requires --best)" << endl
	    << "  --strandfix        attempt to fix strand biases" << endl
	    << "Output:" << endl
	    << "  -S/--sam           write hits in SAM format" << endl
	    << "  --concise          write hits in concise format" << endl
	    << "  -t/--time          print wall-clock time taken by search phases" << endl
	    << "  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)" << endl
	    << "  --quiet            print nothing but the alignments" << endl
	    << "  --refout           write alignments to files refXXXXX.map, 1 map per reference" << endl
	    << "  --refidx           refer to ref. seqs by 0-based index rather than name" << endl
	    << "  --al <fname>       write aligned reads/pairs to file(s) <fname>" << endl
	    << "  --un <fname>       write unaligned reads/pairs to file(s) <fname>" << endl
	    << "  --max <fname>      write reads/pairs over -m limit to file(s) <fname>" << endl
	    << "  --sam-nohead       supppress header lines (starting with @) for SAM output" << endl
	    << "  --sam-nosq         supppress @SQ header lines for SAM output" << endl
	    << "  --fullref          write entire ref name (default: only up to 1st space)" << endl
	    << "Performance:" << endl
#ifdef BOWTIE_PTHREADS
	    << "  -p/--threads <int> number of alignment threads to launch (default: 1)" << endl
#endif
	    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
#ifdef BOWTIE_MM
	    << "  --mm               use memory-mapped I/O for index; many 'bowtie's can share" << endl
#endif
#ifdef BOWTIE_SHARED_MEM
	    << "  --shmem            use shared mem for index; many 'bowtie's can share" << endl
#endif
	    << "Other:" << endl
	    << "  --seed <int>       seed for random number generator" << endl
	    << "  --verbose          verbose output (for debugging)" << endl
	    << "  --version          print version information and quit" << endl
	    << "  -h/--help          print detailed description of tool and its options" << endl
	    << "  --usage            print this usage message" << endl
	    ;
}

/**
 * Print a detailed usage message to the provided output stream.
 *
 * Manual text converted to C++ string with something like:
 * cat MANUAL  | head -1015 | tail -940 | sed -e 's/\"/\\\"/g' | \
 *   sed -e 's/^/"/' | sed -e 's/$/\\n"/'
 */
static void printLongUsage(ostream& out) {
	out << "\n"
	" Using the 'bowtie' Aligner\n"
	" --------------------------\n"
	"\n"
	" The 'bowtie' aligner takes an index and a set of reads as input and\n"
	" outputs a list of alignments.  Alignments are selected according to a\n"
	" combination of the -v/-n/-e/-l options (plus the -I/-X/--fr/--rf/--ff\n"
	" options for paired-end alignment), which define which alignments are\n"
	" legal, and the -k/-a/-m/--best/--strata options which define which and\n"
	" how many legal alignments should be reported.\n"
	"\n"
	" By default, Bowtie enforces an alignment policy equivalent to Maq's\n"
	" quality-aware policy (http://maq.sf.net) (-n 2 -l 28 -e 70), but it\n"
	" can also be made to enforce an end-to-end k-difference policy\n"
	" equivalent to SOAP's (http://soap.genomics.org.cn/) (-v 2).\n"
	"\n"
	" Bowtie is designed to be very fast for read sets where a) many of the\n"
	" reads have at least one good, valid alignment, b) many of the reads\n"
	" are relatively high-quality, c) the number of alignments reported per\n"
	" read is small (close to 1).  These criteria are generally satisfied in\n"
	" the context of modern short-read analyses such as RNA-seq, ChIP-seq,\n"
	" other types of -seq, and especially mammalian genotyping (e.g. the\n"
	" 1000 Genomes Project).  You may observe longer running times in other\n"
	" research contexts.  If you find Bowtie's performance to be\n"
	" disappointingly slow, please try the hints described in the \"High\n"
	" Performance Tips\" section below.  If Bowtie continues to be too slow,\n"
	" please contact us and tell us the nature of your research application\n"
	" and the parameters you are using to run Bowtie.  We are eager to hear\n"
	" your feedback.\n"
	" \n"
	" A result of Bowtie's indexing strategy is that alignments involving\n"
	" one or more ambiguous reference characters ('N', '-', 'R', 'Y', etc.)\n"
	" are considered invalid by Bowtie, regardless of the alignment policy.\n"
	" This is true only for ambiguous characters in the reference;\n"
	" alignments involving ambiguous characters in the read are legal,\n"
	" subject to the alignment policy.\n"
	" \n"
	" Also, alignments that \"fall off\" the reference sequence are not\n"
	" considered legal by Bowtie, though some such alignments will become\n"
	" legal once gapped alignment is implemented.\n"
	"\n"
	" The process by which bowtie chooses an alignment to report is\n"
	" randomized in order to avoid \"mapping bias\" - the phenomenon whereby\n"
	" an aligner systematically fails to report a particular class of good\n"
	" alignments, causing spurious \"holes\" in the comparative assembly.\n"
	" Whenever bowtie reports a subset of the valid alignments that exist,\n"
	" it makes an effort to sample them randomly.  This randomness flows\n"
	" from a simple seeded pseudo-random number generator and is\n"
	" \"deterministic\" in the sense that Bowtie will always produce the same\n"
	" results for the same read when run with the same initial \"seed\" value\n"
	" (see documentation for --seed option).\n"
	" \n"
	" In the default mode, bowtie can exhibit strand bias.  Strand bias\n"
	" occurs when input reference and reads are such that (a) some reads\n"
	" align equally well to sites on the forward and reverse strands of the\n"
	" reference, and (b) the number of such sites on one strand is different\n"
	" from the number on the other strand.  When this happens for a given\n"
	" read, bowtie effectively chooses one strand or the other with 50%\n"
	" probability, then reports a randomly-selected alignment for that read\n"
	" from among the sites on the selected strand.  This tends to overassign\n"
	" alignments to the sites on the strand with fewer sites and underassign\n"
	" to sites on the strand with more sites.  The effect is mitigated,\n"
	" though it may not be eliminated, when reads are longer or when paired-\n"
	" end reads are used.  Running Bowtie in --best mode eliminates strand\n"
	" bias by forcing Bowtie to select one strand or the other with a\n"
	" probability that is proportional to the number of best sites on the\n"
	" strand. \n"
	"\n"
	" Gapped alignments are not currently supported, but we do plan to\n"
	" implement this in the future.  Alignment in ABI \"color space\" is also\n"
	" not currently supported.\n"
	"\n"
	"  Maq-like Policy\n"
	"  ---------------\n"
	"\n"
	"  When the -n option is specified (and it is by default), Bowtie\n"
	"  determines which alignments are valid according to the following\n"
	"  policy, which is equivalent to Maq's default policy:\n"
	"\n"
	"  1. Alignments may have no more than N mismatches in the first L\n"
	"     bases on the high-quality end of the read.\n"
	"\n"
	"  2. The sum of the quality values at all mismatched positions may not\n"
	"     exceed E (where each position has a quality value on a phred-like\n"
	"     scale of 0 up to about 40).\n"
	"\n"
	"  The N, L and E parameters are configured using Bowtie's -n, -l and\n"
	"  -e options.  The -n (Maq-like) option is mutually exclusive with the\n"
	"  -v (end-to-end k-difference) option.\n"
	" \n"
	"  If there are many possible alignments satisfying these criteria,\n"
	"  Bowtie will prefer to report alignments with fewer mismatches and\n"
	"  where the sum from criterion 2 is smaller.  However, Bowtie does not\n"
	"  guarantee that the reported alignment(s) are \"best\" in terms of the\n"
	"  number of mismatches (i.e. the alignment \"stratum\") or in terms of\n"
	"  the quality values at the mismatched positions unless the --best\n"
	"  option is specified.  Bowtie is about 1 to 2.5 times slower when\n"
	"  --best is specified.\n"
	"\n"
	"  Note that Maq internally rounds base qualities to the nearest 10 and\n"
	"  rounds qualities greater than 30 to 30.  To maintain compatibility\n"
	"  with Maq, Bowtie does the same.  Rounding can be suppressed with the\n"
	"  --nomaqround option.\n"
	" \n"
	"  Bowtie is not fully sensitive in -n 2 and -n 3 modes by default.  In\n"
	"  these modes Bowtie imposes a \"backtracking limit\" to limit effort\n"
	"  spent trying to find valid alignments for low-quality reads that are\n"
	"  unlikely to have any.  This may cause bowtie to miss some legal 2-\n"
	"  and 3-mismatch alignments.  The limit is set to a reasonable default\n"
	"  (125 without --best, 800 with --best), but the user may decrease or\n"
	"  increase the limit using the --maxbts and/or -y options.  -y mode is\n"
	"  slow but guarantees full sensitivity.\n"
	" \n"
	"  End-to-end k-difference Policy\n"
	"  ------------------------------\n"
	"  \n"
	"  The policy has one criterion: Alignments may have no more than V\n"
	"  mismatches.  Quality values are ignored.  The number of mismatches\n"
	"  permitted is configurable with the -v option.  The -v (end-to-end)\n"
	"  option is mutually exclusive with the -n (Maq-like) option.\n"
	"\n"
	"  If there are many possible alignments satisfying this criterion,\n"
	"  Bowtie will prefer to report alignments with fewer mismatches.\n"
	"  However, for reads where the \"best\" alignment has one or more\n"
	"  mismatches, Bowtie does not guarantee that the reported alignment(s)\n"
	"  will be best unless the --best option is specified.  Bowtie is\n"
	"  typically about 1 to 2.5 times slower when --best is specified.\n"
	"  \n"
	"  Reporting Modes\n"
	"  ---------------\n"
	"  \n"
	"  With the -k, -a, -m, --best and --strata options, Bowtie gives the\n"
	"  user a great deal of flexibility in selecting which alignments get\n"
	"  reported.  Here we give a few examples that demonstrate a few ways\n"
	"  they can be combined to achieve a desired result.  All examples are\n"
	"  using the e_coli index that comes packaged with Bowtie.  Alignment\n"
	"  summary output is elided.\n"
	"  \n"
	"    Example 1: -a\n"
	"    -------------\n"
	"    \n"
	"    $ ./bowtie -a -v 2 e_coli --concise -c ATGCATCATGCGCCAT\n"
	"    1-:<0,148810,2>\n"
	"    1-:<0,2852852,1>\n"
	"    1-:<0,4930433,2>\n"
	"    1-:<0,905664,2>\n"
	"    1+:<0,1093035,2>\n"
	"  \n"
	"    Specifying -a instructs bowtie to report all valid alignments,\n"
	"    subject to the alignment policy: -v 2.  In this case, bowtie finds\n"
	"    5 inexact hits in the E. coli genome; 1 hit (the 2nd one listed)\n"
	"    has 1 mismatch and 4 hits have 2 mismatches.  Note that they are\n"
	"    not necessarily listed in best-to-worst order.\n"
	"    \n"
	"    Example 2: -k 3\n"
	"    ---------------\n"
	"    \n"
	"    $ ./bowtie -k 3 -v 2 e_coli --concise -c ATGCATCATGCGCCAT\n"
	"    1-:<0,148810,2>\n"
	"    1-:<0,2852852,1>\n"
	"    1-:<0,4930433,2>\n"
	"  \n"
	"    Specifying -k 3 instructs bowtie to report up to 3 valid\n"
	"    alignments.  In this case, a total of 5 valid alignments exist (see\n"
	"    Example 1); bowtie reports 3 out of those 5.  -k can be set to any\n"
	"    integer greater than 0.\n"
	"\n"
	"    Example 3: -k 6\n"
	"    ---------------\n"
	"    \n"
	"    $ ./bowtie -k 6 -v 2 e_coli --concise -c ATGCATCATGCGCCAT\n"
	"    1-:<0,148810,2>\n"
	"    1-:<0,2852852,1>\n"
	"    1-:<0,4930433,2>\n"
	"    1-:<0,905664,2>\n"
	"    1+:<0,1093035,2>\n"
	"  \n"
	"    Specifying -k 6 instructs bowtie to report up to 6 valid\n"
	"    alignments.  In this case, a total of 5 valid alignments exist, so\n"
	"    bowtie reports all 5.\n"
	"    \n"
	"    Example 4: default (-k 1)\n"
	"    -------------------------\n"
	"    \n"
	"    $ ./bowtie -v 2 e_coli --concise -c ATGCATCATGCGCCAT\n"
	"    1-:<0,148810,2>\n"
	"    \n"
	"    Leaving the reporting options at their defaults causes Bowtie to\n"
	"    report the first valid alignment it encounters.  Because --best was\n"
	"    not specified, we are not guaranteed that bowtie will report the\n"
	"    best alignment, and in this case it does not (the 1-mismatch\n"
	"    alignment from the previous example would have been better).  The\n"
	"    default reporting mode is equivalent to -k 1.\n"
	"\n"
	"    Example 5: -a --best\n"
	"    --------------------\n"
	"    \n"
	"    $ ./bowtie -a --best -v 2 e_coli --concise -c ATGCATCATGCGCCAT\n"
	"    1-:<0,2852852,1>\n"
	"    1+:<0,1093035,2>\n"
	"    1-:<0,905664,2>\n"
	"    1-:<0,148810,2>\n"
	"    1-:<0,4930433,2>\n"
	"  \n"
	"    Specifying -a --best results in the same alignments being printed\n"
	"    as if just -a had been specified, but they are guaranteed to be\n"
	"    reported in best-to-worst order.\n"
	"\n"
	"    Example 6: -a --best --strata\n"
	"    -----------------------------\n"
	"    \n"
	"    $ ./bowtie -a --best --strata -v 2 e_coli --concise \\\n"
	"               -c ATGCATCATGCGCCAT\n"
	"    1-:<0,2852852,1>\n"
	"\n"
	"    Specifying --strata in addition to -a and --best causes Bowtie to\n"
	"    report only those alignments in the best alignment \"stratum\".  The\n"
	"    alignments in the best stratum are those having the least number of\n"
	"    mismatches (or mismatches just in the \"seed\" portion of the\n"
	"    alignment in the case of -n mode).  Note that if --strata is\n"
	"    specified, --best must also be specified.\n"
	"\n"
	"    Example 7: -a -m 3\n"
	"    ------------------\n"
	"    \n"
	"    $ ./bowtie -a -m 3 -v 2 e_coli --concise -c ATGCATCATGCGCCAT\n"
	"    No alignments\n"
	"  \n"
	"    Specifying -m 3 instructs bowtie to refrain from reporting any\n"
	"    alignments for reads having more than 3 reportable alignments.  The\n"
	"    -m option is useful when the user would like to guarantee that\n"
	"    reported alignments are \"unique\", for some definition of unique.\n"
	"    \n"
	"    Example 1 showed that the read has 5 reportable alignments when -a\n"
	"    and -v 2 are specified, so the -m 3 limit causes bowtie to output\n"
	"    no alignments.\n"
	"\n"
	"    Example 8: -a -m 5\n"
	"    ------------------\n"
	"    \n"
	"    $ ./bowtie -a -m 5 -v 2 e_coli --concise -c ATGCATCATGCGCCAT\n"
	"    1-:<0,148810,2>\n"
	"    1-:<0,2852852,1>\n"
	"    1-:<0,4930433,2>\n"
	"    1-:<0,905664,2>\n"
	"    1+:<0,1093035,2>\n"
	"  \n"
	"    Specifying -m 5 instructs bowtie to refrain from reporting any\n"
	"    alignments for reads having more than 5 reportable alignments.\n"
	"    Since the read has exactly 5 reportable alignments, the -m 5 limit\n"
	"    allows bowtie to print them as usual. \n"
	"\n"
	"    Example 9: -a -m 3 --best --strata\n"
	"    ----------------------------------\n"
	"    \n"
	"    $ ./bowtie -a -m 3 --best --strata -v 2 e_coli --concise \\n"
	"               -c ATGCATCATGCGCCAT\n"
	"    1-:<0,2852852,1>\n"
	"  \n"
	"    Specifying -m 3 instructs bowtie to refrain from reporting any\n"
	"    alignments for reads having more than 3 reportable alignments.\n"
	"    As we saw in Example 6, the read has only 1 reportable alignment\n"
	"    when -a, --best and --strata are specified, so the -m 3 limit\n"
	"    allows bowtie to print that alignment as usual.\n"
	"    \n"
	"    Intuitively, the -m option, when combined with the --best and\n"
	"    --strata options, guarntees a principled, though somewhat weaker\n"
	"    form of \"uniqueness.\"  A stronger form of uniqueness is enforced\n"
	"    when -m is specified but --best --strata are not.\n"
	"\n"
	"  Paired-end Alignment\n"
	"  --------------------\n"
	"  \n"
	"  Bowtie can align paired-end reads when paired read files are\n"
	"  specified using the -1 and -2 options (for pairs of raw, FASTA, or\n"
	"  FASTQ read files), or using the --12 option (for Tab-delimited read\n"
	"  files).  A valid paired-end alignment satisfies the following\n"
	"  criteria:\n"
	"  \n"
	"  1. Both mates have a valid alignment according to the alignment\n"
	"     policy specified by the -v/-n/-e/-l options.\n"
	"  2. The relative orientation and position of the mates satisfy the\n"
	"     constraints given by the -I/-X/--fr/--rf/--ff options. \n"
	"  \n"
	"  Policies governing which paired-end alignments are reported for a\n"
	"  given read are specified using the -k, -a and -m options as usual.\n"
	"  The --strata and --best options do not apply in paired-end mode.\n"
	"  \n"
	"  A paired-end alignment is reported as a pair of mate alignments, both\n"
	"  on a separate line, where the alignment for each mate is formatted\n"
	"  the same as an unpaired (singleton) alignment.  The alignment for the\n"
	"  mate that occurs closest to the beginning of the reference sequence\n"
	"  (the \"upstream\" mate) is always printed before the alignment for the\n"
	"  downstream mate.  Reads files containing paired-end reads will\n"
	"  sometimes name the reads according to whether they are the #1 or #2\n"
	"  mates by appending a \"/1\" or \"/2\" suffix to the read name.  If no\n"
	"  such suffix is present in Bowtie's input, the suffix will be added\n"
	"  when Bowtie prints read names in alignments.\n"
	"  \n"
	"  Finding a valid paired-end alignment where both mates align to\n"
	"  repetitive regions of the reference can be very time-consuming.  By\n"
	"  default, Bowtie avoids much of this cost by imposing a limit on the\n"
	"  number of \"tries\" it makes to match an alignment for one mate with a\n"
	"  nearby alignment for the other.  The default limit is 100.  This\n"
	"  causes Bowtie to miss some valid paired-end alignments where both\n"
	"  mates lie in repetitive regions, but the user may use the --pairtries\n"
	"  or -y options to increase Bowtie's sensitivity as desired.\n"
	"  \n"
	"  Paired-end alignments where one mate's alignment is entirely\n"
	"  contained within the other's are considered invalid.\n"
	" \n"
	"  Because Bowtie uses an in-memory representation of the original\n"
	"  reference string when finding paired-end alignments, its memory\n"
	"  footprint is larger when aligning paired-end reads.  For example, the\n"
	"  human index has a memory footprint of about 2.2 GB in single-end mode\n"
	"  and 2.9 GB in paired-end mode.\n"
	"  \n"
	"  High Performance Tips\n"
	"  ---------------------\n"
	"\n"
	"  Tip 1: Use 64-bit bowtie if possible\n"
	"\n"
	"  The 64-bit version of Bowtie is substantially faster (usually more\n"
	"  than 50% faster) than the 32-bit version, due to Bowtie's use of\n"
	"  64-bit arithmetic when searching both in the index and in the\n"
	"  reference.  If possible, download the 64-bit binaries for Bowtie and\n"
	"  run them on a 64-bit machine.  If you are building Bowtie from\n"
	"  sources, you may need to pass the -m64 option to g++ to compile the\n"
	"  64-bit version; you can do this by supplying argument BITS=64 to the\n"
	"  'make' command; e.g.: 'make BITS=64 bowtie'.  To determine whether\n"
	"  your version of bowtie is 64-bit or 32-bit, run 'bowtie --version'.\n"
	"  \n"
	"  Tip 2: If your computer has multiple processors/cores, try -p\n"
	"   \n"
	"  The -p <int> option causes Bowtie to launch <int> parallel search\n"
	"  threads.  Each thread runs on a different processor/core and all\n"
	"  threads find alignments in parallel, increasing alignment throughput\n"
	"  by approximately a multiple of <int>.\n"
	"  \n"
	"  Tip 3: If reporting many alignments per read, try tweaking\n"
	"         'bowtie-build --offrate'\n"
	"   \n"
	"  If you are using the -k, -a or -m options and Bowtie is reporting\n"
	"  many alignments per read (an average of more than about 10 per read)\n"
	"  and you have some physical memory to spare, then consider building\n"
	"  an index with a denser SA sample.\n"
	"  \n"
	"  To build an index with a denser SA sample, specify a smaller\n"
	"  --offrate value when running 'bowtie-build'.  A denser SA sample\n"
	"  leads to a larger index, but is also particularly effective at\n"
	"  speeding up alignment when then number of alignments reported per\n"
	"  read is large.  For example, if the number of alignments per read is\n"
	"  very large, decreasing the index's --offrate by 1 could as much as\n"
	"  double alignment performance, and decreasing by 2 could quadruple\n"
	"  alignment performance, etc.\n"
	"  \n"
	"  On the other hand, decreasing --offrate increases the size of the\n"
	"  Bowtie index, both on disk and in memory when aligning reads.  At the\n"
	"  default --offrate of 5, the SA sample for the human genome occupies\n"
	"  about 375 MB of memory when aligning reads.  Decreasing the --offrate\n"
	"  by 1 doubles the memory taken by the SA sample, and decreasing by 2\n"
	"  quadruples the memory taken, etc.\n"
	"  \n"
	"  Tip 4: If bowtie \"thrashes\", try tweaking 'bowtie --offrate'\n"
	"  \n"
	"  If 'bowtie' is very slow and consistently triggers more than a few\n"
	"  page faults per second (as observed via top or vmstat on Mac/Linux,\n"
	"  or via a tool like Process Explorer on Windows), then try giving\n"
	"  bowtie the --offrate <int> option with a larger <int> value than the\n"
	"  value used when building the index.  For example, bowtie-build's\n"
	"  default --offrate is 5 and all pre-built indexes available from the\n"
	"  Bowtie website are built with --offrate 5; so if bowtie thrashes when\n"
	"  querying such an index, try using 'bowtie --offrate 6'.  If bowtie\n"
	"  still thrashes, try 'bowtie --offrate 7', etc.  A higher --offrate\n"
	"  causes bowtie to use a sparser sample of the suffix-array than is\n"
	"  stored in the index; this saves memory but makes alignment reporting\n"
	"  slower (which is especially slow when using -a or large -k).\n"
	"\n"
	"  Command Line\n"
	"  ------------\n"
	"\n"
	"  The following is a detailed description of the options used to control\n"
	"  the 'bowtie' aligner:\n"
	"\n"
	" Usage:\n"
	" \n"
	"  bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]\n"
	"\n"
	"  <ebwt>             The basename of the index to be searched.  The\n"
	"                     basename is the name of any of the index files up\n"
	"                     to but not including the final .1.ebwt /\n"
	"                     .rev.1.ebwt / etc.  bowtie looks for the specified\n"
	"                     index first in the current directory, then in the\n"
	"                     'indexes' subdirectory under the directory where\n"
	"                     the currently-running 'bowtie' executable is\n"
	"                     located, then looks in the directory specified in\n"
	"                     the BOWTIE_INDEXES environment variable.\n"
	"\n"
	"  <m1>               Comma-separated list of files containing the #1\n"
	"                     mates (filename usually includes \"_1\"), or, if -c\n"
	"                     is specified, the mate sequences themselves.\n"
	"                     E.g., this might be \"flyA_1.fq,flyB_1.fq\", or, if\n"
	"                     -c is given, this might be \"GGTCATCCT,ACGGGTCGT\".\n"
	"                     Sequences specified with this option must\n"
	"                     correspond file-for-file and read-for-read with\n"
	"                     those specified in <m2>.  Reads may be a mix of\n"
	"                     different lengths.  If \"-\" is specified, Bowtie\n"
	"                     will read the #1 mates from stdin.\n"
	"\n"
	"  <m2>               Comma-separated list of files containing the #2\n"
	"                     mates (filename usually includes \"_2\"), or, if -c\n"
	"                     is specified, the mate sequences themselves.\n"
	"                     E.g., this might be \"flyA_2.fq,flyB_2.fq\", or, if\n"
	"                     -c is given, this might be \"GTATGCTG,AATTCAGGCTG\".\n"
	"                     Sequences specified with this option must\n"
	"                     correspond file-for-file and read-for-read with\n"
	"                     those specified in <m1>.  Reads may be a mix of\n"
	"                     different lengths.  If \"-\" is specified, Bowtie\n"
	"                     will read the #2 mates from stdin.\n"
	"\n"
	"  <r>                Comma-separated list of files containing a mix of\n"
	"                     unpaired and paired-end reads in Tab-delimited\n"
	"                     format.  Tab-delimited format is a 1-read-per-line\n"
	"                     format where unpaired reads consist of a read\n"
	"                     name, sequence and quality string each separated\n"
	"                     by tabs.  A paired-end read consists of a read\n"
	"                     name, sequnce of the /1 mate, quality values of\n"
	"                     the /1 mate, sequence of the /2 mate, and quality\n"
	"                     values of the /2 mate separated by tabs.  Quality\n"
	"                     values can be expressed using any of the scales\n"
	"                     supported in FASTQ files.  Reads may be a mix of\n"
	"                     different lengths and paired-end and unpaired\n"
	"                     reads may be intermingled in the same file.  If\n"
	"                     \"-\" is specified, Bowtie will read the Tab-\n"
	"                     delimited reads from stdin.\n"
	"\n"
	"  <s>                A comma-separated list of files containing\n"
	"                     unpaired reads to be aligned, or, if -c is\n"
	"                     specified, the unpaired read sequences themselves.\n"
	"                     E.g., this might be\n"
	"                     \"lane1.fq,lane2.fq,lane3.fq,lane4.fq\", or, if -c\n"
	"                     is specified, this might be \"GGTCATCCT,ACGGGTCGT\".\n"
	"                     Reads may be a mix of different lengths.  If \"-\"\n"
	"                     is specified, Bowtie gets the reads from stdin.\n"
	"\n"
	"  <hit>              File to write alignments to.  By default,\n"
	"                     alignments are written to stdout (the console),\n"
	"                     but a <hits> file must be specified if the\n"
	"                     -b/--binout option is also specified.\n"
	"\n"
	" Options:\n"
	" ========\n"
	"\n"
	"   Input:\n"
	"   ------\n"
	"\n"
	"  -q                 The query input files (specified either as <m1>\n"
	"                     and <m2>, or as <s>) are FASTQ files (usually\n"
	"                     having extension .fq or .fastq).  This is the\n"
	"                     default.  See also: --solexa-quals and\n"
	"                     --integer-quals.\n"
	"\n"
	"  -f                 The query input files (specified either as <m1>\n"
	"                     and <m2>, or as <s>) are FASTA files (usually\n"
	"                     having extension .fa, .mfa, .fna or similar).  All\n"
	"                     quality values are assumed to be 40 on the Phred\n"
	"                     scale.\n"
	"\n"
	"  -r                 The query input files (specified either as <m1>\n"
	"                     and <m2>, or as <s>) are Raw files: one sequence\n"
	"                     per line, without quality values or names.  All\n"
	"                     quality values are assumed to be 40 on the Phred\n"
	"                     scale.\n"
	"\n"
	"  -c                 The query sequences are given on command line.\n"
	"                     I.e. <m1>, <m2> and <singles> are comma-separated\n"
	"                     lists of reads rather than lists of read files.\n"
	"\n"
	"  -s/--skip <int>    Skip (i.e. do not align) the first <int> reads or\n"
	"                     pairs in the input.\n"
	"\n"
	"  -u/--qupto <int>   Only align the first <int> reads or read pairs\n"
	"                     from the input (after the -s/--skip reads or pairs\n"
	"                     have been skipped).  Default: no limit.\n"
	"\n"
	"  -5/--trim5 <int>   Trim <int> bases from high-quality (left) end of\n"
	"                     each read before alignment (default: 0).\n"
	"\n"
	"  -3/--trim3 <int>   Trim <int> bases from low-quality (right) end of\n"
	"                     each read before alignment (default: 0).\n"
	"\n"
	"  --phred33-quals    Input qualities are ASCII chars equal to the Phred\n"
	"                     quality plus 33.  Default: on.\n"
	"\n"
	"  --phred64-quals    Input qualities are ASCII chars equal to the Phred\n"
	"                     quality plus 64.  Default: off.\n"
	"\n"
	"  --solexa-quals     Convert input qualities from solexa-scaled (which\n"
	"                     can be negative) to phred-scaled (which can't).\n"
	"                     The formula for conversion is phred-qual =\n"
	"                     10 * log(1 + 10 ** (solexa-qual/10.0)) / log(10).\n"
	"                     This is usually the right option for use with\n"
	"                     (unconverted) reads emitted by GA Pipeline\n"
	"                     versions prior to 1.3.  Default: off.\n"
	"\n"
	"  --solexa1.3-quals  Same as --phred64-quals.  This is usually the\n"
	"                     right option for use with (unconverted) reads\n"
	"                     emitted by GA Pipeline version 1.3 or later. \n"
	"                     Default: off.\n"
	"\n"
	"  --integer-quals    Quality values are represented in the read input\n"
	"                     file as space-separated ASCII integers, e.g.,\n"
	"                     \"40 40 30 40...\", rather than ASCII characters,\n"
	"                     e.g., \"II?I...\".  Integers are treated as being on\n"
	"                     the Phred scale unless --solexa-quals is also\n"
	"                     specified.  Default: off.\n"
	"\n"
	"   Alignment:\n"
	"   ----------\n"
	"\n"
	"  -n/--seedmms <int> The maximum number of mismatches permitted in the\n"
	"                     \"seed\", which is the first 28 base pairs of the\n"
	"                     read by default (see -l/--seedlen).  This may be\n"
	"                     0, 1, 2 or 3 and the default is 2.  This option is\n"
	"                     mutually exclusive with the -v option.\n"
	" \n"
	"  -e/--maqerr <int>  The maximum permitted total of quality values at\n"
	"                     mismatched read positions.  This total is also\n"
	"                     called the \"quality-weighted hamming distance\" or\n"
	"                     \"Q-distance.\"  This is analogous to the -e option\n"
	"                     for \"maq map\".  The default is 70.  Note that,\n"
	"                     like Maq, Bowtie rounds quality values to the\n"
	"                     nearest 10 and saturates at 30.\n"
	"  \n"
	"  -l/--seedlen <int> The \"seed length\"; i.e., the number of bases on\n"
	"                     the high-quality end of the read to which the -n\n"
	"                     ceiling applies.  The lowest permitted value for\n"
	"                     this parameter is 5.  The default is 28.\n"
	"                     Depending on the length of the reference genome,\n"
	"                     Bowtie can become extremely slow as this parameter\n"
	"                     is adjusted downward.\n"
	"\n"
	"  --nomaqround       Maq accepts quality values in the Phred scale, but\n"
	"                     internally rounds quality values to the nearest 10\n"
	"                     saturating at 30.  By default, Bowtie imitates\n"
	"                     this behavior.  Use --nomaqround to prevent this\n"
	"                     type of rounding in Bowtie.\n"
	"\n"
	"  -v <int>           Forego the Maq-like alignment policy and use a\n"
	"                     SOAP-like alignment policy.  I.e., report end-to-\n"
	"                     end alignments with at most <int> mismatches.  If\n"
	"                     -v is specified, base quality values and the -e\n"
	"                     and -l options are ignored.  -v is mutually\n"
	"                     exclusive with -n.\n"
	"\n"
	"  -I/--minins <int>  The minimum insert size for valid paired-end\n"
	"                     alignments.  E.g. if -I 60 is specified and a\n"
	"                     paired-end alignment consists of two 20-bp\n"
	"                     alignments in the appropriate orientation with a\n"
	"                     20-bp gap between them, that alignment is\n"
	"                     considered valid (as long as -X is also\n"
	"                     satisfied).  A 19-bp gap would not be valid in\n"
	"                     that case.  If trimming options -3 or -5 are also\n"
	"                     used, the -I constraint is applied with respect to\n"
	"                     the untrimmed mates, not the trimmed mates.\n"
	"                     Default: 0.\n"
	"\n"
	"  -X/--maxins <int>  The maximum insert size for valid paired-end\n"
	"                     alignments.  E.g. if -X 100 is specified and a\n"
	"                     paired-end alignment consists of two 20-bp\n"
	"                     alignments in the proper orientation with a 60-bp\n"
	"                     gap between them, that alignment is considered\n"
	"                     valid (as long as -I is also satisfied).  A 61-bp\n"
	"                     gap would not be valid in that case.  If trimming\n"
	"                     options -3 or -5 are also used, the -X constraint\n"
	"                     is applied with respect to the untrimmed mates,\n"
	"                     not the trimmed mates.  Default: 250.\n"
	"\n"
	"  --fr/--rf/--ff     The upstream/downstream mate orientations for a\n"
	"                     valid paired-end alignment against the forward\n"
	"                     reference strand.  E.g., if --fr is specified and\n"
	"                     there is a candidate paired-end alignment where\n"
	"                     mate1 appears upstream of the reverse complement\n"
	"                     of mate2 and the insert length constraints are\n"
	"                     met, that alignment is valid.  Also, if mate2\n"
	"                     appears upstream of the reverse complement of\n"
	"                     mate1 and all other constraints are met, that too\n"
	"                     is valid.  --rf likewise requires that an upstream\n"
	"                     mate1 be reverse-complemented and a downstream\n"
	"                     mate2 be forward-oriented.  --ff requires both an\n"
	"                     upstream mate1 and a downstream mate2 to be\n"
	"                     forward-oriented.  Default: --fr (appropriate for\n"
	"                     the Illumina short insert library).\n"
	"\n"
	"  --nofw/--norc      If --nofw is specified, Bowtie will not attempt to\n"
	"                     align against the forward reference strand.  If\n"
	"                     --norc is specified, Bowtie will not attempt to\n"
	"                     align against the reverse-complement reference\n"
	"                     strand.  For paired-end reads using --fr or --rf\n"
	"                     modes, --nofw and --norc apply to the forward and\n"
	"                     reverse-complement pair orientations.  I.e.\n"
	"                     specifying --nofw and --fr will only find reads in\n"
	"                     the R/F orientation where mate 2 occurs upstream\n"
	"                     of mate 1 with respect to the forward reference\n"
	"                     strand.\n"
	"\n"
	"  --maxbts           The maximum number of backtracks permitted when\n"
	"                     aligning a read in -n 2 or -n 3 mode (default:\n"
	"                     125 without --best, 800 with --best).  A\n"
	"                     \"backtrack\" is the introduction of a speculative\n"
	"                     substitution into the alignment.  Without this\n"
	"                     limit, the default parameters will sometimes\n"
	"                     require that 'bowtie' try 100s or 1,000s of\n"
	"                     backtracks to align a read, especially if the read\n"
	"                     has many low-quality bases and/or has no valid\n"
	"                     alignments, slowing bowtie down significantly.\n"
	"                     However, this limit may cause some valid\n"
	"                     alignments to be missed.  Higher limits yield\n"
	"                     greater sensitivity at the expensive of longer\n"
	"                     running times.  See also: -y/--tryhard.\n"
	"\n"
	"  --pairtries <int>  For paired-end alignment, this is the maximum\n"
	"                     number of attempts Bowtie will make to match an\n"
	"                     alignment for one mate up with an alignment for\n"
	"                     the opposite mate.  Most paired-end alignments\n"
	"                     require only a few such attempts, but pairs where\n"
	"                     both mates occur in highly repetitive regions of\n"
	"                     the reference can require significantly more.\n"
	"                     Setting this to a higher number allows Bowtie to\n"
	"                     find more paired-end alignments for repetitive\n"
	"                     pairs at the expense of speed.  The default is\n"
	"                     100.  See also: -y/--tryhard.\n"
	"\n"
	"  -y/--tryhard       Try as hard as possible to find valid alignments\n"
	"                     when they exist, including paired-end alignments.\n"
	"                     This is equivalent to specifying very high values\n"
	"                     for the --maxbts and --pairtries options.  This\n"
	"                     mode is generally MUCH SLOWER than the default\n"
	"                     settings, but can be useful for certain research\n"
	"                     problems.  This mode is slower when (a) the\n"
	"                     reference is very repetitive, (b) the reads are\n"
	"                     low quality, or (c) not many reads have valid\n"
	"                     alignments.\n"
	"\n"
	"  --chunkmbs <int>   The number of megabytes of memory a given thread\n"
	"                     is given to store path descriptors in --best mode.\n"
	"                     Best-first search must keep track of many paths at\n"
	"                     once to ensure it is always extending the path\n"
	"                     with the lowest cumulative cost.  Bowtie tries to\n"
	"                     minimize the memory impact of the descriptors, but\n"
	"                     they can still grow very large in some cases.  If\n"
	"                     you receive an error message saying that chunk\n"
	"                     memory has been exhausted in --best mode, try\n"
	"                     adjusting this parameter up to dedicate more\n"
	"                     memory to the descriptors.  Default: 32.\n"
	"\n"
	"   Reporting:\n"
	"   ----------\n"
	"\n"
	"  -k <int>           Report up to <int> valid alignments per read or\n"
	"                     pair (default: 1).  Validity of alignments is\n"
	"                     determined by the alignment policy (combined\n"
	"                     effects of -n, -v, -l, and -e).  If more than one\n"
	"                     valid alignment exists and the --best and --strata\n"
	"                     options are specified, then only those alignments\n"
	"                     belonging to the best alignment \"stratum\" (i.e.\n"
	"                     those with the fewest mismatches) will be\n"
	"                     reported.  Bowtie is designed to be very fast for\n"
	"                     small -k but bowtie can become significantly\n"
	"                     slower as -k increases.  If you would like to use\n"
	"                     Bowtie for larger values of -k, consider building\n"
	"                     an index with a denser suffix-array sample, i.e.\n"
	"                     specify a smaller '--offrate' when invoking\n"
	"                     'bowtie-build' for the relevant index (see\n"
	"                     Performance Tips section for details).\n"
	"\n"
	"  -a/--all           Report all valid alignments per read or pair\n"
	"                     (default: off).  Validity of alignments is\n"
	"                     determined by the alignment policy (combined\n"
	"                     effects of -n, -v, -l, and -e).  If more than one\n"
	"                     valid alignment exists and the --best and --strata\n"
	"                     options are specified, then only those alignments\n"
	"                     belonging to the best alignment \"stratum\" (i.e.\n"
	"                     those with the fewest mismatches) will be\n"
	"                     reported.  Bowtie is designed to be very\n"
	"                     fast for small -k but bowtie can become\n"
	"                     significantly slower if -a/--all is specified.  If\n"
	"                     you would like to use Bowtie with -a, consider\n"
	"                     building an index with a denser suffix-array\n"
	"                     sample, i.e. specify a smaller '--offrate' when\n"
	"                     invoking 'bowtie-build' for the relevant index\n"
	"                     (see Performance Tips section for details).\n"
	"\n"
	"  -m <int>           Suppress all alignments for a particular read or\n"
	"                     pair if more than <int> reportable alignments\n"
	"                     exist for it.  Reportable alignments are those\n"
	"                     that would be reported given the -n, -v, -l, -e,\n"
	"                     -k, -a, --best, and --strata options.  Default:\n"
	"                     no limit.  Bowtie is designed to be very fast for\n"
	"                     small -m but bowtie can become significantly\n"
	"                     slower for larger values of -m.    If you would\n"
	"                     like to use Bowtie for larger values of -k,\n"
	"                     consider building an index with a denser suffix-\n"
	"                     array sample, i.e. specify a smaller '--offrate'\n"
	"                     when invoking 'bowtie-build' for the relevant\n"
	"                     index (see Performance Tips section for details).\n"
	"\n"
	"  --best             Make Bowtie guarantee that reported singleton\n"
	"                     alignments are \"best\" in terms of stratum (i.e.\n"
	"                     number of mismatches, or mismatches in the seed in\n"
	"                     the case of -n mode) and in terms of the quality\n"
	"                     values at the mismatched position(s).  Stratum\n"
	"                     always trumps quality; e.g. a 1-mismatch alignment\n"
	"                     where the mismatched position has Phred quality 40\n"
	"                     is preferred over a 2-mismatch alignment where the\n"
	"                     mismatched positions both have Phred quality 10.\n"
	"                     When --best is not specified, Bowtie may report\n"
	"                     alignments that are sub-optimal in terms of\n"
	"                     stratum and/or quality (though an effort is made\n"
	"                     to report the best alignment).  --best mode also\n"
	"                     removes all strand bias.  Note that --best does\n"
	"                     not affect which alignments are considered \"valid\"\n"
	"                     by Bowtie, only which valid alignments are\n"
	"                     reported by Bowtie.  When --best is specified and\n"
	"                     multiple hits are allowed (via -k or -a), the\n"
	"                     alignments for a given read are guaranteed to\n"
	"                     appear in best-to-worst order in Bowtie's output.\n"
	"                     Bowtie is about 1-2.5 times slower when --best is\n"
	"                     specified.\n"
	"\n"
	"  --strata           If many valid alignments exist and are reportable\n"
	"                     (e.g. are not disallowed via the -k option) and\n"
	"                     they fall into more than one alignment \"stratum\",\n"
	"                     report only those alignments that fall into the\n"
	"                     best stratum.  By default, Bowtie reports all\n"
	"                     reportable alignments regardless of whether they\n"
	"                     fall into multiple strata.  When --strata is\n"
	"                     specified, --best must also be specified. \n"
	"\n"
	"   Output:\n"
	"   -------\n"
	"\n"
	"  -S/--sam           Print alignments in SAM format.  See the SAM\n"
	"                     format spec (http://samtools.sf.net/) and the \"SAM\n"
	"                     output\" section of the manual for details.  To\n"
	"                     suppress all SAM headers, use --sam-nohead.  To\n"
	"                     suppress just the @SQ headers (e.g. if the\n"
	"                     alignment is against a very large number of\n"
	"                     reference sequences), use --sam-nosq.  Bowtie\n"
	"                     does not write BAM files directly, but SAM output\n"
	"                     can be converted to BAM on the fly by piping\n"
	"                     Bowtie's output to 'samtools view'.  -S/--sam is\n"
	"                     not compatible with --refout.\n"
	"\n"
	"  --concise          Print alignments in a concise format. Each line\n"
	"                     has format 'read_idx{-|+}:<ref_idx,ref_off,mms>',\n"
	"                     where read_idx is the index of the read mapped,\n"
	"                     {-|+} is the orientation of the read, ref_idx is\n"
	"                     the index of the reference sequence aligned to,\n"
	"                     ref_off is the offset into the reference sequence,\n"
	"                     and mms is the number of mismatches in the\n"
	"                     alignment.  Each alignment appears on a separate\n"
	"                     line.\n"
	"\n"
	"  -t/--time          Print the amount of wall-clock time taken by each\n"
	"                     phase.\n"
	"\n"
	"  -B/--offbase <int> When outputting alignments, number the first base\n"
	"                     of a reference sequence as <int>.  Default: 0.\n"
	"                     (Default is likely to change to 1 in Bowtie 1.0.)\n"
	"\n"
	"  --quiet            Print nothing besides alignments.\n"
	"\n"
	"  --refout           Write alignments to a set of files named\n"
	"                     refXXXXX.map, where XXXXX is the 0-padded index of\n"
	"                     the reference sequence aligned to.  This can be a\n"
	"                     useful way to break up work for downstream\n"
	"                     analyses when dealing with, for example, large\n"
	"                     numbers of reads aligned to the assembled human\n"
	"                     genome.  If <hits> is also specified, it will be\n"
	"                     ignored.\n"
	"\n"
	"  --refidx           When a reference sequence is referred to in a\n"
	"                     reported alignment, refer to it by 0-based index\n"
	"                     (its offset into the list of references that were\n"
	"                     indexed) rather than by name.\n"
	"\n"
	"  --al <filename>    Write all reads for which at least one alignment\n"
	"                     was reported to a file with name <filename>.\n"
	"                     Written reads will appear as they did in the\n"
	"                     input, without any of the trimming or translation\n"
	"                     of quality values that may have taken place within\n"
	"                     Bowtie.  Paired-end reads will be written to two\n"
	"                     parallel files with \"_1\" and \"_2\" inserted in the\n"
	"                     filename, e.g., if <filename> is aligned.fq, the\n"
	"                     #1 and #2 mates that fail to align will be written\n"
	"                     to aligned_1.fq and aligned_2.fq respectively.\n"
	"\n"
	"  --un <filename>    Write all reads that could not be aligned to a\n"
	"                     file with name <filename>.  Written reads will\n"
	"                     appear as they did in the input, without any of\n"
	"                     the trimming or translation of quality values that\n"
	"                     may have taken place within Bowtie.  Paired-end\n"
	"                     reads will be written to two parallel files with\n"
	"                     \"_1\" and \"_2\" inserted in the filename, e.g., if\n"
	"                     <filename> is unaligned.fq, the #1 and #2 mates\n"
	"                     that fail to align will be written to\n"
	"                     unaligned_1.fq and unaligned_2.fq respectively.\n"
	"                     Unless --max is also specified, reads with a\n"
	"                     number of valid alignments exceeding the limit set\n"
	"                     with the -m option are also written to <filename>.\n"
	"\n"
	"  --max <filename>   Write all reads with a number of valid alignments\n"
	"                     exceeding the limit set with the -m option to a\n"
	"                     file with name <filename>.  Written reads will\n"
	"                     appear as they did in the input, without any of\n"
	"                     the trimming or translation of quality values that\n"
	"                     may have taken place within Bowtie.  Paired-end\n"
	"                     reads will be written to two parallel files with\n"
	"                     \"_1\" and \"_2\" inserted in the filename, e.g., if\n"
	"                     <filename> is max.fq, the #1 and #2 mates\n"
	"                     that fail to align will be written to\n"
	"                     max_1.fq and max_2.fq respectively.  These reads\n"
	"                     are not written to the file specified with --un.\n"
	"\n"
	"  --sam-nohead       Suppress header lines (starting with @) when\n"
	"                     output is SAM.\n"
	"\n"
	"  --sam-nosq         Suppress @SQ header lines when output is SAM.\n"
	"\n"
	"  --fullref          Print the full refernce sequence name, including\n"
	"                     whitespace, in alignment output.  By default\n"
	"                     bowtie prints everything up to but not including\n"
	"                     the first whitespace.\n"
	"   Performance:\n"
	"   ------------\n"
	"\n"
	"  -p/--threads <int> Launch <int> parallel search threads (default: 1).\n"
	"                     Threads will run on separate processors/cores and\n"
	"                     synchronize when grabbing reads and outputting\n"
	"                     alignments.  Searching for alignments is highly\n"
	"                     parallel, and speedup is fairly close to linear.\n"
	"                     This option is only available if bowtie is linked\n"
	"                     with the pthreads library (i.e. if\n"
	"                     BOWTIE_PTHREADS=0 is not specified at build time).\n"
	"\n"
	"  -o/--offrate <int> Override the offrate of the index with <int>.  If\n"
	"                     <int> is greater than the offrate used to build\n"
	"                     the index, then some row markings are discarded\n"
	"                     when the index is read into memory.  This reduces\n"
	"                     the memory footprint of the aligner but requires\n"
	"                     more time to calculate text offsets.  <int> must\n"
	"                     be greater than the value used to build the index.\n"
	"\n"
	"  --mm               Use memory-mapped I/O to load the index, rather\n"
	"                     than normal POSIX/C file I/O.  Memory-mapping the\n"
	"                     index allows many concurrent bowtie processes on\n"
	"                     the same computer to share the same memory image\n"
	"                     of the index (i.e. you pay the memory overhead\n"
	"                     just once).  This facilitates memory-efficient\n"
	"                     parallelization of Bowtie in situations where\n"
	"                     using -p is not desirable.\n"
	"\n"
	"  --shmem            Use shared memory to load the index, rather than\n"
	"                     normal POSIX/C file I/O.  Using shared memory\n"
	"                     allows many concurrent bowtie processes on\n"
	"                     the same computer to share the same memory image\n"
	"                     of the index (i.e. you pay the memory overhead\n"
	"                     just once).  This facilitates memory-efficient\n"
	"                     parallelization of Bowtie in situations where\n"
	"                     using -p is not desirable.  Unlike --mm, --shmem\n"
	"                     installs the index into shared memory permanently,\n"
	"                     or until the user deletes the shared memory chunks\n"
	"                     manually.  See your operating system documentation\n"
	"                     for details on how to manually list and remove\n"
	"                     shared memory chunks (on Linux and Mac OS X, these\n"
	"                     commands are ipcs and ipcrm).  You may also need\n"
	"                     to increase your OS's maximum shared-memory chunk\n"
	"                     size to accomodate larger indexes; see your OS\n"
	"                     documentation.\n"
	"\n"
	"   Other:\n"
	"   ------\n"
	"\n"
	"  --seed <int>       Use <int> as the seed for pseudo-random number\n"
	"                     generator.\n"
	"\n"
	"  --verbose          Print verbose output (for debugging).\n"
	"\n"
	"  --version          Print version information and quit.\n"
	"\n"
	"  -h/--help          Print detailed description of tool and its options\n"
	"                     (from MANUAL).\n"
	"\n"
	"  Default output\n"
	"  --------------\n"
	"\n"
	"  The 'bowtie' aligner outputs each alignment on a separate line.  Each\n"
	"  line is a collection of 8 fields separated by tabs; from left to\n"
	"  right, the fields are:\n"
	"\n"
	"   1. Name of read that aligned\n"
	"\n"
	"   2. Orientation of read in the alignment, '-' for reverse complement,\n"
	"      '+' otherwise\n"
	"\n"
	"   3. Name of reference sequence where alignment occurs, or ordinal ID\n"
	"      if no name was provided\n"
	"\n"
	"   4. 0-based offset into the forward reference strand where leftmost\n"
	"      character of the alignment occurs\n"
	"\n"
	"   5. Read sequence (reverse-complemented if orientation is '-')\n"
	"\n"
	"   6. ASCII-encoded read qualities (reversed if orientation is '-').\n"
	"      The encoded quality values are on the Phred scale and the\n"
	"      encoding is ASCII-offset by 33 (ASCII char '!'). \n"
	"\n"
	"   7. Number of other instances where the same read aligns against the\n"
	"      same reference characters as were aligned against in this\n"
	"      alignment.  This is *not* the number of other places the read\n"
	"      aligns with the same number of mismatches.  The number in this\n"
	"      column is generally not a good proxy for that number (e.g., the\n"
	"      number in this column may be '0' while the number of other\n"
	"      alignments with the same number of mismatches might be large).\n"
	"      This column was previously described as \"Reserved\".\n"
	"\n"
	"   8. Comma-separated list of mismatch descriptors.  If there are no\n"
	"      mismatches in the alignment, this field is empty.  A single\n"
	"      descriptor has the format offset:reference-base>read-base.  The\n"
	"      offset is expressed as a 0-based offset from the high-quality\n"
	"      (5') end of the read. \n"
	"\n"
	"  SAM output\n"
	"  ----------\n"
	"\n"
	"  Following is a brief description of the SAM format as output by\n"
	"  Bowtie when the -S/--sam option is specified.  For more details, see\n"
	"  the SAM format specification (http://samtools.sf.net/SAM1.pdf).\n"
	"\n"
	"  When -S/--sam is specified, Bowtie will always print a SAM header\n"
	"  with @HD, @SQ and @PG lines.\n"
	"\n"
	"  Each subsequnt line corresponds to a read or an alignment.  Each line\n"
	"  is a collection of at least 12 fields separated by tabs; from left to\n"
	"  right, the fields are:\n"
	"\n"
	"   1. Name of read that aligned\n"
	"\n"
	"   2. Sum of all applicable flags.  Flags relevant to Bowtie are:\n"
	"\n"
	"        1: The read is one of a pair\n"
	"        2: The alignment is one end of a proper paired-end alignment\n"
	"        4: The read has no reported alignments\n"
	"        8: The read is one of a pair and has no reported alignments\n"
	"       16: The alignment is to the reverse reference strand\n"
	"       32: The other mate in the paired-end alignment is aligned to the\n"
	"           reverse reference strand\n"
	"       64: The read is the first (#1) mate in a pair\n"
	"      128: The read is the second (#2) mate in a pair\n"
	"\n"
	"      Thus, an unpaired read that aligns to the reverse reference\n"
	"      strand will have flag 16.  A paired-end read that aligns and is\n"
	"      the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1).\n"
	"\n"
	"   3. Name of reference sequence where alignment occurs, or ordinal ID\n"
	"      if no name was provided\n"
	"\n"
	"   4. 1-based offset into the forward reference strand where leftmost\n"
	"      character of the alignment occurs\n"
	"\n"
	"   5. Mapping quality (always 255)\n"
	"\n"
	"   6. CIGAR string representation of alignment\n"
	"\n"
	"   7. Name of reference sequence where mate's alignment occurs.  Set to\n"
	"      \"=\" if the mate's reference sequence is the same as this\n"
	"      alignment's, or \"*\" if there is no mate.\n"
	"\n"
	"   8. 1-based offset into the forward reference strand where leftmost\n"
	"      character of the mate's alignment occurs.  Offset is 0 if there\n"
	"      is no mate.\n"
	"\n"
	"   9. Inferred insert size.  Size is negative if the mate's alignment\n"
	"      occurs upstream of this alignment.  Size is 0 if there is no\n"
	"      mate.\n"
	"\n"
	"  10. Read sequence (reverse-complemented if aligned to the reverse\n"
	"      strand)\n"
	"\n"
	"  11. ASCII-encoded read qualities (reversed if aligned to the reverse\n"
	"      strand).  The encoded quality values are on the Phred scale and\n"
	"      the encoding is ASCII-offset by 33 (ASCII char '!'). \n"
	"\n"
	" 12+. Optional fields.  Fields are tab-separated.  For descriptions of\n"
	"      all possible optional fields, see the SAM format specification.\n"
	"      Bowtie outputs one or more of these optional fields for each\n"
	"      alignment, depending on the type of the alignment:\n"
	"\n"
	"      NM:i:<N> : Aligned read has an edit distance of <N>\n"
	"      MD:Z:<S> : For aligned reads, <S> is a string representation of\n"
	"                 the mismatched reference bases in the alignment.  See\n"
	"                 SAM format specification for details.\n"
	"      XA:i:<N> : Aligned read belongs to stratum <N>\n"
	"      XM:i:<N> : For a read with no reported alignments, <N> is 0 if\n"
	"                 the read had no alignments, or 1 if the read had\n"
	"                 alignments that were suppressed by the -m option.\n"
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
 * Parse a T string 'str'.
 */
template<typename T>
T parse(const char *s) {
	T tmp;
	stringstream ss(s);
	ss >> tmp;
	return tmp;
}

/**
 * Parse a pair of Ts from a string, 'str', delimited with 'delim'.
 */
template<typename T>
pair<T, T> parsePair(const char *str, char delim) {
	string s(str);
	vector<string> ss;
	tokenize(s, delim, ss);
	pair<T, T> ret;
	ret.first = parse<T>(ss[0].c_str());
	ret.second = parse<T>(ss[1].c_str());
	return ret;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	if(startVerbose) { cerr << "Parsing options: "; logTime(cerr, true); }
	do {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		switch (next_option) {
			case '1': tokenize(optarg, ",", mates1); break;
			case '2': tokenize(optarg, ",", mates2); break;
			case ARG_ONETWO: tokenize(optarg, ",", mates12); format = TAB_MATE; break;
			case 'f': format = FASTA; break;
			case 'F': {
				format = FASTA_CONT;
				pair<size_t, size_t> p = parsePair<size_t>(optarg, ',');
				fastaContLen = p.first;
				fastaContFreq = p.second;
				break;
			}
			case 'q': format = FASTQ; break;
			case 'r': format = RAW; break;
			case 'c': format = CMDLINE; break;
			case ARG_CHAININ: format = INPUT_CHAIN; break;
			case 'I':
				minInsert = (uint32_t)parseInt(0, "-I arg must be positive");
				break;
			case 'X':
				maxInsert = (uint32_t)parseInt(1, "-X arg must be at least 1");
				break;
			case 's':
				skipReads = (uint32_t)parseInt(0, "-s arg must be positive");
				break;
			case ARG_FF: mate1fw = true;  mate2fw = true;  break;
			case ARG_RF: mate1fw = false; mate2fw = true;  break;
			case ARG_FR: mate1fw = true;  mate2fw = false; break;
			case ARG_RANDOM_READS: format = RANDOM; break;
			case ARG_RANDOM_READS_NOSYNC:
				format = RANDOM;
				randReadsNoSync = true;
				break;
			case ARG_RANGE: rangeMode = true; break;
			case ARG_CONCISE: outType = OUTPUT_CONCISE; break;
			case ARG_CHAINOUT: outType = OUTPUT_CHAIN; break;
			case 'b': outType = OUTPUT_BINARY; break;
			case 'S': outType = OUTPUT_SAM; break;
			case ARG_REFOUT: refOut = true; break;
			case ARG_NOOUT: outType = OUTPUT_NONE; break;
			case ARG_REFMAP: refMapFile = optarg; break;
			case ARG_ANNOTMAP: annotMapFile = optarg; break;
			case ARG_USE_SPINLOCK: useSpinlock = false; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_MM: {
#ifdef BOWTIE_MM
				useMm = true;
				break;
#else
				cerr << "Memory-mapped I/O mode is disabled because bowtie was not compiled with" << endl
				     << "BOWTIE_MM defined.  Memory-mapped I/O is not supported under Windows.  If you" << endl
				     << "would like to use memory-mapped I/O on a platform that supports it, please" << endl
				     << "refrain from specifying BOWTIE_MM=0 when compiling Bowtie." << endl;
				throw 1;
#endif
			}
			case ARG_MMSWEEP: mmSweep = true; break;
			case ARG_HADOOPOUT: hadoopOut = true; break;
			case ARG_AL: dumpAlBase = optarg; break;
			case ARG_UN: dumpUnalBase = optarg; break;
			case ARG_MAXDUMP: dumpMaxBase = optarg; break;
			case ARG_SOLEXA_QUALS: solexaQuals = true; break;
			case ARG_integerQuals: integerQuals = true; break;
			case ARG_PHRED64: phred64Quals = true; break;
			case ARG_PHRED33: solexaQuals = false; phred64Quals = false; break;
			case ARG_FORGIVE_INPUT: forgiveInput = true; break;
			case ARG_NOMAQROUND: noMaqRound = true; break;
			case 'z': fullIndex = false; break;
			case ARG_REFIDX: noRefNames = true; break;
			case ARG_STATEFUL: stateful = true; break;
			case ARG_FUZZY: fuzzy = true; break;
			case ARG_REPORTSE: reportSe = true; break;
			case ARG_FULLREF: fullRef = true; break;
			case ARG_PREFETCH_WIDTH:
				prefetchWidth = parseInt(1, "--prewidth must be at least 1");
				break;
			case 'B':
				offBase = parseInt(-999999, "-B/--offbase cannot be a large negative number");
				break;
			case ARG_SEED:
				seed = parseInt(0, "--seed arg must be at least 0");
				break;
			case 'u':
				qUpto = (uint32_t)parseInt(1, "-u/--qupto arg must be at least 1");
				break;
			case 'k':
				khits = (uint32_t)parseInt(1, "-k arg must be at least 1");
				break;
			case 'm':
				mhits = (uint32_t)parseInt(1, "-m arg must be at least 1");
				break;
			case 'x':
				mixedThresh = (uint32_t)parseInt(0, "-x arg must be at least 0");
				break;
			case ARG_MIXED_ATTEMPTS:
				mixedAttemptLim = (uint32_t)parseInt(1, "--mixatt arg must be at least 1");
				break;
			case ARG_CACHE_LIM:
				cacheLimit = (uint32_t)parseInt(1, "--cachelim arg must be at least 1");
				break;
			case ARG_CACHE_SZ:
				cacheSize = (uint32_t)parseInt(1, "--cachesz arg must be at least 1");
				cacheSize *= (1024 * 1024); // convert from MB to B
				break;
			case ARG_NO_RECONCILE:
				dontReconcileMates = true;
				break;
			case 'p':
#ifndef BOWTIE_PTHREADS
				cerr << "-p/--threads is disabled because bowtie was not compiled with pthreads support" << endl;
				throw 1;
#endif
				nthreads = parseInt(1, "-p/--threads arg must be at least 1");
				break;
			case ARG_FILEPAR:
#ifndef BOWTIE_PTHREADS
				cerr << "--filepar is disabled because bowtie was not compiled with pthreads support" << endl;
				throw 1;
#endif
				fileParallel = true;
				break;
			case 'v':
				maqLike = 0;
				mismatches = parseInt(0, "-v arg must be at least 0");
				if(mismatches > 3) {
					cerr << "-v arg must be at most 3" << endl;
					throw 1;
				}
				break;
			case '3': trim3 = parseInt(0, "-3/--trim3 arg must be at least 0"); break;
			case '5': trim5 = parseInt(0, "-5/--trim5 arg must be at least 0"); break;
			case 'o': offRate = parseInt(1, "-o/--offrate arg must be at least 1"); break;
			case ARG_ISARATE: isaRate = parseInt(0, "--isarate arg must be at least 0"); break;
			case 'e': qualThresh = parseInt(1, "-e/--err arg must be at least 1"); break;
			case 'n': seedMms = parseInt(0, "-n/--seedmms arg must be at least 0"); maqLike = 1; break;
			case 'l': seedLen = parseInt(5, "-l/--seedlen arg must be at least 5"); break;
			case 'h': printLongUsage(cout); throw 0; break;
			case ARG_USAGE: printUsage(cout); throw 0; break;
			case 'a': allHits = true; break;
			case 'y': tryHard = true; break;
			case ARG_RECAL: recal = true; break;
			case ARG_CHUNKMBS: chunkPoolMegabytes = parseInt(1, "--chunkmbs arg must be at least 1"); break;
			case ARG_CHUNKSZ: chunkSz = parseInt(1, "--chunksz arg must be at least 1"); break;
			case ARG_CHUNKVERBOSE: chunkVerbose = true; break;
			case ARG_BETTER: stateful = true; better = true; oldBest = false; break;
			case ARG_OLDBEST: oldBest = true; stateful = false; break;
			case ARG_BEST: stateful = true; useV1 = false; oldBest = false; break;
			case ARG_STRATA: strata = true; break;
			case ARG_VERBOSE: verbose = true; break;
			case ARG_STARTVERBOSE: startVerbose = true; break;
			case ARG_QUIET: quiet = true; break;
			case ARG_SANITY: sanityCheck = true; break;
			case 't': timing = true; break;
			case ARG_NO_FW: nofw = true; break;
			case ARG_NO_RC: norc = true; break;
			case ARG_STATS: stats = true; break;
			case ARG_PEV2: useV1 = false; break;
			case ARG_SAM_NOHEAD: samNoHead = true; break;
			case ARG_SAM_NOSQ: samNoSQ = true; break;
			case ARG_MAXBTS: {
				maxBts  = parseInt(0, "--maxbts must be positive");
				maxBtsBetter = maxBts;
				break;
			}
			case ARG_DUMP_PATS: patDumpfile = optarg; break;
			case ARG_STRAND_FIX: strandFix = true; break;
			case ARG_RANDOMIZE_QUALS: randomizeQuals = true; break;
			case ARG_PARTITION: partitionSz = parseInt(0, "--partition must be positive"); break;
			case ARG_ORIG:
				if(optarg == NULL || strlen(optarg) == 0) {
					cerr << "--orig arg must be followed by a string" << endl;
					printUsage(cerr);
					throw 1;
				}
				origString = optarg;
				break;
			case -1: break; /* Done with options. */
			case 0:
				if (long_options[option_index].flag != 0)
					break;
			default:
				printUsage(cerr);
				throw 1;
		}
	} while(next_option != -1);
	bool paired = mates1.size() > 0 || mates2.size() > 0 || mates12.size() > 0;
	if(rangeMode) {
		// Tell the Ebwt loader to ignore the suffix-array portion of
		// the index.  We don't need it because the user isn't asking
		// for bowtie to report reference positions (just matrix
		// ranges).
		offRate = 32;
	}
	if(maqLike) {
		revcomp = true;
	} else if(mismatches == 3 && fullIndex) {
		// Much faster than normal 3-mismatch mode
		stateful = true;
	}
	if(mates1.size() != mates2.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
		     << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
		     << "sequences must be specified with -1 and -2." << endl;
		throw 1;
	}
	// Check for duplicate mate input files
	if(format != CMDLINE) {
		for(size_t i = 0; i < mates1.size(); i++) {
			for(size_t j = 0; j < mates2.size(); j++) {
				if(mates1[i] == mates2[j] && !quiet) {
					cerr << "Warning: Same mate file \"" << mates1[i] << "\" appears as argument to both -1 and -2" << endl;
				}
			}
		}
	}
	// Check for all the parameter combinations that aren't compatible
	// with -z/--phased mode.
	if(!fullIndex) {
		bool error = false;
		if(khits > 1) {
			cerr << "When -z/--phased is used, -k X for X > 1 is unavailable" << endl;
			error = true;
		}
		if(mhits != 0xffffffff) {
			cerr << "When -z/--phased is used, -m is unavailable" << endl;
			error = true;
		}
		if(oldBest) {
			cerr << "When -z/--phased is used, --oldbest is unavailable" << endl;
			error = true;
		}
		if(stateful && better) {
			cerr << "When -z/--phased is used, --better is unavailable" << endl;
			error = true;
		}
		else if(stateful) {
			cerr << "When -z/--phased is used, --best is unavailable" << endl;
			error = true;
		}
		if(allHits && strata) {
			cerr << "When -a/--all and -z/--phased are used, --strata cannot also be used." << endl
			     << "Stratified all-hits search cannot be combined with phased search." << endl;
			error = true;
		}
		if(paired) {
			cerr << "When -z/--phased is used, paired-end mode is unavailable" << endl;
			error = true;
		}
		if(!dumpAlBase.empty() ||
		   !dumpAlFaBase.empty() ||
		   !dumpAlFqBase.empty())
		{
			cerr << "When -z/--phased is used, the --al option is unavailable" << endl;
			error = true;
		}
		if(error) throw 1;
	}
	if(tryHard) {
		// Increase backtracking limit to huge number
		maxBts = maxBtsBetter = INT_MAX;
		// Increase number of paired-end scan attempts to huge number
		mixedAttemptLim = UINT_MAX;
	}
	if(fullIndex && strata && !stateful && !oldBest) {
		cerr << "--strata must be combined with --best" << endl;
		throw 1;
	}
	if(strata && !allHits && khits == 1 && mhits == 0xffffffff) {
		cerr << "--strata has no effect unless combined with -k, -m or -a" << endl;
		throw 1;
	}
	if(fuzzy && (!stateful && !paired)) {
		cerr << "--fuzzy must be combined with --best or paired-end alignment" << endl;
		throw 1;
	}
	// If both -s and -u are used, we need to adjust qUpto accordingly
	// since it uses patid to know if we've reached the -u limit (and
	// patids are all shifted up by skipReads characters)
	if(qUpto + skipReads > qUpto) {
		qUpto += skipReads;
	}
	if(useShmem && useMm && !quiet) {
		cerr << "Warning: --shmem overrides --mm..." << endl;
		useMm = false;
	}
	if(format == INPUT_CHAIN) {
		bool error = false;
		if(!stateful) {
			cerr << "Error: --chainin must be combined with --best; aborting..." << endl;
			error = true;
		}
		if(paired) {
			cerr << "Error: --chainin cannot be combined with paired-end alignment; aborting..." << endl;
			error = true;
		}
		if(error) throw 1;
	}

	if(outType == OUTPUT_CHAIN) {
		bool error = false;
		if(refOut) {
			cerr << "Error: --chainout is not compatible with --refout; aborting..." << endl;
			error = true;
		}
		if(!stateful) {
			cerr << "Error: --chainout must be combined with --best; aborting..." << endl;
			error = true;
		}
		if(paired) {
			cerr << "Error: --chainout cannot be combined with paired-end alignment; aborting..." << endl;
			error = true;
		}
		if(error) throw 1;
	}
	if(outType == OUTPUT_SAM && refOut) {
		cerr << "Error: --refout cannot be combined with -S/--sam" << endl;
		throw 1;
	}
	if(mate1fw && trim5 > 0) {
		if(verbose) {
			cerr << "Adjusting -I/-X down by " << trim5
				 << " because mate1 is FW & trim5 is " << trim5 << endl;
		}
		maxInsert = max<int>(0, (int)maxInsert - trim5);
		minInsert = max<int>(0, (int)minInsert - trim5);
	}
	if(mate2fw && trim3 > 0) {
		if(verbose) {
			cerr << "Adjusting -I/-X down by " << trim3
				 << " because mate2 is FW & trim3 is " << trim3 << endl;
		}
		maxInsert = max<int>(0, (int)maxInsert - trim3);
		minInsert = max<int>(0, (int)minInsert - trim3);
	}
	if(!mate1fw && trim3 > 0) {
		if(verbose) {
			cerr << "Adjusting -I/-X down by " << trim3
				 << " because mate1 is RC & trim3 is " << trim3 << endl;
		}
		maxInsert = max<int>(0, (int)maxInsert - trim3);
		minInsert = max<int>(0, (int)minInsert - trim3);
	}
	if(!mate2fw && trim5 > 0) {
		if(verbose) {
			cerr << "Adjusting -I/-X down by " << trim5
				 << " because mate2 is RC & trim5 is " << trim5 << endl;
		}
		maxInsert = max<int>(0, (int)maxInsert - trim5);
		minInsert = max<int>(0, (int)minInsert - trim5);
	}
}

static const char *argv0 = NULL;

#define FINISH_READ(p) \
	/* Don't do finishRead if the read isn't legit or if the read was skipped by the doneMask */ \
	if(!p->empty()) { \
		sink->finishRead(*p, true, !skipped); \
	} \
	skipped = false;

static inline void finishReadWithHitmask(PatternSourcePerThread* p,
                                         HitSinkPerThread* sink,
                                         SyncBitset& hitMask,
                                         bool r,
                                         bool& skipped)
{
	/* Don't do finishRead if the read isn't legit */
	if(!p->empty()) {
		/* r = whether to consider reporting the read as unaligned */
		bool reportUnAl = r;
		if(reportUnAl) {
			/* If the done-mask already shows the read as done, */
			/* then we already reported the unaligned read and */
			/* should refrain from re-reporting*/
			reportUnAl = !skipped;
			if(reportUnAl) {
				/* If there hasn't been a hit reported, then report */
				/* read as unaligned */
				reportUnAl = !hitMask.test(p->patid());
			}
		}
		if(sink->finishRead(*p, true, reportUnAl) > 0) {
			/* We reported a hit for the read, so we set the */
			/* appropriate bit in the hitMask to prevent it from */
			/* being reported as unaligned. */
			if(!reportUnAl && sink->dumpsReads()) {
				hitMask.setOver(p->patid());
			}
		}
	}
	skipped = false;
}

/// Macro for getting the next read, possibly aborting depending on
/// whether the result is empty or the patid exceeds the limit, and
/// marshaling the read into convenient variables.
#define GET_READ(p) \
	p->nextReadPair(); \
	if(p->empty() || p->patid() >= qUpto) { \
		p->bufa().clearAll(); \
		break; \
	} \
	assert(!empty(p->bufa().patFw)); \
	String<Dna5>& patFw  = p->bufa().patFw;  \
	patFw.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<Dna5>& patRc  = p->bufa().patRc;  \
	patRc.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& qual = p->bufa().qual; \
	qual.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& qualRev = p->bufa().qualRev; \
	qualRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<Dna5>& patFwRev  = p->bufa().patFwRev;  \
	patFwRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<Dna5>& patRcRev  = p->bufa().patRcRev;  \
	patRcRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& name   = p->bufa().name;   \
	name.data_begin += 0; /* suppress "unused" compiler warning */ \
	uint32_t      patid  = p->patid();       \
	params.setPatId(patid);

/// Macro for getting the forward oriented version of next read,
/// possibly aborting depending on whether the result is empty or the
/// patid exceeds the limit, and marshaling the read into convenient
/// variables.
#define GET_READ_FW(p) \
	p->nextReadPair(); \
	if(p->empty() || p->patid() >= qUpto) { \
		p->bufa().clearAll(); \
		break; \
	} \
	params.setPatId(p->patid()); \
	assert(!empty(p->bufa().patFw)); \
	String<Dna5>& patFw  = p->bufa().patFw;  \
	patFw.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& qual = p->bufa().qual; \
	qual.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<Dna5>& patFwRev  = p->bufa().patFwRev;  \
	patFwRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& qualRev = p->bufa().qualRev; \
	qualRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& name   = p->bufa().name;   \
	name.data_begin += 0; /* suppress "unused" compiler warning */ \
	uint32_t      patid  = p->patid();

#ifdef BOWTIE_PTHREADS
#define WORKER_EXIT() \
	patsrcFact->destroy(patsrc); \
	delete patsrcFact; \
	sinkFact->destroy(sink); \
	delete sinkFact; \
	if(tid > 0) { \
		pthread_exit(NULL); \
	} \
	return NULL;
#else
#define WORKER_EXIT() \
	patsrcFact->destroy(patsrc); \
	delete patsrcFact; \
	sinkFact->destroy(sink); \
	delete sinkFact; \
	return NULL;
#endif

#ifdef CHUD_PROFILING
#define CHUD_START() chudStartRemotePerfMonitor("Bowtie");
#define CHUD_STOP()  chudStopRemotePerfMonitor();
#else
#define CHUD_START()
#define CHUD_STOP()
#endif

/// Create a PatternSourcePerThread for the current thread according
/// to the global params and return a pointer to it
static PatternSourcePerThreadFactory*
createPatsrcFactory(PairedPatternSource& _patsrc, int tid) {
	PatternSourcePerThreadFactory *patsrcFact;
	if(randReadsNoSync) {
		patsrcFact = new RandomPatternSourcePerThreadFactory(numRandomReads, lenRandomReads, nthreads, tid);
	} else {
		patsrcFact = new WrappedPatternSourcePerThreadFactory(_patsrc);
	}
	assert(patsrcFact != NULL);
	return patsrcFact;
}

/**
 * Allocate a HitSinkPerThreadFactory on the heap according to the
 * global params and return a pointer to it.
 */
static HitSinkPerThreadFactory*
createSinkFactory(HitSink& _sink, bool sanity) {
	HitSinkPerThreadFactory *sink = NULL;
	if(format == INPUT_CHAIN) {
		assert(stateful);
		sink = new ChainingHitSinkPerThreadFactory(_sink, allHits ? 0xffffffff : khits, mhits, sanity, strata);
	} else if(!strata) {
		// Unstratified
		if(!allHits) {
			if(oldBest) {
				// First N best, spanning strata
				sink = new NBestHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			} else {
				// First N good; "good" inherently ignores strata
				sink = new NGoodHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			}
		} else {
			// All hits, spanning strata
			sink = new AllHitSinkPerThreadFactory(_sink, mhits, sanity);
		}
	} else {
		// Stratified
		assert(oldBest || stateful);
		if(!allHits) {
			if(oldBest) {
				// First N best, not spanning strata
				sink = new NBestStratHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			} else {
				assert(stateful);
				// Buffer best hits, assuming they're arriving in best-
				// to-worst order
				sink = new NBestFirstStratHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			}
		} else {
			if(oldBest) {
				// All hits, not spanning strata
				sink = new AllStratHitSinkPerThreadFactory(_sink, mhits, sanity);
			} else {
				assert(stateful);
				// Buffer best hits, assuming they're arriving in best-
				// to-worst order
				sink = new NBestFirstStratHitSinkPerThreadFactory(_sink, 0xffffffff/2, mhits, sanity);
			}
		}
	}
	assert(sink != NULL);
	return sink;
}

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static PairedPatternSource*   exactSearch_patsrc;
static HitSink*               exactSearch_sink;
static Ebwt<String<Dna> >*    exactSearch_ebwt;
static vector<String<Dna5> >* exactSearch_os;
static BitPairReference*      exactSearch_refs;
static void *exactSearchWorker(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource& _patsrc = *exactSearch_patsrc;
	HitSink& _sink               = *exactSearch_sink;
	Ebwt<String<Dna> >& ebwt     = *exactSearch_ebwt;
	vector<String<Dna5> >& os    = *exactSearch_os;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	// Per-thread initialization
	PatternSourcePerThreadFactory *patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread *patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);
	HitSinkPerThread* sink = sinkFact->create();
	EbwtSearchParams<String<Dna> > params(
	        *sink,      // HitSink
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        true,       // read is forward
	        true,       // index is forward
	        rangeMode); // range mode
	GreedyDFSRangeSource bt(
	        &ebwt, params,
	        0xffffffff,     // qualThresh
	        0xffffffff,    // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	bool skipped = false;
	while(true) {
		FINISH_READ(patsrc);
		GET_READ(patsrc);
		#include "search_exact.c"
	}
	FINISH_READ(patsrc);
	WORKER_EXIT();
}

/**
 * A statefulness-aware worker driver.  Uses UnpairedExactAlignerV1.
 */
static void *exactSearchWorkerStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource& _patsrc = *exactSearch_patsrc;
	HitSink& _sink               = *exactSearch_sink;
	Ebwt<String<Dna> >& ebwt     = *exactSearch_ebwt;
	vector<String<Dna5> >& os    = *exactSearch_os;
	BitPairReference* refs       =  exactSearch_refs;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);

	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);
	UnpairedExactAlignerV1Factory alSEfact(
			ebwt,
			NULL,
			!nofw,
			!norc,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			os,
			!noMaqRound,
			!better,
			strandFix,
			rangeMode,
			verbose,
			quiet,
			seed);
	PairedExactAlignerV1Factory alPEfact(
			ebwt,
			NULL,
			!nofw,
			!norc,
			useV1,
			_sink,
			*sinkFact,
			mate1fw,
			mate2fw,
			minInsert,
			maxInsert,
			dontReconcileMates,
			mhits,       // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs, os,
			reportSe,
			!noMaqRound,
			strandFix,
			!better,
			rangeMode,
			verbose,
			quiet,
			seed);
	{
		MixedMultiAligner multi(
				prefetchWidth,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) pthread_exit(NULL);
#endif
	return NULL;
}

#define SET_A_FW(bt, p, params) \
	bt.setQuery(&p->bufa().patFw, &p->bufa().qual, &p->bufa().name); \
	params.setFw(true);
#define SET_A_RC(bt, p, params) \
	bt.setQuery(&p->bufa().patRc, &p->bufa().qualRev, &p->bufa().name); \
	params.setFw(false);
#define SET_B_FW(bt, p, params) \
	bt.setQuery(&p->bufb().patFw, &p->bufb().qual, &p->bufb().name); \
	params.setFw(true);
#define SET_B_RC(bt, p, params) \
	bt.setQuery(&p->bufb().patRc, &p->bufb().qualRev, &p->bufb().name); \
	params.setFw(false);

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static void exactSearch(PairedPatternSource& _patsrc,
                        HitSink& _sink,
                        Ebwt<String<Dna> >& ebwt,
                        vector<String<Dna5> >& os,
                        bool paired = false)
{
	exactSearch_patsrc = &_patsrc;
	exactSearch_sink   = &_sink;
	exactSearch_ebwt   = &ebwt;
	exactSearch_os     = &os;

	assert(!ebwt.isInMemory());
	{
		// Load the rest of (vast majority of) the backward Ebwt into
		// memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwt.loadIntoMemory(!noRefNames, startVerbose);
	}

	BitPairReference *refs = NULL;
	if((mates1.size() > 0 || mates12.size() > 0) && mixedThresh < 0xffffffff) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, sanityCheck, NULL, &os, false, useMm, useShmem, mmSweep, verbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	exactSearch_refs   = refs;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	int numAdditionalThreads = nthreads-1;
	pthread_t *threads = new pthread_t[numAdditionalThreads];
	int *tids = new int[numAdditionalThreads];
#endif
	CHUD_START();
	{
		Timer _t(cerr, "Time for 0-mismatch search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < numAdditionalThreads; i++) {
			tids[i] = i+1;
			if(stateful)
				pthread_create(&threads[i], &pt_attr, exactSearchWorkerStateful, (void *)&tids[i]);
			else
				pthread_create(&threads[i], &pt_attr, exactSearchWorker, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		if(stateful) exactSearchWorkerStateful((void*)&tmp);
		else         exactSearchWorker((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < numAdditionalThreads; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
	delete[] tids;
#endif
	if(refs != NULL) delete refs;
}

/**
 * Search through a pair of Ebwt indexes, one for the forward direction
 * and one for the backward direction, for exact end-to-end hits and 1-
 * mismatch end-to-end hits.  In my experience, this is slightly faster
 * than Maq (default) mode with the -n 1 option.
 *
 * Forward Ebwt (ebwtFw) is already loaded into memory and backward
 * Ebwt (ebwtBw) is not loaded into memory.
 */
static PairedPatternSource*           mismatchSearch_patsrc;
static HitSink*                       mismatchSearch_sink;
static Ebwt<String<Dna> >*            mismatchSearch_ebwtFw;
static Ebwt<String<Dna> >*            mismatchSearch_ebwtBw;
static vector<String<Dna5> >*         mismatchSearch_os;
static SyncBitset*                    mismatchSearch_doneMask;
static SyncBitset*                    mismatchSearch_hitMask;
static BitPairReference*              mismatchSearch_refs;

static void* mismatchSearchWorkerPhase1(void *vp){
	int tid = *((int*)vp);
	PairedPatternSource&   _patsrc       = *mismatchSearch_patsrc;
	HitSink&               _sink         = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw        = *mismatchSearch_ebwtFw;
	vector<String<Dna5> >& os            = *mismatchSearch_os;
	SyncBitset&            doneMask      = *mismatchSearch_doneMask;
	SyncBitset&            hitMask       = *mismatchSearch_hitMask;
    bool sanity = sanityCheck && !os.empty() && !rangeMode;
    PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread* patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);
	HitSinkPerThread* sink = sinkFact->create();
	EbwtSearchParams<String<Dna> > params(
			*sink,      // HitSinkPerThread
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        false,      // read is forward
	        true,       // index is forward
	        rangeMode); // range mode
	GreedyDFSRangeSource bt(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        0xffffffff,    // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
		GET_READ(patsrc);
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_1mm_phase1.c"
		#undef DONEMASK_SET
	} // End read loop
	finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
    WORKER_EXIT();
}

static void* mismatchSearchWorkerPhase2(void *vp){
	int tid = *((int*)vp);
	PairedPatternSource&   _patsrc      = *mismatchSearch_patsrc;
	HitSink&               _sink        = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtBw       = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os           = *mismatchSearch_os;
	SyncBitset&            doneMask     = *mismatchSearch_doneMask;
	SyncBitset&            hitMask      = *mismatchSearch_hitMask;
    // Per-thread initialization
    bool sanity = sanityCheck && !os.empty() && !rangeMode;
    PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread* patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);
	HitSinkPerThread* sink = sinkFact->create();
	EbwtSearchParams<String<Dna> > params(
			*sink,      // HitSinkPerThread
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        true,       // read is forward
	        false,      // index is mirror index
	        rangeMode); // range mode
	GreedyDFSRangeSource bt(
			&ebwtBw, params,
	        0xffffffff,     // qualThresh
	        0xffffffff,    // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, true, skipped);
		GET_READ(patsrc);
		if(doneMask.test(patid)) { skipped = true; continue; }
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		#include "search_1mm_phase2.c"
	} // End read loop
	finishReadWithHitmask(patsrc, sink, hitMask, true, skipped);
    WORKER_EXIT();
}

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static void mismatchSearch(PairedPatternSource& _patsrc,
                           HitSink& _sink,
                           Ebwt<String<Dna> >& ebwtFw,
                           Ebwt<String<Dna> >& ebwtBw,
                           vector<String<Dna5> >& os)
{
	uint32_t numQs = ((qUpto == 0xffffffff) ? 16 * 1024 * 1024 : qUpto);
	if(fullIndex) numQs = 0;
	SyncBitset doneMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the read mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");
	// No need to keep track of which reads are aligned because the
	// user hasn't requested an unaligned-read dump
	if(!_sink.dumpsReads()) numQs = 0;
	SyncBitset hitMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the hit mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");

	mismatchSearch_patsrc       = &_patsrc;
	mismatchSearch_sink         = &_sink;
	mismatchSearch_ebwtFw       = &ebwtFw;
	mismatchSearch_ebwtBw       = &ebwtBw;
	mismatchSearch_doneMask     = &doneMask;
	mismatchSearch_hitMask      = &hitMask;
	mismatchSearch_os           = &os;

	assert(!ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());

	{
		// Load the rest of (vast majority of) the backward Ebwt into
		// memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(!noRefNames, startVerbose);
	}

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
	int *tids = new int[nthreads-1];
#endif
    CHUD_START();
	// Phase 1
    {
		Timer _t(cerr, "Time for 1-mismatch Phase 1 of 2: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerPhase1, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		mismatchSearchWorkerPhase1((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
    }

	// Release most of the memory associated with the forward Ebwt
    ebwtFw.evictFromMemory();
	{
		// Load the rest of (vast majority of) the backward Ebwt into
		// memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(!noRefNames, startVerbose);
	}
    _patsrc.reset();          // reset pattern source to 1st pattern
	// Sanity-check the restored version of the Ebwt
	if(sanityCheck && !os.empty()) {
		ebwtBw.checkOrigs(os, true);
	}

	// Phase 2
	{
		Timer _t(cerr, "Time for 1-mismatch Phase 2 of 2: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerPhase2, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		mismatchSearchWorkerPhase2((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
	delete[] tids;
#endif
}

/**
 * A statefulness-aware worker driver.  Uses Unpaired/Paired1mmAlignerV1.
 */
static void *mismatchSearchWorkerFullStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource&   _patsrc = *mismatchSearch_patsrc;
	HitSink&               _sink   = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw  = *mismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw  = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os      = *mismatchSearch_os;
	BitPairReference*      refs    =  mismatchSearch_refs;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
    PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);
	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);

	Unpaired1mmAlignerV1Factory alSEfact(
			ebwtFw,
			&ebwtBw,
			!nofw,
			!norc,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			os,
			!noMaqRound,
			!better,
			strandFix,
			rangeMode,
			verbose,
			quiet,
			seed);
	Paired1mmAlignerV1Factory alPEfact(
			ebwtFw,
			&ebwtBw,
			!nofw,
			!norc,
			useV1,
			_sink,
			*sinkFact,
			mate1fw,
			mate2fw,
			minInsert,
			maxInsert,
			dontReconcileMates,
			mhits,     // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs, os,
			reportSe,
			!noMaqRound,
			!better,
			strandFix,
			rangeMode,
			verbose,
			quiet,
			seed);
	{
		MixedMultiAligner multi(
				prefetchWidth,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) pthread_exit(NULL);
#endif
	return NULL;
}

static void* mismatchSearchWorkerFull(void *vp){
	int tid = *((int*)vp);
	PairedPatternSource&   _patsrc      = *mismatchSearch_patsrc;
	HitSink&               _sink        = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw       = *mismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw       = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os           = *mismatchSearch_os;
    // Per-thread initialization
    bool sanity = sanityCheck && !os.empty() && !rangeMode;
    PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread* patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);
	HitSinkPerThread* sink = sinkFact->create();
	EbwtSearchParams<String<Dna> > params(
			*sink,      // HitSinkPerThread
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        true,       // read is forward
	        false,      // index is mirror index
	        rangeMode); // range mode
	GreedyDFSRangeSource bt(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        0xffffffff,     // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	bool skipped = false;
	while(true) {
		FINISH_READ(patsrc);
		GET_READ(patsrc);
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		#define DONEMASK_SET(p)
		#include "search_1mm_phase1.c"
		#include "search_1mm_phase2.c"
		#undef DONEMASK_SET
	} // End read loop
	FINISH_READ(patsrc);
    WORKER_EXIT();
}

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static void mismatchSearchFull(PairedPatternSource& _patsrc,
                               HitSink& _sink,
                               Ebwt<String<Dna> >& ebwtFw,
                               Ebwt<String<Dna> >& ebwtBw,
                               vector<String<Dna5> >& os)
{
	mismatchSearch_patsrc       = &_patsrc;
	mismatchSearch_sink         = &_sink;
	mismatchSearch_ebwtFw       = &ebwtFw;
	mismatchSearch_ebwtBw       = &ebwtBw;
	mismatchSearch_doneMask     = NULL;
	mismatchSearch_hitMask      = NULL;
	mismatchSearch_os           = &os;

	assert(!ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(!noRefNames, startVerbose);
	}
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(!noRefNames, startVerbose);
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	if((mates1.size() > 0 || mates12.size() > 0) && mixedThresh < 0xffffffff) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, sanityCheck, NULL, &os, false, useMm, useShmem, mmSweep, verbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	mismatchSearch_refs = refs;

#ifdef BOWTIE_PTHREADS
	// Allocate structures for threads
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
	int *tids = new int[nthreads-1];
#endif
    CHUD_START();
    {
		Timer _t(cerr, "Time for 1-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			if(stateful)
				pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerFullStateful, (void *)&tids[i]);
			else
				pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerFull, (void *)&tids[i]);
		}
#endif
		// Go to town
		int tmp = 0;
		if(stateful) mismatchSearchWorkerFullStateful((void*)&tmp);
		else         mismatchSearchWorkerFull((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
    }
#ifdef BOWTIE_PTHREADS
	delete[] threads;
	delete[] tids;
#endif
	if(refs != NULL) delete refs;
}

#define SWITCH_TO_FW_INDEX() { \
	/* Evict the mirror index from memory if necessary */ \
	if(ebwtBw.isInMemory()) ebwtBw.evictFromMemory(); \
	assert(!ebwtBw.isInMemory()); \
	/* Load the forward index into memory if necessary */ \
	if(!ebwtFw.isInMemory()) { \
		Timer _t(cerr, "Time loading forward index: ", timing); \
		ebwtFw.loadIntoMemory(!noRefNames, startVerbose); \
	} \
	assert(ebwtFw.isInMemory()); \
	_patsrc.reset(); /* rewind pattern source to first pattern */ \
}

#define SWITCH_TO_BW_INDEX() { \
	/* Evict the forward index from memory if necessary */ \
	if(ebwtFw.isInMemory()) ebwtFw.evictFromMemory(); \
	assert(!ebwtFw.isInMemory()); \
	/* Load the forward index into memory if necessary */ \
	if(!ebwtBw.isInMemory()) { \
		Timer _t(cerr, "Time loading mirror index: ", timing); \
		ebwtBw.loadIntoMemory(!noRefNames, startVerbose); \
	} \
	assert(ebwtBw.isInMemory()); \
	_patsrc.reset(); /* rewind pattern source to first pattern */ \
}

#define ASSERT_NO_HITS_FW(ebwtfw) \
	if(sanityCheck && os.size() > 0) { \
		vector<Hit> hits; \
		uint32_t threeRevOff = (seedMms <= 3) ? s : 0; \
		uint32_t twoRevOff   = (seedMms <= 2) ? s : 0; \
		uint32_t oneRevOff   = (seedMms <= 1) ? s : 0; \
		uint32_t unrevOff    = (seedMms == 0) ? s : 0; \
		::naiveOracle( \
		        os, \
				patFw, \
				plen, \
		        qual, \
		        name, \
		        patid, \
		        hits, \
		        qualCutoff, \
		        unrevOff, \
		        oneRevOff, \
		        twoRevOff, \
		        threeRevOff, \
		        true,        /* fw */ \
		        ebwtfw,      /* ebwtFw */ \
		        0,           /* iham */ \
		        NULL,        /* muts */ \
		        !noMaqRound, /* maqRound */ \
		        false,       /* halfAndHalf */ \
		        true,        /* reportExacts */ \
		        ebwtfw);     /* invert */ \
		if(hits.size() > 0) { \
			/* Print offending hit obtained by oracle */ \
			::printHit( \
				os, \
				hits[0], \
				patFw, \
				plen, \
			    unrevOff, \
			    oneRevOff, \
			    twoRevOff, \
			    threeRevOff, \
			    ebwtfw);  /* ebwtFw */ \
		} \
		assert_eq(0, hits.size()); \
	}

#define ASSERT_NO_HITS_RC(ebwtfw) \
	if(sanityCheck && os.size() > 0) { \
		vector<Hit> hits; \
		uint32_t threeRevOff = (seedMms <= 3) ? s : 0; \
		uint32_t twoRevOff   = (seedMms <= 2) ? s : 0; \
		uint32_t oneRevOff   = (seedMms <= 1) ? s : 0; \
		uint32_t unrevOff    = (seedMms == 0) ? s : 0; \
		::naiveOracle( \
		        os, \
				patRc, \
				plen, \
		        qualRev, \
		        name, \
		        patid, \
		        hits, \
		        qualCutoff, \
		        unrevOff, \
		        oneRevOff, \
		        twoRevOff, \
		        threeRevOff, \
		        false,       /* fw */ \
		        ebwtfw,      /* ebwtFw */ \
		        0,           /* iham */ \
		        NULL,        /* muts */ \
		        !noMaqRound, /* maqRound */ \
		        false,       /* halfAndHalf */ \
		        true,        /* reportExacts */ \
		        !ebwtfw);    /* invert */ \
		if(hits.size() > 0) { \
			/* Print offending hit obtained by oracle */ \
			::printHit( \
				os, \
				hits[0], \
				patRc, \
				plen, \
			    unrevOff, \
			    oneRevOff, \
			    twoRevOff, \
			    threeRevOff, \
			    ebwtfw);  /* ebwtFw */ \
		} \
		assert_eq(0, hits.size()); \
	}

static PairedPatternSource*           twoOrThreeMismatchSearch_patsrc;
static HitSink*                       twoOrThreeMismatchSearch_sink;
static Ebwt<String<Dna> >*            twoOrThreeMismatchSearch_ebwtFw;
static Ebwt<String<Dna> >*            twoOrThreeMismatchSearch_ebwtBw;
static vector<String<Dna5> >*         twoOrThreeMismatchSearch_os;
static SyncBitset*                    twoOrThreeMismatchSearch_doneMask;
static SyncBitset*                    twoOrThreeMismatchSearch_hitMask;
static bool                           twoOrThreeMismatchSearch_two;
static BitPairReference*              twoOrThreeMismatchSearch_refs;

#define TWOTHREE_WORKER_SETUP() \
	int tid = *((int*)vp); \
	PairedPatternSource&           _patsrc  = *twoOrThreeMismatchSearch_patsrc;   \
	HitSink&                       _sink    = *twoOrThreeMismatchSearch_sink;     \
	vector<String<Dna5> >&         os       = *twoOrThreeMismatchSearch_os;       \
	bool                           two      = twoOrThreeMismatchSearch_two; \
    PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid); \
	PatternSourcePerThread* patsrc = patsrcFact->create(); \
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, false); \
	HitSinkPerThread* sink = sinkFact->create(); \
	/* Per-thread initialization */ \
	EbwtSearchParams<String<Dna> > params( \
			*sink,       /* HitSink */ \
	        os,          /* reference sequences */ \
	        revcomp,     /* forward AND reverse complement? */ \
	        true,        /* read is forward */ \
	        true,        /* index is forward */ \
	        rangeMode);  /* range mode (irrelevant here) */

static void* twoOrThreeMismatchSearchWorkerPhase1(void *vp) {
	TWOTHREE_WORKER_SETUP();
	SyncBitset& doneMask = *twoOrThreeMismatchSearch_doneMask;
	SyncBitset& hitMask  = *twoOrThreeMismatchSearch_hitMask;
	Ebwt<String<Dna> >& ebwtFw = *twoOrThreeMismatchSearch_ebwtFw;
	GreedyDFSRangeSource btr1(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	bool skipped = false;
    while(true) { // Read read-in loop
    	finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
		GET_READ(patsrc);
		// If requested, check that this read has the same length
		// as all the previous ones
		size_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_23mm_phase1.c"
		#undef DONEMASK_SET
    }
	finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
    // Threads join at end of Phase 1
	WORKER_EXIT();
}

static void* twoOrThreeMismatchSearchWorkerPhase2(void *vp) {
	TWOTHREE_WORKER_SETUP();
	SyncBitset& doneMask = *twoOrThreeMismatchSearch_doneMask;
	SyncBitset& hitMask = *twoOrThreeMismatchSearch_hitMask;
	Ebwt<String<Dna> >& ebwtBw = *twoOrThreeMismatchSearch_ebwtBw;
	GreedyDFSRangeSource bt2(
	        &ebwtBw, params,
	        0xffffffff,     // qualThresh
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (no)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
		GET_READ(patsrc);
		if(doneMask.test(patid)) continue;
		size_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_23mm_phase2.c"
		#undef DONEMASK_SET
	}
	finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
	WORKER_EXIT();
}

static void* twoOrThreeMismatchSearchWorkerPhase3(void *vp) {
	TWOTHREE_WORKER_SETUP();
	SyncBitset& doneMask = *twoOrThreeMismatchSearch_doneMask;
	SyncBitset& hitMask = *twoOrThreeMismatchSearch_hitMask;
	Ebwt<String<Dna> >& ebwtFw   = *twoOrThreeMismatchSearch_ebwtFw;
	// GreedyDFSRangeSource to search for seedlings for case 4F
	GreedyDFSRangeSource bt3(
	        &ebwtFw, params,
	        0xffffffff,     // qualThresh (none)
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	GreedyDFSRangeSource bthh3(
	        &ebwtFw, params,
	        0xffffffff,     // qualThresh
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false,          // considerQuals
	        true);          // halfAndHalf
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, true, skipped);
		GET_READ(patsrc);
		if(doneMask.testUnsync(patid)) { skipped = true; continue; }
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_23mm_phase3.c"
		#undef DONEMASK_SET
	}
	finishReadWithHitmask(patsrc, sink, hitMask, true, skipped);
	WORKER_EXIT();
}

template<typename TStr>
static void twoOrThreeMismatchSearch(
        PairedPatternSource& _patsrc,   /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of mirror text
        vector<String<Dna5> >& os,      /// text strings, if available (empty otherwise)
        bool two = true)                /// true -> 2, false -> 3
{
	// Global initialization
	assert(revcomp);
	assert(!fullIndex);
	assert(!ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());

	uint32_t numQs = ((qUpto == 0xffffffff) ? 16 * 1024 * 1024 : qUpto);
	SyncBitset doneMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the read mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");
	// No need to keep track of which reads are aligned because the
	// user hasn't requested an unaligned-read dump
	if(!_sink.dumpsReads()) numQs = 0;
	SyncBitset hitMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the hit mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");

	twoOrThreeMismatchSearch_patsrc   = &_patsrc;
	twoOrThreeMismatchSearch_sink     = &_sink;
	twoOrThreeMismatchSearch_ebwtFw   = &ebwtFw;
	twoOrThreeMismatchSearch_ebwtBw   = &ebwtBw;
	twoOrThreeMismatchSearch_os       = &os;
	twoOrThreeMismatchSearch_doneMask = &doneMask;
	twoOrThreeMismatchSearch_hitMask  = &hitMask;
	twoOrThreeMismatchSearch_two      = two;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
	int *tids = new int[nthreads-1];
#endif

	// Load forward index
	SWITCH_TO_FW_INDEX();
	{ // Phase 1
		Timer _t(cerr, "End-to-end 2/3-mismatch Phase 1 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerPhase1, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		twoOrThreeMismatchSearchWorkerPhase1((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
    }
	// Unload forward index and load mirror index
	SWITCH_TO_BW_INDEX();
	{ // Phase 2
		Timer _t(cerr, "End-to-end 2/3-mismatch Phase 2 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerPhase2, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		twoOrThreeMismatchSearchWorkerPhase2((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
	SWITCH_TO_FW_INDEX();
	{ // Phase 3
		Timer _t(cerr, "End-to-end 2/3-mismatch Phase 3 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerPhase3, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		twoOrThreeMismatchSearchWorkerPhase3((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
	delete[] tids;
#endif
	return;
}

/**
 * A statefulness-aware worker driver.  Uses UnpairedExactAlignerV1.
 */
static void *twoOrThreeMismatchSearchWorkerStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource&   _patsrc = *twoOrThreeMismatchSearch_patsrc;
	HitSink&               _sink   = *twoOrThreeMismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw  = *twoOrThreeMismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw  = *twoOrThreeMismatchSearch_ebwtBw;
	vector<String<Dna5> >& os      = *twoOrThreeMismatchSearch_os;
	BitPairReference*      refs    =  twoOrThreeMismatchSearch_refs;
	static bool            two     =  twoOrThreeMismatchSearch_two;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);

	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);
	Unpaired23mmAlignerV1Factory alSEfact(
			ebwtFw,
			&ebwtBw,
			two,
			!nofw,
			!norc,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			os,
			!noMaqRound,
			!better,
			strandFix,
			rangeMode,
			verbose,
			quiet,
			seed);
	Paired23mmAlignerV1Factory alPEfact(
			ebwtFw,
			&ebwtBw,
			!nofw,
			!norc,
			useV1,
			two,
			_sink,
			*sinkFact,
			mate1fw,
			mate2fw,
			minInsert,
			maxInsert,
			dontReconcileMates,
			mhits,       // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs, os,
			reportSe,
			!noMaqRound,
			!better,
			strandFix,
			rangeMode,
			verbose,
			quiet,
			seed);
	{
		MixedMultiAligner multi(
				prefetchWidth,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) pthread_exit(NULL);
#endif
	return NULL;
}

static void* twoOrThreeMismatchSearchWorkerFull(void *vp) {
	TWOTHREE_WORKER_SETUP();
	Ebwt<String<Dna> >& ebwtFw = *twoOrThreeMismatchSearch_ebwtFw;
	Ebwt<String<Dna> >& ebwtBw = *twoOrThreeMismatchSearch_ebwtBw;
	GreedyDFSRangeSource btr1(
	        &ebwtFw, params,
	        0xffffffff,     // qualThresh
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	GreedyDFSRangeSource bt2(
	        &ebwtBw, params,
	        0xffffffff,     // qualThresh
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (no)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	GreedyDFSRangeSource bt3(
	        &ebwtFw, params,
	        0xffffffff,     // qualThresh (none)
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false);         // considerQuals
	GreedyDFSRangeSource bthh3(
	        &ebwtFw, params,
	        0xffffffff,     // qualThresh
	        // Do not impose maximums in 2/3-mismatch mode
	        0xffffffff,     // max backtracks (no limit)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        &os,
	        false,          // considerQuals
	        true);          // halfAndHalf
	bool skipped = false;
	while(true) { // Read read-in loop
		FINISH_READ(patsrc);
		GET_READ(patsrc);
		patid += 0; // kill unused variable warning
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		#define DONEMASK_SET(p)
		#include "search_23mm_phase1.c"
		#include "search_23mm_phase2.c"
		#include "search_23mm_phase3.c"
		#undef DONEMASK_SET
	}
	FINISH_READ(patsrc);
	// Threads join at end of Phase 1
	WORKER_EXIT();
}

template<typename TStr>
static void twoOrThreeMismatchSearchFull(
		PairedPatternSource& _patsrc,   /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of mirror text
        vector<String<Dna5> >& os,      /// text strings, if available (empty otherwise)
        bool two = true)                /// true -> 2, false -> 3
{
	// Global initialization
	assert(revcomp);
	assert(fullIndex);
	assert(!ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(!noRefNames, startVerbose);
	}
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(!noRefNames, startVerbose);
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	if((mates1.size() > 0 || mates12.size() > 0) && mixedThresh < 0xffffffff) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, sanityCheck, NULL, &os, false, useMm, useShmem, mmSweep, verbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	twoOrThreeMismatchSearch_refs     = refs;
	twoOrThreeMismatchSearch_patsrc   = &_patsrc;
	twoOrThreeMismatchSearch_sink     = &_sink;
	twoOrThreeMismatchSearch_ebwtFw   = &ebwtFw;
	twoOrThreeMismatchSearch_ebwtBw   = &ebwtBw;
	twoOrThreeMismatchSearch_os       = &os;
	twoOrThreeMismatchSearch_doneMask = NULL;
	twoOrThreeMismatchSearch_hitMask  = NULL;
	twoOrThreeMismatchSearch_two      = two;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
	int *tids = new int[nthreads-1];
#endif
    CHUD_START();
    {
		Timer _t(cerr, "End-to-end 2/3-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			if(stateful)
				pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerStateful, (void *)&tids[i]);
			else
				pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerFull, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		if(stateful) twoOrThreeMismatchSearchWorkerStateful((void*)&tmp);
		else         twoOrThreeMismatchSearchWorkerFull((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
    }
#ifdef BOWTIE_PTHREADS
	delete[] threads;
	delete[] tids;
#endif
	if(refs != NULL) delete refs;
	return;
}

static PairedPatternSource*     seededQualSearch_patsrc;
static HitSink*                 seededQualSearch_sink;
static Ebwt<String<Dna> >*      seededQualSearch_ebwtFw;
static Ebwt<String<Dna> >*      seededQualSearch_ebwtBw;
static vector<String<Dna5> >*   seededQualSearch_os;
static SyncBitset*              seededQualSearch_doneMask;
static SyncBitset*              seededQualSearch_hitMask;
static PartialAlignmentManager* seededQualSearch_pamFw;
static PartialAlignmentManager* seededQualSearch_pamRc;
static int                      seededQualSearch_qualCutoff;
static BitPairReference*        seededQualSearch_refs;

#define SEEDEDQUAL_WORKER_SETUP() \
	int tid = *((int*)vp); \
	PairedPatternSource&     _patsrc    = *seededQualSearch_patsrc;    \
	HitSink&                 _sink      = *seededQualSearch_sink;      \
	vector<String<Dna5> >&   os         = *seededQualSearch_os;        \
	int                      qualCutoff = seededQualSearch_qualCutoff; \
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid); \
	PatternSourcePerThread* patsrc = patsrcFact->create(); \
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, false); \
	HitSinkPerThread* sink = sinkFact->create(); \
	/* Per-thread initialization */ \
	EbwtSearchParams<String<Dna> > params( \
	        *sink,       /* HitSink */ \
	        os,          /* reference sequences */ \
	        revcomp,     /* forward AND reverse complement? */ \
	        true,        /* read is forward */ \
	        true,        /* index is forward */ \
	        rangeMode);  /* range mode (irrelevant here) */

static void* seededQualSearchWorkerPhase1(void *vp) {
	SEEDEDQUAL_WORKER_SETUP();
	SyncBitset& doneMask = *seededQualSearch_doneMask;
	SyncBitset& hitMask = *seededQualSearch_hitMask;
	Ebwt<String<Dna> >& ebwtFw = *seededQualSearch_ebwtFw;
	uint32_t s = seedLen;
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	// GreedyDFSRangeSource for finding exact hits for the forward-
	// oriented read
	GreedyDFSRangeSource btf1(
	        &ebwtFw, params,
	        qualCutoff,            // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partials
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,
	        false);                // considerQuals
	GreedyDFSRangeSource bt1(
	        &ebwtFw, params,
	        qualCutoff,            // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partials
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,
	        true,                  // considerQuals
	        false, !noMaqRound);
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
		GET_READ(patsrc);
		size_t plen = length(patFw);
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_seeded_phase1.c"
		#undef DONEMASK_SET
	}
	finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
	WORKER_EXIT();
}

static void* seededQualSearchWorkerPhase2(void *vp) {
	SEEDEDQUAL_WORKER_SETUP();
	SyncBitset& doneMask = *seededQualSearch_doneMask;
	SyncBitset& hitMask  = *seededQualSearch_hitMask;
	uint32_t s = seedLen;
	uint32_t s3 = s >> 1; /* length of 3' half of seed */
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	Ebwt<String<Dna> >& ebwtBw = *seededQualSearch_ebwtBw;
	PartialAlignmentManager* pamRc = seededQualSearch_pamRc;
	// GreedyDFSRangeSource to search for hits for cases 1F, 2F, 3F
	GreedyDFSRangeSource btf2(
	        &ebwtBw, params,
	        qualCutoff,            // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,                     // reportPartials (no)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partial alignment manager
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for partial alignments for case 4R
	GreedyDFSRangeSource btr2(
	        &ebwtBw, params,
	        qualCutoff,            // qualThresh (none)
	        maxBtsBetter,          // max backtracks
	        seedMms,               // report partials (up to seedMms mms)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamRc,                 // partial alignment manager
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, seedMms == 0, skipped);
		GET_READ(patsrc);
		if(doneMask.test(patid)) { skipped = true; continue; }
		size_t plen = length(patFw);
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs3 = (qs >> 1);
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_seeded_phase2.c"
		#undef DONEMASK_SET
	}
	finishReadWithHitmask(patsrc, sink, hitMask, seedMms == 0, skipped);
	WORKER_EXIT();
}

static void* seededQualSearchWorkerPhase3(void *vp) {
	SEEDEDQUAL_WORKER_SETUP();
	SyncBitset& doneMask = *seededQualSearch_doneMask;
	SyncBitset& hitMask  = *seededQualSearch_hitMask;
	uint32_t s = seedLen;
	uint32_t s3 = s >> 1; /* length of 3' half of seed */
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	Ebwt<String<Dna> >& ebwtFw        = *seededQualSearch_ebwtFw;
	PartialAlignmentManager* pamFw    = seededQualSearch_pamFw;
	PartialAlignmentManager* pamRc    = seededQualSearch_pamRc;
	// GreedyDFSRangeSource to search for seedlings for case 4F
	GreedyDFSRangeSource btf3(
	        &ebwtFw, params,
	        qualCutoff,            // qualThresh (none)
	        maxBtsBetter,          // max backtracks
	        seedMms,               // reportPartials (do)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamFw,                 // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for hits for case 4R by extending
	// the partial alignments found in Phase 2
	GreedyDFSRangeSource btr3(
	        &ebwtFw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// The half-and-half GreedyDFSRangeSource
	GreedyDFSRangeSource btr23(
	        &ebwtFw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,
	        true,    // considerQuals
	        true,    // halfAndHalf
	        !noMaqRound);
	vector<PartialAlignment> pals;
	String<QueryMutation> muts;
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
		GET_READ(patsrc);
		size_t plen = length(patFw);
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs3 = (qs >> 1);
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		if(doneMask.test(patid)) continue;
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_seeded_phase3.c"
		#undef DONEMASK_SET
	}
	finishReadWithHitmask(patsrc, sink, hitMask, false, skipped);
	WORKER_EXIT();
}

static void* seededQualSearchWorkerPhase4(void *vp) {
	SEEDEDQUAL_WORKER_SETUP();
	SyncBitset& doneMask = *seededQualSearch_doneMask;
	SyncBitset& hitMask  = *seededQualSearch_hitMask;
	uint32_t s = seedLen;
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	Ebwt<String<Dna> >& ebwtBw = *seededQualSearch_ebwtBw;
	PartialAlignmentManager* pamFw = seededQualSearch_pamFw;
	// GreedyDFSRangeSource to search for hits for case 4F by extending
	// the partial alignments found in Phase 3
	GreedyDFSRangeSource btf4(
	        &ebwtBw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter, // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// Half-and-half GreedyDFSRangeSource for forward read
	GreedyDFSRangeSource btf24(
	        &ebwtBw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter, // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,
	        true,    // considerQuals
	        true,    // halfAndHalf
	        !noMaqRound);
	vector<PartialAlignment> pals;
	String<QueryMutation> muts;
	bool skipped = false;
	while(true) {
		finishReadWithHitmask(patsrc, sink, hitMask, true, skipped);
		GET_READ_FW(patsrc);
		if(doneMask.testUnsync(patid)) { skipped = true; continue; }
		size_t plen = length(patFw);
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_seeded_phase4.c"
		#undef DONEMASK_SET
	}
	finishReadWithHitmask(patsrc, sink, hitMask, true, skipped);
	WORKER_EXIT();
}

static void* seededQualSearchWorkerFull(void *vp) {
	SEEDEDQUAL_WORKER_SETUP();
	Ebwt<String<Dna> >& ebwtFw = *seededQualSearch_ebwtFw;
	Ebwt<String<Dna> >& ebwtBw = *seededQualSearch_ebwtBw;
	PartialAlignmentManager * pamRc = NULL;
	PartialAlignmentManager * pamFw = NULL;
	if(seedMms > 0) {
		pamRc = new PartialAlignmentManager(64);
		pamFw = new PartialAlignmentManager(64);
	}
	vector<PartialAlignment> pals;
	// GreedyDFSRangeSource for finding exact hits for the forward-
	// oriented read
	GreedyDFSRangeSource btf1(
	        &ebwtFw, params,
	        qualCutoff,            // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,
	        false);                // considerQuals
	GreedyDFSRangeSource bt1(
	        &ebwtFw, params,
	        qualCutoff,            // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for hits for cases 1F, 2F, 3F
	GreedyDFSRangeSource btf2(
	        &ebwtBw, params,
	        qualCutoff,            // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,                     // reportPartials (no)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partial alignment manager
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for partial alignments for case 4R
	GreedyDFSRangeSource btr2(
	        &ebwtBw, params,
	        qualCutoff,            // qualThresh (none)
	        maxBtsBetter,          // max backtracks
	        seedMms,               // report partials (up to seedMms mms)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamRc,                 // partial alignment manager
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for seedlings for case 4F
	GreedyDFSRangeSource btf3(
	        &ebwtFw, params,
	        qualCutoff,            // qualThresh (none)
	        maxBtsBetter,          // max backtracks
	        seedMms,               // reportPartials (do)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamFw,                 // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for hits for case 4R by extending
	// the partial alignments found in Phase 2
	GreedyDFSRangeSource btr3(
	        &ebwtFw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// The half-and-half GreedyDFSRangeSource
	GreedyDFSRangeSource btr23(
	        &ebwtFw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,
	        true,    // considerQuals
	        true,    // halfAndHalf
	        !noMaqRound);
	// GreedyDFSRangeSource to search for hits for case 4F by extending
	// the partial alignments found in Phase 3
	GreedyDFSRangeSource btf4(
	        &ebwtBw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// Half-and-half GreedyDFSRangeSource for forward read
	GreedyDFSRangeSource btf24(
	        &ebwtBw, params,
	        qualCutoff, // qualThresh
	        maxBtsBetter,          // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
	        NULL,    // mutations
	        verbose, // verbose
	        &os,
	        true,    // considerQuals
	        true,    // halfAndHalf
	        !noMaqRound);
	String<QueryMutation> muts;
	bool skipped = false;
	while(true) {
		FINISH_READ(patsrc);
		GET_READ(patsrc);
		size_t plen = length(patFw);
		uint32_t s = seedLen;
		uint32_t s3 = (s >> 1); /* length of 3' half of seed */
		uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs3 = qs >> 1;
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		#define DONEMASK_SET(p)
		#include "search_seeded_phase1.c"
		#include "search_seeded_phase2.c"
		#include "search_seeded_phase3.c"
		#include "search_seeded_phase4.c"
		#undef DONEMASK_SET
	}
	FINISH_READ(patsrc);
	if(seedMms > 0) {
		delete pamRc;
		delete pamFw;
	}
	WORKER_EXIT();
}

static void* seededQualSearchWorkerFullStateful(void *vp) {
	int tid = *((int*)vp);
	PairedPatternSource&     _patsrc    = *seededQualSearch_patsrc;
	HitSink&                 _sink      = *seededQualSearch_sink;
	Ebwt<String<Dna> >&      ebwtFw     = *seededQualSearch_ebwtFw;
	Ebwt<String<Dna> >&      ebwtBw     = *seededQualSearch_ebwtBw;
	vector<String<Dna5> >&   os         = *seededQualSearch_os;
	int                      qualCutoff = seededQualSearch_qualCutoff;
	BitPairReference*        refs       = seededQualSearch_refs;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);
	ChunkPool *pool = new ChunkPool(chunkSz * 1024, chunkPoolMegabytes * 1024 * 1024, chunkVerbose);

	AlignerMetrics *metrics = NULL;
	if(stats) {
		metrics = new AlignerMetrics();
	}
	UnpairedSeedAlignerFactory alSEfact(
			ebwtFw,
			&ebwtBw,
			!nofw,
			!norc,
			seedMms,
			seedLen,
			qualCutoff,
			maxBts,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			os,
			!noMaqRound,
			!better,
			strandFix,
			rangeMode,
			verbose,
			quiet,
			seed,
			metrics);
	PairedSeedAlignerFactory alPEfact(
			ebwtFw,
			&ebwtBw,
			useV1,
			!nofw,
			!norc,
			seedMms,
			seedLen,
			qualCutoff,
			maxBts,
			_sink,
			*sinkFact,
			mate1fw,
			mate2fw,
			minInsert,
			maxInsert,
			dontReconcileMates,
			mhits,       // for symCeiling
			mixedThresh,
			mixedAttemptLim,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			pool,
			refs,
			os,
			reportSe,
			!noMaqRound,
			!better,
			strandFix,
			rangeMode,
			verbose,
			quiet,
			seed);
	{
		MixedMultiAligner multi(
				prefetchWidth,
				qUpto,
				alSEfact,
				alPEfact,
				*patsrcFact);
		// Run that mother
		multi.run();
		// MultiAligner must be destroyed before patsrcFact
	}
	if(metrics != NULL) {
		metrics->printSummary();
		delete metrics;
	}

	delete patsrcFact;
	delete sinkFact;
	delete pool;
#ifdef BOWTIE_PTHREADS
	if(tid > 0) {
		pthread_exit(NULL);
	}
#endif
	return NULL;
}

/**
 * Search for a good alignments for each read using criteria that
 * correspond somewhat faithfully to Maq's.  Search is aided by a pair
 * of Ebwt indexes, one for the original references, and one for the
 * transpose of the references.  Neither index should be loaded upon
 * entry to this function.
 *
 * Like Maq, we treat the first 24 base pairs of the read (those
 * closest to the 5' end) differently from the remainder of the read.
 * We call the first 24 base pairs the "seed."
 */
template<typename TStr>
static void seededQualCutoffSearch(
        int seedLen,                  /// length of seed (not a maq option)
        int qualCutoff,               /// maximum sum of mismatch qualities
                                      /// like maq map's -e option
                                      /// default: 70
        int seedMms,                  /// max # mismatches allowed in seed
                                      /// (like maq map's -n option)
                                      /// Can only be 1 or 2, default: 1
        PairedPatternSource& _patsrc, /// pattern source
        HitSink& _sink,               /// hit sink
        Ebwt<TStr>& ebwtFw,           /// index of original text
        Ebwt<TStr>& ebwtBw,           /// index of mirror text
        vector<String<Dna5> >& os)    /// text strings, if available (empty otherwise)
{
	// Global intialization
	assert(revcomp);
	assert(!fullIndex);
	assert_leq(seedMms, 3);
	uint32_t numQs = ((qUpto == 0xffffffff) ? 16 * 1024 * 1024 : qUpto);
	SyncBitset doneMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the read mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");
	// No need to keep track of which reads are aligned because the
	// user hasn't requested an unaligned-read dump
	if(!_sink.dumpsReads()) numQs = 0;
	SyncBitset hitMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the hit mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");

	seededQualSearch_patsrc   = &_patsrc;
	seededQualSearch_sink     = &_sink;
	seededQualSearch_ebwtFw   = &ebwtFw;
	seededQualSearch_ebwtBw   = &ebwtBw;
	seededQualSearch_os       = &os;
	seededQualSearch_doneMask = &doneMask;
	seededQualSearch_hitMask  = &hitMask;
	seededQualSearch_pamFw    = NULL;
	seededQualSearch_pamRc    = NULL;
	seededQualSearch_qualCutoff = qualCutoff;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
	int *tids = new int[nthreads-1];
#endif

	SWITCH_TO_FW_INDEX();
	{
		// Phase 1: Consider cases 1R and 2R
		const char * msg = "Seeded quality search Phase 1 of 4: ";
		if(seedMms == 0) {
			msg = "Seeded quality search Phase 1 of 2: ";
		}
		Timer _t(cerr, msg, timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase1, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		seededQualSearchWorkerPhase1((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
	// Unload forward index and load mirror index
	SWITCH_TO_BW_INDEX();
	PartialAlignmentManager *pamRc = NULL;
	try {
		if(seedMms > 0) pamRc = new PartialAlignmentManager();
	} catch(bad_alloc& ba) {
		cerr << "Could not reserve space for PartialAlignmentManager" << endl;
		cerr << "Please subdivide the read set and invoke bowtie separately for each subdivision" << endl;
		throw 1;
	}
	seededQualSearch_pamRc = pamRc;
	{
		// Phase 2: Consider cases 1F, 2F and 3F and generate seedlings
		// for case 4R
		const char * msg = "Seeded quality search Phase 2 of 4: ";
		if(seedMms == 0) {
			msg = "Seeded quality search Phase 2 of 2: ";
		}
		Timer _t(cerr, msg, timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase2, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		seededQualSearchWorkerPhase2((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
	if(seedMms == 0) {
		// If we're not allowing any mismatches in the seed, then there
		// is no need to continue to phases 3 and 4
		assert(pamRc == NULL);
		return;
	}
	// Unload mirror index and load forward index
	SWITCH_TO_FW_INDEX();
	PartialAlignmentManager *pamFw = NULL;
	try {
		if(seedMms > 0) pamFw = new PartialAlignmentManager();
	} catch(bad_alloc& ba) {
		cerr << "Could not reserve space for PartialAlignmentManager" << endl;
		cerr << "Please subdivide the read set and invoke bowtie separately for each subdivision" << endl;
		throw 1;
	}
	seededQualSearch_pamFw = pamFw;
	{
		// Phase 3: Consider cases 3R and 4R and generate seedlings for
		// case 4F
		Timer _t(cerr, "Seeded quality search Phase 3 of 4: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase3, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		seededQualSearchWorkerPhase3((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
	// Some with the reverse-complement partial alignments
	if(pamRc != NULL) {
		delete pamRc;
		pamRc = NULL;
		seededQualSearch_pamRc = NULL;
	}
	// Unload forward index and load mirror index
	SWITCH_TO_BW_INDEX();
	{
		// Phase 4: Consider case 4F
		Timer _t(cerr, "Seeded quality search Phase 4 of 4: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase4, (void *)&tids[i]);
		}
#endif
		int tmp = 0;
		seededQualSearchWorkerPhase4((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
	if(pamFw != NULL) {
		delete pamFw;
		pamFw = NULL;
		seededQualSearch_pamFw = NULL;
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
	delete[] tids;
#endif
}

/**
 * Search for a good alignments for each read using criteria that
 * correspond somewhat faithfully to Maq's.  Search is aided by a pair
 * of Ebwt indexes, one for the original references, and one for the
 * transpose of the references.  Neither index should be loaded upon
 * entry to this function.
 *
 * Like Maq, we treat the first 24 base pairs of the read (those
 * closest to the 5' end) differently from the remainder of the read.
 * We call the first 24 base pairs the "seed."
 */
template<typename TStr>
static void seededQualCutoffSearchFull(
        int seedLen,                    /// length of seed (not a maq option)
        int qualCutoff,                 /// maximum sum of mismatch qualities
                                        /// like maq map's -e option
                                        /// default: 70
        int seedMms,                    /// max # mismatches allowed in seed
                                        /// (like maq map's -n option)
                                        /// Can only be 1 or 2, default: 1
        PairedPatternSource& _patsrc,   /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of mirror text
        vector<String<Dna5> >& os)      /// text strings, if available (empty otherwise)
{
	// Global intialization
	assert(fullIndex);
	assert(revcomp);
	assert_leq(seedMms, 3);

	seededQualSearch_patsrc   = &_patsrc;
	seededQualSearch_sink     = &_sink;
	seededQualSearch_ebwtFw   = &ebwtFw;
	seededQualSearch_ebwtBw   = &ebwtBw;
	seededQualSearch_os       = &os;
	seededQualSearch_doneMask = NULL;
	seededQualSearch_hitMask  = NULL;
	seededQualSearch_pamFw    = NULL;
	seededQualSearch_pamRc    = NULL;
	seededQualSearch_qualCutoff = qualCutoff;

	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	if((mates1.size() > 0 || mates12.size() > 0) && mixedThresh < 0xffffffff) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, sanityCheck, NULL, &os, false, useMm, useShmem, mmSweep, verbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	seededQualSearch_refs = refs;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
	int *tids = new int[nthreads-1];
#endif

	SWITCH_TO_FW_INDEX();
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(!noRefNames, startVerbose);
	}
	CHUD_START();
	{
		// Phase 1: Consider cases 1R and 2R
		Timer _t(cerr, "Seeded quality full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			tids[i] = i+1;
			if(stateful) pthread_create(&threads[i], &pt_attr,
			                            seededQualSearchWorkerFullStateful,
			                            (void*)&tids[i]);
			else         pthread_create(&threads[i], &pt_attr,
			                            seededQualSearchWorkerFull,
			                            (void*)&tids[i]);
		}
#endif
		int tmp = 0;
		if(stateful) seededQualSearchWorkerFullStateful((void*)&tmp);
		else         seededQualSearchWorkerFull((void*)&tmp);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				throw 1;
			}
		}
#endif
	}
	if(refs != NULL) {
		delete refs;
	}
	ebwtBw.evictFromMemory();
#ifdef BOWTIE_PTHREADS
	delete[] threads;
	delete[] tids;
#endif
}

/**
 * Return a new dynamically allocated PatternSource for the given
 * format, using the given list of strings as the filenames to read
 * from or as the sequences themselves (i.e. if -c was used).
 */
static PatternSource*
patsrcFromStrings(int format, const vector<string>& qs) {
	switch(format) {
		case FASTA:
			return new FastaPatternSource (qs, randomizeQuals,
			                               useSpinlock,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               forgiveInput,
			                               skipReads);
		case FASTA_CONT:
			return new FastaContinuousPatternSource (
			                               qs, fastaContLen,
			                               fastaContFreq,
			                               useSpinlock,
			                               patDumpfile, verbose,
			                               skipReads);
		case RAW:
			return new RawPatternSource   (qs, randomizeQuals,
			                               useSpinlock,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               skipReads);
		case FASTQ:
			return new FastqPatternSource (qs, randomizeQuals,
			                               useSpinlock,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               forgiveInput,
			                               solexaQuals, phred64Quals,
			                               integerQuals, fuzzy,
			                               skipReads);
		case TAB_MATE:
			return new TabbedPatternSource(qs, randomizeQuals,
			                               useSpinlock,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               forgiveInput,
			                               solexaQuals, phred64Quals,
			                               integerQuals, skipReads);
		case CMDLINE:
			return new VectorPatternSource(qs, randomizeQuals,
			                               useSpinlock,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               skipReads);
		case INPUT_CHAIN:
			return new ChainPatternSource (qs, useSpinlock, patDumpfile,
			                               verbose, skipReads);
		case RANDOM:
			return new RandomPatternSource(2000000, lenRandomReads,
			                               useSpinlock, patDumpfile,
			                               verbose, seed);
		default: {
			cerr << "Internal error; bad patsrc format: " << format << endl;
			throw 1;
		}
	}
}

#define PASS_DUMP_FILES \
	dumpAlFaBase, dumpAlFqBase, dumpAlBase, \
	dumpUnalFaBase, dumpUnalFqBase, dumpUnalBase, \
	dumpMaxFaBase, dumpMaxFqBase, dumpMaxBase

static string argstr;

template<typename TStr>
static void driver(const char * type,
                   const string& ebwtFileBase,
                   const string& query,
                   const vector<string>& queries,
                   const string& outfile)
{
	if(verbose || startVerbose)  {
		cerr << "Entered driver(): "; logTime(cerr, true);
	}
	// Vector of the reference sequences; used for sanity-checking
	vector<String<Dna5> > os;
	// Read reference sequences from the command-line or from a FASTA file
	if(!origString.empty()) {
		// Determine if it's a file by looking at whether it has a FASTA-like
		// extension
		size_t len = origString.length();
		if((len >= 6 && origString.substr(len-6) == ".fasta") ||
		   (len >= 4 && origString.substr(len-4) == ".mfa")   ||
		   (len >= 4 && origString.substr(len-4) == ".fas")   ||
		   (len >= 4 && origString.substr(len-4) == ".fna")   ||
		   (len >= 3 && origString.substr(len-3) == ".fa"))
		{
			// Read fasta file
			vector<string> origFiles;
			tokenize(origString, ",", origFiles);
			readSequenceFiles<String<Dna5>, Fasta>(origFiles, os);
		} else {
			// Read sequence
			readSequenceString(origString, os);
		}
	}
	// Adjust
	adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, verbose);

	vector<PatternSource*> patsrcs_a;
	vector<PatternSource*> patsrcs_b;
	vector<PatternSource*> patsrcs_ab;

	// If there were any first-mates specified, we will operate in
	// stateful mode
	bool paired = mates1.size() > 0 || mates12.size() > 0;
	if(paired) stateful = true;

	// Create list of pattern sources for paired reads appearing
	// interleaved in a single file
	if(verbose || startVerbose) {
		cerr << "Creating paired-end patsrcs: "; logTime(cerr, true);
	}
	for(size_t i = 0; i < mates12.size(); i++) {
		if(mates12[i] == "-" && !fullIndex) {
			cerr << "Input file \"-\" is not compatible with -z/--phased" << endl;
			throw 1;
		}
		const vector<string>* qs = &mates12;
		vector<string> tmp;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmp;
			tmp.push_back(mates12[i]);
			assert_eq(1, tmp.size());
		}
		patsrcs_ab.push_back(patsrcFromStrings(format, *qs));
		if(!fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < mates1.size(); i++) {
		if(mates1[i] == "-" && !fullIndex) {
			cerr << "Input file \"-\" is not compatible with -z/--phased" << endl;
			throw 1;
		}
		const vector<string>* qs = &mates1;
		vector<string> tmp;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmp;
			tmp.push_back(mates1[i]);
			assert_eq(1, tmp.size());
		}
		patsrcs_a.push_back(patsrcFromStrings(format, *qs));
		if(!fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < mates2.size(); i++) {
		if(mates2[i] == "-" && !fullIndex) {
			cerr << "Input file \"-\" is not compatible with -z/--phased" << endl;
			throw 1;
		}
		const vector<string>* qs = &mates2;
		vector<string> tmp;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmp;
			tmp.push_back(mates2[i]);
			assert_eq(1, tmp.size());
		}
		patsrcs_b.push_back(patsrcFromStrings(format, *qs));
		if(!fileParallel) {
			break;
		}
	}
	// All mates/mate files must be paired
	assert_eq(patsrcs_a.size(), patsrcs_b.size());

	// Create list of pattern sources for the unpaired reads
	if(verbose || startVerbose) {
		cerr << "Creating single-end patsrcs: "; logTime(cerr, true);
	}
	for(size_t i = 0; i < queries.size(); i++) {
		if(queries[i] == "-" && !fullIndex) {
			cerr << "Input file \"-\" is not compatible with -z/--phased" << endl;
			throw 1;
		}
		const vector<string>* qs = &queries;
		PatternSource* patsrc = NULL;
		vector<string> tmp;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmp;
			tmp.push_back(queries[i]);
			assert_eq(1, tmp.size());
		}
		patsrc = patsrcFromStrings(format, *qs);
		assert(patsrc != NULL);
		patsrcs_a.push_back(patsrc);
		patsrcs_b.push_back(NULL);
		if(!fileParallel) {
			break;
		}
	}

	if(verbose || startVerbose) {
		cerr << "Creating PatternSource: "; logTime(cerr, true);
	}
	PairedPatternSource *patsrc = NULL;
	if(mates12.size() > 0) {
		patsrc = new PairedSoloPatternSource(patsrcs_ab);
	} else {
		patsrc = new PairedDualPatternSource(patsrcs_a, patsrcs_b);
	}

	// Open hit output file
	if(verbose || startVerbose) {
		cerr << "Opening hit output file: "; logTime(cerr, true);
	}
	OutFileBuf *fout;
	if(!outfile.empty()) {
		if(refOut) {
			fout = NULL;
			if(!quiet) {
				cerr << "Warning: ignoring alignment output file " << outfile << " because --refout was specified" << endl;
			}
		} else {
			fout = new OutFileBuf(outfile.c_str(), outType == OUTPUT_BINARY);
		}
	} else {
		if(outType == OUTPUT_BINARY && !refOut) {
			cerr << "Error: Must specify an output file when output mode is binary" << endl;
			throw 1;
		}
		else if(outType == OUTPUT_CHAIN) {
			cerr << "Error: Must specify an output file when output mode is --chain" << endl;
			throw 1;
		}
		fout = new OutFileBuf();
	}
	ReferenceMap* rmap = NULL;
	if(refMapFile != NULL) {
		if(verbose || startVerbose) {
			cerr << "About to load in a reference map file with name "
			     << refMapFile << ": "; logTime(cerr, true);
		}
		rmap = new ReferenceMap(refMapFile, !noRefNames);
	}
	AnnotationMap* amap = NULL;
	if(annotMapFile != NULL) {
		if(verbose || startVerbose) {
			cerr << "About to load in an annotation map file with name "
			     << annotMapFile << ": "; logTime(cerr, true);
		}
		amap = new AnnotationMap(annotMapFile);
	}
	// Initialize Ebwt object and read in header
	if(verbose || startVerbose) {
		cerr << "About to initialize fw Ebwt: "; logTime(cerr, true);
	}
	Ebwt<TStr> ebwt(adjustedEbwtFileBase,
	                true,     // index is for the forward direction
	                /* overriding: */ offRate,
	                /* overriding: */ isaRate,
	                useMm,    // whether to use memory-mapped files
	                useShmem, // whether to use shared memory
	                mmSweep,  // sweep memory-mapped files
	                !noRefNames, // load names?
	                rmap,     // reference map, or NULL if none is needed
	                verbose, // whether to be talkative
	                startVerbose, // talkative during initialization
	                false /*passMemExc*/,
	                sanityCheck);
	Ebwt<TStr>* ebwtBw = NULL;
	// We need the mirror index if mismatches are allowed
	if(mismatches > 0 || maqLike) {
		if(verbose || startVerbose) {
			cerr << "About to initialize rev Ebwt: "; logTime(cerr, true);
		}
		ebwtBw = new Ebwt<TStr>(adjustedEbwtFileBase + ".rev",
		                        false, // index is for the reverse direction
		                        /* overriding: */ offRate,
		                        /* overriding: */ isaRate,
		                        useMm,    // whether to use memory-mapped files
		                        useShmem, // whether to use shared memory
		                        mmSweep,  // sweep memory-mapped files
		                        !noRefNames, // load names?
		                        rmap,     // reference map, or NULL if none is needed
		                        verbose,  // whether to be talkative
		                        startVerbose, // talkative during initialization
		                        false /*passMemExc*/,
		                        sanityCheck);
	}
	if(sanityCheck && !os.empty()) {
		// Sanity check number of patterns and pattern lengths in Ebwt
		// against original strings
		assert_eq(os.size(), ebwt.nPat());
		for(size_t i = 0; i < os.size(); i++) {
			assert_eq(length(os[i]), ebwt.plen()[i]);
		}
	}
	// Sanity-check the restored version of the Ebwt
	if(sanityCheck && !os.empty()) {
		ebwt.loadIntoMemory(!noRefNames, startVerbose);
		ebwt.checkOrigs(os, false);
		ebwt.evictFromMemory();
	}
	{
		Timer _t(cerr, "Time searching: ", timing);
		if(verbose || startVerbose) {
			cerr << "Creating HitSink: "; logTime(cerr, true);
		}
		// Set up hit sink; if sanityCheck && !os.empty() is true,
		// then instruct the sink to "retain" hits in a vector in
		// memory so that we can easily sanity check them later on
		HitSink *sink;
		RecalTable *table = NULL;
		if(recal) {
			table = new RecalTable(recalMaxCycle, recalMaxQual, recalQualShift);
		}
		vector<string>* refnames = &ebwt.refnames();
		if(noRefNames) refnames = NULL;
		switch(outType) {
			case OUTPUT_FULL:
				if(refOut) {
					sink = new VerboseHitSink(
							ebwt.nPat(), offBase, rmap, amap,
							fullRef, PASS_DUMP_FILES,
							format == TAB_MATE,
							table, refnames, partitionSz);
				} else {
					sink = new VerboseHitSink(
							fout, offBase, rmap, amap,
							fullRef, PASS_DUMP_FILES,
							format == TAB_MATE,
							table, refnames, partitionSz);
				}
				break;
			case OUTPUT_SAM:
				if(refOut) {
					throw 1;
				} else {
					SAMHitSink *sam = new SAMHitSink(
							fout, 1, rmap, amap,
							fullRef, PASS_DUMP_FILES,
							format == TAB_MATE,
							table, refnames);
					if(!samNoHead) {
						vector<string> refnames;
						if(!samNoSQ) {
							readEbwtRefnames(adjustedEbwtFileBase, refnames);
						}
						sam->appendHeaders(
								sam->out(0), ebwt.nPat(),
								refnames, samNoSQ, rmap,
								ebwt.plen(), fullRef,
								argstr.c_str());
					}
					sink = sam;
				}
				break;
			case OUTPUT_CONCISE:
				if(refOut) {
					sink = new ConciseHitSink(
							ebwt.nPat(), offBase,
							PASS_DUMP_FILES,
							format == TAB_MATE,
							table, refnames, reportOpps);
				} else {
					sink = new ConciseHitSink(
							fout, offBase,
							PASS_DUMP_FILES,
							format == TAB_MATE,
							table, refnames, reportOpps);
				}
				break;
			case OUTPUT_BINARY:
				if(refOut) {
					sink = new BinaryHitSink(
							ebwt.nPat(), offBase,
							PASS_DUMP_FILES,
							format == TAB_MATE,
							table, refnames);
				} else {
					sink = new BinaryHitSink(
							fout, offBase,
							PASS_DUMP_FILES,
							format == TAB_MATE,
							table, refnames);
				}
				break;
			case OUTPUT_CHAIN:
				assert(stateful);
				sink = new ChainingHitSink(
						fout, strata, amap,
						PASS_DUMP_FILES,
						true,
						table, refnames);
				break;
			case OUTPUT_NONE:
				sink = new StubHitSink();
				break;
			default:
				cerr << "Invalid output type: " << outType << endl;
				throw 1;
		}
    	if(verbose || startVerbose) {
    		cerr << "Dispatching to search driver: "; logTime(cerr, true);
    	}
		if(maqLike) {
			if(!fullIndex) {
				seededQualCutoffSearch(seedLen,
				                       qualThresh,
				                       seedMms,
				                       *patsrc,
				                       *sink,
				                       ebwt,    // forward index
				                       *ebwtBw, // mirror index (not optional)
				                       os);     // references, if available
			} else {
				seededQualCutoffSearchFull(seedLen,
				                           qualThresh,
				                           seedMms,
				                           *patsrc,
				                           *sink,
				                           ebwt,    // forward index
				                           *ebwtBw, // mirror index (not optional)
				                           os);     // references, if available
			}
		}
		else if(mismatches > 0) {
			if(mismatches == 1) {
				if(!fullIndex) {
					mismatchSearch(*patsrc, *sink, ebwt, *ebwtBw, os);
				} else {
					assert(ebwtBw != NULL);
					mismatchSearchFull(*patsrc, *sink, ebwt, *ebwtBw, os);
				}
			} else if(mismatches == 2 || mismatches == 3) {
				if(!fullIndex) {
					twoOrThreeMismatchSearch(*patsrc, *sink, ebwt, *ebwtBw, os, mismatches == 2);
				} else {
					twoOrThreeMismatchSearchFull(*patsrc, *sink, ebwt, *ebwtBw, os, mismatches == 2);
				}
			} else {
				cerr << "Error: " << mismatches << " is not a supported number of mismatches" << endl;
				throw 1;
			}
		} else {
			// Search without mismatches
			// Note that --fast doesn't make a difference here because
			// we're only loading half of the index anyway
			exactSearch(*patsrc, *sink, ebwt, os, paired);
		}
		// Evict any loaded indexes from memory
		if(ebwt.isInMemory()) {
			ebwt.evictFromMemory();
		}
		if(ebwtBw != NULL) {
			delete ebwtBw;
		}
		if(!quiet) {
			sink->finish(hadoopOut); // end the hits section of the hit file
		}
		for(size_t i = 0; i < patsrcs_a.size(); i++) {
			assert(patsrcs_a[i] != NULL);
			delete patsrcs_a[i];
		}
		for(size_t i = 0; i < patsrcs_b.size(); i++) {
			if(patsrcs_b[i] != NULL) {
				delete patsrcs_b[i];
			}
		}
		for(size_t i = 0; i < patsrcs_ab.size(); i++) {
			if(patsrcs_ab[i] != NULL) {
				delete patsrcs_ab[i];
			}
		}
		delete patsrc;
		delete sink;
		delete amap;
		delete rmap;
		if(fout != NULL) delete fout;
	}
}

// C++ name mangling is disabled for the bowtie() function to make it
// easier to use Bowtie as a library.
extern "C" {

/**
 * Main bowtie entry function.  Parses argc/argv style command-line
 * options, sets global configuration variables, and calls the driver()
 * function.
 */
int bowtie(int argc, const char **argv) {
	try {
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();
		for(int i = 0; i < argc; i++) {
			argstr += argv[i];
			if(i < argc-1) argstr += " ";
		}
		string ebwtFile;  // read serialized Ebwt from this file
		string query;   // read query string(s) from this file
		vector<string> queries;
		string outfile; // write query results to this file
		if(startVerbose) { cerr << "Entered main(): "; logTime(cerr, true); }
		parseOptions(argc, argv);
		argv0 = argv[0];
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
	#ifdef CHUD_PROFILING
		chudInitialize();
		chudAcquireRemoteAccess();
	#endif
		{
			Timer _t(cerr, "Overall time: ", timing);
			if(startVerbose) {
				cerr << "Parsing index and read arguments: "; logTime(cerr, true);
			}

			// Get index basename
			if(optind >= argc) {
				cerr << "No index, query, or output file specified!" << endl;
				printUsage(cerr);
				return 1;
			}
			ebwtFile = argv[optind++];

			// Get query filename
			if(optind >= argc) {
				if(mates1.size() > 0 || mates12.size() > 0) {
					query = "";
				} else {
					cerr << "No query or output file specified!" << endl;
					printUsage(cerr);
					return 1;
				}
			} else if (mates1.size() == 0 && mates12.size() == 0) {
				query = argv[optind++];
				// Tokenize the list of query files
				tokenize(query, ",", queries);
				if(queries.size() < 1) {
					cerr << "Tokenized query file list was empty!" << endl;
					printUsage(cerr);
					return 1;
				}
			}

			// Get output filename
			if(optind < argc) {
				outfile = argv[optind++];
			}

			// Extra parametesr?
			if(optind < argc) {
				cerr << "Extra parameter(s) specified: ";
				for(int i = optind; i < argc; i++) {
					cerr << "\"" << argv[i] << "\"";
					if(i < argc-1) cerr << ", ";
				}
				cerr << endl;
				if(mates1.size() > 0) {
					cerr << "Note that if <mates> files are specified using -1/-2, a <singles> file cannot" << endl
						 << "also be specified.  Please run bowtie separately for mates and singles." << endl;
				}
				throw 1;
			}

			// Optionally summarize
			if(verbose) {
				cout << "Input ebwt file: \"" << ebwtFile << "\"" << endl;
				cout << "Query inputs (DNA, " << file_format_names[format] << "):" << endl;
				for(size_t i = 0; i < queries.size(); i++) {
					cout << "  " << queries[i] << endl;
				}
				cout << "Output file: \"" << outfile << "\"" << endl;
				cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
				cout << "Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
			#ifdef NDEBUG
				cout << "Assertions: disabled" << endl;
			#else
				cout << "Assertions: enabled" << endl;
			#endif
			}
			if(ipause) {
				cout << "Press key to continue..." << endl;
				getchar();
			}
			driver<String<Dna, Alloc<> > >("DNA", ebwtFile, query, queries, outfile);
			CHUD_STOP();
		}
#ifdef CHUD_PROFILING
		chudReleaseRemoteAccess();
#endif
		return 0;
	} catch(exception& e) {
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
} // bowtie()
} // extern "C"
