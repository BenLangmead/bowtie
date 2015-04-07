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
#include "ebwt_search.h"
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
bool quiet;        // print nothing but the alignments
static int sanityCheck;   // enable expensive sanity checks
static int format;        // default read format is FASTQ
static string origString; // reference text, or filename(s)
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
static bool noRefNames;       // true -> print reference indexes; not names
static string dumpAlBase;     // basename of same-format files to dump aligned reads to
static string dumpUnalBase;   // basename of same-format files to dump unaligned reads to
static string dumpMaxBase;    // basename of same-format files to dump reads with more than -m valid alignments to
static uint32_t khits;  // number of hits per read; >1 is much slower
static uint32_t mhits;  // don't report any hits if there are > mhits
static bool better;     // true -> guarantee alignments from best possible stratum
static bool strata;     // true -> don't stop at stratum boundaries
static bool refOut;     // if true, alignments go to per-ref files
static int partitionSz; // output a partitioning key in first field
static bool noMaqRound; // true -> don't round quals to nearest 10 like maq
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
static bool mateFwSet;         // true -> user set --ff/--fr/--rf
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
static bool samNoQnameTrunc; // don't truncate QNAME field at first whitespace
static bool samNoHead; // don't print any header lines in SAM output
static bool samNoSQ;   // don't print @SQ header lines
bool color;     // true -> inputs are colorspace
bool colorExEnds; // true -> nucleotides on either end of decoded cspace alignment should be excluded
static string rgs; // SAM outputs for @RG header line
int snpPhred; // probability of SNP, for scoring colorspace alignments
static Bitset suppressOuts(64); // output fields to suppress
static bool sampleMax; // whether to report a random alignment when maxed-out via -m/-M
static int defaultMapq; // default mapping quality to print in SAM mode
bool colorSeq; // true -> show colorspace alignments as colors, not decoded bases
bool colorQual; // true -> show colorspace qualities as original quals, not decoded quals
static bool printCost; // true -> print stratum and cost
bool showSeed;
static vector<string> qualities;
static vector<string> qualities1;
static vector<string> qualities2;
static string wrapper; // Type of wrapper script
bool gAllowMateContainment;
bool gReportColorPrimer;
MUTEX_T gLock;

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
	noRefNames				= false; // true -> print reference indexes; not names
	dumpAlBase				= "";    // basename of same-format files to dump aligned reads to
	dumpUnalBase			= "";    // basename of same-format files to dump unaligned reads to
	dumpMaxBase				= "";    // basename of same-format files to dump reads with more than -m valid alignments to
	khits					= 1;     // number of hits per read; >1 is much slower
	mhits					= 0xffffffff; // don't report any hits if there are > mhits
	better					= false; // true -> guarantee alignments from best possible stratum
	strata					= false; // true -> don't stop at stratum boundaries
	refOut					= false; // if true, alignments go to per-ref files
	partitionSz				= 0;     // output a partitioning key in first field
	noMaqRound				= false; // true -> don't round quals to nearest 10 like maq
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
	mateFwSet				= false; // true -> user set mate1fw/mate2fw with --ff/--fr/--rf
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
	chunkPoolMegabytes		= 64;    // max MB to dedicate to best-first search frames per thread
	chunkSz					= 256;   // size of single chunk disbursed by ChunkPool (in KB)
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
	samNoQnameTrunc         = false; // don't truncate at first whitespace?
	samNoHead				= false; // don't print any header lines in SAM output
	samNoSQ					= false; // don't print @SQ header lines
	color					= false; // don't align in colorspace by default
	colorExEnds				= true;  // true -> nucleotides on either end of decoded cspace alignment should be excluded
	rgs						= "";    // SAM outputs for @RG header line
	snpPhred				= 30;    // probability of SNP, for scoring colorspace alignments
	suppressOuts.clear();            // output fields to suppress
	sampleMax				= false;
	defaultMapq				= 255;
	colorSeq				= false; // true -> show colorspace alignments as colors, not decoded bases
	colorQual				= false; // true -> show colorspace qualities as original quals, not decoded quals
	printCost				= false; // true -> print cost and stratum
	showSeed				= false; // true -> print per-read pseudo-random seed
	qualities.clear();
	qualities1.clear();
	qualities2.clear();
	wrapper.clear();
	gAllowMateContainment	= false; // true -> alignments where one mate lies inside the other are valid
	gReportColorPrimer		= false; // true -> print flag with trimmed color primer and downstream color
}

// mating constraints

static const char *short_options = "fF:qbzhcu:rv:s:at3:5:o:e:n:l:w:p:k:m:M:1:2:I:X:x:B:ySCQ:";

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
	ARG_SANITY,
	ARG_OLDBEST,
	ARG_BETTER,
	ARG_BEST,
	ARG_REFOUT,
	ARG_ISARATE,
	ARG_PARTITION,
	ARG_integerQuals,
	ARG_NOMAQROUND,
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
	ARG_REFMAP,
	ARG_ANNOTMAP,
	ARG_REPORTSE,
	ARG_HADOOPOUT,
	ARG_FUZZY,
	ARG_FULLREF,
	ARG_USAGE,
	ARG_SNPPHRED,
	ARG_SNPFRAC,
	ARG_SAM_NO_QNAME_TRUNC,
	ARG_SAM_NOHEAD,
	ARG_SAM_NOSQ,
	ARG_SAM_RG,
	ARG_SUPPRESS_FIELDS,
	ARG_DEFAULT_MAPQ,
	ARG_COLOR_SEQ,
	ARG_COLOR_QUAL,
	ARG_COST,
	ARG_COLOR_KEEP_ENDS,
	ARG_SHOWSEED,
	ARG_QUALS1,
	ARG_QUALS2,
	ARG_ALLOW_CONTAIN,
	ARG_COLOR_PRIMER,
	ARG_WRAPPER
};

static struct option long_options[] = {
	{(char*)"verbose",      no_argument,       0,            ARG_VERBOSE},
	{(char*)"startverbose", no_argument,       0,            ARG_STARTVERBOSE},
	{(char*)"quiet",        no_argument,       0,            ARG_QUIET},
	{(char*)"sanity",       no_argument,       0,            ARG_SANITY},
	{(char*)"pause",        no_argument,       &ipause,      1},
	{(char*)"orig",         required_argument, 0,            ARG_ORIG},
	{(char*)"all",          no_argument,       0,            'a'},
	{(char*)"concise",      no_argument,       0,            ARG_CONCISE},
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
	{(char*)"quals",        required_argument, 0,            'Q'},
	{(char*)"Q1",           required_argument, 0,            ARG_QUALS1},
	{(char*)"Q2",           required_argument, 0,            ARG_QUALS2},
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
	{(char*)"refout",       no_argument,       0,            ARG_REFOUT},
	{(char*)"partition",    required_argument, 0,            ARG_PARTITION},
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
	{(char*)"refmap",       required_argument, 0,            ARG_REFMAP},
	{(char*)"annotmap",     required_argument, 0,            ARG_ANNOTMAP},
	{(char*)"reportse",     no_argument,       0,            ARG_REPORTSE},
	{(char*)"hadoopout",    no_argument,       0,            ARG_HADOOPOUT},
	{(char*)"fuzzy",        no_argument,       0,            ARG_FUZZY},
	{(char*)"fullref",      no_argument,       0,            ARG_FULLREF},
	{(char*)"usage",        no_argument,       0,            ARG_USAGE},
	{(char*)"sam",          no_argument,       0,            'S'},
	{(char*)"sam-no-qname-trunc", no_argument, 0,            ARG_SAM_NO_QNAME_TRUNC},
	{(char*)"sam-nohead",   no_argument,       0,            ARG_SAM_NOHEAD},
	{(char*)"sam-nosq",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"sam-noSQ",     no_argument,       0,            ARG_SAM_NOSQ},
	{(char*)"color",        no_argument,       0,            'C'},
	{(char*)"sam-RG",       required_argument, 0,            ARG_SAM_RG},
	{(char*)"snpphred",     required_argument, 0,            ARG_SNPPHRED},
	{(char*)"snpfrac",      required_argument, 0,            ARG_SNPFRAC},
	{(char*)"suppress",     required_argument, 0,            ARG_SUPPRESS_FIELDS},
	{(char*)"mapq",         required_argument, 0,            ARG_DEFAULT_MAPQ},
	{(char*)"col-cseq",     no_argument,       0,            ARG_COLOR_SEQ},
	{(char*)"col-cqual",    no_argument,       0,            ARG_COLOR_QUAL},
	{(char*)"col-keepends", no_argument,       0,            ARG_COLOR_KEEP_ENDS},
	{(char*)"cost",         no_argument,       0,            ARG_COST},
	{(char*)"showseed",     no_argument,       0,            ARG_SHOWSEED},
	{(char*)"allow-contain",no_argument,       0,            ARG_ALLOW_CONTAIN},
	{(char*)"col-primer",   no_argument,       0,            ARG_COLOR_PRIMER},
	{(char*)"wrapper",      required_argument, 0,            ARG_WRAPPER},
	{(char*)0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
#ifdef BOWTIE_64BIT_INDEX
	string tool_name = "bowtie-build-l";
#else
	string tool_name = "bowtie-build-s";
#endif
	if(wrapper == "basic-0") {
		tool_name = "bowtie";
	}

	out << "Usage: " << endl
        << tool_name << " [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]" << endl
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
	    << "  -C                 reads and index are in colorspace" << endl
	    << "  -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C" << endl
	    << "  --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively" << endl
	    << "  -s/--skip <int>    skip the first <int> reads/pairs in the input" << endl
	    << "  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads" << endl
		<< "  --phred33-quals    input quals are Phred+33 (default)" << endl
		<< "  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)" << endl
		<< "  --solexa-quals     input quals are from GA Pipeline ver. < 1.3" << endl
		<< "  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3" << endl
		<< "  --integer-quals    qualities are given as space-separated integers (not ASCII)" << endl;
		if(wrapper == "basic-0") {
		out << "  --large-index      force usage of a 'large' index, even if a small one is present" << endl;
		}
		out << "Alignment:" << endl
	    << "  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities" << endl
	    << "    or" << endl
	    << "  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)" << endl
	    << "  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)" << endl
	    << "  -l/--seedlen <int> seed length for -n (default: 28)" << endl
		<< "  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)" << endl
	    << "  -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)" << endl
	    << "  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)" << endl
	    << "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)" << endl
	    << "  --nofw/--norc      do not align to forward/reverse-complement reference strand" << endl
	    << "  --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)" << endl
	    << "  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)" << endl
	    << "  -y/--tryhard       try hard to find valid alignments, at the expense of speed" << endl
	    << "  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)" << endl
	    << "Reporting:" << endl
	    << "  -k <int>           report up to <int> good alignments per read (default: 1)" << endl
	    << "  -a/--all           report all alignments per read (much slower than low -k)" << endl
	    << "  -m <int>           suppress all alignments if > <int> exist (def: no limit)" << endl
	    << "  -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best" << endl
	    << "  --best             hits guaranteed best stratum; ties broken by quality" << endl
	    << "  --strata           hits in sub-optimal strata aren't reported (requires --best)" << endl
	    << "Output:" << endl
	    << "  -t/--time          print wall-clock time taken by search phases" << endl
	    << "  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)" << endl
	    << "  --quiet            print nothing but the alignments" << endl
	    << "  --refout           write alignments to files refXXXXX.map, 1 map per reference" << endl
	    << "  --refidx           refer to ref. seqs by 0-based index rather than name" << endl
	    << "  --al <fname>       write aligned reads/pairs to file(s) <fname>" << endl
	    << "  --un <fname>       write unaligned reads/pairs to file(s) <fname>" << endl
	    << "  --max <fname>      write reads/pairs over -m limit to file(s) <fname>" << endl
	    << "  --suppress <cols>  suppresses given columns (comma-delim'ed) in default output" << endl
	    << "  --fullref          write entire ref name (default: only up to 1st space)" << endl
	    << "Colorspace:" << endl
	    << "  --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)" << endl
	    << "     or" << endl
	    << "  --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred" << endl
	    << "  --col-cseq         print aligned colorspace seqs as colors, not decoded bases" << endl
	    << "  --col-cqual        print original colorspace quals, not decoded quals" << endl
	    << "  --col-keepends     keep nucleotides at extreme ends of decoded alignment" << endl
	    << "SAM:" << endl
	    << "  -S/--sam           write hits in SAM format" << endl
	    << "  --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments" << endl
	    << "  --sam-nohead       supppress header lines (starting with @) for SAM output" << endl
	    << "  --sam-nosq         supppress @SQ header lines for SAM output" << endl
	    << "  --sam-RG <text>    add <text> (usually \"lab=value\") to @RG line of SAM header" << endl
	    << "Performance:" << endl
	    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
	    << "  -p/--threads <int> number of alignment threads to launch (default: 1)" << endl
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
	    << "  -h/--help          print this usage message" << endl
	    ;
	if(wrapper.empty()) {
		cerr << endl
		     << "*** Warning ***" << endl
			 << tool_name << " was run directly.  It is recommended that you run the wrapper script 'bowtie' instead." << endl
			 << endl;
	}
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, int upper, const char *errmsg, const char *arg) {
	long l;
	char *endPtr= NULL;
	l = strtol(arg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower || l > upper) {
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
 * Parse from optarg by default.
 */
static int parseInt(int lower, const char *errmsg) {
	return parseInt(lower, INT_MAX, errmsg, optarg);
}

/**
 * Upper is INT_MAX by default.
 */
static int parseInt(int lower, const char *errmsg, const char *arg) {
	return parseInt(lower, INT_MAX, errmsg, arg);
}

/**
 * Upper is INT_MAX, parse from optarg by default.
 */
static int parseInt(int lower, int upper, const char *errmsg) {
	return parseInt(lower, upper, errmsg, optarg);
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
			case ARG_WRAPPER: wrapper = optarg; break;
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
			case 'C': color = true; break;
			case 'I':
				minInsert = (uint32_t)parseInt(0, "-I arg must be positive");
				break;
			case 'X':
				maxInsert = (uint32_t)parseInt(1, "-X arg must be at least 1");
				break;
			case 's':
				skipReads = (uint32_t)parseInt(0, "-s arg must be positive");
				break;
			case ARG_FF: mate1fw = true;  mate2fw = true;  mateFwSet = true; break;
			case ARG_RF: mate1fw = false; mate2fw = true;  mateFwSet = true; break;
			case ARG_FR: mate1fw = true;  mate2fw = false; mateFwSet = true; break;
			case ARG_RANDOM_READS: format = RANDOM; break;
			case ARG_RANDOM_READS_NOSYNC:
				format = RANDOM;
				randReadsNoSync = true;
				break;
			case ARG_RANGE: rangeMode = true; break;
			case ARG_CONCISE: outType = OUTPUT_CONCISE; break;
			case 'S': outType = OUTPUT_SAM; break;
			case ARG_REFOUT: refOut = true; break;
			case ARG_NOOUT: outType = OUTPUT_NONE; break;
			case ARG_REFMAP: refMapFile = optarg; break;
			case ARG_ANNOTMAP: annotMapFile = optarg; break;
			case ARG_SHMEM: useShmem = true; break;
			case ARG_COLOR_SEQ: colorSeq = true; break;
			case ARG_COLOR_QUAL: colorQual = true; break;
			case ARG_SHOWSEED: showSeed = true; break;
			case ARG_ALLOW_CONTAIN: gAllowMateContainment = true; break;
			case ARG_COLOR_PRIMER: gReportColorPrimer = true; break;
			case ARG_SUPPRESS_FIELDS: {
				vector<string> supp;
				tokenize(optarg, ",", supp);
				for(size_t i = 0; i < supp.size(); i++) {
					int ii = parseInt(1, "--suppress arg must be at least 1", supp[i].c_str());
					suppressOuts.set(ii-1);
				}
				break;
			}
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
			case ARG_NOMAQROUND: noMaqRound = true; break;
			case ARG_COLOR_KEEP_ENDS: colorExEnds = false; break;
			case ARG_SNPPHRED: snpPhred = parseInt(0, "--snpphred must be at least 0"); break;
			case ARG_SNPFRAC: {
				double p = parse<double>(optarg);
				if(p <= 0.0) {
					cerr << "Error: --snpfrac parameter must be > 0.0" << endl;
					throw 1;
				}
				p = (log10(p) * -10);
				snpPhred = (int)(p + 0.5);
				if(snpPhred < 10)
				cout << "snpPhred: " << snpPhred << endl;
				break;
			}
			case 'z': {
				cerr << "Error: -z/--phased mode is no longer supported" << endl;
				throw 1;
			}
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
			case 'Q':
				tokenize(optarg, ",", qualities);
				integerQuals = true;
				break;
			case ARG_QUALS1:
				tokenize(optarg, ",", qualities1);
				integerQuals = true;
				break;
			case ARG_QUALS2:
				tokenize(optarg, ",", qualities2);
				integerQuals = true;
				break;
			case 'M':
				sampleMax = true;
			case 'm':
				mhits = (uint32_t)parseInt(1, "-m arg must be at least 1");
				break;
			case 'x':
				mixedThresh = (uint32_t)parseInt(0, "-x arg must be at least 0");
				break;
			case ARG_MIXED_ATTEMPTS:
				mixedAttemptLim = (uint32_t)parseInt(1, "--pairtries arg must be at least 1");
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
				nthreads = parseInt(1, "-p/--threads arg must be at least 1");
				break;
			case ARG_FILEPAR:
				fileParallel = true;
				break;
			case 'v':
				maqLike = 0;
				mismatches = parseInt(0, 3, "-v arg must be at least 0 and at most 3");
				break;
			case '3': trim3 = parseInt(0, "-3/--trim3 arg must be at least 0"); break;
			case '5': trim5 = parseInt(0, "-5/--trim5 arg must be at least 0"); break;
			case 'o': offRate = parseInt(1, "-o/--offrate arg must be at least 1"); break;
			case ARG_ISARATE: isaRate = parseInt(0, "--isarate arg must be at least 0"); break;
			case 'e': qualThresh = parseInt(1, "-e/--err arg must be at least 1"); break;
			case 'n': seedMms = parseInt(0, 3, "-n/--seedmms arg must be at least 0 and at most 3"); maqLike = 1; break;
			case 'l': seedLen = parseInt(5, "-l/--seedlen arg must be at least 5"); break;
			case 'h': printUsage(cout); throw 0; break;
			case ARG_USAGE: printUsage(cout); throw 0; break;
			case 'a': allHits = true; break;
			case 'y': tryHard = true; break;
			case ARG_RECAL: recal = true; break;
			case ARG_CHUNKMBS: chunkPoolMegabytes = parseInt(1, "--chunkmbs arg must be at least 1"); break;
			case ARG_CHUNKSZ: chunkSz = parseInt(1, "--chunksz arg must be at least 1"); break;
			case ARG_CHUNKVERBOSE: chunkVerbose = true; break;
			case ARG_BETTER: stateful = true; better = true; break;
			case ARG_BEST: stateful = true; useV1 = false; break;
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
			case ARG_SAM_NO_QNAME_TRUNC: samNoQnameTrunc = true; break;
			case ARG_SAM_NOHEAD: samNoHead = true; break;
			case ARG_SAM_NOSQ: samNoSQ = true; break;
			case ARG_SAM_RG: {
				if(!rgs.empty()) rgs += '\t';
				rgs += optarg;
				break;
			}
			case ARG_COST: printCost = true; break;
			case ARG_DEFAULT_MAPQ:
				defaultMapq = parseInt(0, "--mapq must be positive");
				break;
			case ARG_MAXBTS: {
				maxBts  = parseInt(0, "--maxbts must be positive");
				maxBtsBetter = maxBts;
				break;
			}
			case ARG_DUMP_PATS: patDumpfile = optarg; break;
			case ARG_STRAND_FIX: strandFix = true; break;
			case ARG_RANDOMIZE_QUALS: randomizeQuals = true; break;
			case ARG_PARTITION: partitionSz = parse<int>(optarg); break;
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
	if(!maqLike && mismatches == 3) {
		// Much faster than normal 3-mismatch mode
		stateful = true;
	}
	if(mates1.size() != mates2.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
		     << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
		     << "sequences must be specified with -1 and -2." << endl;
		throw 1;
	}
	if(qualities.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with -Q but -f was not" << endl
		     << "enabled.  -Q works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities.size() && !color) {
		cerr << "Error: one or more quality files were specified with -Q but -C was not" << endl
		     << "enabled.  -Q works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q1 but -f was not" << endl
		     << "enabled.  --Q1 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() && !color) {
		cerr << "Error: one or more quality files were specified with --Q1 but -C was not" << endl
		     << "enabled.  --Q1 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities2.size() && format != FASTA) {
		cerr << "Error: one or more quality files were specified with --Q2 but -f was not" << endl
		     << "enabled.  --Q2 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities2.size() && !color) {
		cerr << "Error: one or more quality files were specified with --Q2 but -C was not" << endl
		     << "enabled.  --Q2 works only in combination with -f and -C." << endl;
		throw 1;
	}
	if(qualities1.size() > 0 && mates1.size() != qualities1.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << qualities1.size() << endl
		     << "quality files were specified with --Q1.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -1 and --Q1." << endl;
		throw 1;
	}
	if(qualities2.size() > 0 && mates2.size() != qualities2.size()) {
		cerr << "Error: " << mates2.size() << " mate files/sequences were specified with -2, but " << qualities2.size() << endl
		     << "quality files were specified with --Q2.  The same number of mate and quality" << endl
		     << "files must sequences must be specified with -2 and --Q2." << endl;
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
	if(tryHard) {
		// Increase backtracking limit to huge number
		maxBts = maxBtsBetter = INT_MAX;
		// Increase number of paired-end scan attempts to huge number
		mixedAttemptLim = UINT_MAX;
	}
	if(!stateful && sampleMax) {
		if(!quiet) {
			cerr << "Warning: -M was specified w/o --best; automatically enabling --best" << endl;
		}
		stateful = true;
	}
	if(strata && !stateful) {
		cerr << "--strata must be combined with --best" << endl;
		throw 1;
	}
	if(strata && !allHits && khits == 1 && mhits == 0xffffffff) {
		cerr << "--strata has no effect unless combined with -m, -a, or -k N where N > 1" << endl;
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
	if(snpPhred <= 10 && color && !quiet) {
		cerr << "Warning: the colorspace SNP penalty (--snpphred) is very low: " << snpPhred << endl;
	}
	if(outType == OUTPUT_SAM && refOut) {
		cerr << "Error: --refout cannot be combined with -S/--sam" << endl;
		throw 1;
	}
	if(!mateFwSet) {
		if(color) {
			// Set colorspace default (--ff)
			mate1fw = true;
			mate2fw = true;
		} else {
			// Set nucleotide space default (--fr)
			mate1fw = true;
			mate2fw = false;
		}
	}
	if(outType != OUTPUT_FULL && suppressOuts.count() > 0 && !quiet) {
		cerr << "Warning: Ignoring --suppress because output type is not default." << endl;
		cerr << "         --suppress is only available for the default output type." << endl;
		suppressOuts.clear();
	}
}

static const char *argv0 = NULL;

#define FINISH_READ(p) \
	/* Don't do finishRead if the read isn't legit or if the read was skipped by the doneMask */ \
	if(!p->empty()) { \
		sink->finishRead(*p, true, !skipped); \
	} \
	skipped = false;

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

#define WORKER_EXIT() \
	patsrcFact->destroy(patsrc); \
	delete patsrcFact; \
	sinkFact->destroy(sink); \
	delete sinkFact; \
	return;

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
createSinkFactory(HitSink& _sink) {
	HitSinkPerThreadFactory *sink = NULL;
	if(!strata) {
		// Unstratified
		if(!allHits) {
			// First N good; "good" inherently ignores strata
			sink = new NGoodHitSinkPerThreadFactory(_sink, khits, mhits);
		} else {
			// All hits, spanning strata
			sink = new AllHitSinkPerThreadFactory(_sink, mhits);
		}
	} else {
		// Stratified
		assert(stateful);
		if(!allHits) {
			assert(stateful);
			// Buffer best hits, assuming they're arriving in best-
			// to-worst order
			sink = new NBestFirstStratHitSinkPerThreadFactory(_sink, khits, mhits);
		} else {
			assert(stateful);
			// Buffer best hits, assuming they're arriving in best-
			// to-worst order
			sink = new NBestFirstStratHitSinkPerThreadFactory(_sink, 0xffffffff/2, mhits);
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
#ifdef WITH_TBB
void exactSearchWorker::operator()() {
#else
static void exactSearchWorker(void *vp) {
	int tid = *((int*)vp);
#endif
	PairedPatternSource& _patsrc = *exactSearch_patsrc;
	HitSink& _sink               = *exactSearch_sink;
	Ebwt<String<Dna> >& ebwt     = *exactSearch_ebwt;
	vector<String<Dna5> >& os    = *exactSearch_os;
	const BitPairReference* refs =  exactSearch_refs;

	// Per-thread initialization
	PatternSourcePerThreadFactory *patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread *patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
	HitSinkPerThread* sink = sinkFact->create();
	EbwtSearchParams<String<Dna> > params(
	        *sink,      // HitSink
	        os,         // reference sequences
	        true,       // read is forward
	        true);       // index is forward
	GreedyDFSRangeSource bt(
	        &ebwt, params,
	        refs,           // reference sequence (for colorspace)
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
		#include "search_exact.c"
	}
	FINISH_READ(patsrc);
	WORKER_EXIT();
}

/**
 * A statefulness-aware worker driver.  Uses UnpairedExactAlignerV1.
 */
#ifdef WITH_TBB
void exactSearchWorkerStateful::operator()() {
#else
static void exactSearchWorkerStateful(void *vp) {
	int tid = *((int*)vp);
#endif
	PairedPatternSource& _patsrc = *exactSearch_patsrc;
	HitSink& _sink               = *exactSearch_sink;
	Ebwt<String<Dna> >& ebwt     = *exactSearch_ebwt;
	vector<String<Dna5> >& os    = *exactSearch_os;
	BitPairReference* refs       =  exactSearch_refs;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);

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
			refs,
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
			color,
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
	return;
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
                        vector<String<Dna5> >& os)
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
		ebwt.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose);
	}

	BitPairReference *refs = NULL;
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	if(color || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, color, sanityCheck, NULL, &os, false, true, useMm, useShmem, mmSweep, verbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	exactSearch_refs   = refs;
#ifdef WITH_TBB
	tbb::task_group tbb_grp;
#else
	AutoArray<tthread::thread*> threads(nthreads+1);
	AutoArray<int> tids(nthreads+1);
#endif
	CHUD_START();
	{
		Timer _t(cerr, "Time for 0-mismatch search: ", timing);

		for(int i = 1; i <= nthreads; i++) {
#ifdef WITH_TBB
			if(stateful) {
				tbb_grp.run(exactSearchWorkerStateful(i));
			} else {
				tbb_grp.run(exactSearchWorker(i));
			}
		}
		tbb_grp.wait();
#else
			tids[i] = i;
			if(stateful) {
                                threads[i] = new tthread::thread(exactSearchWorkerStateful, (void*)&tids[i]);
			} else {
                                threads[i] = new tthread::thread(exactSearchWorker, (void*)&tids[i]);
			}
		}

		for(int i = 1; i <= nthreads; i++)
                    threads[i]->join();
#endif
	}
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

/**
 * A statefulness-aware worker driver.  Uses Unpaired/Paired1mmAlignerV1.
 */
#ifdef WITH_TBB
void mismatchSearchWorkerFullStateful::operator()() {
#else
static void mismatchSearchWorkerFullStateful(void *vp) {
	int tid = *((int*)vp);
#endif
	PairedPatternSource&   _patsrc = *mismatchSearch_patsrc;
	HitSink&               _sink   = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw  = *mismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw  = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os      = *mismatchSearch_os;
	BitPairReference*      refs    =  mismatchSearch_refs;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
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
			refs,
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
			color,
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
	return;
}
#ifdef WITH_TBB
void mismatchSearchWorkerFull::operator()(){
#else
static void mismatchSearchWorkerFull(void *vp){
	int tid = *((int*)vp);
#endif
	PairedPatternSource&   _patsrc   = *mismatchSearch_patsrc;
	HitSink&               _sink     = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw    = *mismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw    = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os        = *mismatchSearch_os;
	const BitPairReference* refs     =  mismatchSearch_refs;

	// Per-thread initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread* patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
	HitSinkPerThread* sink = sinkFact->create();
	EbwtSearchParams<String<Dna> > params(
	        *sink,      // HitSinkPerThread
	        os,         // reference sequences
	        true,       // read is forward
	        false);     // index is mirror index
	GreedyDFSRangeSource bt(
	        &ebwtFw, params,
	        refs,           // reference sequence (for colorspace)
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
		ebwtFw.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose);
	}
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose);
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	if(color || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, color, sanityCheck, NULL, &os, false, true, useMm, useShmem, mmSweep, verbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	mismatchSearch_refs = refs;

#ifdef WITH_TBB
	tbb::task_group tbb_grp;
#else
	AutoArray<tthread::thread*> threads(nthreads+1);
	AutoArray<int> tids(nthreads+1);
#endif

    CHUD_START();
    {
		Timer _t(cerr, "Time for 1-mismatch full-index search: ", timing);

		for(int i = 1; i <= nthreads; i++) {
#ifdef WITH_TBB
			if(stateful) {
				tbb_grp.run(mismatchSearchWorkerFullStateful(i));
			} else {
				tbb_grp.run(mismatchSearchWorkerFull(i));
			}
		}
		tbb_grp.wait();
#else
			tids[i] = i;
			if(stateful) {
                                threads[i] = new tthread::thread(mismatchSearchWorkerFullStateful, (void*)&tids[i]);
			} else {
                                threads[i] = new tthread::thread(mismatchSearchWorkerFull, (void*)&tids[i]);
			}
		}

		for(int i = 1; i <= nthreads; i++)
                    threads[i]->join();
#endif
    }
	if(refs != NULL) delete refs;
}

#define SWITCH_TO_FW_INDEX() { \
	/* Evict the mirror index from memory if necessary */ \
	if(ebwtBw.isInMemory()) ebwtBw.evictFromMemory(); \
	assert(!ebwtBw.isInMemory()); \
	/* Load the forward index into memory if necessary */ \
	if(!ebwtFw.isInMemory()) { \
		Timer _t(cerr, "Time loading forward index: ", timing); \
		ebwtFw.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose); \
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
		ebwtBw.loadIntoMemory(color ? 1 : 0, !noRefNames, startVerbose); \
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


/**
 * A statefulness-aware worker driver.  Uses UnpairedExactAlignerV1.
 */
#ifdef WITH_TBB
void twoOrThreeMismatchSearchWorkerStateful::operator()(){
#else
static void twoOrThreeMismatchSearchWorkerStateful(void *vp) {
	int tid = *((int*)vp);
#endif
	PairedPatternSource&   _patsrc = *twoOrThreeMismatchSearch_patsrc;
	HitSink&               _sink   = *twoOrThreeMismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw  = *twoOrThreeMismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw  = *twoOrThreeMismatchSearch_ebwtBw;
	vector<String<Dna5> >& os      = *twoOrThreeMismatchSearch_os;
	BitPairReference*      refs    =  twoOrThreeMismatchSearch_refs;
	static bool            two     =  twoOrThreeMismatchSearch_two;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);

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
			refs,
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
			color,
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
	return;
}
#ifdef WITH_TBB
void twoOrThreeMismatchSearchWorkerFull::operator()(){
#else
static void twoOrThreeMismatchSearchWorkerFull(void *vp) {
	int tid = *((int*)vp);
#endif
	PairedPatternSource&           _patsrc  = *twoOrThreeMismatchSearch_patsrc;
	HitSink&                       _sink    = *twoOrThreeMismatchSearch_sink;
	vector<String<Dna5> >&         os       = *twoOrThreeMismatchSearch_os;
	bool                           two      = twoOrThreeMismatchSearch_two;
    PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread* patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
	HitSinkPerThread* sink = sinkFact->create();
	/* Per-thread initialization */
	EbwtSearchParams<String<Dna> > params(
			*sink,       /* HitSink */
	        os,          /* reference sequences */
	        true,        /* read is forward */
	        true);       /* index is forward */
	Ebwt<String<Dna> >& ebwtFw = *twoOrThreeMismatchSearch_ebwtFw;
	Ebwt<String<Dna> >& ebwtBw = *twoOrThreeMismatchSearch_ebwtBw;
	const BitPairReference* refs = twoOrThreeMismatchSearch_refs;
	GreedyDFSRangeSource btr1(
	        &ebwtFw, params,
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	assert(!ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose);
	}
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose);
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	if(color || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, color, sanityCheck, NULL, &os, false, true, useMm, useShmem, mmSweep, verbose, startVerbose);
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

#ifdef WITH_TBB
	tbb::task_group tbb_grp;
#else
	AutoArray<tthread::thread*> threads(nthreads+1);
	AutoArray<int> tids(nthreads+1);
#endif

        CHUD_START();
    {
		Timer _t(cerr, "End-to-end 2/3-mismatch full-index search: ", timing);
		for(int i = 1; i <= nthreads; i++) {
#ifdef WITH_TBB
			if(stateful) {
				tbb_grp.run(twoOrThreeMismatchSearchWorkerStateful(i));
			} else {
				tbb_grp.run(twoOrThreeMismatchSearchWorkerFull(i));
			}
		}
		tbb_grp.wait();
#else
			tids[i] = i;
			if(stateful) {
                                threads[i] = new tthread::thread(twoOrThreeMismatchSearchWorkerStateful, (void*)&tids[i]);
			} else {
                                threads[i] = new tthread::thread(twoOrThreeMismatchSearchWorkerFull, (void*)&tids[i]);
			}
		}

		for(int i = 1; i <= nthreads; i++)
                    threads[i]->join();
#endif
    }
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

#ifdef WITH_TBB
void seededQualSearchWorkerFull::operator()(){
#else
static void seededQualSearchWorkerFull(void *vp) {
	int tid = *((int*)vp);
#endif
	PairedPatternSource&     _patsrc    = *seededQualSearch_patsrc;
	HitSink&                 _sink      = *seededQualSearch_sink;
	vector<String<Dna5> >&   os         = *seededQualSearch_os;
	int                      qualCutoff = seededQualSearch_qualCutoff;
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	PatternSourcePerThread* patsrc = patsrcFact->create();
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
	HitSinkPerThread* sink = sinkFact->create();
	/* Per-thread initialization */
	EbwtSearchParams<String<Dna> > params(
	        *sink,       /* HitSink */
	        os,          /* reference sequences */
	        true,        /* read is forward */
	        true);       /* index is forward */
	Ebwt<String<Dna> >& ebwtFw = *seededQualSearch_ebwtFw;
	Ebwt<String<Dna> >& ebwtBw = *seededQualSearch_ebwtBw;
	PartialAlignmentManager * pamRc = NULL;
	PartialAlignmentManager * pamFw = NULL;
	if(seedMms > 0) {
		pamRc = new PartialAlignmentManager(64);
		pamFw = new PartialAlignmentManager(64);
	}
	vector<PartialAlignment> pals;
	const BitPairReference* refs = seededQualSearch_refs;
	// GreedyDFSRangeSource for finding exact hits for the forward-
	// oriented read
	GreedyDFSRangeSource btf1(
	        &ebwtFw, params,
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
	        refs,           // reference sequence (for colorspace)
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
		uint32_t plen = (uint32_t)length(patFw);
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
#ifdef WITH_TBB
void seededQualSearchWorkerFullStateful::operator()(){
#else
static void seededQualSearchWorkerFullStateful(void *vp) {
	int tid = *((int*)vp);
#endif
	PairedPatternSource&     _patsrc    = *seededQualSearch_patsrc;
	HitSink&                 _sink      = *seededQualSearch_sink;
	Ebwt<String<Dna> >&      ebwtFw     = *seededQualSearch_ebwtFw;
	Ebwt<String<Dna> >&      ebwtBw     = *seededQualSearch_ebwtBw;
	vector<String<Dna5> >&   os         = *seededQualSearch_os;
	int                      qualCutoff = seededQualSearch_qualCutoff;
	BitPairReference*        refs       = seededQualSearch_refs;

	// Global initialization
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, tid);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink);
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
			refs,
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
			color,
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
	return;
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
	bool pair = mates1.size() > 0 || mates12.size() > 0;
	if(color || (pair && mixedThresh < 0xffffffff)) {
		Timer _t(cerr, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, color, sanityCheck, NULL, &os, false, true, useMm, useShmem, mmSweep, verbose, startVerbose);
		if(!refs->loaded()) throw 1;
	}
	seededQualSearch_refs = refs;

#ifdef WITH_TBB
	tbb::task_group tbb_grp;
#else
	AutoArray<tthread::thread*> threads(nthreads+1);
	AutoArray<int> tids(nthreads+1);
#endif

	SWITCH_TO_FW_INDEX();
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cerr, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose);
	}
	CHUD_START();
	{
		// Phase 1: Consider cases 1R and 2R
		Timer _t(cerr, "Seeded quality full-index search: ", timing);

		for(int i = 1; i <= nthreads; i++) {
#ifdef WITH_TBB
			if(stateful) {
				tbb_grp.run(seededQualSearchWorkerFullStateful(i));
			} else {
				tbb_grp.run(seededQualSearchWorkerFull(i));
			}
		}
		tbb_grp.wait();
#else
			tids[i] = i;
			if(stateful) {
                                threads[i] = new tthread::thread(seededQualSearchWorkerFullStateful, (void*)&tids[i]);
			} else {
                                threads[i] = new tthread::thread(seededQualSearchWorkerFull, (void*)&tids[i]);
			}
		}

		for(int i = 1; i <= nthreads; i++)
                    threads[i]->join();
#endif
	}
	if(refs != NULL) {
		delete refs;
	}
	ebwtBw.evictFromMemory();
}

/**
 * Return a new dynamically allocated PatternSource for the given
 * format, using the given list of strings as the filenames to read
 * from or as the sequences themselves (i.e. if -c was used).
 */
static PatternSource*
patsrcFromStrings(int format,
                  const vector<string>& reads,
                  const vector<string>* quals)
{
	switch(format) {
		case FASTA:
			return new FastaPatternSource (seed, reads, quals, color,
			                               randomizeQuals,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               solexaQuals, phred64Quals,
			                               integerQuals,
			                               skipReads);
		case FASTA_CONT:
			return new FastaContinuousPatternSource (
			                               seed, reads, fastaContLen,
			                               fastaContFreq,
			                               patDumpfile, verbose,
			                               skipReads);
		case RAW:
			return new RawPatternSource   (seed, reads, color,
			                               randomizeQuals,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               skipReads);
		case FASTQ:
			return new FastqPatternSource (seed, reads, color,
			                               randomizeQuals,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               solexaQuals, phred64Quals,
			                               integerQuals, fuzzy,
			                               skipReads);
		case TAB_MATE:
			return new TabbedPatternSource(seed, reads, color,
			                               randomizeQuals,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               skipReads);
		case CMDLINE:
			return new VectorPatternSource(seed, reads, color,
			                               randomizeQuals,
			                               patDumpfile, verbose,
			                               trim3, trim5,
			                               skipReads);
		case RANDOM:
			return new RandomPatternSource(seed, 2000000, lenRandomReads,
			                               patDumpfile,
			                               verbose);
		default: {
			cerr << "Internal error; bad patsrc format: " << format << endl;
			throw 1;
		}
	}
}

#define PASS_DUMP_FILES dumpAlBase, dumpUnalBase, dumpMaxBase

static string argstr;

template<typename TStr>
static void driver(const char * type,
                   const string& ebwtFileBase,
                   const string& query,
                   const vector<string>& queries,
                   const vector<string>& qualities,
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
		const vector<string>* qs = &mates12;
		vector<string> tmp;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmp;
			tmp.push_back(mates12[i]);
			assert_eq(1, tmp.size());
		}
		patsrcs_ab.push_back(patsrcFromStrings(format, *qs, NULL));
		if(!fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < mates1.size(); i++) {
		const vector<string>* qs = &mates1;
		const vector<string>* quals = &qualities1;
		vector<string> tmpSeq;
		vector<string> tmpQual;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(mates1[i]);
			quals = &tmpSeq;
			tmpQual.push_back(qualities1[i]);
			assert_eq(1, tmpSeq.size());
		}
		if(quals->empty()) quals = NULL;
		patsrcs_a.push_back(patsrcFromStrings(format, *qs, quals));
		if(!fileParallel) {
			break;
		}
	}

	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < mates2.size(); i++) {
		const vector<string>* qs = &mates2;
		const vector<string>* quals = &qualities2;
		vector<string> tmpSeq;
		vector<string> tmpQual;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(mates2[i]);
			quals = &tmpQual;
			tmpQual.push_back(qualities2[i]);
			assert_eq(1, tmpSeq.size());
		}
		if(quals->empty()) quals = NULL;
		patsrcs_b.push_back(patsrcFromStrings(format, *qs, quals));
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
		const vector<string>* qs = &queries;
		const vector<string>* quals = &qualities;
		PatternSource* patsrc = NULL;
		vector<string> tmpSeq;
		vector<string> tmpQual;
		if(fileParallel) {
			// Feed query files one to each PatternSource
			qs = &tmpSeq;
			tmpSeq.push_back(queries[i]);
			quals = &tmpQual;
			tmpQual.push_back(qualities[i]);
			assert_eq(1, tmpSeq.size());
		}
		if(quals->empty()) quals = NULL;
		patsrc = patsrcFromStrings(format, *qs, quals);
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
		patsrc = new PairedSoloPatternSource(patsrcs_ab, seed);
	} else {
		patsrc = new PairedDualPatternSource(patsrcs_a, patsrcs_b, seed);
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
			fout = new OutFileBuf(outfile.c_str(), false);
		}
	} else {
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
	                color,  // index is colorspace
	                -1,     // don't care about entireReverse
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
		ebwtBw = new Ebwt<TStr>(
			adjustedEbwtFileBase + ".rev",
			color,  // index is colorspace
			-1,     // don't care about entireReverse
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
	if(!os.empty()) {
		for(size_t i = 0; i < os.size(); i++) {
			size_t olen = seqan::length(os[i]);
			int longestStretch = 0;
			int curStretch = 0;
			for(size_t j = 0; j < olen; j++) {
				if((int)os[i][j] < 4) {
					curStretch++;
					if(curStretch > longestStretch) longestStretch = curStretch;
				} else {
					curStretch = 0;
				}
			}
			if(longestStretch < (color ? 2 : 1)) {
				os.erase(os.begin() + i);
				i--;
			}
		}
	}
	if(sanityCheck && !os.empty()) {
		// Sanity check number of patterns and pattern lengths in Ebwt
		// against original strings
		assert_eq(os.size(), ebwt.nPat());
		for(size_t i = 0; i < os.size(); i++) {
			assert_eq(length(os[i]), ebwt.plen()[i] + (color ? 1 : 0));
		}
		ebwt.loadIntoMemory(color ? 1 : 0, -1, !noRefNames, startVerbose);
		ebwt.checkOrigs(os, color, false);
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
							ebwt.nPat(), offBase,
							colorSeq, colorQual, printCost,
							suppressOuts, rmap, amap,
							fullRef, PASS_DUMP_FILES,
							format == TAB_MATE, sampleMax,
							table, refnames, partitionSz);
				} else {
					sink = new VerboseHitSink(
							fout, offBase,
							colorSeq, colorQual, printCost,
							suppressOuts, rmap, amap,
							fullRef, PASS_DUMP_FILES,
							format == TAB_MATE, sampleMax,
							table, refnames, partitionSz);
				}
				break;
			case OUTPUT_SAM:
				if(refOut) {
					throw 1;
				} else {
					SAMHitSink *sam = new SAMHitSink(
							fout, 1, rmap, amap,
							fullRef, samNoQnameTrunc, defaultMapq,
							PASS_DUMP_FILES,
							format == TAB_MATE, sampleMax,
							table, refnames);
					if(!samNoHead) {
						vector<string> refnames;
						if(!samNoSQ) {
							readEbwtRefnames(adjustedEbwtFileBase, refnames);
						}
						sam->appendHeaders(
								sam->out(0), ebwt.nPat(),
								refnames, color, samNoSQ, rmap,
								ebwt.plen(), fullRef,
								samNoQnameTrunc,
								argstr.c_str(),
								rgs.empty() ? NULL : rgs.c_str());
					}
					sink = sam;
				}
				break;
			case OUTPUT_CONCISE:
				if(refOut) {
					sink = new ConciseHitSink(
							ebwt.nPat(), offBase,
							PASS_DUMP_FILES,
							format == TAB_MATE,  sampleMax,
							table, refnames, reportOpps);
				} else {
					sink = new ConciseHitSink(
							fout, offBase,
							PASS_DUMP_FILES,
							format == TAB_MATE,  sampleMax,
							table, refnames, reportOpps);
				}
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
			seededQualCutoffSearchFull(seedLen,
									   qualThresh,
									   seedMms,
									   *patsrc,
									   *sink,
									   ebwt,    // forward index
									   *ebwtBw, // mirror index (not optional)
									   os);     // references, if available
		}
		else if(mismatches > 0) {
			if(mismatches == 1) {
				assert(ebwtBw != NULL);
				mismatchSearchFull(*patsrc, *sink, ebwt, *ebwtBw, os);
			} else if(mismatches == 2 || mismatches == 3) {
				twoOrThreeMismatchSearchFull(*patsrc, *sink, ebwt, *ebwtBw, os, mismatches == 2);
			} else {
				cerr << "Error: " << mismatches << " is not a supported number of mismatches" << endl;
				throw 1;
			}
		} else {
			// Search without mismatches
			// Note that --fast doesn't make a difference here because
			// we're only loading half of the index anyway
			exactSearch(*patsrc, *sink, ebwt, os);
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
				cout << "Quality inputs:" << endl;
				for(size_t i = 0; i < qualities.size(); i++) {
					cout << "  " << qualities[i] << endl;
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
			driver<String<Dna, Alloc<> > >("DNA", ebwtFile, query, queries, qualities, outfile);
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
