#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>
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
#include "aligner.h"
#include "aligner_0mm.h"
#include "aligner_1mm.h"
#include "aligner_23mm.h"
#ifdef CHUD_PROFILING
#include <CHUD/CHUD.h>
#endif

using namespace std;
using namespace seqan;

static vector<string> mates1; // mated reads (first mate)
static vector<string> mates2; // mated reads (second mate)
static string adjustedEbwtFileBase = "";
static int verbose				= 0; // be talkative
static bool quiet				= false; // print nothing but the alignments
static int sanityCheck			= 0;  // enable expensive sanity checks
static int format				= FASTQ; // default read format is FASTQ
static string origString		= ""; // reference text, or filename(s)
static int revcomp				= 1; // search for reverse complements?
static int seed					= 0; // srandom() seed
static int timing				= 0; // whether to report basic timing data
static bool allHits				= false; // for multihits, report just one
static bool rangeMode			= false; // report BWT ranges instead of ref locs
static int showVersion			= 0; // just print version and quit?
static int ipause				= 0; // pause before maching?
static uint32_t qUpto			= 0xffffffff; // max # of queries to read
static int skipSearch			= 0; // abort before searching
static int qSameLen				= 0; // abort before searching
static int trim5				= 0; // amount to trim from 5' end
static int trim3				= 0; // amount to trim from 3' end
static int reportOpps			= 0; // whether to report # of other mappings
static int offRate				= -1; // keep default offRate
static int isaRate				= -1; // keep default isaRate
static int mismatches			= 0; // allow 0 mismatches by default
static char *patDumpfile		= NULL; // filename to dump patterns to
static bool solexa_quals		= false; //quality strings are solexa qualities, instead of phred
static bool integer_quals		= false; //quality strings are space-separated strings of integers, instead of ASCII
static int maqLike				= 1; // do maq-like searching
static int seedLen              = 28; // seed length (changed in Maq 0.6.4 from 24)
static int seedMms              = 2;  // # mismatches allowed in seed (maq's -n)
static int qualThresh           = 70; // max qual-weighted hamming dist (maq's -e)
static int maxBts               = 125; // max # backtracks allowed in half-and-half mode
static int maxBts0              = 28; // max # backtracks allowed in half-and-half mode
static int maxBts1              = 16; // max # backtracks allowed in half-and-half mode
static int maxBts2              = 12; // max # backtracks allowed in half-and-half mode
static int nsPolicy             = NS_TO_NS; // policy for handling no-confidence bases
static int nthreads             = 1;     // number of pthreads operating concurrently
static output_types outType		= FULL;  // style of output
static bool randReadsNoSync     = false; // true -> generate reads from per-thread random source
static int numRandomReads       = 50000000; // # random reads (see Random*PatternSource in pat.h)
static int lenRandomReads       = 35;    // len of random reads (see Random*PatternSource in pat.h)
static bool fullIndex           = true;  // load halves one at a time and proceed in phases
static bool noRefNames          = false; // true -> print reference indexes; not names
static ofstream *dumpNoHits     = NULL;  // file to dump non-hitting reads to (for performance study)
static ofstream *dumpHHHits     = NULL;  // file to dump half-and-half hits to (for performance study)
static string dumpUnalFaBase    = "";    // basename of FASTA files to dump unaligned reads to
static string dumpUnalFqBase    = "";    // basename of FASTQ files to dump unaligned reads to
static string dumpMaxFaBase     = "";    // basename of FASTA files to dump reads with more than -m valid alignments to
static string dumpMaxFqBase     = "";    // basename of FASTQ files to dump reads with more than -m valid alignments to
static uint32_t khits           = 1;     // number of hits per read; >1 is much slower
static uint32_t mhits           = 0xffffffff; // don't report any hits if there are > mhits
static bool onlyBest			= false; // true -> guarantee alignments from best possible stratum
static bool spanStrata			= false; // true -> don't stop at stratum boundaries
static bool refOut				= false; // if true, alignments go to per-ref files
static bool seedAndExtend		= false; // use seed-and-extend aligner; for metagenomics recruitment
static int partitionSz          = 0;     // output a partitioning key in first field
static bool noMaqRound          = false; // true -> don't round quals to nearest 10 like maq
static bool forgiveInput        = false; // let read input be a little wrong w/o complaining or dying
static bool useSpinlock         = true;  // false -> don't use of spinlocks even if they're #defines
static bool fileParallel        = false; // separate threads read separate input files in parallel
static bool useShmem            = false; // use shared memory to hold _ebwt[] and _offs[] arrays?
static bool stateful            = false; // use stateful aligners
static uint32_t prefetchWidth   = 1;     // number of reads to process in parallel w/ --stateful
static uint32_t minInsert       = 0;     // minimum insert size (Maq = 0, SOAP = 400)
static uint32_t maxInsert       = 250;   // maximum insert size (Maq = 250, SOAP = 600)
static bool mate1fw             = true;  // -1 mate aligns in fw orientation on fw strand
static bool mate2fw             = false; // -2 mate aligns in rc orientation on fw strand
static uint32_t mixedThresh     = 4;     // threshold for when to switch to paired-end mixed mode (see aligner.h)
static uint32_t mixedAttemptLim = 100;   // number of attempts to make in "mixed mode" before giving up on orientation
static bool dontReconcileMates  = true;  // suppress pairwise all-versus-all way of resolving mates
static uint32_t cacheLimit      = 5;     // ranges w/ size > limit will be cached
static uint32_t cacheSize       = 0;     // # words per range cache
static int offBase              = 0;     // offsets are 0-based by default, but configurable
static bool tryHard             = false; // set very high maxBts, mixedAttemptLim
static uint32_t skipReads       = 0;     // # reads/read pairs to skip
static bool nofw                = false; // don't align fw orientation of read
static bool norc                = false; // don't align rc orientation of read
// mating constraints

static const char *short_options = "fqbzh?cu:rv:s:at3:5:o:e:n:l:w:p:k:m:1:2:I:X:x:B:y";

enum {
	ARG_ORIG = 256,
	ARG_SEED,
	ARG_DUMP_PATS,
	ARG_RANGE,
	ARG_CONCISE,
	ARG_SOLEXA_QUALS,
	ARG_MAXBTS,
	ARG_MAXBTS0,
	ARG_MAXBTS1,
	ARG_MAXBTS2,
	ARG_VERBOSE,
	ARG_QUIET,
	ARG_RANDOM_READS,
	ARG_RANDOM_READS_NOSYNC,
	ARG_NOOUT,
	ARG_FAST,
	ARG_UNFA,
	ARG_UNFQ,
	ARG_MAXFA,
	ARG_MAXFQ,
	ARG_REFIDX,
	ARG_DUMP_NOHIT,
	ARG_DUMP_HHHIT,
	ARG_SANITY,
	ARG_BEST,
	ARG_SPANSTRATA,
	ARG_REFOUT,
	ARG_ISARATE,
	ARG_SEED_EXTEND,
	ARG_PARTITION,
	ARG_INTEGER_QUALS,
	ARG_FORGIVE_INPUT,
	ARG_NOMAQROUND,
	ARG_USE_SPINLOCK,
	ARG_FILEPAR,
	ARG_SHARED_MEM,
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
	ARG_SKIP
};

static struct option long_options[] = {
	{"verbose",      no_argument,       0,            ARG_VERBOSE},
	{"quiet",        no_argument,       0,            ARG_QUIET},
	{"sanity",       no_argument,       0,            ARG_SANITY},
	{"exact",        no_argument,       0,            '0'},
	{"1mm",          no_argument,       0,            '1'},
	{"2mm",          no_argument,       0,            '2'},
	{"pause",        no_argument,       &ipause,      1},
	{"orig",         required_argument, 0,            ARG_ORIG},
	{"all",          no_argument,       0,            'a'},
	{"concise",      no_argument,       0,            ARG_CONCISE},
	{"binout",       no_argument,       0,            'b'},
	{"noout",        no_argument,       0,            ARG_NOOUT},
	{"solexa-quals", no_argument,       0,            ARG_SOLEXA_QUALS},
	{"integer-quals",no_argument,       0,            ARG_INTEGER_QUALS},
	{"time",         no_argument,       0,            't'},
	{"trim3",        required_argument, 0,            '3'},
	{"trim5",        required_argument, 0,            '5'},
	{"seed",         required_argument, 0,            ARG_SEED},
	{"qupto",        required_argument, 0,            'u'},
	{"unfa",         required_argument, 0,            ARG_UNFA},
	{"unfq",         required_argument, 0,            ARG_UNFQ},
	{"maxfa",        required_argument, 0,            ARG_MAXFA},
	{"maxfq",        required_argument, 0,            ARG_MAXFQ},
	{"offrate",      required_argument, 0,            'o'},
	{"isarate",      required_argument, 0,            ARG_ISARATE},
	{"skipsearch",   no_argument,       &skipSearch,  1},
	{"qsamelen",     no_argument,       &qSameLen,    1},
	{"reportopps",   no_argument,       &reportOpps,  1},
	{"version",      no_argument,       &showVersion, 1},
	{"ntoa",         no_argument,       &nsPolicy,    NS_TO_AS},
	{"dumppats",     required_argument, 0,            ARG_DUMP_PATS},
	{"revcomp",      no_argument,       0,            'r'},
	{"maqerr",       required_argument, 0,            'e'},
	{"seedlen",      required_argument, 0,            'l'},
	{"seedmms",      required_argument, 0,            'n'},
	{"filepar",      no_argument,       0,            ARG_FILEPAR},
	{"help",         no_argument,       0,            'h'},
	{"threads",      required_argument, 0,            'p'},
	{"khits",        required_argument, 0,            'k'},
	{"mhits",        required_argument, 0,            'm'},
	{"minins",       required_argument, 0,            'I'},
	{"maxins",       required_argument, 0,            'X'},
	{"best",         no_argument,       0,            ARG_BEST},
	{"nostrata",     no_argument,       0,            ARG_SPANSTRATA},
	{"nomaqround",   no_argument,       0,            ARG_NOMAQROUND},
	{"refidx",       no_argument,       0,            ARG_REFIDX},
	{"range",        no_argument,       0,            ARG_RANGE},
	{"maxbts",       required_argument, 0,            ARG_MAXBTS},
	{"maxbts0",      required_argument, 0,            ARG_MAXBTS0},
	{"maxbts1",      required_argument, 0,            ARG_MAXBTS1},
	{"maxbts2",      required_argument, 0,            ARG_MAXBTS2},
	{"randread",     no_argument,       0,            ARG_RANDOM_READS},
	{"randreadnosync", no_argument,     0,            ARG_RANDOM_READS_NOSYNC},
	{"phased",       no_argument,       0,            'z'},
	{"dumpnohit",    no_argument,       0,            ARG_DUMP_NOHIT},
	{"dumphhhit",    no_argument,       0,            ARG_DUMP_HHHIT},
	{"refout",       no_argument,       0,            ARG_REFOUT},
	{"seedextend",   no_argument,       0,            ARG_SEED_EXTEND},
	{"partition",    required_argument, 0,            ARG_PARTITION},
	{"forgive",      no_argument,       0,            ARG_FORGIVE_INPUT},
	{"nospin",       no_argument,       0,            ARG_USE_SPINLOCK},
	{"sharedmem",    no_argument,       0,            ARG_SHARED_MEM},
	{"stateful",     no_argument,       0,            ARG_STATEFUL},
	{"prewidth",     required_argument, 0,            ARG_PREFETCH_WIDTH},
	{"ff",           no_argument,       0,            ARG_FF},
	{"fr",           no_argument,       0,            ARG_FR},
	{"rf",           no_argument,       0,            ARG_RF},
	{"mixthresh",    required_argument, 0,            'x'},
	{"pairtries",    required_argument, 0,            ARG_MIXED_ATTEMPTS},
	{"noreconcile",  no_argument,       0,            ARG_NO_RECONCILE},
	{"cachelim",     required_argument, 0,            ARG_CACHE_LIM},
	{"cachesz",      required_argument, 0,            ARG_CACHE_SZ},
	{"nofw",         no_argument,       0,            ARG_NO_FW},
	{"norc",         no_argument,       0,            ARG_NO_RC},
	{"offbase",      required_argument, 0,            'B'},
	{"tryhard",      no_argument,       0,            'y'},
	{"skip",         required_argument, 0,            's'},
	{0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: bowtie [options]* <ebwt> {-1 <mates1> -2 <mates2> | <singles>} [<hits>]" << endl
	    << "  <ebwt>    base filename for index files (minus trailing .1.ebwt/.2.ebwt/etc.)" << endl
	    << "  <mates1>  comma-separated list of files containing upstream mates (or the" << endl
	    << "            sequences themselves, if -c is set) paired with mates in <mates2>" << endl
	    << "  <mates2>  comma-separated list of files containing downstream mates (or" << endl
	    << "            sequences themselves if -c is set) paired with mates in <mates1>" << endl
	    << "  <singles> comma-separated list of files containing unpaired reads, or the" << endl
	    << "            sequences themselves, if -c is set" << endl
	    << "  <hits>    file to write hits to (default: stdout)" << endl
	    << "Input:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    << "  -r                 query input files are raw one-sequence-per-line" << endl
	    << "  -c                 query sequences given on cmd line (as <mates>, <singles>)" << endl
	    << "  -s/--skip <int>    skip the first <int> reads/pairs in the input" << endl
	    << "  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads" << endl
		<< "  --solexa-quals     convert quals from solexa (can be < 0) to phred (can't)" << endl
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
	    << "  --maxbts <int>     max number of backtracks allowed for -n 2/3 (default: 125)" << endl
	    << "  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)" << endl
	    << "  -y/--tryhard       try hard to find valid alignments, at the expense of speed" << endl
	    << "Reporting:" << endl
	    << "  -k <int>           report up to <int> good alignments per read (default: 1)" << endl
	    << "  -a/--all           report all alignments per read (much slower than low -k)" << endl
	    << "  -m <int>           suppress all alignments if > <int> exist (def.: no limit)" << endl
	    << "  --best             guarantee reported alignments are at best possible stratum" << endl
	    << "  --nostrata         if reporting >1 alignment, don't quit at stratum boundaries" << endl
	    << "Output:" << endl
	    << "  --concise          write hits in concise format" << endl
	    << "  -b/--binout        write hits in binary format (<hits> argument not optional)" << endl
	    << "  -t/--time          print wall-clock time taken by search phases" << endl
	    << "  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)" << endl
	    << "  --quiet            print nothing but the alignments" << endl
	    << "  --refout           write alignments to files refXXXXX.map, 1 map per reference" << endl
	    << "  --refidx           refer to ref. seqs by 0-based index rather than name" << endl
	    << "  --unfa <fname>     write unaligned reads to FASTA file(s) <fname>*.fa" << endl
	    << "  --unfq <fname>     write unaligned reads to FASTQ file(s) <fname>*.fq" << endl
	    << "  --maxfa <fname>    write reads exceeding -m limit to FASTA file(s) <fname>*.fa" << endl
	    << "  --maxfq <fname>    write reads exceeding -m limit to FASTQ file(s) <fname>*.fq" << endl
	    << "Performance:" << endl
#ifdef BOWTIE_PTHREADS
	    << "  -p/--threads <int> number of alignment threads to launch (default: 1)" << endl
#endif
		<< "  -z/--phased        alternate between index halves; slower, but uses 1/2 mem" << endl
	    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
	    << "Other:" << endl
	    << "  --seed <int>       seed for random number generator" << endl
	    << "  --verbose          verbose output (for debugging)" << endl
	    << "  --version          print version information and quit" << endl
	    << "  -h/--help          print detailed description of tool and its options" << endl
	    ;
}

/**
 * Print a detailed usage message to the provided output stream.
 *
 * Manual text converted to C++ string with something like:
 * cat MANUAL  | head -415 | tail -340 | sed -e 's/\"/\\\"/g' | \
 *   sed -e 's/^/"/' | sed -e 's/$/\\n"/'
 */
static void printLongUsage(ostream& out) {
	out <<
	"\n"
	" Using the 'bowtie' Aligner\n"
	" --------------------------\n"
	"\n"
	" The 'bowtie' aligner takes an index and a set of reads as input and\n"
	" outputs a list of alignments.  Alignments are selected according to a\n"
	" combination of the -v/-n/-e/-l options, which define which alignments\n"
	" are legal, and the -k/-a/-m/--best/--nostrata options which define\n"
	" which and how many legal alignments should be reported.\n"
	"\n"
	" By default, Bowtie enforces an alignment policy equivalent to Maq's\n"
	" quality-aware policy (http://maq.sf.net) (-n 2 -l 28 -e 70), but it\n"
	" can also be made to enforce an end-to-end k-difference policy\n"
	" equivalent to SOAP's (http://soap.genomics.org.cn/) (-v 2).\n"
	"\n"
	" The process by which bowtie chooses an alignment to report is\n"
	" randomized in order to avoid \"mapping bias\" - the phenomenon whereby\n"
	" an aligner systematically fails to report a particular class of good\n"
	" alignments, causing spurious \"holes\" in the comparative assembly.\n"
	" Whenever bowtie reports a subset of the valid alignments that exist,\n"
	" it makes an effort to sample them randomly.  This randomness flows\n"
	" from a simple seeded pseudo-random number generator and is\n"
	" \"deterministic\" in the sense that Bowtie will always produce the same\n"
	" results when run with the same initial \"seed\" value (see documentation\n"
	" for --seed option).\n"
	"\n"
	" Gapped alignments are not currently supported, but we do plan to\n"
	" implement this in the future.  Alignment in ABI \"color space\" is also\n"
	" not currently supported.\n"
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
	"  -e options.\n"
	" \n"
	"  If there are many possible alignments that satisfy both criteria,\n"
	"  Bowtie will make an effort to give preference to alignments with\n"
	"  where the sum from criterion 2 is smaller.  Bowtie does not guarantee\n"
	"  that it will report the minimum-sum alignment.\n"
	"  \n"
	"  Note that Maq internally rounds base qualities to the nearest 10 and\n"
	"  truncates qualities greater than 30 to 30.  To maintain compatibility\n"
	"  with Maq, Bowtie does the same.\n"
	" \n"
	"  Also note that bowtie is not fully sensitive in -n 2 and -n 3\n"
	"  modes.  In those modes, bowtie imposes a \"backtracking limit\" to\n"
	"  limit the amount of effort spent trying to find valid alignments for\n"
	"  low-quality reads that are unlikely have any.  Since the limit is\n"
	"  arbitrary, it may cause bowtie to miss some legal 2- and 3-mismatch\n"
	"  alignments.  We have set the limit to what we consider a reasonable\n"
	"  default (125), but the user may decrease or increase the limit using\n"
	"  the --maxbts option.  Setting the limit to a very large number\n"
	"  (>10000) guarantees full sensitivity.\n"
	" \n"
	"  End-to-end k-difference Policy\n"
	"  ------------------------------\n"
	"  \n"
	"  The policy has one criterion: Alignments may have no more than V\n"
	"  mismatches.  Quality values are ignored.  The number of mismatches\n"
	"  permitted is configurable with the -v option.\n"
	"  \n"
	"  Paired-end Alignment\n"
	"  --------------------\n"
	"  \n"
	"  As of version 0.9.9, Bowtie contains preliminary support for\n"
	"  paired-end alignment.  Support is described as \"preliminary\" because\n"
	"  paired-end alignment only works under -v 0 or -v 1 alignment policies\n"
	"  and because performance has not been thoroughly optimized.  Also,\n"
	"  there is currently no way to get Bowtie to report a mix of paired-end\n"
	"  alignments and singleton alignments (for unpaired mates) in one run;\n"
	"  if a mix of alignments is desired, the user must first perform a\n"
	"  paired-end run using the --unfa or --unfq option to write unaligned\n"
	"  pairs to files, then run Bowtie again in single-ended mode using\n"
	"  those files as input.  Future versions of Bowtie will expand paired-\n"
	"  end support and, hopefully, improve performance as well.\n"
	"  \n"
	"  A valid paired-end alignment satisfies the following criteria:\n"
	"  \n"
	"  1. Both mates have a valid alignment according to the alignment\n"
	"     policy specified by the -v/-n/-e/-l options.\n"
	"  2. The relative orientation and position of the mates satisfy the\n"
	"     constraints given by the -I/-X/--lr/--rl/--ll options. \n"
	"  \n"
	"  Policies governing which paired-end alignments are reported for a\n"
	"  given read are specified using the -k, -a and -m options as usual.\n"
	"  The --nostrata and --best options do not apply in paired-end mode.\n"
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
	"  or -y options to increase Bowtie's sensitivity as desired.  \n"
	"  \n"
	"  High Performance Tips\n"
	"  ---------------------\n"
	"  \n"
	"  Tip 1: If your computer has multiple processors/cores, try -p\n"
	"   \n"
	"  The -p <int> option causes Bowtie to launch <int> parallel search\n"
	"  threads.  Each thread runs on a different processor/core and all\n"
	"  threads find alignments in parallel, increasing alignment throughput\n"
	"  by approximately a multiple of <int>.\n"
	"  \n"
	"  Tip 2: If reporting many alignments per read, try tweaking\n"
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
	"  Tip 3: If bowtie \"thrashes\", try tweaking 'bowtie --offrate'\n"
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
	"  \n"
	"  A note from firsthand experience: I have a MacBook Pro with 2 GB of\n"
	"  RAM and I noticed that 'bowtie -z' would thrash and run very slowly\n"
	"  when aligning reads against the pre-built human genome available from\n"
	"  the Bowtie website.  Using 'bowtie -z --offrate 6' prevented the\n"
	"  thrashing and allowed bowtie to run much faster.\n"
	"\n"
	"  Command Line\n"
	"  ------------\n"
	"\n"
	"  The following is a detailed description of the options used to control\n"
	"  the 'bowtie' aligner:\n"
	"\n"
	" Usage: bowtie [options]* <ebwt> {-1 <mates1> -2 <mates2> | <singles>} [<hits>]\n"
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
	"  <mates1>           Comma-separated list of files containing the #1\n"
	"                     mates (filename usually includes \"_1\"), or, if -c\n"
	"                     is specified, the mate sequences themselves.\n"
	"                     E.g., this might be \"flyA_1.fq,flyB_1.fq\", or, if\n"
	"                     -c is given, this might be \"GGTCATCCT,ACGGGTCGT\".\n"
	"                     Sequences specified with this option must\n"
	"                     correspond file-for-file and read-for-read with\n"
	"                     those specified in <mates2>.  Reads may be a mix\n"
	"                     of different lengths.  Note that as of version\n"
	"                     0.9.9, bowtie can perform paired-end alignment\n"
	"                     only when the alignment policy is -v 0 or -v 1.\n"
	"                     Paired-end support will be expanded in future\n"
	"                     versions. \n"
	"\n"
	"  <mates2>           Comma-separated list of files containing the #2\n"
	"                     mates (filename usually includes \"_2\"), or, if -c\n"
	"                     is specified, the mate sequences themselves.\n"
	"                     E.g., this might be \"flyA_2.fq,flyB_2.fq\", or, if\n"
	"                     -c is given, this might be \"GTATGCTG,AATTCAGGCTG\".\n"
	"                     Sequences specified with this option must\n"
	"                     correspond file-for-file and read-for-read with\n"
	"                     those specified in <mates1>.  Reads may be a mix\n"
	"                     of different lengths.  Note that as of version\n"
	"                     0.9.9, bowtie can perform paired-end alignment\n"
	"                     only when the alignment policy is -v 0 or -v 1.\n"
	"                     Paired-end support will be expanded in future\n"
	"                     versions. \n"
	"\n"
	"  <singles>          A comma-separated list of files containing the\n"
	"                     reads to be aligned, or, if -c is specified, the\n"
	"                     sequences themselves. E.g., this might be\n"
	"                     \"lane1.fq,lane2.fq,lane3.fq,lane4.fq\", or, if -c\n"
	"                     is specified, this might be \"GGTCATCCT,ACGGGTCGT\".\n"
	"                     Reads may be a mix of different lengths.\n"
	"\n"
	"  <hits>             File to write alignments to.  By default,\n"
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
	"  -q                 The query input files (specified as <mates1>,\n"
	"                     <mates2> and <singles>) are FASTQ files (usually\n"
	"                     having extension .fq or .fastq).  This is the\n"
	"                     default.  See also: --solexa-quals and\n"
	"                     --integer-quals.\n"
	"\n"
	"  -f                 The query input files (specified as <mates1>,\n"
	"                     <mates2> and <singles>) are FASTA files (usually\n"
	"                     having extension .fa, .mfa, .fna or similar).  All\n"
	"                     quality values are assumed to be 40 on the Phred\n"
	"                     scale.\n"
	"\n"
	"  -r                 The query input files (specified as <mates1>,\n"
	"                     <mates2> and <singles>) are Raw files: one\n"
	"                     sequence per line, without quality values or\n"
	"                     names.  All quality values are assumed to be 40 on\n"
	"                     the Phred scale.\n"
	"\n"
	"  -c                 The query sequences are given on command line.\n"
	"                     I.e. <mates1>, <mates2> and <singles> are comma-\n"
	"                     separated lists of reads rather than a list of\n"
	"                     read files.\n"
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
	"  --solexa-quals     Convert FASTQ qualities from solexa-scaled (which\n"
	"                     can be negative) to phred-scaled (which can't).\n"
	"                     The formula for conversion is phred-qual =\n"
	"                     10 * log(1 + 10 ** (solexa-qual/10.0)) / log(10).\n"
	"                     Used with -q.  Default: off.\n"
	"\n"
	"  --integer-quals    Quality values are represented in the FASTQ file\n"
	"                     as space-separated ASCII integers, e.g.,\n"
	"                     \"40 40 30 40...\", rather than ASCII characters,\n"
	"                     e.g., \"II?I...\".  Used with -q.  Default: off.\n"
	"\n"
	"   Alignment:\n"
	"   ----------\n"
	"\n"
	"  -n/--seedmms <int> The maximum number of mismatches permitted in the\n"
	"                     \"seed\", which is the first 28 base pairs of the\n"
	"                     read by default (see -l/--seedlen).  This may be\n"
	"                     0, 1, 2 or 3 and the default is 2.\n"
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
	"                     ceiling applies.  The default is 28.\n"
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
	"                     -v is specified, base quality values and the -e,\n"
	"                     -l and -n options are ignored.\n"
	"\n"
	"  -I/--minins <int>  The minimum insert size for valid paired-end\n"
	"                     alignments.  E.g. if -I 60 is specified and a\n"
	"                     paired-end alignment consists of two 20-bp\n"
	"                     alignments in the appropriate orientation with a\n"
	"                     20-bp gap between them, that alignment is\n"
	"                     considered valid (as long as -X is also\n"
	"                     satisfied).  A 19-bp gap would not be valid in\n"
	"                     that case.  Default: 0.\n"
	"\n"
	"  -X/--maxins <int>  The maximum insert size for valid paired-end\n"
	"                     alignments.  E.g. if -X 100 is specified and a\n"
	"                     paired-end alignment consists of two 20-bp\n"
	"                     alignments in the proper orientation with a 60-bp\n"
	"                     gap between them, that alignment is considered\n"
	"                     valid (as long as -I is also satisfied).  A 61-bp\n"
	"                     gap would not be valid in that case.  Default:\n"
	"                     250.\n"
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
	"                     mate2 be forward-oriented.  --ll requires both an\n"
	"                     upstream mate1 and a downstream mate2 to be\n"
	"                     forward-oriented.  Default: --fr (appropriate for\n"
	"                     the Illumina short insert library). \n"
	"\n"
	"  --maxbts           The maximum number of backtracks permitted when\n"
	"                     aligning a read in -n 2 or -n 3 mode (default:\n"
	"                     125).  A \"backtrack\" is the introduction of a\n"
	"                     speculative substitution into the alignment.\n"
	"                     Without this limit, the default parameters will\n"
	"                     sometimes require that 'bowtie' try 100s or 1,000s\n"
	"                     of backtracks to align a read, especially if the\n"
	"                     read has many low-quality bases and/or has no\n"
	"                     valid alignments, slowing bowtie down\n"
	"                     significantly.  The drawback of having a limit is\n"
	"                     that some valid alignments may be missed.  Higher\n"
	"                     limits yield greater sensitivity at the expensive\n"
	"                     of longer running times.  See also: -y/--tryhard.\n"
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
	"   Reporting:\n"
	"   ----------\n"
	"\n"
	"  -k <int>           Report up to <int> valid alignments per read or\n"
	"                     pair (default: 1).  Validity of alignments is\n"
	"                     determined by the alignment policy (combined\n"
	"                     effects of -n, -v, -l, and -e).  If many\n"
	"                     alignments are reported, they may be subject to\n"
	"                     stratification; see --best, --nostrata.  Bowtie is\n"
	"                     designed to be very fast for small -k but BOWTIE\n"
	"                     CAN BECOME VERY SLOW AS -k INCREASES.  If you\n"
	"                     would like to use Bowtie for larger values of -k,\n"
	"                     consider building an index with a denser suffix-\n"
	"                     array sample, i.e. specify a smaller '--offrate'\n"
	"                     when invoking 'bowtie-build' for the relevant\n"
	"                     index.  This will increase the memory footprint of\n"
	"                     the index, but makes alignment much faster for\n"
	"                     large -k.\n"
	"\n"
	"  -a/--all           Report all valid alignments per read or pair\n"
	"                     (default: off).  Validity of alignments is\n"
	"                     determined by the alignment policy (combined\n"
	"                     effects of -n, -v, -l, and -e).  Reported\n"
	"                     alignments may be subject to stratification; see\n"
	"                     --best, --nostrata.  Bowtie is designed to be very\n"
	"                     fast for small -k; BOWTIE CAN BECOME VERY SLOW\n"
	"                     IF -a/--all IS SPECIFIED.  If you would like to\n"
	"                     use Bowtie with -a, consider building an index\n"
	"                     with a denser suffix-array sample, i.e. specify a\n"
	"                     smaller '--offrate' when invoking 'bowtie-build'\n"
	"                     for the relevant index.  This will increase the\n"
	"                     memory footprint of the index, but makes alignment\n"
	"                     much faster in -a mode.\n"
	"\n"
	"  -m <int>           Suppress all alignments for a particular read or\n"
	"                     pair if more than <int> reportable alignments\n"
	"                     exist for it.  Reportable alignments are those\n"
	"                     that would be reported given the -n, -v, -l, -e,\n"
	"                     -k, -a, --best, and --nostrata options.  Default:\n"
	"                     no limit.  Bowtie is designed to be very fast for\n"
	"                     small -m but BOWTIE CAN BECOME VERY SLOW AS -m\n"
	"                     INCREASES.\n"
	"\n"
	"  --best             Reported singleton alignments must belong to the\n"
	"                     best possible alignment \"stratum\" (default: off).\n"
	"                     A stratum is a category defined by the number of\n"
	"                     mismatches present in the alignment (for -n, the\n"
	"                     number of mismatches present in the seed region of\n"
	"                     the alignment).  E.g., if --best is not specified,\n"
	"                     Bowtie may sometimes report an alignment with 2\n"
	"                     mismatches in the seed even though there exists an\n"
	"                     unreported alignment with 1 mismatch in the seed.\n"
	"                     BOWTIE IS ABOUT 3-5 TIMES SLOWER WHEN --best IS\n"
	"                     SPECIFIED.\n"
	"\n"
	"  --nostrata         If many valid singleton alignments exist and are\n"
	"                     reportable (according to the --best and -k\n"
	"                     options) and they fall into more than one\n"
	"                     alignment \"stratum\", report all of them.  By\n"
	"                     default, Bowtie only reports those alignments that\n"
	"                     fall into the best stratum for which alignments\n"
	"                     were found, i.e., the one with fewest mismatches.\n"
	"                     BOWTIE CAN BECOME VERY SLOW WHEN --nostrata IS\n"
	"                     COMBINED WITH -k OR -a.\n"
	"\n"
	"   Output:\n"
	"   -------\n"
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
	"  -b/--binout        Output alignments in a concise binary format.  If\n"
	"                     this is specified, <hit_outfile> must also be\n"
	"                     specified.\n"
	"\n"
	"  -t/--time          Print the amount of wall-clock time taken by each\n"
	"                     search phase and index turnover.\n"
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
	"  --unfa <filename>  Write all reads that fail to align to a FASTA file\n"
	"                     with name <filename>.  Paired-end reads will be\n"
	"                     written to two parallel files with \"_1\" and \"_2\"\n"
	"                     inserted in the filename, e.g., if <filename> is\n"
	"                     unaligned.fa, the #1 and #2 mates that fail to\n"
	"                     align will be written to unaligned_1.fa and\n"
	"                     unaligned_2.fa respectively.  Unless --maxfa is\n"
	"                     also specified, reads with a number of valid\n"
	"                     alignments exceeding the limit set with the -m\n"
	"                     option are also written to <filename>.\n"
	"\n"
	"  --unfq <filename>  Write all reads that fail to align to a FASTQ file\n"
	"                     with name <filename>.  Paired-end reads will be\n"
	"                     written to two parallel files with \"_1\" and \"_2\"\n"
	"                     inserted in the filename, e.g., if <filename> is\n"
	"                     unaligned.fq, the #1 and #2 mates that fail to\n"
	"                     align will be written to unaligned_1.fq and\n"
	"                     unaligned_2.fq respectively.  Unless --maxfq is\n"
	"                     also specified, reads with a number of valid\n"
	"                     alignments exceeding the limit set with the -m\n"
	"                     option are also written to <filename>.\n"
	"\n"
	"  --maxfa <filename> Write all reads with a number of valid alignments\n"
	"                     exceeding the limit set with the -m option to a\n"
	"                     FASTA file with given name.  Paired-end reads that\n"
	"                     exceed the -m limit are written to parallel files\n"
	"                     with \"_1\" and \"_2\" inserted into the filename, in\n"
	"                     the same way as with --unfa.  These reads are not\n"
	"                     written to the file specified with --unfa.\n"
	"\n"
	"  --maxfq <filename> Write all reads with a number of valid alignments\n"
	"                     exceeding the limit set with the -m option to a\n"
	"                     FASTQ file with given name.  Paired-end reads that\n"
	"                     exceed the -m limit are written to parallel files\n"
	"                     with \"_1\" and \"_2\" inserted into the filename, in\n"
	"                     the same way as with --unfq.  These reads are not\n"
	"                     written to the file specified with --unfq.\n"
	"\n"
	"   Performance:\n"
	"   ------------\n"
	"\n"
	"  -p/--threads <int> Launch <int> parallel search threads (default: 1).\n"
	"                     Threads will run on separate processors/cores and\n"
	"                     synchronize when grabbing reads and outputting\n"
	"                     alignments.  Searching for alignments is almost\n"
	"                     totally parallel, and speedup is close to linear.\n"
	"                     Speedup suffers somewhat in -z mode.  This option\n"
	"                     is only available if bowtie is linked with the\n"
	"                     pthreads library (i.e. if BOWTIE_PTHREADS=0 is not\n"
	"                     specified at build time).\n"
	"\n"
	"  -z/--phased        Alternate between using the forward and mirror\n"
	"                     indexes in a series of phases such that only one\n"
	"                     \"half\" of the index is resident in memory at one\n"
	"                     time.  This uses about half the amount of memory\n"
	"                     as the default (which keeps both forward and\n"
	"                     mirror indexes resident in memory at once), but is\n"
	"                     somewhat slower, scales worse (see -p), and is\n"
	"                     incompatible with use of --best or -k greater than\n"
	"                     1.\n"
	"\n"
	"  -o/--offrate <int> Override the offrate of the index with <int>.  If\n"
	"                     <int> is greater than the offrate used to build\n"
	"                     the index, then some row markings are discarded\n"
	"                     when the index is read into memory.  This reduces\n"
	"                     the memory footprint of the aligner but requires\n"
	"                     more time to calculate text offsets.  <int> must\n"
	"                     be greater than the value used to build the index.\n"
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
	"   4. 0-based offset into the reference sequence where leftmost\n"
	"      character of the alignment occurs\n"
	"\n"
	"   5. Read sequence (reverse-complemented if orientation is '-')\n"
	"\n"
	"   6. Read qualities (reversed if orientation is '-')\n"
	"\n"
	"   7. Reserved\n"
	"\n"
	"   8. Comma-separated list of mismatch descriptors.  If there are no\n"
	"      mismatches in the alignment, this field is empty.  A single\n"
	"      descriptor has the format offset:reference-base>read-base.  The\n"
	"      offset is expressed as a 0-based offset from the high-quality\n"
	"      (5') end of the read. \n"
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
			case '1': tokenize(optarg, ",", mates1); break;
			case '2': tokenize(optarg, ",", mates2); break;
	   		case 'f': format = FASTA; break;
	   		case 'q': format = FASTQ; break;
	   		case 'r': format = RAW; break;
	   		case 'c': format = CMDLINE; break;
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
	   		case ARG_CONCISE: outType = CONCISE; break;
	   		case 'b': outType = BINARY; break;
	   		case ARG_REFOUT: refOut = true; break;
	   		case ARG_SEED_EXTEND: seedAndExtend = true; break;
	   		case ARG_NOOUT: outType = NONE; break;
	   		case ARG_USE_SPINLOCK: useSpinlock = false; break;
	   		case ARG_SHARED_MEM: {
#ifdef BOWTIE_SHARED_MEM
	   			useShmem = true; break;
#else
	   			cerr << "Shared-memory mode is disabled because bowtie was not compiled with the" << endl
	   			     << "BOWTIE_SHARED_MEM option." << endl;
	   			exit(1);
#endif
	   		}
	   		case ARG_DUMP_NOHIT: dumpNoHits = new ofstream(".nohits.dump"); break;
	   		case ARG_DUMP_HHHIT: dumpHHHits = new ofstream(".hhhits.dump"); break;
	   		case ARG_UNFA: dumpUnalFaBase = optarg; break;
	   		case ARG_UNFQ: dumpUnalFqBase = optarg; break;
	   		case ARG_MAXFA: dumpMaxFaBase = optarg; break;
	   		case ARG_MAXFQ: dumpMaxFqBase = optarg; break;
			case ARG_SOLEXA_QUALS: solexa_quals = true; break;
			case ARG_INTEGER_QUALS: integer_quals = true; break;
			case ARG_FORGIVE_INPUT: forgiveInput = true; break;
			case ARG_NOMAQROUND: noMaqRound = true; break;
			case 'z': fullIndex = false; break;
			case ARG_REFIDX: noRefNames = true; break;
			case ARG_STATEFUL: stateful = true; break;
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
	   			exit(1);
#endif
	   			nthreads = parseInt(1, "-p/--threads arg must be at least 1");
	   			break;
	   		case ARG_FILEPAR:
#ifndef BOWTIE_PTHREADS
	   			cerr << "--filepar is disabled because bowtie was not compiled with pthreads support" << endl;
	   			exit(1);
#endif
	   			fileParallel = true;
	   			break;
	   		case 'v':
	   			maqLike = 0;
	   			mismatches = parseInt(0, "-v arg must be at least 0");
	   			if(mismatches > 3) {
	   				cerr << "-v arg must be at most 3" << endl;
	   				exit(1);
	   			}
	   			break;
	   		case '3': trim3 = parseInt(0, "-3/--trim3 arg must be at least 0"); break;
	   		case '5': trim5 = parseInt(0, "-5/--trim5 arg must be at least 0"); break;
	   		case 'o': offRate = parseInt(1, "-o/--offrate arg must be at least 1"); break;
	   		case ARG_ISARATE: isaRate = parseInt(0, "--isarate arg must be at least 0"); break;
	   		case 'e': qualThresh = parseInt(1, "-e/--err arg must be at least 1"); break;
	   		case 'n': seedMms = parseInt(0, "-n/--seedmms arg must be at least 0"); maqLike = 1; break;
	   		case 'l': seedLen = parseInt(20, "-l/--seedlen arg must be at least 20"); break;
	   		case 'h': printLongUsage(cout); exit(0); break;
	   		case '?': printUsage(cerr); exit(1); break;
	   		case 'a': allHits = true; break;
	   		case 'y': tryHard = true; break;
	   		case ARG_BEST: onlyBest = true; break;
	   		case ARG_SPANSTRATA: spanStrata = true; break;
	   		case ARG_VERBOSE: verbose = true; break;
	   		case ARG_QUIET: quiet = true; break;
	   		case ARG_SANITY: sanityCheck = true; break;
	   		case 't': timing = true; break;
	   		case ARG_NO_FW: nofw = true; break;
	   		case ARG_NO_RC: norc = true; break;
			case ARG_MAXBTS:  maxBts  = parseInt(0, "--maxbts must be positive");  break;
			case ARG_MAXBTS0: maxBts0 = parseInt(0, "--maxbts0 must be positive"); break;
			case ARG_MAXBTS1: maxBts1 = parseInt(0, "--maxbts1 must be positive"); break;
			case ARG_MAXBTS2: maxBts2 = parseInt(0, "--maxbts2 must be positive"); break;
	   		case ARG_DUMP_PATS: patDumpfile = optarg; break;
	   		case ARG_PARTITION: partitionSz = parseInt(0, "--partition must be positive"); break;
	   		case ARG_ORIG:
   				if(optarg == NULL || strlen(optarg) == 0) {
   					cerr << "--orig arg must be followed by a string" << endl;
   					printUsage(cerr);
   					exit(1);
   				}
   				origString = optarg;
	   			break;

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
	if(rangeMode) {
		// Tell the Ebwt loader to ignore the suffix-array portion of
		// the index.  We don't need it because the user isn't asking
		// for bowtie to report reference positions (just matrix
		// ranges).
		offRate = 32;
	}
	if(maqLike) {
		revcomp = true;
	}
	if(mates1.size() != mates2.size()) {
		cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
		     << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
		     << "sequences must be specified with -1 and -2." << endl;
		exit(1);
	}
	if(format != CMDLINE) {
		for(size_t i = 0; i < mates1.size(); i++) {
			for(size_t j = 0; j < mates2.size(); j++) {
				if(mates1[i] == mates2[j]) {
					cerr << "Warning: Same mate file \"" << mates1[i] << "\" appears as argument to both -1 and -2" << endl;
				}
			}
		}
	}
	if(mates1.size() > 0 || mates2.size() > 0) {
		spanStrata = true;
	}
	if((mates1.size() > 0 || mates2.size() > 0) && maqLike) {
		cerr << "Paired-end mode is not yet compatible with -n mode; use -v 0 or -v 1 instead." << endl;
		exit(1);
	}
	if(!fullIndex) {
		bool error = false;
		if(khits > 1) {
			cerr << "When -z/--phased is used, -k X for X > 1 is unavailable" << endl;
			error = true;
		}
		if(onlyBest) {
			cerr << "When -z/--phased is used, --best is unavailable" << endl;
			error = true;
		}
		if(allHits && !spanStrata) {
			cerr << "When -a/--all and -z/--phased are used, --nostrata must also be used." << endl
			     << "Stratified all-hits search cannot be combined with phased search." << endl;
			error = true;
		}
		if(mates1.size() > 0 || mates2.size() > 0) {
			cerr << "When -z/--phased is used, paired-end mode is unavailable" << endl;
			error = true;
		}
		if(error) exit(1);
	}
	//if(!maqLike && onlyBest) {
	//	onlyBest = false;
	//	stateful = true;
	//}
	if(tryHard) {
		// Increase backtracking limit to huge number
		maxBts = INT_MAX;
		// Increase number of paired-end scan attempts to huge number
		mixedAttemptLim = UINT_MAX;
	}
	// If both -s and -u are used, we need to adjust qUpto accordingly
	// since it uses patid to know if we've reached the -u limit (and
	// patids are all shifted up by skipReads characters)
	if(qUpto + skipReads > qUpto) {
		qUpto += skipReads;
	}
}

static char *argv0 = NULL;

#define FINISH_READ(p) \
	/* Don't do finishRead if the read isn't legit or if the read was skipped by the doneMask */ \
	if(!p->empty()) { \
		sink->finishRead(*p, !skipped); \
	} \
	skipped = false;

#define FINISH_READ_WITH_HITMASK(p, r) \
	/* Don't do finishRead if the read isn't legit */ \
	if(!p->empty()) { \
		/* r = whether to consider reporting the read as unaligned */ \
		bool reportUnAl = r; \
		if(reportUnAl) { \
			/* If the done-mask already shows the read as done, */ \
			/* then we already reported the unaligned read and */ \
			/* should refrain from re-reporting*/ \
			reportUnAl = !skipped; \
			if(reportUnAl) { \
				/* If there hasn't been a hit reported, then report */ \
				/* read as unaligned */ \
				reportUnAl = !hitMask.test(p->patid()); \
			} \
		} \
		if(sink->finishRead(*p, reportUnAl) > 0) { \
			/* We reported a hit for the read, so we set the */ \
			/* appropriate bit in the hitMask to prevent it from */ \
			/* being reported as unaligned. */ \
			if(!reportUnAl && sink->dumpsReads()) { \
				hitMask.setOver(p->patid()); \
			} \
		} \
	} \
	skipped = false;

/// Macro for getting the next read, possibly aborting depending on
/// whether the result is empty or the patid exceeds the limit, and
/// marshaling the read into convenient variables.
#define GET_READ(p) \
	p->nextReadPair(); \
	if(p->empty() || p->patid() >= qUpto) break; \
	assert(!empty(p->bufa().patFw)); \
	String<Dna5>& patFw  = p->bufa().patFw;  \
	String<Dna5>& patRc  = p->bufa().patRc;  \
	String<char>& qualFw = p->bufa().qualFw; \
	String<char>& qualRc = p->bufa().qualRc; \
	String<char>& name   = p->bufa().name;   \
	uint32_t      patid  = p->patid();       \
	params.setPatId(patid); \
	if(lastLen == 0) lastLen = seqan::length(patFw); \
	if(qSameLen && seqan::length(patFw) != lastLen) { \
		throw runtime_error("All reads must be the same length"); \
	}

/// Macro for getting the forward oriented version of next read,
/// possibly aborting depending on whether the result is empty or the
/// patid exceeds the limit, and marshaling the read into convenient
/// variables.
#define GET_READ_FW(p) \
	p->nextReadPair(); \
	if(p->empty() || p->patid() >= qUpto) break; \
	params.setPatId(p->patid()); \
	assert(!empty(p->bufa().patFw)); \
	String<Dna5>& patFw  = p->bufa().patFw;  \
	String<char>& qualFw = p->bufa().qualFw; \
	String<char>& name   = p->bufa().name;   \
	uint32_t      patid  = p->patid();     \
	if(lastLen == 0) lastLen = seqan::length(patFw); \
	if(qSameLen && seqan::length(patFw) != lastLen) { \
		throw runtime_error("All reads must be the same length"); \
	}

#ifdef BOWTIE_PTHREADS
#define WORKER_EXIT() \
	if((long)vp != 0L) { \
    	pthread_exit(NULL); \
    } \
	delete patsrc; \
	delete sink; \
    return NULL;
#else
#define WORKER_EXIT() \
	delete patsrc; \
	delete sink; \
	return NULL;
#endif


/// Create a PatternSourcePerThread for the current thread according
/// to the global params and return a pointer to it
static PatternSourcePerThreadFactory*
createPatsrcFactory(PairedPatternSource& _patsrc, int tid) {
	PatternSourcePerThreadFactory *patsrcFact;
	if(randReadsNoSync) {
		patsrcFact = new RandomPatternSourcePerThreadFactory(numRandomReads, lenRandomReads, nthreads, tid, false);
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
    if(spanStrata) {
		if(!allHits) {
			if(onlyBest) {
				// First N best, spanning strata
				sink = new FirstNBestHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			} else {
				// First N good; "good" inherently ignores strata
				sink = new FirstNGoodHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			}
		} else {
			// All hits, spanning strata
			sink = new AllHitSinkPerThreadFactory(_sink, mhits, sanity);
		}
    } else {
		if(!allHits) {
			if(onlyBest) {
				// First N best, not spanning strata
				sink = new FirstNBestStratifiedHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			} else {
				// First N good; "good" inherently ignores strata
				sink = new FirstNGoodHitSinkPerThreadFactory(_sink, khits, mhits, sanity);
			}
		} else {
			// All hits, not spanning strata
			sink = new AllStratifiedHitSinkPerThreadFactory(_sink, mhits, sanity);
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
	PairedPatternSource& _patsrc = *exactSearch_patsrc;
	HitSink& _sink               = *exactSearch_sink;
	Ebwt<String<Dna> >& ebwt     = *exactSearch_ebwt;
	vector<String<Dna5> >& os    = *exactSearch_os;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	// Per-thread initialization
	uint32_t lastLen = 0;
	PatternSourcePerThread *patsrc = createPatsrcFactory(_patsrc, (int)(long)vp)->create();
	HitSinkPerThread* sink = createSinkFactory(_sink, sanity)->create();
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
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed,           // seed
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
	PairedPatternSource& _patsrc = *exactSearch_patsrc;
	HitSink& _sink               = *exactSearch_sink;
	Ebwt<String<Dna> >& ebwt     = *exactSearch_ebwt;
	vector<String<Dna5> >& os    = *exactSearch_os;
	BitPairReference* refs       =  exactSearch_refs;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, (int)(long)vp);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);

	UnpairedExactAlignerV1Factory alSEfact(
			ebwt,
			NULL,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			os,
			rangeMode,
			verbose,
			seed);
	PairedExactAlignerV1Factory alPEfact(
			ebwt,
			NULL,
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
			refs, os,
			rangeMode,
			verbose,
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
#ifdef BOWTIE_PTHREADS
	if((long)vp != 0L) pthread_exit(NULL);
#endif
	return NULL;
}

#define SET_A_FW(bt, p, params) \
	bt.setQuery(&p->bufa().patFw, &p->bufa().qualFw, &p->bufa().name); \
	params.setFw(true);
#define SET_A_RC(bt, p, params) \
	bt.setQuery(&p->bufa().patRc, &p->bufa().qualRc, &p->bufa().name); \
	params.setFw(false);
#define SET_B_FW(bt, p, params) \
	bt.setQuery(&p->bufb().patFw, &p->bufb().qualFw, &p->bufb().name); \
	params.setFw(true);
#define SET_B_RC(bt, p, params) \
	bt.setQuery(&p->bufb().patRc, &p->bufb().qualRc, &p->bufb().name); \
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
		Timer _t(cout, "Time loading forward index: ", timing);
		ebwt.loadIntoMemory();
	}

	BitPairReference *refs = NULL;
	if(mates1.size() > 0 && mixedThresh < 0xffffffff) {
		Timer _t(cout, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, sanityCheck, NULL, &os);
		if(!refs->loaded()) exit(1);
	}
	exactSearch_refs   = refs;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	int numAdditionalThreads = nthreads-1;
	pthread_t *threads = new pthread_t[numAdditionalThreads];
#endif

	{
		Timer _t(cout, "Time for 0-mismatch search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			if(stateful)
				pthread_create(&threads[i], &pt_attr, exactSearchWorkerStateful, (void *)(long)(i+1));
			else
				pthread_create(&threads[i], &pt_attr, exactSearchWorker, (void *)(long)(i+1));
		}
#endif
		if(stateful) exactSearchWorkerStateful((void*)0L);
		else         exactSearchWorker((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < numAdditionalThreads; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
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
	PairedPatternSource&   _patsrc       = *mismatchSearch_patsrc;
	HitSink&               _sink         = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw        = *mismatchSearch_ebwtFw;
	vector<String<Dna5> >& os            = *mismatchSearch_os;
	SyncBitset&            doneMask      = *mismatchSearch_doneMask;
	SyncBitset&            hitMask       = *mismatchSearch_hitMask;
    bool sanity = sanityCheck && !os.empty() && !rangeMode;
	uint32_t lastLen = 0; // for checking if all reads have same length
	PatternSourcePerThread* patsrc = createPatsrcFactory(_patsrc, (int)(long)vp)->create();
	HitSinkPerThread* sink = createSinkFactory(_sink, sanity)->create();
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
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed,           // seed
	        &os,
	        false);         // considerQuals
	bool skipped = false;
	while(true) {
		FINISH_READ_WITH_HITMASK(patsrc, false);
		GET_READ(patsrc);
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_1mm_phase1.c"
		#undef DONEMASK_SET
	} // End read loop
	FINISH_READ_WITH_HITMASK(patsrc, false);
    WORKER_EXIT();
}

static void* mismatchSearchWorkerPhase2(void *vp){
	PairedPatternSource&   _patsrc      = *mismatchSearch_patsrc;
	HitSink&               _sink        = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtBw       = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os           = *mismatchSearch_os;
	SyncBitset&            doneMask     = *mismatchSearch_doneMask;
	SyncBitset&            hitMask      = *mismatchSearch_hitMask;
    // Per-thread initialization
    bool sanity = sanityCheck && !os.empty() && !rangeMode;
	uint32_t lastLen = 0; // for checking if all reads have same length
	PatternSourcePerThread* patsrc = createPatsrcFactory(_patsrc, (int)(long)vp)->create();
	HitSinkPerThread* sink = createSinkFactory(_sink, sanity)->create();
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
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed+1,         // seed
	        &os,
	        false);         // considerQuals
	bool skipped = false;
	while(true) {
		FINISH_READ_WITH_HITMASK(patsrc, true);
		GET_READ(patsrc);
		if(doneMask.test(patid)) { skipped = true; continue; }
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		#include "search_1mm_phase2.c"
	} // End read loop
	FINISH_READ_WITH_HITMASK(patsrc, true);
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
		Timer _t(cout, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory();
	}

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

    _patsrc.setReverse(false); // don't reverse patterns

	// Phase 1
    {
		Timer _t(cout, "Time for 1-mismatch Phase 1 of 2: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerPhase1, (void *)(long)(i+1));
		}
#endif
		mismatchSearchWorkerPhase1((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
    }

	// Release most of the memory associated with the forward Ebwt
    ebwtFw.evictFromMemory();
	{
		// Load the rest of (vast majority of) the backward Ebwt into
		// memory
		Timer _t(cout, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory();
	}
    _patsrc.reset();          // reset pattern source to 1st pattern
    _patsrc.setReverse(true); // reverse patterns
	// Sanity-check the restored version of the Ebwt
	if(sanityCheck && !os.empty()) {
		ebwtBw.checkOrigs(os, true);
	}

	// Phase 2
	{
		Timer _t(cout, "Time for 1-mismatch Phase 2 of 2: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerPhase2, (void *)(long)(i+1));
		}
#endif
		mismatchSearchWorkerPhase2((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
}

/**
 * A statefulness-aware worker driver.  Uses Unpaired/Paired1mmAlignerV1.
 */
static void *mismatchSearchWorkerFullStateful(void *vp) {
	PairedPatternSource&   _patsrc = *mismatchSearch_patsrc;
	HitSink&               _sink   = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw  = *mismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw  = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os      = *mismatchSearch_os;
	BitPairReference*      refs    =  mismatchSearch_refs;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, (int)(long)vp);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);

	Unpaired1mmAlignerV1Factory alSEfact(
			ebwtFw,
			&ebwtBw,
			_sink,
			*sinkFact,
			NULL, //&cacheFw,
			NULL, //&cacheBw,
			cacheLimit,
			os,
			rangeMode,
			verbose,
			seed);
	Paired1mmAlignerV1Factory alPEfact(
			ebwtFw,
			&ebwtBw,
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
			refs, os,
			rangeMode,
			verbose,
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
#ifdef BOWTIE_PTHREADS
	if((long)vp != 0L) pthread_exit(NULL);
#endif
	return NULL;
}

static void* mismatchSearchWorkerFull(void *vp){
	PairedPatternSource&   _patsrc      = *mismatchSearch_patsrc;
	HitSink&               _sink        = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw       = *mismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw       = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os           = *mismatchSearch_os;
    // Per-thread initialization
    bool sanity = sanityCheck && !os.empty() && !rangeMode;
	uint32_t lastLen = 0; // for checking if all reads have same length
	PatternSourcePerThread* patsrc = createPatsrcFactory(_patsrc, (int)(long)vp)->create();
	HitSinkPerThread* sink = createSinkFactory(_sink, sanity)->create();
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
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed,           // seed
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
		patsrc->reverseRead();
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
		Timer _t(cout, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory();
	}
	{
		// Load the other half of the index into memory
		Timer _t(cout, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory();
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	if(mates1.size() > 0 && mixedThresh < 0xffffffff) {
		Timer _t(cout, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, sanityCheck, NULL, &os);
		if(!refs->loaded()) exit(1);
	}
	mismatchSearch_refs = refs;

#ifdef BOWTIE_PTHREADS
	// Allocate structures for threads
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

    _patsrc.setReverse(false); // don't reverse patterns
    {
		Timer _t(cout, "Time for 1-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			if(stateful)
				pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerFullStateful, (void *)(long)(i+1));
			else
				pthread_create(&threads[i], &pt_attr, mismatchSearchWorkerFull, (void *)(long)(i+1));
		}
#endif
		// Go to town
		if(stateful) mismatchSearchWorkerFullStateful((void*)0L);
		else         mismatchSearchWorkerFull((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
    }
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
	if(refs != NULL) delete refs;
}

#define SWITCH_TO_FW_INDEX() { \
	/* Evict the mirror index from memory if necessary */ \
	if(ebwtBw.isInMemory()) ebwtBw.evictFromMemory(); \
	assert(!ebwtBw.isInMemory()); \
	/* Load the forward index into memory if necessary */ \
	if(!ebwtFw.isInMemory()) { \
		Timer _t(cout, "Time loading forward index: ", timing); \
		ebwtFw.loadIntoMemory(); \
	} \
	assert(ebwtFw.isInMemory()); \
	_patsrc.reset(); /* rewind pattern source to first pattern */ \
	_patsrc.setReverse(false); /* tell pattern source not to reverse patterns */ \
}

#define SWITCH_TO_BW_INDEX() { \
	/* Evict the forward index from memory if necessary */ \
	if(ebwtFw.isInMemory()) ebwtFw.evictFromMemory(); \
	assert(!ebwtFw.isInMemory()); \
	/* Load the forward index into memory if necessary */ \
	if(!ebwtBw.isInMemory()) { \
		Timer _t(cout, "Time loading mirror index: ", timing); \
		ebwtBw.loadIntoMemory(); \
	} \
	assert(ebwtBw.isInMemory()); \
	_patsrc.reset(); /* rewind pattern source to first pattern */ \
	_patsrc.setReverse(true); /* tell pattern source to reverse patterns */ \
}

#define ASSERT_NO_HITS_FW(ebwtfw) \
	if(sanityCheck && os.size() > 0) { \
		vector<Hit> hits; \
		vector<int> strata; \
		uint32_t threeRevOff = (seedMms <= 3) ? s : 0; \
		uint32_t twoRevOff   = (seedMms <= 2) ? s : 0; \
		uint32_t oneRevOff   = (seedMms <= 1) ? s : 0; \
		uint32_t unrevOff    = (seedMms == 0) ? s : 0; \
		::naiveOracle( \
		        os, \
				patFw, \
				plen, \
		        qualFw, \
		        name, \
		        patid, \
		        hits, \
		        strata, \
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
		vector<int> strata; \
		uint32_t threeRevOff = (seedMms <= 3) ? s : 0; \
		uint32_t twoRevOff   = (seedMms <= 2) ? s : 0; \
		uint32_t oneRevOff   = (seedMms <= 1) ? s : 0; \
		uint32_t unrevOff    = (seedMms == 0) ? s : 0; \
		::naiveOracle( \
		        os, \
				patRc, \
				plen, \
		        qualRc, \
		        name, \
		        patid, \
		        hits, \
		        strata, \
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
	PairedPatternSource&           _patsrc  = *twoOrThreeMismatchSearch_patsrc;   \
	HitSink&                       _sink    = *twoOrThreeMismatchSearch_sink;     \
	vector<String<Dna5> >&         os       = *twoOrThreeMismatchSearch_os;       \
	bool                           two      = twoOrThreeMismatchSearch_two; \
	uint32_t lastLen = 0; \
	PatternSourcePerThread* patsrc = createPatsrcFactory(_patsrc, (int)(long)vp)->create(); \
	HitSinkPerThread* sink = createSinkFactory(_sink, false)->create(); \
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
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed,           // seed
	        &os,
	        false);         // considerQuals
	bool skipped = false;
    while(true) { // Read read-in loop
    	FINISH_READ_WITH_HITMASK(patsrc, false);
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
	FINISH_READ_WITH_HITMASK(patsrc, false);
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
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (no)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+1,         // seed
		    &os,
		    false);         // considerQuals
	bool skipped = false;
    while(true) {
    	FINISH_READ_WITH_HITMASK(patsrc, false);
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
	FINISH_READ_WITH_HITMASK(patsrc, false);
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
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+3,         // seed
		    &os,
		    false);         // considerQuals
	GreedyDFSRangeSource bthh3(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+5,         // seed
		    &os,
		    false,          // considerQuals
		    true);          // halfAndHalf
	bool skipped = false;
    while(true) {
    	FINISH_READ_WITH_HITMASK(patsrc, true);
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
	FINISH_READ_WITH_HITMASK(patsrc, true);
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
#endif

	// Load forward index
	SWITCH_TO_FW_INDEX();
    { // Phase 1
		Timer _t(cout, "End-to-end 2/3-mismatch Phase 1 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerPhase1, (void *)(long)(i+1));
		}
#endif
		twoOrThreeMismatchSearchWorkerPhase1((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
    }
	// Unload forward index and load mirror index
	SWITCH_TO_BW_INDEX();
	{ // Phase 2
		Timer _t(cout, "End-to-end 2/3-mismatch Phase 2 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerPhase2, (void *)(long)(i+1));
		}
#endif
		twoOrThreeMismatchSearchWorkerPhase2((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
	}
	SWITCH_TO_FW_INDEX();
	{ // Phase 3
		Timer _t(cout, "End-to-end 2/3-mismatch Phase 3 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerPhase3, (void *)(long)(i+1));
		}
#endif
		twoOrThreeMismatchSearchWorkerPhase3((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
	return;
}

/**
 * A statefulness-aware worker driver.  Uses UnpairedExactAlignerV1.
 */
static void *twoOrThreeMismatchSearchWorkerStateful(void *vp) {
	PairedPatternSource&   _patsrc = *twoOrThreeMismatchSearch_patsrc;
	HitSink&               _sink   = *twoOrThreeMismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw  = *twoOrThreeMismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw  = *twoOrThreeMismatchSearch_ebwtBw;
	vector<String<Dna5> >& os      = *twoOrThreeMismatchSearch_os;
	BitPairReference*      refs    =  twoOrThreeMismatchSearch_refs;
	static bool            two     =  twoOrThreeMismatchSearch_two;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	PatternSourcePerThreadFactory* patsrcFact = createPatsrcFactory(_patsrc, (int)(long)vp);
	HitSinkPerThreadFactory* sinkFact = createSinkFactory(_sink, sanity);

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
			os,
			rangeMode,
			verbose,
			seed);
	Paired23mmAlignerV1Factory alPEfact(
			ebwtFw,
			&ebwtBw,
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
			refs, os,
			rangeMode,
			verbose,
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
#ifdef BOWTIE_PTHREADS
	if((long)vp != 0L) pthread_exit(NULL);
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
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed,           // seed
	        &os,
	        false);         // considerQuals
	GreedyDFSRangeSource bt2(
			&ebwtBw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (no)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+1,         // seed
		    &os,
		    false);         // considerQuals
	GreedyDFSRangeSource bt3(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh (none)
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+3,         // seed
		    &os,
		    false);         // considerQuals
	GreedyDFSRangeSource bthh3(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        rangeMode,      // reportRanges
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+5,         // seed
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
		patsrc->reverseRead();
		#include "search_23mm_phase2.c"
		patsrc->reverseRead();
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
		Timer _t(cout, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory();
	}
	{
		// Load the other half of the index into memory
		Timer _t(cout, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory();
	}
	// Create range caches, which are shared among all aligners
	BitPairReference *refs = NULL;
	if(mates1.size() > 0 && mixedThresh < 0xffffffff) {
		Timer _t(cout, "Time loading reference: ", timing);
		refs = new BitPairReference(adjustedEbwtFileBase, sanityCheck, NULL, &os);
		if(!refs->loaded()) exit(1);
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
#endif

    {
		Timer _t(cout, "End-to-end 2/3-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			if(stateful)
				pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerStateful, (void *)(long)(i+1));
			else
				pthread_create(&threads[i], &pt_attr, twoOrThreeMismatchSearchWorkerFull, (void *)(long)(i+1));
		}
#endif
		if(stateful) twoOrThreeMismatchSearchWorkerStateful((void*)0L);
		else         twoOrThreeMismatchSearchWorkerFull((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
    }
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
	return;
}

static PairedPatternSource*           seededQualSearch_patsrc;
static HitSink*                       seededQualSearch_sink;
static Ebwt<String<Dna> >*            seededQualSearch_ebwtFw;
static Ebwt<String<Dna> >*            seededQualSearch_ebwtBw;
static vector<String<Dna5> >*         seededQualSearch_os;
static SyncBitset*                    seededQualSearch_doneMask;
static SyncBitset*                    seededQualSearch_hitMask;
static PartialAlignmentManager*       seededQualSearch_pamFw;
static PartialAlignmentManager*       seededQualSearch_pamRc;
static int                            seededQualSearch_qualCutoff;

#define SEEDEDQUAL_WORKER_SETUP() \
	PairedPatternSource&     _patsrc    = *seededQualSearch_patsrc;    \
	HitSink&                 _sink      = *seededQualSearch_sink;      \
	vector<String<Dna5> >&   os         = *seededQualSearch_os;        \
	int                      qualCutoff = seededQualSearch_qualCutoff; \
	uint32_t lastLen = 0; \
	PatternSourcePerThread* patsrc = createPatsrcFactory(_patsrc, (int)(long)vp)->create(); \
	HitSinkPerThread* sink = createSinkFactory(_sink, false)->create(); \
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
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partials
	        NULL,                  // mutations
	        verbose,               // verbose
	        seed,                  // seed
	        &os,
	        false);                // considerQuals
	GreedyDFSRangeSource bt1(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partials
	        NULL,                  // mutations
	        verbose,               // verbose
	        seed,                  // seed
	        &os,
	        true,                  // considerQuals
	        false, !noMaqRound);
	bool skipped = false;
    while(true) {
    	FINISH_READ_WITH_HITMASK(patsrc, false);
    	GET_READ(patsrc);
		size_t plen = length(patFw);
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_seeded_phase1.c"
		#undef DONEMASK_SET
    }
	FINISH_READ_WITH_HITMASK(patsrc, false);
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
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (no)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+1,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for partial alignments for case 4R
	GreedyDFSRangeSource btr2(
			&ebwtBw, params,
	        qualCutoff,            // qualThresh (none)
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // report partials (up to seedMms mms)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamRc,                 // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+2,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	bool skipped = false;
    while(true) {
    	FINISH_READ_WITH_HITMASK(patsrc, seedMms == 0);
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
	FINISH_READ_WITH_HITMASK(patsrc, seedMms == 0);
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
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // reportPartials (do)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamFw,                 // seedlings
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+3,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for hits for case 4R by extending
	// the partial alignments found in Phase 2
	GreedyDFSRangeSource btr3(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
		    seed+4,  // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// The half-and-half GreedyDFSRangeSource
	GreedyDFSRangeSource btr23(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
		    seed+5,  // seed
		    &os,
		    true,    // considerQuals
		    true,    // halfAndHalf
		    !noMaqRound);
	vector<PartialAlignment> pals;
	String<QueryMutation> muts;
	bool skipped = false;
    while(true) {
    	FINISH_READ_WITH_HITMASK(patsrc, false);
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
	FINISH_READ_WITH_HITMASK(patsrc, false);
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
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
	        seed+6,  // seed
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// Half-and-half GreedyDFSRangeSource for forward read
	GreedyDFSRangeSource btf24(
			&ebwtBw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
	        seed+7,  // seed
	        &os,
	        true,    // considerQuals
	        true,    // halfAndHalf
	        !noMaqRound);
	vector<PartialAlignment> pals;
	String<QueryMutation> muts;
	bool skipped = false;
    while(true) {
    	FINISH_READ_WITH_HITMASK(patsrc, true);
		GET_READ_FW(patsrc);
		if(doneMask.testUnsync(patid)) { skipped = true; continue; }
		size_t plen = length(patFw);
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		#define DONEMASK_SET(p) doneMask.set(p)
		#include "search_seeded_phase4.c"
		#undef DONEMASK_SET
    }
	FINISH_READ_WITH_HITMASK(patsrc, true);
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
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        seed,                  // seed
	        &os,
	        false);                // considerQuals
	GreedyDFSRangeSource bt1(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        seed,                  // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for hits for cases 1F, 2F, 3F
	GreedyDFSRangeSource btf2(
			&ebwtBw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (no)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        NULL,                  // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+1,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for partial alignments for case 4R
	GreedyDFSRangeSource btr2(
			&ebwtBw, params,
	        qualCutoff,            // qualThresh (none)
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // report partials (up to seedMms mms)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamRc,                 // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+2,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for seedlings for case 4F
	GreedyDFSRangeSource btf3(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh (none)
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // reportPartials (do)
	        true,                  // reportExacts
	        rangeMode,             // reportRanges
	        pamFw,                 // seedlings
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+3,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// GreedyDFSRangeSource to search for hits for case 4R by extending
	// the partial alignments found in Phase 2
	GreedyDFSRangeSource btr3(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
		    seed+4,  // seed
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// The half-and-half GreedyDFSRangeSource
	GreedyDFSRangeSource btr23(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
		    seed+5,  // seed
		    &os,
		    true,    // considerQuals
		    true,    // halfAndHalf
		    !noMaqRound);
	// GreedyDFSRangeSource to search for hits for case 4F by extending
	// the partial alignments found in Phase 3
	GreedyDFSRangeSource btf4(
			&ebwtBw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),  // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
	        seed+6,  // seed
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// Half-and-half GreedyDFSRangeSource for forward read
	GreedyDFSRangeSource btf24(
			&ebwtBw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),  // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        rangeMode,// reportRanges
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
	        seed+7,  // seed
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
		patsrc->reverseRead();
		#include "search_seeded_phase2.c"
		patsrc->reverseRead();
		#include "search_seeded_phase3.c"
		patsrc->reverseRead();
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
#endif

	SWITCH_TO_FW_INDEX();
	{
		// Phase 1: Consider cases 1R and 2R
		const char * msg = "Seeded quality search Phase 1 of 4: ";
		if(seedMms == 0) {
			msg = "Seeded quality search Phase 1 of 2: ";
		}
		Timer _t(cout, msg, timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase1, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase1((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
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
		exit(1);
	}
	seededQualSearch_pamRc = pamRc;
	{
		// Phase 2: Consider cases 1F, 2F and 3F and generate seedlings
		// for case 4R
		const char * msg = "Seeded quality search Phase 2 of 4: ";
		if(seedMms == 0) {
			msg = "Seeded quality search Phase 2 of 2: ";
		}
		Timer _t(cout, msg, timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase2, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase2((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
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
		exit(1);
	}
	seededQualSearch_pamFw = pamFw;
	{
		// Phase 3: Consider cases 3R and 4R and generate seedlings for
		// case 4F
		Timer _t(cout, "Seeded quality search Phase 3 of 4: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase3, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase3((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
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
		Timer _t(cout, "Seeded quality search Phase 4 of 4: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerPhase4, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase4((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
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

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

	SWITCH_TO_FW_INDEX();
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cout, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory();
	}
	{
		// Phase 1: Consider cases 1R and 2R
		Timer _t(cout, "Seeded quality full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pt_attr, seededQualSearchWorkerFull, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerFull((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			int ret;
			if((ret = pthread_join(threads[i], NULL)) != 0) {
				cerr << "Error: pthread_join returned non-zero status: " << ret << endl;
				exit(1);
			}
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
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
			return new FastaPatternSource (qs, false, useSpinlock,
			                               patDumpfile, trim3, trim5,
			                               nsPolicy, forgiveInput,
			                               skipReads);
		case RAW:
			return new RawPatternSource   (qs, false, useSpinlock,
			                               patDumpfile, trim3, trim5,
			                               nsPolicy, skipReads);
		case FASTQ:
			return new FastqPatternSource (qs, false, useSpinlock,
			                               patDumpfile, trim3, trim5,
			                               nsPolicy, forgiveInput,
			                               solexa_quals,
			                               integer_quals, skipReads);
		case CMDLINE:
			return new VectorPatternSource(qs, false, useSpinlock,
			                               patDumpfile, trim3,
			                               trim5, nsPolicy, skipReads);
		case RANDOM:
			return new RandomPatternSource(2000000, lenRandomReads,
			                               useSpinlock, patDumpfile,
			                               seed);
		default: {
			cerr << "Internal error; bad patsrc format: " << format << endl;
			exit(1);
		}
	}
}

template<typename TStr>
static void driver(const char * type,
                   const string& ebwtFileBase,
                   const string& query,
                   const vector<string>& queries,
                   const string& outfile)
{
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

	// If there were any first-mates specified, we will operate in
	// stateful mode
	bool paired = mates1.size() > 0;
	if(paired) stateful = true;
	// Create list of pattern sources for paired reads
	for(size_t i = 0; i < mates1.size(); i++) {
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
	for(size_t i = 0; i < queries.size(); i++) {
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
	assert_gt(patsrcs_a.size(), 0);
	PairedPatternSource patsrc(patsrcs_a, patsrcs_b);

	if(skipSearch) return;
	// Open hit output file
	OutFileBuf *fout;
	if(!outfile.empty()) {
		if(refOut) {
			fout = NULL;
			cerr << "Warning: ignoring alignment output file " << outfile << " because --refout was specified" << endl;
		} else {
			fout = new OutFileBuf(outfile.c_str(), outType == BINARY);
		}
	} else {
		if(outType == BINARY && !refOut) {
			cerr << "Error: Must specify an output file when output mode is binary" << endl;
			exit(1);
		}
		fout = new OutFileBuf();
	}
	// Initialize Ebwt object and read in header
    Ebwt<TStr> ebwt(adjustedEbwtFileBase,
                    true,     // index is for the forward direction
                    /* overriding: */ offRate,
                    /* overriding: */ isaRate,
                    useShmem, // whether to use shared memory
                    verbose,  // whether to be talkative
                    false /*passMemExc*/,
                    sanityCheck);
    Ebwt<TStr>* ebwtBw = NULL;
    // We need the mirror index if mismatches are allowed
    if(mismatches > 0 || maqLike) {
    	ebwtBw = new Ebwt<TStr>(adjustedEbwtFileBase + ".rev",
    	                        false, // index is for the reverse direction
    	                        /* overriding: */ offRate,
    	                        /* overriding: */ isaRate,
    	                        useShmem, // whether to use shared memory
    	                        verbose,  // whether to be talkative
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
		ebwt.loadIntoMemory();
		ebwt.checkOrigs(os, false);
		ebwt.evictFromMemory();
	}
	{
		Timer _t(cout, "Time searching: ", timing);
		// Set up hit sink; if sanityCheck && !os.empty() is true,
		// then instruct the sink to "retain" hits in a vector in
		// memory so that we can easily sanity check them later on
		HitSink *sink;
		vector<string>* refnames = &ebwt.refnames();
		if(noRefNames) refnames = NULL;
		switch(outType) {
			case FULL:
				if(refOut) {
					sink = new VerboseHitSink(ebwt.nPat(), offBase, dumpUnalFaBase, dumpUnalFqBase, dumpMaxFaBase, dumpMaxFqBase, refnames, partitionSz);
				} else {
					sink = new VerboseHitSink(fout, offBase, dumpUnalFaBase, dumpUnalFqBase, dumpMaxFaBase, dumpMaxFqBase, refnames, partitionSz);
				}
				break;
			case CONCISE:
				if(refOut) {
					sink = new ConciseHitSink(ebwt.nPat(), offBase, dumpUnalFaBase, dumpUnalFqBase, dumpMaxFaBase, dumpMaxFqBase, reportOpps, refnames);
				} else {
					sink = new ConciseHitSink(fout, offBase, dumpUnalFaBase, dumpUnalFqBase, dumpMaxFaBase, dumpMaxFqBase, reportOpps, refnames);
				}
				break;
			case BINARY:
				if(refOut) {
					sink = new BinaryHitSink(ebwt.nPat(), offBase, dumpUnalFaBase, dumpUnalFqBase, dumpMaxFaBase, dumpMaxFqBase, refnames);
				} else {
					sink = new BinaryHitSink(fout, offBase, dumpUnalFaBase, dumpUnalFqBase, dumpMaxFaBase, dumpMaxFqBase, refnames);
				}
				break;
			case NONE:
				sink = new StubHitSink();
				break;
			default:
				cerr << "Invalid output type: " << outType << endl;
				exit(1);
		}
#ifdef CHUD_PROFILING
		chudStartRemotePerfMonitor("Bowtie");
#endif
		if(maqLike) {
			if(!fullIndex) {
				seededQualCutoffSearch(seedLen,
									   qualThresh,
									   seedMms,
									   patsrc,
									   *sink,
									   ebwt,    // forward index
									   *ebwtBw, // mirror index (not optional)
									   os);     // references, if available
			} else {
				seededQualCutoffSearchFull(seedLen,
				                           qualThresh,
				                           seedMms,
				                           patsrc,
				                           *sink,
				                           ebwt,    // forward index
				                           *ebwtBw, // mirror index (not optional)
				                           os);     // references, if available
			}
		}
		else if(mismatches > 0) {
			if(mismatches == 1) {
				if(!fullIndex) {
					mismatchSearch(patsrc, *sink, ebwt, *ebwtBw, os);
				} else {
					assert(ebwtBw != NULL);
					mismatchSearchFull(patsrc, *sink, ebwt, *ebwtBw, os);
				}
			} else if(mismatches == 2 || mismatches == 3) {
				if(!fullIndex) {
					twoOrThreeMismatchSearch(patsrc, *sink, ebwt, *ebwtBw, os, mismatches == 2);
				} else {
					twoOrThreeMismatchSearchFull(patsrc, *sink, ebwt, *ebwtBw, os, mismatches == 2);
				}
			} else {
				cerr << "Error: " << mismatches << " is not a supported number of mismatches" << endl;
				exit(1);
			}
		} else {
			// Search without mismatches
			// Note that --fast doesn't make a difference here because
			// we're only loading half of the index anyway
			exactSearch(patsrc, *sink, ebwt, os, paired);
		}
#ifdef CHUD_PROFILING
		chudStopRemotePerfMonitor();
#endif
		// Evict any loaded indexes from memory
		if(ebwt.isInMemory()) {
			ebwt.evictFromMemory();
		}
		if(ebwtBw != NULL) {
			delete ebwtBw;
		}
		if(!quiet) {
			sink->finish(); // end the hits section of the hit file
		}
		if(dumpHHHits != NULL) dumpHHHits->close();
		if(dumpNoHits != NULL) dumpNoHits->close();
		for(size_t i = 0; i < patsrcs_a.size(); i++) {
			assert(patsrcs_a[i] != NULL);
			delete patsrcs_a[i];
		}
		for(size_t i = 0; i < patsrcs_a.size(); i++) {
			if(patsrcs_b[i] != NULL) {
				delete patsrcs_b[i];
			}
		}
		delete sink;
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
		cout << "Source hash: " << EBWT_SEARCH_HASH << endl;
		return 0;
	}
#ifdef CHUD_PROFILING
	chudInitialize();
	chudAcquireRemoteAccess();
#endif
	{
		Timer _t(cout, "Overall time: ", timing);

		// Get index basename
		if(optind >= argc) {
			cerr << "No index, query, or output file specified!" << endl;
			printUsage(cerr);
			return 1;
		}
		ebwtFile = argv[optind++];

		// Get query filename
		if(optind >= argc) {
			if(mates1.size() > 0) {
				query = "";
			} else {
				cerr << "No query or output file specified!" << endl;
				printUsage(cerr);
				return 1;
			}
		} else if (mates1.size() == 0) {
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
			exit(1);
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
	}
#ifdef CHUD_PROFILING
	chudReleaseRemoteAccess();
#endif
#ifdef BOWTIE_PTHREADS
	pthread_exit(NULL);
#else
	return 0;
#endif
}
