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

using namespace std;
using namespace seqan;

static int verbose				= 0; // be talkative
static bool quiet				= false; // print nothing but the alignments
static int sanityCheck			= 0;  // enable expensive sanity checks
static int format				= FASTQ; // default read format is FASTQ
static string origString		= ""; // reference text, or filename(s)
static int revcomp				= 1; // search for reverse complements?
static int seed					= 0; // srandom() seed
static int timing				= 0; // whether to report basic timing data
static bool allHits				= false; // for multihits, report just one
static bool arrowMode			= false; // report SA arrows instead of locs
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
static int maxNs                = 999999; // max # Ns allowed in read
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
static ofstream *dumpUnalignFa  = NULL;  // FASTA file to dump aligned reads to
static ofstream *dumpUnalignFq  = NULL;  // FASTQ file to dump unaligned reads to
static uint32_t khits           = 1;     // number of hits per read; >1 is much slower
static uint32_t mhits           = 0xffffffff; // don't report any hits if there are > mhits
static bool onlyBest			= false; // true -> guarantee alignments from best possible stratum
static bool spanStrata			= false; // true -> don't stop at stratum boundaries
static bool refOut				= false; // if true, alignments go to per-ref files
static bool seedAndExtend		= false; // use seed-and-extend aligner; for metagenomics recruitment
static int partitionSz          = 0;     // output a partitioning key in first field
static bool noMaqRound          = false;

static const char *short_options = "fqbzh?cu:rv:sat3:5:o:e:n:l:w:p:k:m:";

enum {
	ARG_ORIG = 256,
	ARG_SEED,
	ARG_DUMP_PATS,
	ARG_ARROW,
	ARG_CONCISE,
	ARG_SOLEXA_QUALS,
	ARG_MAXBTS,
	ARG_MAXBTS0,
	ARG_MAXBTS1,
	ARG_MAXBTS2,
	ARG_VERBOSE,
	ARG_QUIET,
	ARG_MAXNS,
	ARG_RANDOM_READS,
	ARG_RANDOM_READS_NOSYNC,
	ARG_NOOUT,
	ARG_FAST,
	ARG_UNFA,
	ARG_UNFQ,
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
	ARG_NOMAQROUND
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
	{"offrate",      required_argument, 0,            'o'},
	{"isarate",      required_argument, 0,            ARG_ISARATE},
	{"skipsearch",   no_argument,       &skipSearch,  1},
	{"qsamelen",     no_argument,       &qSameLen,    1},
	{"reportopps",   no_argument,       &reportOpps,  1},
	{"version",      no_argument,       &showVersion, 1},
	{"maq",          no_argument,       &maqLike,     1},
	{"ntoa",         no_argument,       &nsPolicy,    NS_TO_AS},
	{"dumppats",     required_argument, 0,            ARG_DUMP_PATS},
	{"revcomp",      no_argument,       0,            'r'},
	{"maqerr",       required_argument, 0,            'e'},
	{"seedlen",      required_argument, 0,            'l'},
	{"seedmms",      required_argument, 0,            'n'},
	{"help",         no_argument,       0,            'h'},
	{"threads",      required_argument, 0,            'p'},
	{"khits",        required_argument, 0,            'k'},
	{"mhits",        required_argument, 0,            'm'},
	{"best",         no_argument,       0,            ARG_BEST},
	{"nostrata",     no_argument,       0,            ARG_SPANSTRATA},
	{"nomaqround",   no_argument,       0,            ARG_NOMAQROUND},
	{"refidx",       no_argument,       0,            ARG_REFIDX},
	{"arrows",       no_argument,       0,            ARG_ARROW},
	{"maxbts",       required_argument, 0,            ARG_MAXBTS},
	{"maxbts0",      required_argument, 0,            ARG_MAXBTS0},
	{"maxbts1",      required_argument, 0,            ARG_MAXBTS1},
	{"maxbts2",      required_argument, 0,            ARG_MAXBTS2},
	{"maxns",        required_argument, 0,            ARG_MAXNS},
	{"randread",     no_argument,       0,            ARG_RANDOM_READS},
	{"randreadnosync", no_argument,     0,            ARG_RANDOM_READS_NOSYNC},
	{"phased",       no_argument,       0,            'z'},
	{"dumpnohit",    no_argument,       0,            ARG_DUMP_NOHIT},
	{"dumphhhit",    no_argument,       0,            ARG_DUMP_HHHIT},
	{"refout",       no_argument,       0,            ARG_REFOUT},
	{"seedextend",   no_argument,       0,            ARG_SEED_EXTEND},
	{"partition",    required_argument, 0,            ARG_PARTITION},
	{0, 0, 0, 0} // terminator
};

/**
 * Print a summary usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: bowtie [options]* <ebwt_base> <query_in> [<hit_outfile>]" << endl
	    << "  <ebwt_base>        ebwt filename minus trailing .1.ebwt/.2.ebwt" << endl
	    << "  <query_in>         comma-separated list of files containing query reads" << endl
	    << "                     (or the sequences themselves, if -c is specified)" << endl
	    << "  <hit_outfile>      file to write hits to (default: stdout)" << endl
	    << "Options:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    << "  -r                 query input files are raw one-sequence-per-line" << endl
	    //<< "  -m                 query input files are Maq .bfq" << endl
	    //<< "  -x                 query input files are Solexa _seq.txt" << endl
	    << "  -c                 query sequences given on command line (as <query_in>)" << endl
	    << "  -e/--maqerr <int>  max sum of mismatch quals (rounds like maq; default: 70)" << endl
	    << "  -l/--seedlen <int> seed length (default: 28)" << endl
	    << "  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: 2)" << endl
	    << "  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities" << endl
	    << "  -k <int>           report up to <int> good alignments per read (default: 1)" << endl
	    << "  -a/--all           report all alignments per read (much slower than low -k)" << endl
	    << "  -m <int>           suppress all alignments if > <int> exist (def.: no limit)" << endl
	    << "  --best             guarantee reported alignments are at best possible stratum" << endl
	    << "  --nostrata         if reporting >1 alignment, don't quit at stratum boundaries" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads" << endl
#ifdef BOWTIE_PTHREADS
	    << "  -p/--threads <int> number of search threads to launch (default: 1)" << endl
#endif
	    << "  -u/--qupto <int>   stop after the first <int> reads" << endl
	    << "  --unfa <filename>  write unaligned reads to FASTA file with name <filename>" << endl
	    << "  --unfq <filename>  write unaligned reads to FASTQ file with name <filename>" << endl
	    << "  -t/--time          print wall-clock time taken by search phases" << endl
		<< "  -z/--phased        alternate between index halves; slower, but uses 1/2 mem" << endl
		<< "  --solexa-quals     convert quals from solexa (can be < 0) to phred (can't)" << endl
		<< "  --integer-quals    qualities are given as space-separated integers (not ASCII)" << endl
		<< "  --nomaqround       disable Maq-like quality rounding (to nearest 10 <= 30)" << endl
		<< "  --ntoa             Ns in reads become As; default: Ns match nothing" << endl
	    //<< "  --sanity           enable sanity checks (increases runtime and mem usage!)" << endl
	    //<< "  --orig <str>       specify original string (for sanity-checking)" << endl
	    //<< "  --qsamelen         die with error if queries don't all have the same length" << endl
	    //<< "  --reportopps       report # of other potential mapping targets for each hit" << endl
	    //<< "  --arrows           report hits as top/bottom offsets into SA" << endl
	    //<< "  --randomReads      generate random reads; ignore -q/-f/-r and <query_in>" << endl
	    << "  --concise          write hits in concise format" << endl
	    << "  -b/--binout        write hits in binary format (<hit_outfile> not optional)" << endl
	    << "  --refout           write alignments to files refXXXXX.map, 1 map per reference" << endl
	    << "  --refidx           refer to ref. seqs by 0-based index rather than name" << endl
	    << "  --maxbts <int>     max number of backtracks allowed for -n 2/3 (default: 125)" << endl
	    << "  --maxns <int>      skip reads w/ >n no-confidence bases (default: no limit)" << endl
	    //<< "  --dumppats <file>  dump all patterns read to a file" << endl
	    << "  -o/--offrate <int> override offrate of index; must be >= index's offrate" << endl
	    << "  --seed <int>       seed for random number generator" << endl
	    << "  --verbose          verbose output (for debugging)" << endl
	    << "  --quiet            print nothing but the alignments" << endl
	    << "  -h/--help          print detailed description of tool and its options" << endl
	    << "  --version          print version information and quit" << endl
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
	" are legal, and the -k/-a/--best/--nostrata options which define which\n"
	" and how many legal alignments should be reported.\n"
	"\n"
	" Bowtie is designed to be very fast for read sets where a) many of the\n"
	" reads have at least one good, valid alignment, b) many of the reads\n"
	" are relatively high-quality, c) the number of alignments reported per\n"
	" read is small (close to 1).  These criteria are generally satisfied in\n"
	" the context of mammalian resequencing projects, but you may observe\n"
	" longer running times in other contexts.\n"
	"\n"
	" By default, Bowtie enforces a policy that is equivalent to Maq's\n"
	" quality-aware policy (http://maq.sf.net) (-n 2 -l 28 -e 70), but it\n"
	" can also be made to enforce an end-to-end k-difference policy\n"
	" equivalent to SOAP's (http://soap.genomics.org.cn/) (-v 2).\n"
	"\n"
	" The process by which bowtie chooses an alignment to report is\n"
	" randomized in order to avoid \"mapping bias\" - the phenomenon whereby\n"
	" an aligner systematically fails to report a particular class of good\n"
	" alignments, causing spurious \"holes\" in the comparative assembly.\n"
	" Whenever bowtie reports a subset of the valid alignments that exist,\n"
	" it makes an effort to sample them randomly.  Some bias may still\n"
	" exist.\n"
	"\n"
	" Indels and paired-end alignment are not currently supported.\n"
	" Alignment in ABI \"color space\" is also not currently supported.\n"
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
	"  permitted is configurable with the -V option.\n"
	"  \n"
	"  Command Line\n"
	"  ------------\n"
	"\n"
	"  The following is a detailed description of the options used to control\n"
	"  the 'bowtie' aligner:\n"
	"\n"
	" Usage: bowtie [options]* <ebwt_base> <query_in> [<hit_outfile>]\n"
	"\n"
	"  <ebwt_base>        The basename of the index to be searched.  The\n"
	"                     basename is the name of any of the four index\n"
	"                     files up to but not including the first period.\n"
	"                     bowtie first looks in the current directory for\n"
	"                     the index files, then looks in the 'indexes'\n"
	"                     subdirectory under the directory where the\n"
	"                     currently-running 'bowtie' executable is located,\n"
	"                     then looks in the directory specified in the\n"
	"                     BOWTIE_INDEXES environment variable.\n"
	"\n"
	"  <query_in>         A comma-separated list of files containing the\n"
	"                     reads to be aligned, or, if -c is specified, the\n"
	"                     sequences themselves. E.g., this might be\n"
	"                     \"lane1.fq,lane2.fq,lane3.fq,lane4.fq\", or, if -c\n"
	"                     is specified, this might be \"GGTCATCCT,ACGGGTCGT\"\n"
	"\n"
	"  <hit_outfile>      File to write alignments to.  By default,\n"
	"                     alignments are written to stdout (the console).\n"
	"\n"
	" Options:\n"
	" \n"
	"  -q                 The query input files (specified as <query_in>)\n"
	"                     are FASTQ files (usually having extension .fq or\n"
	"                     .fastq).  This is the default.  See also:\n"
	"                     --solexa-quals.\n"
	"\n"
	"  -f                 The query input files (specified as <query_in>)\n"
	"                     are FASTA files (usually having extension .fa,\n"
	"                     .mfa, .fna or similar).  All quality values are\n"
	"                     assumed to be 40.\n"
	"\n"
	"  -r                 The query input files (specified as <query_in>)\n"
	"                     are Raw files: one sequence per line, without\n"
	"                     quality values or names.\n"
	"\n"
	"  -c                 The query sequences are given on command line.\n"
	"                     I.e. <query_in> is a comma-separated list of\n"
	"                     reads rather than a list of read files.\n"
	"\n"
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
	"  -n/--seedmms <int> The maximum number of mismatches permitted in the\n"
	"                     seed.  This may be 0, 1, 2 or 3 and the default is\n"
	"                     2.\n"
	"\n"
	"  -v <int>           Forego the Maq-like alignment policy and use a\n"
	"                     SOAP-like alignment policy.  I.e., report end-to-\n"
	"                     end alignments with at most <int> mismatches.  If\n"
	"                     -v is specified, quality values and the -e, -l and\n"
	"                     -n options are ignored.\n"
	"\n"
	"  -k <int>           Report up to <int> valid alignments per read\n"
	"                     (default: 1).  Validity of alignments is\n"
	"                     determined by the alignment policy (combined\n"
	"                     effects of -n, -v, -l, and -e).  If many\n"
	"                     alignments are reported, they may be subject to\n"
	"                     stratification; see --best, --nostrata.  Bowtie is\n"
	"                     designed to be very fast for small -k but BOWTIE\n"
	"                     CAN BECOME VERY SLOW AS -k INCREASES.\n"
	"\n"
	"  -a/--all           Report all valid alignments per read (default:\n"
	"                     off).  Validity of alignments is determined by the\n"
	"                     alignment policy (combined effects of -n, -v, -l,\n"
	"                     and -e).  Reported alignments may be subject to\n"
	"                     stratification; see --best, --nostrata.  Bowtie is\n"
	"                     designed to be very fast for small -k; BOWTIE CAN\n"
	"                     CAN BECOME VERY SLOW IF -a/--all IS SPECIFIED.\n"
	"\n"
	"  -m <int>           Suppress all alignments for a particular read if\n"
	"                     more than <int> reportable alignments exist for\n"
	"                     it.  Reportable alignments are those that would be\n"
	"                     reported given the -n, -v, -l, -e, -k, -a, --best,\n"
	"                     and --nostrata options.  Default: no limit.\n"
	"\n"
	"  --best             Reported alignments must belong to the best\n"
	"                     possible alignment \"stratum\" (default: off).  A\n"
	"                     stratum is a category defined by the number of\n"
	"                     mismatches present in the alignment (for -n, the\n"
	"                     number of mismatches present in the seed region of\n"
	"                     the alignment).  E.g., if --best is not specified,\n"
	"                     Bowtie may sometimes report an alignment with 2\n"
	"                     mismatches in the seed even though there exists an\n"
	"                     unreported alignment with 1 mismatch in the seed.\n"
	"                     BOWTIE IS ABOUT 3-5 TIMES SLOWER WHEN --best IS\n"
	"                     SPECIFIED.\n"
	"\n"
	"  --nostrata         If many valid alignments exist and are reportable\n"
	"                     (according to the --best and -k options) and they\n"
	"                     fall into various alignment \"strata\", report all\n"
	"                     of them.  By default, Bowtie only reports those\n"
	"                     alignments that fall into the best stratum, i.e.,\n"
	"                     the one with fewest mismatches.  BOWTIE CAN BECOME\n"
	"                     VERY SLOW WHEN --nostrata IS COMBINED WITH -k OR\n"
	"                     -a. \n"
	"\n"
	"  -5/--trim5 <int>   Trim <int> bases from high-quality (left) end of\n"
	"                     each read before alignment (default: 0).\n"
	"\n"
	"  -3/--trim3 <int>   Trim <int> bases from low-quality (right) end of\n"
	"                     each read before alignment (default: 0).\n"
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
	"  -u/--qupto <int>   Only align the first <int> reads from the\n"
	"                     specified read set.  Default: no limit.\n"
	"\n"
	"  --unfa <filename>  Write all reads that fail to align to a FASTA file\n"
	"                     with name <filename>.\n"
	"\n"
	"  --unfq <filename>  Write all reads that fail to align to a FASTQ file\n"
	"                     with name <filename>.\n"
	"\n"
	"  -t/--time          Print the amount of wall-clock time taken by each\n"
	"                     search phase and index turnover.\n"
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
	"  --nomaqround       Maq accepts quality values in the Phred scale, but\n"
	"                     internally rounds quality values to the nearest 10\n"
	"                     saturating at 30.  By default, Bowtie imitates\n"
	"                     this behavior.  Use --nomaqround to cancel\n"
	"                     rounding in Bowtie.\n"
	"\n"
	"  --ntoa             No-confidence bases in reads (usually 'N' or '.')\n"
	"                     are converted to As before alignment.  By default,\n"
	"                     no-confidence bases do not match any base. \n"
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
	"  -b/--binout        Outout alignments in a concise binary format.  If\n"
	"                     this is specified, <hit_outfile> must also be\n"
	"                     specified. \n"
	"\n"
	"  --refout           Write alignments to a set of files named\n"
	"                     refXXXXX.map, where XXXXX is the 0-padded index of\n"
	"                     the reference sequence aligned to.  This can be a\n"
	"                     useful way to break up work for downstream\n"
	"                     analyses when dealing with, for example, large\n"
	"                     numbers of reads aligned to the assembled human\n"
	"                     genome.  If <hit_outfile> is also specified, it\n"
	"                     will be ignored.\n"
	"\n"
	"  --refidx           When a reference sequence is referred to in a\n"
	"                     reported alignment, refer to it by 0-based index\n"
	"                     (its offset into the list of references that were\n"
	"                     indexed) rather than by name.\n"
	"\n"
	"  --maxbts           The maximum number of backtracks permitted when\n"
	"                     aligning a read in -n 2 or -n 3 mode (default:\n"
	"                     125).  A \"backtrack\" is the introduction of a\n"
	"                     speculative substitution into the alignment.\n"
	"                     Without this limit, the default paramters will\n"
	"                     sometimes require that 'bowtie' try 100s or 1,000s\n"
	"                     of backtracks to align a read, especially if the\n"
	"                     read has many low-quality bases and/or has no\n"
	"                     valid alignments, slowing bowtie down\n"
	"                     significantly.  The drawback of having a limit is\n"
	"                     that some valid alignments may be missed.  Higher\n"
	"                     limits yield greater sensitivity at the expensive\n"
	"                     of longer running times.  \n"
	"\n"
	"  --maxns <int>      Skip reads with more than <int> N's (no-confidence\n"
	"                     bases) in their sequence.  Default: no limit.\n"
	"\n"
	"  -o/--offrate <int> Override the offrate of the index with <int>.  If\n"
	"                     <int> is greater than the offrate used to build\n"
	"                     the index, then some row markings are discarded\n"
	"                     when the index is read into memory.  This reduces\n"
	"                     the memory footprint of the aligner but requires\n"
	"                     more time to calculate text offsets.  <int> must\n"
	"                     be greater than the value used to build the index.\n"
	"\n"
	"  --seed <int>       Use <int> as the seed for pseudo-random number\n"
	"                     generator.\n"
	"\n"
	"  --verbose          Print verbose output (for debugging).\n"
	"\n"
	"  -h/--help          Print detailed description of tool and its options\n"
	"                     (from MANUAL).\n"
	"\n"
	"  --version          Print version information and quit.\n"
	"\n"
	"  Output\n"
	"  ------\n"
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
	   		case 'f': format = FASTA; break;
	   		case 'q': format = FASTQ; break;
	   		case 'r': format = RAW; break;
	   		case 'c': format = CMDLINE; break;
	   		case ARG_RANDOM_READS: format = RANDOM; break;
	   		case ARG_RANDOM_READS_NOSYNC:
	   			format = RANDOM;
	   			randReadsNoSync = true;
	   			break;
	   		case ARG_ARROW: arrowMode = true; break;
	   		case ARG_CONCISE: outType = CONCISE; break;
	   		case 'b': outType = BINARY; break;
	   		case ARG_REFOUT: refOut = true; break;
	   		case ARG_SEED_EXTEND: seedAndExtend = true; break;
	   		case ARG_NOOUT: outType = NONE; break;
	   		case ARG_DUMP_NOHIT: dumpNoHits = new ofstream(".nohits.dump"); break;
	   		case ARG_DUMP_HHHIT: dumpHHHits = new ofstream(".hhhits.dump"); break;
	   		case ARG_UNFA: {
	   			dumpUnalignFa = new ofstream(optarg, ios_base::out); break;
	   			if(!dumpUnalignFa->good()) {
	   				cerr << "Could not open unaligned FASTA file " << optarg << " for appending" << endl;
	   				exit(1);
	   			}
	   			break;
	   		}
	   		case ARG_UNFQ: {
	   			dumpUnalignFq = new ofstream(optarg, ios_base::out); break;
	   			if(!dumpUnalignFq->good()) {
	   				cerr << "Could not open unaligned FASTQ file " << optarg << " for appending" << endl;
	   				exit(1);
	   			}
	   			break;
	   		}
			case ARG_SOLEXA_QUALS: solexa_quals = true; break;
			case ARG_INTEGER_QUALS: integer_quals = true; break;
			case ARG_NOMAQROUND: noMaqRound = true; break;
			case 'z': fullIndex = false; break;
			case ARG_REFIDX: noRefNames = true; break;
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
	   		case 'p':
#ifndef BOWTIE_PTHREADS
	   			cerr << "-p/--threads is disabled because bowtie was not compiled with pthreads support" << endl;
	   			exit(1);
#endif
	   			nthreads = parseInt(1, "-p/--threads arg must be at least 1");
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
	   		case 'n': seedMms = parseInt(0, "-n/--seedmms arg must be at least 0"); break;
	   		case 'l': seedLen = parseInt(20, "-l/--seedlen arg must be at least 20"); break;
	   		case 'h': printLongUsage(cout); exit(0); break;
	   		case '?': printUsage(cerr); exit(1); break;
	   		case ARG_MAXNS: maxNs = parseInt(0, "--maxns arg must be at least 0"); break;
	   		case 'a': allHits = true; break;
	   		case ARG_BEST: onlyBest = true; break;
	   		case ARG_SPANSTRATA: spanStrata = true; break;
	   		case ARG_VERBOSE: verbose = true; break;
	   		case ARG_QUIET: quiet = true; break;
	   		case ARG_SANITY: sanityCheck = true; break;
	   		case 't': timing = true; break;
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
	if(maqLike) {
		revcomp = true;
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
		if(error) exit(1);
	}
}

static char *argv0 = NULL;

#define FINISH_READ(p) \
	/* Don't do finishRead if the read isn't legit or if the read was skipped by the doneMask */ \
	if(!p->empty()) { \
		sink->finishRead(*p, skipped ? NULL : dumpUnalignFa, skipped ? NULL : dumpUnalignFq); \
	} \
	skipped = false;

#define FINISH_READ_WITH_HITMASK(p, r) \
	/* Don't do finishRead if the read isn't legit or if the read was skipped by the doneMask */ \
	if(!p->empty()) { \
		bool reportUnAl = r; \
		if(reportUnAl) { \
			reportUnAl = !skipped; \
			if(reportUnAl) { \
				reportUnAl = !hitMask.test(p->patid()); \
			} \
		} \
		if(sink->finishRead(*p, reportUnAl ? dumpUnalignFa : NULL, reportUnAl ? dumpUnalignFq : NULL) > 0) { \
			if(!reportUnAl && (dumpUnalignFa != NULL || dumpUnalignFq != NULL)) \
			hitMask.setOver(p->patid()); \
		} \
	} \
	skipped = false;

/// Macro for getting the next read, possibly aborting depending on
/// whether the result is empty or the patid exceeds the limit, and
/// marshaling the read into convenient variables.
#define GET_READ(p) \
	p->nextRead(); \
	if(p->empty() || p->patid() >= qUpto) { break; } \
	assert(!empty(p->patFw())); \
	String<Dna5>& patFw  = p->patFw();  \
	String<Dna5>& patRc  = p->patRc();  \
	String<char>& qualFw = p->qualFw(); \
	String<char>& qualRc = p->qualRc(); \
	String<char>& name   = p->name(); \
	uint32_t      patid  = p->patid(); \
	params.setPatId(patid); \
	if(lastLen == 0) lastLen = length(patFw); \
	if(qSameLen && length(patFw) != lastLen) { \
		throw runtime_error("All reads must be the same length"); \
	}

/// Macro for getting the forward oriented version of next read,
/// possibly aborting depending on whether the result is empty or the
/// patid exceeds the limit, and marshaling the read into convenient
/// variables.
#define GET_READ_FW(p) \
	p->nextRead(); \
	if(p->empty() || p->patid() >= qUpto) break; \
	params.setPatId(p->patid()); \
	assert(!empty(p->patFw())); \
	String<Dna5>& patFw  = p->patFw();  \
	String<char>& qualFw = p->qualFw(); \
	String<char>& name   = p->name(); \
	uint32_t      patid  = p->patid(); \
	if(lastLen == 0) lastLen = length(patFw); \
	if(qSameLen && length(patFw) != lastLen) { \
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
static PatternSourcePerThread* createPatSrc(PatternSource& _patsrc, int tid) {
	PatternSourcePerThread *patsrc;
	if(randReadsNoSync) {
		patsrc = new RandomPatternSourcePerThread(numRandomReads, lenRandomReads, nthreads, tid, false);
	} else {
		patsrc = new WrappedPatternSourcePerThread(_patsrc);
	}
    assert(patsrc != NULL);
    return patsrc;
}

/// Create a HitSinkPerThread according to the global params and return
/// a pointer to it
static HitSinkPerThread* createSink(HitSink& _sink, bool sanity) {
    HitSinkPerThread *sink = NULL;
    if(spanStrata) {
		if(!allHits) {
			if(onlyBest) {
				// First N best, spanning strata
				sink = new FirstNBestHitSinkPerThread(_sink, khits, mhits, sanity);
			} else {
				// First N good; "good" inherently ignores strata
				sink = new FirstNGoodHitSinkPerThread(_sink, khits, mhits, sanity);
			}
		} else {
			// All hits, spanning strata
			sink = new AllHitSinkPerThread(_sink, mhits, sanity);
		}
    } else {
		if(!allHits) {
			if(onlyBest) {
				// First N best, not spanning strata
				sink = new FirstNBestStratifiedHitSinkPerThread(_sink, khits, mhits, sanity);
			} else {
				// First N good; "good" inherently ignores strata
				sink = new FirstNGoodHitSinkPerThread(_sink, khits, mhits, sanity);
			}
		} else {
			// All hits, not spanning strata
			sink = new AllStratifiedHitSinkPerThread(_sink, mhits, sanity);
		}
    }
    assert(sink != NULL);
    return sink;
}

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static PatternSource*                 exactSearch_patsrc;
static HitSink*                       exactSearch_sink;
static Ebwt<String<Dna> >*            exactSearch_ebwt;
static vector<String<Dna5> >*         exactSearch_os;
static void *exactSearchWorker(void *vp) {
	PatternSource& _patsrc               = *exactSearch_patsrc;
	HitSink& _sink                       = *exactSearch_sink;
	Ebwt<String<Dna> >& ebwt             = *exactSearch_ebwt;
	vector<String<Dna5> >& os            = *exactSearch_os;

	// Global initialization
	bool sanity = sanityCheck && !os.empty();
	// Per-thread initialization
	uint32_t lastLen = 0;
	PatternSourcePerThread *patsrc = createPatSrc(_patsrc, (int)(long)vp);
	HitSinkPerThread* sink = createSink(_sink, sanity);
	EbwtSearchParams<String<Dna> > params(
			*sink,      // HitSink
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        true,       // read is forward
	        true,       // index is forward
	        arrowMode); // arrow mode
	BacktrackManager<String<Dna> > bt(
			&ebwt, params,
	        0xffffffff,     // qualThresh
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
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
    	// Process forward-oriented read
    	bt.setOffs(0, 0, plen, plen, plen, plen);
    	bt.setQuery(&patFw, &qualFw, &name);
    	bool hit = bt.backtrack();
	    // If the forward direction matched exactly, ignore the
	    // reverse complement
	    if(hit) {
	    	continue;
	    }
	    if(!revcomp) continue;
	    // Process reverse-complement read
		params.setFw(false);
		bt.setQuery(&patRc, &qualRc, &name);
	    bt.backtrack();
		params.setFw(true);
    }
	FINISH_READ(patsrc);
    WORKER_EXIT();
}

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
static void exactSearch(PatternSource& _patsrc,
                        HitSink& _sink,
                        Ebwt<String<Dna> >& ebwt,
                        vector<String<Dna5> >& os)
{
	exactSearch_patsrc = &_patsrc;
	exactSearch_sink   = &_sink;
	exactSearch_ebwt   = &ebwt;
	exactSearch_os     = &os;
#ifdef BOWTIE_PTHREADS
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];

	{
		Timer _t(cout, "Time for 0-mismatch search: ", timing);
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, exactSearchWorker, (void *)(long)(i+1));
		}
#endif
		exactSearchWorker((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
	}
#endif
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
static PatternSource*                 mismatchSearch_patsrc;
static HitSink*                       mismatchSearch_sink;
static Ebwt<String<Dna> >*            mismatchSearch_ebwtFw;
static Ebwt<String<Dna> >*            mismatchSearch_ebwtBw;
static vector<String<Dna5> >*         mismatchSearch_os;
static SyncBitset*                    mismatchSearch_doneMask;
static SyncBitset*                    mismatchSearch_hitMask;

static void* mismatchSearchWorkerPhase1(void *vp){
	PatternSource&         _patsrc       = *mismatchSearch_patsrc;
	HitSink&               _sink         = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw        = *mismatchSearch_ebwtFw;
	vector<String<Dna5> >& os            = *mismatchSearch_os;
	SyncBitset&            doneMask      = *mismatchSearch_doneMask;
	SyncBitset&            hitMask       = *mismatchSearch_hitMask;
    bool sanity = sanityCheck && !os.empty() && !arrowMode;
	uint32_t lastLen = 0; // for checking if all reads have same length
	PatternSourcePerThread* patsrc = createPatSrc(_patsrc, (int)(long)vp);
	HitSinkPerThread* sink = createSink(_sink, sanity);
	EbwtSearchParams<String<Dna> > params(
			*sink,      // HitSinkPerThread
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        false,      // read is forward
	        true,       // index is forward
	        arrowMode); // arrow mode
	BacktrackManager<String<Dna> > bt(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
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
	PatternSource&         _patsrc      = *mismatchSearch_patsrc;
	HitSink&               _sink        = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtBw       = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os           = *mismatchSearch_os;
	SyncBitset&            doneMask     = *mismatchSearch_doneMask;
	SyncBitset&            hitMask      = *mismatchSearch_hitMask;
    // Per-thread initialization
    bool sanity = sanityCheck && !os.empty() && !arrowMode;
	uint32_t lastLen = 0; // for checking if all reads have same length
	PatternSourcePerThread* patsrc = createPatSrc(_patsrc, (int)(long)vp);
	HitSinkPerThread* sink = createSink(_sink, sanity);
	EbwtSearchParams<String<Dna> > params(
			*sink,      // HitSinkPerThread
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        true,       // read is forward
	        false,      // index is mirror index
	        arrowMode); // arrow mode
	BacktrackManager<String<Dna> > bt(
			&ebwtBw, params,
	        0xffffffff,     // qualThresh
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed,           // seed
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
static void mismatchSearch(PatternSource& _patsrc,
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
	// doneMask is sufficient to know whether a read has an alignment
	if(!allHits && khits == 1) numQs = 0;
	// No need to keep track of which reads are aligned because the
	// user hasn't requested an unaligned-read dump
	if(dumpUnalignFa == NULL && dumpUnalignFq == NULL) numQs = 0;
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

	assert(ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

    _patsrc.setReverse(false); // don't reverse patterns

	// Phase 1
    {
		Timer _t(cout, "Time for 1-mismatch Phase 1 of 2: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, mismatchSearchWorkerPhase1, (void *)(long)(i+1));
		}
#endif
		mismatchSearchWorkerPhase1((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
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
			pthread_create(&threads[i], &pthread_custom_attr, mismatchSearchWorkerPhase2, (void *)(long)(i+1));
		}
#endif
		mismatchSearchWorkerPhase2((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
}

static void* mismatchSearchWorkerFull(void *vp){
	PatternSource&         _patsrc      = *mismatchSearch_patsrc;
	HitSink&               _sink        = *mismatchSearch_sink;
	Ebwt<String<Dna> >&    ebwtFw       = *mismatchSearch_ebwtFw;
	Ebwt<String<Dna> >&    ebwtBw       = *mismatchSearch_ebwtBw;
	vector<String<Dna5> >& os           = *mismatchSearch_os;
    // Per-thread initialization
    bool sanity = sanityCheck && !os.empty() && !arrowMode;
	uint32_t lastLen = 0; // for checking if all reads have same length
	PatternSourcePerThread* patsrc = createPatSrc(_patsrc, (int)(long)vp);
	HitSinkPerThread* sink = createSink(_sink, sanity);
	EbwtSearchParams<String<Dna> > params(
			*sink,      // HitSinkPerThread
	        os,         // reference sequences
	        revcomp,    // forward AND reverse complement?
	        true,       // read is forward
	        false,      // index is mirror index
	        arrowMode); // arrow mode
	BacktrackManager<String<Dna> > bt(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        BacktrackLimits(), // max backtracks (no max)
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
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
static void mismatchSearchFull(PatternSource& _patsrc,
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

	assert(ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cout, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory();
	}

#ifdef BOWTIE_PTHREADS
	// Allocate structures for threads
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

    _patsrc.setReverse(false); // don't reverse patterns
    {
		Timer _t(cout, "Time for 1-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, mismatchSearchWorkerFull, (void *)(long)(i+1));
		}
#endif
		// Go to town
		mismatchSearchWorkerFull((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
    }
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
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
		BacktrackManager<String<Dna> >::naiveOracle( \
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
			BacktrackManager<String<Dna> >::printHit( \
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
		BacktrackManager<String<Dna> >::naiveOracle( \
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
			BacktrackManager<String<Dna> >::printHit( \
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

static PatternSource*                 twoOrThreeMismatchSearch_patsrc;
static HitSink*                       twoOrThreeMismatchSearch_sink;
static Ebwt<String<Dna> >*            twoOrThreeMismatchSearch_ebwtFw;
static Ebwt<String<Dna> >*            twoOrThreeMismatchSearch_ebwtBw;
static vector<String<Dna5> >*         twoOrThreeMismatchSearch_os;
static SyncBitset*                    twoOrThreeMismatchSearch_doneMask;
static SyncBitset*                    twoOrThreeMismatchSearch_hitMask;
static bool                           twoOrThreeMismatchSearch_two;

#define TWOTHREE_WORKER_SETUP() \
	PatternSource&                 _patsrc  = *twoOrThreeMismatchSearch_patsrc;   \
	HitSink&                       _sink    = *twoOrThreeMismatchSearch_sink;     \
	vector<String<Dna5> >&         os       = *twoOrThreeMismatchSearch_os;       \
	bool                           two      = twoOrThreeMismatchSearch_two; \
	uint32_t lastLen = 0; \
	PatternSourcePerThread* patsrc = createPatSrc(_patsrc, (int)(long)vp); \
	HitSinkPerThread* sink = createSink(_sink, false); \
	/* Per-thread initialization */ \
	EbwtSearchParams<String<Dna> > params( \
			*sink,       /* HitSink */ \
	        os,          /* reference sequences */ \
	        revcomp,     /* forward AND reverse complement? */ \
	        true,        /* read is forward */ \
	        true,        /* index is forward */ \
	        arrowMode);  /* arrow mode (irrelevant here) */

static void* twoOrThreeMismatchSearchWorkerPhase1(void *vp) {
	TWOTHREE_WORKER_SETUP();
	SyncBitset& doneMask = *twoOrThreeMismatchSearch_doneMask;
	SyncBitset& hitMask = *twoOrThreeMismatchSearch_hitMask;
	Ebwt<String<Dna> >& ebwtFw = *twoOrThreeMismatchSearch_ebwtFw;
	BacktrackManager<String<Dna> > btr1(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
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
	BacktrackManager<String<Dna> > bt2(
			&ebwtBw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (no)
	        true,           // reportExacts
	        false,          // reportArrows
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
	// BacktrackManager to search for seedlings for case 4F
	BacktrackManager<String<Dna> > bt3(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh (none)
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+3,         // seed
		    &os,
		    false);         // considerQuals
	BacktrackManager<String<Dna> > bthh3(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
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
        PatternSource& _patsrc,         /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of mirror text
        vector<String<Dna5> >& os,      /// text strings, if available (empty otherwise)
        bool two = true)                /// true -> 2, false -> 3
{
	// Global initialization
	assert(revcomp);
	assert(!fullIndex);
	assert(ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());

	uint32_t numQs = ((qUpto == 0xffffffff) ? 16 * 1024 * 1024 : qUpto);
	SyncBitset doneMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the read mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");
	// doneMask is sufficient to know whether a read has an alignment
	if(!allHits && khits == 1) numQs = 0;
	// No need to keep track of which reads are aligned because the
	// user hasn't requested an unaligned-read dump
	if(dumpUnalignFa == NULL && dumpUnalignFq == NULL) numQs = 0;
	SyncBitset hitMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the hit mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");

	uint32_t numPats = 0;

	twoOrThreeMismatchSearch_patsrc   = &_patsrc;
	twoOrThreeMismatchSearch_sink     = &_sink;
	twoOrThreeMismatchSearch_ebwtFw   = &ebwtFw;
	twoOrThreeMismatchSearch_ebwtBw   = &ebwtBw;
	twoOrThreeMismatchSearch_os       = &os;
	twoOrThreeMismatchSearch_doneMask = &doneMask;
	twoOrThreeMismatchSearch_hitMask  = &hitMask;
	twoOrThreeMismatchSearch_two      = two;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

    { // Phase 1
		Timer _t(cout, "End-to-end 2/3-mismatch Phase 1 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, twoOrThreeMismatchSearchWorkerPhase1, (void *)(long)(i+1));
		}
#endif
		twoOrThreeMismatchSearchWorkerPhase1((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 1
	    numPats = _patsrc.patid();
    }
	// Unload forward index and load mirror index
	SWITCH_TO_BW_INDEX();
	{ // Phase 2
		Timer _t(cout, "End-to-end 2/3-mismatch Phase 2 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, twoOrThreeMismatchSearchWorkerPhase2, (void *)(long)(i+1));
		}
#endif
		twoOrThreeMismatchSearchWorkerPhase2((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 2
	    assert_eq(numPats, _patsrc.patid());
	}
	SWITCH_TO_FW_INDEX();
	{ // Phase 3
		Timer _t(cout, "End-to-end 2/3-mismatch Phase 3 of 3: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, twoOrThreeMismatchSearchWorkerPhase3, (void *)(long)(i+1));
		}
#endif
		twoOrThreeMismatchSearchWorkerPhase3((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 3
	    assert_eq(numPats, _patsrc.patid());
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
	return;
}

static void* twoOrThreeMismatchSearchWorkerFull(void *vp) {
	TWOTHREE_WORKER_SETUP();
	Ebwt<String<Dna> >& ebwtFw = *twoOrThreeMismatchSearch_ebwtFw;
	Ebwt<String<Dna> >& ebwtBw = *twoOrThreeMismatchSearch_ebwtBw;
	BacktrackManager<String<Dna> > btr1(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
	        NULL,           // seedlings
	        NULL,           // mutations
	        verbose,        // verbose
	        seed,           // seed
	        &os,
	        false);         // considerQuals
	BacktrackManager<String<Dna> > bt2(
			&ebwtBw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (no)
	        true,           // reportExacts
	        false,          // reportArrows
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+1,         // seed
		    &os,
		    false);         // considerQuals
	BacktrackManager<String<Dna> > bt3(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh (none)
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
	        NULL,           // seedlings
		    NULL,           // mutations
	        verbose,        // verbose
		    seed+3,         // seed
		    &os,
		    false);         // considerQuals
	BacktrackManager<String<Dna> > bthh3(
			&ebwtFw, params,
	        0xffffffff,     // qualThresh
	        //BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),
	        // Do not impose maximums in 2/3-mismatch mode
	        BacktrackLimits(), // max backtracks
	        0,              // reportPartials (don't)
	        true,           // reportExacts
	        false,          // reportArrows
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
        PatternSource& _patsrc,         /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of mirror text
        vector<String<Dna5> >& os,      /// text strings, if available (empty otherwise)
        bool two = true)                /// true -> 2, false -> 3
{
	// Global initialization
	assert(revcomp);
	assert(fullIndex);
	assert(ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	{
		// Load the other half of the index into memory
		Timer _t(cout, "Time loading mirror index: ", timing);
		ebwtBw.loadIntoMemory();
	}
	twoOrThreeMismatchSearch_patsrc   = &_patsrc;
	twoOrThreeMismatchSearch_sink     = &_sink;
	twoOrThreeMismatchSearch_ebwtFw   = &ebwtFw;
	twoOrThreeMismatchSearch_ebwtBw   = &ebwtBw;
	twoOrThreeMismatchSearch_os       = &os;
	twoOrThreeMismatchSearch_doneMask = NULL;
	twoOrThreeMismatchSearch_hitMask  = NULL;
	twoOrThreeMismatchSearch_two      = two;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

    {
		Timer _t(cout, "End-to-end 2/3-mismatch full-index search: ", timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, twoOrThreeMismatchSearchWorkerFull, (void *)(long)(i+1));
		}
#endif
		twoOrThreeMismatchSearchWorkerFull((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
    }
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
	return;
}

static PatternSource*                 seededQualSearch_patsrc;
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
	PatternSource&                 _patsrc  = *seededQualSearch_patsrc;   \
	HitSink&                       _sink    = *seededQualSearch_sink;     \
	vector<String<Dna5> >&         os       = *seededQualSearch_os;       \
	int                          qualCutoff = seededQualSearch_qualCutoff; \
	uint32_t lastLen = 0; \
	PatternSourcePerThread* patsrc = createPatSrc(_patsrc, (int)(long)vp); \
	HitSinkPerThread* sink = createSink(_sink, false); \
	/* Per-thread initialization */ \
	EbwtSearchParams<String<Dna> > params( \
			*sink,       /* HitSink */ \
	        os,          /* reference sequences */ \
	        revcomp,     /* forward AND reverse complement? */ \
	        true,        /* read is forward */ \
	        true,        /* index is forward */ \
	        arrowMode);  /* arrow mode (irrelevant here) */

static void* seededQualSearchWorkerPhase1(void *vp) {
	SEEDEDQUAL_WORKER_SETUP();
	SyncBitset& doneMask = *seededQualSearch_doneMask;
	SyncBitset& hitMask = *seededQualSearch_hitMask;
	Ebwt<String<Dna> >& ebwtFw = *seededQualSearch_ebwtFw;
	uint32_t s = seedLen;
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	// BacktrackManager for finding exact hits for the forward-
	// oriented read
	BacktrackManager<String<Dna> > btf1(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        NULL,                  // partials
	        NULL,                  // mutations
	        verbose,               // verbose
	        seed,                  // seed
	        &os,
	        false);                // considerQuals
	BacktrackManager<String<Dna> > bt1(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        false,                 // reportArrows
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
	SyncBitset& hitMask = *seededQualSearch_hitMask;
	uint32_t s = seedLen;
	uint32_t s3 = s >> 1; /* length of 3' half of seed */
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	Ebwt<String<Dna> >& ebwtBw = *seededQualSearch_ebwtBw;
	PartialAlignmentManager* pamRc = seededQualSearch_pamRc;
	// BacktrackManager to search for hits for cases 1F, 2F, 3F
	BacktrackManager<String<Dna> > btf2(
			&ebwtBw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (no)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        NULL,                  // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+1,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// BacktrackManager to search for partial alignments for case 4R
	BacktrackManager<String<Dna> > btr2(
			&ebwtBw, params,
	        qualCutoff,            // qualThresh (none)
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // report partials (up to seedMms mms)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        pamRc,                 // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+2,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
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
		#include "search_seeded_phase2.c"
		#undef DONEMASK_SET
    }
	FINISH_READ_WITH_HITMASK(patsrc, false);
	WORKER_EXIT();
}

static void* seededQualSearchWorkerPhase3(void *vp) {
	SEEDEDQUAL_WORKER_SETUP();
	SyncBitset& doneMask = *seededQualSearch_doneMask;
	SyncBitset& hitMask = *seededQualSearch_hitMask;
	uint32_t s = seedLen;
	uint32_t s3 = s >> 1; /* length of 3' half of seed */
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	Ebwt<String<Dna> >& ebwtFw        = *seededQualSearch_ebwtFw;
	PartialAlignmentManager* pamFw    = seededQualSearch_pamFw;
	PartialAlignmentManager* pamRc    = seededQualSearch_pamRc;
	// BacktrackManager to search for seedlings for case 4F
	BacktrackManager<String<Dna> > btf3(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh (none)
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // reportPartials (do)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        pamFw,                 // seedlings
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+3,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// BacktrackManager to search for hits for case 4R by extending
	// the partial alignments found in Phase 2
	BacktrackManager<String<Dna> > btr3(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
		    seed+4,  // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// The half-and-half BacktrackManager
	BacktrackManager<String<Dna> > btr23(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
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
	SyncBitset& hitMask = *seededQualSearch_hitMask;
	uint32_t s = seedLen;
	uint32_t s5 = (s >> 1) + (s & 1); /* length of 5' half of seed */
	Ebwt<String<Dna> >& ebwtBw = *seededQualSearch_ebwtBw;
	PartialAlignmentManager* pamFw = seededQualSearch_pamFw;
	// BacktrackManager to search for hits for case 4F by extending
	// the partial alignments found in Phase 3
	BacktrackManager<String<Dna> > btf4(
			&ebwtBw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
	        seed+6,  // seed
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// Half-and-half BacktrackManager for forward read
	BacktrackManager<String<Dna> > btf24(
			&ebwtBw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
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
		size_t plen = length(patFw);
		uint32_t qs = min<uint32_t>(plen, s);
		uint32_t qs5 = (qs >> 1) + (qs & 1);
		if(doneMask.test(patid)) { skipped = true; continue; }
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
	// BacktrackManager for finding exact hits for the forward-
	// oriented read
	BacktrackManager<String<Dna> > btf1(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        NULL,                  // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        seed,                  // seed
	        &os,
	        false);                // considerQuals
	BacktrackManager<String<Dna> > bt1(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (don't)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        NULL,                  // seedlings
	        NULL,                  // mutations
	        verbose,               // verbose
	        seed,                  // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// BacktrackManager to search for hits for cases 1F, 2F, 3F
	BacktrackManager<String<Dna> > btf2(
			&ebwtBw, params,
	        qualCutoff,            // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,                     // reportPartials (no)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        NULL,                  // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+1,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// BacktrackManager to search for partial alignments for case 4R
	BacktrackManager<String<Dna> > btr2(
			&ebwtBw, params,
	        qualCutoff,            // qualThresh (none)
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // report partials (up to seedMms mms)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        pamRc,                 // partial alignment manager
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+2,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// BacktrackManager to search for seedlings for case 4F
	BacktrackManager<String<Dna> > btf3(
			&ebwtFw, params,
	        qualCutoff,            // qualThresh (none)
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        seedMms,               // reportPartials (do)
	        true,                  // reportExacts
	        false,                 // reportArrows
	        pamFw,                 // seedlings
		    NULL,                  // mutations
	        verbose,               // verbose
		    seed+3,                // seed
	        &os,                   // reference sequences
	        true,                  // considerQuals
	        false, !noMaqRound);
	// BacktrackManager to search for hits for case 4R by extending
	// the partial alignments found in Phase 2
	BacktrackManager<String<Dna> > btr3(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
		    seed+4,  // seed
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// The half-and-half BacktrackManager
	BacktrackManager<String<Dna> > btr23(
			&ebwtFw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
		    seed+5,  // seed
		    &os,
		    true,    // considerQuals
		    true,    // halfAndHalf
		    !noMaqRound);
	// BacktrackManager to search for hits for case 4F by extending
	// the partial alignments found in Phase 3
	BacktrackManager<String<Dna> > btf4(
			&ebwtBw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),  // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
	        seed+6,  // seed
	        &os,     // reference sequences
	        true,    // considerQuals
	        false, !noMaqRound);
	// Half-and-half BacktrackManager for forward read
	BacktrackManager<String<Dna> > btf24(
			&ebwtBw, params,
	        qualCutoff, // qualThresh
	        BacktrackLimits(maxBts, maxBts0, maxBts1, maxBts2),  // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        false,   // reportArrows
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
		int seedLen,                    /// length of seed (not a maq option)
        int qualCutoff,                 /// maximum sum of mismatch qualities
                                        /// like maq map's -e option
                                        /// default: 70
        int seedMms,                    /// max # mismatches allowed in seed
                                        /// (like maq map's -n option)
                                        /// Can only be 1 or 2, default: 1
        PatternSource& _patsrc,         /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of mirror text
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
	// doneMask is sufficient to know whether a read has an alignment
	if(!allHits && khits == 1) numQs = 0;
	// No need to keep track of which reads are aligned because the
	// user hasn't requested an unaligned-read dump
	if(dumpUnalignFa == NULL && dumpUnalignFq == NULL) numQs = 0;
	SyncBitset hitMask(numQs,
		// Error message for if an allocation fails
		"Could not allocate enough memory for the hit mask; please subdivide reads and\n"
		"run bowtie separately on each subset.\n");
	uint32_t numPats;

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
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
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
			pthread_create(&threads[i], &pthread_custom_attr, seededQualSearchWorkerPhase1, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase1((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 1
	    numPats = _patsrc.patid();
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
			pthread_create(&threads[i], &pthread_custom_attr, seededQualSearchWorkerPhase2, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase2((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 2
	    assert_eq(numPats, _patsrc.patid());
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
			pthread_create(&threads[i], &pthread_custom_attr, seededQualSearchWorkerPhase3, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase3((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 3
	    assert_eq(numPats, _patsrc.patid());
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
			pthread_create(&threads[i], &pthread_custom_attr, seededQualSearchWorkerPhase4, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerPhase4((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 4
	    assert_eq(numPats, _patsrc.patid());
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
        PatternSource& _patsrc,         /// pattern source
        HitSink& _sink,                 /// hit sink
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of mirror text
        vector<String<Dna5> >& os)    /// text strings, if available (empty otherwise)
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
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

	SWITCH_TO_FW_INDEX();
	assert(ebwtFw.isInMemory());
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
			pthread_create(&threads[i], &pthread_custom_attr, seededQualSearchWorkerFull, (void *)(long)(i+1));
		}
#endif
		seededQualSearchWorkerFull((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	}
#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
}

static PatternSource*                 seedAndSWExtendSearch_patsrc;
static HitSink*                       seedAndSWExtendSearch_sink;
static Ebwt<String<Dna> >*            seedAndSWExtendSearch_ebwtFw;
static vector<String<Dna5> >*         seedAndSWExtendSearch_os;

#define SEED_SW_EXTEND_WORKER_SETUP() \
	PatternSource&         _patsrc  = *seedAndSWExtendSearch_patsrc;   \
	HitSink&               _sink    = *seedAndSWExtendSearch_sink;     \
	vector<String<Dna5> >& os       = *seedAndSWExtendSearch_os;       \
	uint32_t lastLen = 0; \
	PatternSourcePerThread* patsrc = createPatSrc(_patsrc, (int)(long)vp); \
	allHits = false; khits = 0; onlyBest = false; \
	HitSinkPerThread* sink = new AllHitSinkPerThread(_sink, true); \
	/* Per-thread initialization */ \
	EbwtSearchParams<String<Dna> > params( \
			*sink,       /* HitSink */ \
	        os,          /* reference sequences */ \
	        revcomp,     /* forward AND reverse complement? */ \
	        true,        /* read is forward */ \
	        true,        /* index is forward */ \
	        true);       /* arrow mode (irrelevant here) */

static void* seedAndSWExtendSearchWorkerPhase1(void *vp) {
	SEED_SW_EXTEND_WORKER_SETUP();
	Ebwt<String<Dna> >& ebwtFw = *seedAndSWExtendSearch_ebwtFw;

	// Half-and-half BacktrackManager for forward read
	BacktrackManager<String<Dna> > bt(
			&ebwtFw, params,
	        0xffffffff,        // qualThresh
	        BacktrackLimits(), // max backtracks
	        0,       // reportPartials (don't)
	        true,    // reportExacts
	        true,    // reportArrows
	        NULL,    // seedlings
		    NULL,    // mutations
	        verbose, // verbose
	        seed,    // seed
	        &os,     // orig strings
	        false,   // considerQuals
	        false);  // halfAndHalf
	vector<Hit>& hits = sink->retainedHits();
	bool skipped = false;
    while(true) {
    	FINISH_READ(patsrc);
    	GET_READ(patsrc);
		size_t plen = length(patFw);
		uint32_t s = min<uint32_t>(plen, seedLen);
		patid += 0;
		// First, try exact hits for the forward-oriented read
		params.setFw(true);
		bt.setQuery(&patFw, &qualFw, &name);
		bt.setOffs(0, 0, s, s, s, s);
		bt.backtrack();
		for(size_t i = 0; i < hits.size(); i++) {
			String<Dna5> lhsbuf; String<Dna5> rhsbuf;
			// h.first has 'top' and h.second has 'bot'
			assert_gt(hits[i].h.second, hits[i].h.first);
			for(size_t j = hits[i].h.first; j < hits[i].h.second; j++) {
				ebwtFw.reportReconstruct(patFw, &qualFw, &name,
				                         lhsbuf,           // lbuf
				                         rhsbuf,           // rbuf
				                         NULL,             // mmui32
				                         NULL,             // refcs
				                         0,                // numMms
				                         j,                // i
				                         hits[i].h.first,  // top
				                         hits[i].h.second, // bot
				                         s,                // qlen
				                         0,                // stratum
				                         params,           // params
				                         NULL);            // locus
			}
		}
		hits.clear();
		// First, try exact hits for the forward-oriented read
		params.setFw(false);
		bt.setQuery(&patRc, &qualRc, &name);
		bt.setOffs(0, 0, s, s, s, s);
		bt.backtrack();
		for(size_t i = 0; i < hits.size(); i++) {
			String<Dna5> lhsbuf; String<Dna5> rhsbuf;
			// h.first has 'top' and h.second has 'bot'
			assert_gt(hits[i].h.second, hits[i].h.first);
			for(size_t j = hits[i].h.first; j < hits[i].h.second; j++) {
				ebwtFw.reportReconstruct(patRc, &qualRc, &name,
				                         lhsbuf,           // lbuf
				                         rhsbuf,           // rbuf
				                         NULL,             // mmui32
				                         NULL,             // refcs
				                         0,                // numMms
				                         j,                // i
				                         hits[i].h.first,  // top
				                         hits[i].h.second, // bot
				                         s,                // qlen
				                         0,                // stratum
				                         params,           // params
				                         NULL);            // locus
			}
		}
		hits.clear();
    }
	FINISH_READ(patsrc);
    WORKER_EXIT();
}

template<typename TStr>
static void seedAndSWExtendSearch(
		int seedLen,                /// length of seed (not a maq option)
        int mms,                    /// max # mismatches allowed in seed
                                    /// (like maq map's -n option)
                                    /// Can only be 1 or 2, default: 1
        PatternSource& _patsrc,     /// pattern source
        HitSink& _sink,             /// hit sink
        Ebwt<TStr>& ebwtFw,         /// index of original text
        vector<String<Dna5> >& os)  /// text strings, if available (empty otherwise)
{
	// Global intialization
	assert_eq(0, mms);
	uint32_t numPats;

	seedAndSWExtendSearch_patsrc   = &_patsrc;
	seedAndSWExtendSearch_sink     = &_sink;
	seedAndSWExtendSearch_ebwtFw   = &ebwtFw;
	seedAndSWExtendSearch_os       = &os;

#ifdef BOWTIE_PTHREADS
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	pthread_attr_setdetachstate(&pthread_custom_attr, PTHREAD_CREATE_JOINABLE);
	pthread_t *threads = new pthread_t[nthreads-1];
#endif

	/* Load the forward index into memory if necessary */ \
	if(!ebwtFw.isInMemory()) {
		Timer _t(cout, "Time loading forward index: ", timing);
		ebwtFw.loadIntoMemory();
	}
	assert(ebwtFw.isInMemory());
	_patsrc.reset(); /* rewind pattern source to first pattern */
	_patsrc.setReverse(false); /* tell pattern source not to reverse patterns */
	{
		// Phase 1: Consider cases 1R and 2R
		const char * msg = "Seed-and-SW-extend search Phase 1 of 1: ";
		Timer _t(cout, msg, timing);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_create(&threads[i], &pthread_custom_attr, seedAndSWExtendSearchWorkerPhase1, (void *)(long)(i+1));
		}
#endif
		seedAndSWExtendSearchWorkerPhase1((void*)0L);
#ifdef BOWTIE_PTHREADS
		for(int i = 0; i < nthreads-1; i++) {
			pthread_join(threads[i], NULL);
		}
#endif
	    // Threads join at end of Phase 1
	    numPats = _patsrc.patid();
	}

#ifdef BOWTIE_PTHREADS
	delete[] threads;
#endif
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
	if(sanityCheck && !origString.empty()) {
		// Determine if it's a file by looking at whether it has a FASTA-like
		// extension
		size_t len = origString.length();
		if(len >= 6 && origString.substr(len-6) == ".fasta" ||
		   len >= 4 && origString.substr(len-4) == ".mfa"   ||
		   len >= 4 && origString.substr(len-4) == ".fas"   ||
		   len >= 4 && origString.substr(len-4) == ".fna"   ||
		   len >= 3 && origString.substr(len-3) == ".fa")
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
	string adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, verbose);
	// Seed random number generator
	srand(seed);
	// Create a pattern source for the queries
	PatternSource *patsrc = NULL;
	if(nsPolicy == NS_TO_NS && !maqLike) {
		maxNs = min<int>(maxNs, mismatches);
	}
	switch(format) {
		case FASTA:
			patsrc = new FastaPatternSource (queries, false,
			                                 patDumpfile, trim3, trim5,
			                                 nsPolicy, maxNs);
			break;
		case RAW:
			patsrc = new RawPatternSource   (queries, false,
			                                 patDumpfile, trim3, trim5,
			                                 nsPolicy, maxNs);
			break;
		case FASTQ:
			patsrc = new FastqPatternSource (queries, false,
			                                 patDumpfile, trim3, trim5,
			                                 nsPolicy, solexa_quals,
											 integer_quals, maxNs);
			break;
		case CMDLINE:
			patsrc = new VectorPatternSource(queries, false,
			                                 patDumpfile, trim3,
			                                 trim5, nsPolicy, maxNs);
			break;
		case RANDOM:
			patsrc = new RandomPatternSource(2000000, lenRandomReads, patDumpfile, seed);
			break;
		default: assert(false);
	}
	if(skipSearch) return;
	// Open hit output file
	ostream *fout;
	if(!outfile.empty()) {
		if(refOut) {
			fout = NULL;
			cerr << "Warning: ignoring alignment output file " << outfile << " because --refout was specified" << endl;
		} else {
			fout = new ofstream(outfile.c_str(), ios::binary);
		}
	} else {
		if(outType == BINARY && !refOut) {
			cerr << "Error: Must specify an output file when output mode is binary" << endl;
			exit(1);
		}
		fout = &cout;
	}
	// Initialize Ebwt object and read in header
    Ebwt<TStr> ebwt(adjustedEbwtFileBase,
                    /* overriding: */ offRate,
                    /* overriding: */ isaRate,
                    verbose,
                    false /*passMemExc*/,
                    sanityCheck);
    Ebwt<TStr>* ebwtBw = NULL;
    // We need the mirror index if mismatches are allowed
    if(mismatches > 0 || maqLike) {
    	ebwtBw = new Ebwt<TStr>(adjustedEbwtFileBase + ".rev",
    	                        /* overriding: */ offRate,
    	                        /* overriding: */ isaRate,
    	                        verbose,
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
    // Load rest of (vast majority of) Ebwt into memory
	if(!maqLike) {
		Timer _t(cout, "Time loading forward index: ", timing);
	    ebwt.loadIntoMemory();
	}
	// Sanity-check the restored version of the Ebwt
	if(sanityCheck && !os.empty()) {
		if(maqLike) ebwt.loadIntoMemory();
		ebwt.checkOrigs(os, false);
		if(maqLike) ebwt.evictFromMemory();
	}
    // If sanity-check is enabled and an original text string
    // was specified, sanity-check the Ebwt by confirming that
    // the unpermuted version equals the original.
	// NOTE: Disabled since, with fragments, it's no longer possible to do
	// this straightforwardly with the os vector.  Rather, we need to either
	// split each element of the os vector on Ns, or we need to read the
	// references in differently.  The former seems preferable.
//	if(!maqLike && sanityCheck && !os.empty()) {
//		TStr rs; ebwt.restore(rs);
//		TStr joinedo = Ebwt<TStr>::join(os, ebwt.eh().chunkRate(), seed);
//		assert_leq(length(rs), length(joinedo));
//		assert_geq(length(rs) + ebwt.eh().chunkLen(), length(joinedo));
//		for(size_t i = 0; i < length(rs); i++) {
//			if(rs[i] != joinedo[i]) {
//				cout << "At character " << i << " of " << length(rs) << endl;
//			}
//			assert_eq(rs[i], joinedo[i]);
//		}
//	}
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
					sink = new VerboseHitSink(ebwt.nPat(), refnames, partitionSz);
				} else {
					sink = new VerboseHitSink(*fout, refnames, partitionSz);
				}
				break;
			case CONCISE:
				if(refOut) {
					sink = new ConciseHitSink(ebwt.nPat(), reportOpps, refnames);
				} else {
					sink = new ConciseHitSink(*fout, reportOpps, refnames);
				}
				break;
			case BINARY:
				if(refOut) {
					sink = new BinaryHitSink(ebwt.nPat(), refnames);
				} else {
					sink = new BinaryHitSink(*fout, refnames);
				}
				break;
			case NONE:
				sink = new StubHitSink();
				break;
			default:
				cerr << "Invalid output type: " << outType << endl;
				exit(1);
		}
		if(seedAndExtend) {
			seedAndSWExtendSearch(seedLen,
			                      seedMms,
			                      *patsrc,
			                      *sink,
			                      ebwt,
			                      os);
		}
		else if(maqLike) {
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
				exit(1);
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
			sink->finish(); // end the hits section of the hit file
		}
	    sink->flush();
		if(fout != NULL) {
			((ofstream*)fout)->close();
		}
		if(dumpHHHits != NULL) dumpHHHits->close();
		if(dumpNoHits != NULL) dumpNoHits->close();
		if(dumpUnalignFa != NULL) dumpUnalignFa->close();
		if(dumpUnalignFq != NULL) dumpUnalignFq->close();
		delete patsrc;
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
	Timer _t(cout, "Overall time: ", timing);

	// Get input filename
	if(optind >= argc) {
		cerr << "No input sequence, query, or output file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	ebwtFile = argv[optind++];

	// Get query filename
	if(optind >= argc) {
		cerr << "No query or output file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	query = argv[optind++];

	// Tokenize the list of query files
	tokenize(query, ",", queries);
	if(queries.size() < 1) {
		cerr << "Tokenized query file list was empty!" << endl;
		printUsage(cerr);
		return 1;
	}

	// Get output filename
	if(optind < argc) {
		outfile = argv[optind++];
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
#ifdef BOWTIE_PTHREADS
	pthread_exit(NULL);
#else
	return 0;
#endif
}
