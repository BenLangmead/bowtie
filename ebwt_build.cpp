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
static bool verbose;
static int sanityCheck;
static int format;
static uint32_t bmax;
static uint32_t bmaxMultSqrt;
static uint32_t bmaxDivN;
static int dcv;
static int noDc;
static int entireSA;
static int seed;
static int showVersion;
static bool doubleEbwt;
//   Ebwt parameters
static int32_t lineRate;
static int32_t linesPerSide;
static int32_t offRate;
static int32_t ftabChars;
static int  bigEndian;
static bool nsToAs;
static bool autoMem;
static bool packed;
static bool writeRef;
static bool justRef;
static int reverseType;
bool color;

static void resetOptions() {
	verbose      = true;  // be talkative (default)
	sanityCheck  = 0;     // do slow sanity checks
	format       = FASTA; // input sequence format
	bmax         = 0xffffffff; // max blockwise SA bucket size
	bmaxMultSqrt = 0xffffffff; // same, as multplier of sqrt(n)
	bmaxDivN     = 4;          // same, as divisor of n
	dcv          = 1024;  // bwise SA difference-cover sample sz
	noDc         = 0;     // disable difference-cover sample
	entireSA     = 0;     // 1 = disable blockwise SA
	seed         = 0;     // srandom seed
	showVersion  = 0;     // just print version and quit?
	doubleEbwt   = true;  // build forward and reverse Ebwts
	//   Ebwt parameters
	lineRate     = 6;  // a "line" is 64 bytes
	linesPerSide = 1;  // 1 64-byte line on a side
	offRate      = 5;  // sample 1 out of 32 SA elts
	ftabChars    = 10; // 10 chars in initial lookup table
	bigEndian    = 0;  // little endian
	nsToAs       = false; // convert reference Ns to As prior to indexing
	autoMem      = true;  // automatically adjust memory usage parameters
	packed       = false; //
	writeRef     = true;  // write compact reference to .3.ebwt/.4.ebwt
	justRef      = false; // *just* write compact reference, don't index
	reverseType  = REF_READ_REVERSE_EACH;
	color        = false;
}

// Argument constants for getopts
enum {
	ARG_BMAX = 256,
	ARG_BMAX_MULT,
	ARG_BMAX_DIV,
	ARG_DCV,
	ARG_SEED,
	ARG_CUTOFF,
	ARG_PMAP,
	ARG_NTOA,
	ARG_USAGE,
	ARG_NEW_REVERSE
};

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
	    << "    -C/--color              build a colorspace index" << endl
	    << "    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting" << endl
	    << "    -p/--packed             use packed strings internally; slower, uses less mem" << endl
	    << "    -B                      build both letter- and colorspace indexes" << endl
	    << "    --bmax <int>            max bucket sz for blockwise suffix-array builder" << endl
	    //<< "    --bmaxmultsqrt <int>    max bucket sz as multiple of sqrt(ref len)" << endl
	    << "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)" << endl
	    << "    --dcv <int>             diff-cover period for blockwise (default: 1024)" << endl
	    << "    --nodc                  disable diff-cover (algorithm becomes quadratic)" << endl
	    << "    -r/--noref              don't build .3/.4.ebwt (packed reference) portion" << endl
	    << "    -3/--justref            just build .3/.4.ebwt (packed reference) portion" << endl
	    << "    -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)" << endl
	    << "    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)" << endl
	    << "    --ntoa                  convert Ns in reference to As" << endl
	    //<< "    --big --little          endianness (default: little, this host: "
	    //<< (currentlyBigEndian()? "big":"little") << ")" << endl
	    << "    --seed <int>            seed for random number generator" << endl
	    //<< "    --new-reverse           concatenate then reverse stretches, not vice versa" << endl
	    << "    -q/--quiet              verbose output (for debugging)" << endl
	    << "    -h/--help               print detailed description of tool and its options" << endl
	    << "    --usage                 print this usage message" << endl
	    << "    --version               print version information and quit" << endl
	    ;
}

static const char *short_options = "qraph?nscfl:i:o:t:h:3C";

static struct option long_options[] = {
	{(char*)"quiet",        no_argument,       0,            'q'},
	{(char*)"sanity",       no_argument,       0,            's'},
	{(char*)"packed",       no_argument,       0,            'p'},
	{(char*)"little",       no_argument,       &bigEndian,   0},
	{(char*)"big",          no_argument,       &bigEndian,   1},
	{(char*)"bmax",         required_argument, 0,            ARG_BMAX},
	{(char*)"bmaxmultsqrt", required_argument, 0,            ARG_BMAX_MULT},
	{(char*)"bmaxdivn",     required_argument, 0,            ARG_BMAX_DIV},
	{(char*)"dcv",          required_argument, 0,            ARG_DCV},
	{(char*)"nodc",         no_argument,       &noDc,        1},
	{(char*)"seed",         required_argument, 0,            ARG_SEED},
	{(char*)"entiresa",     no_argument,       &entireSA,    1},
	{(char*)"version",      no_argument,       &showVersion, 1},
	{(char*)"noauto",       no_argument,       0,            'a'},
	{(char*)"noblocks",     required_argument, 0,            'n'},
	{(char*)"linerate",     required_argument, 0,            'l'},
	{(char*)"linesperside", required_argument, 0,            'i'},
	{(char*)"offrate",      required_argument, 0,            'o'},
	{(char*)"ftabchars",    required_argument, 0,            't'},
	{(char*)"help",         no_argument,       0,            'h'},
	{(char*)"ntoa",         no_argument,       0,            ARG_NTOA},
	{(char*)"justref",      no_argument,       0,            '3'},
	{(char*)"noref",        no_argument,       0,            'r'},
	{(char*)"color",        no_argument,       0,            'C'},
	{(char*)"usage",        no_argument,       0,            ARG_USAGE},
	{(char*)"new-reverse",  no_argument,       0,            ARG_NEW_REVERSE},
	{(char*)0, 0, 0, 0} // terminator
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
			throw 1;
		}
		return t;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	throw 1;
	return -1;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, const char **argv) {
	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(
			argc, const_cast<char**>(argv),
			short_options, long_options, &option_index);
		switch (next_option) {
			case 'f': format = FASTA; break;
			case 'c': format = CMDLINE; break;
			case 'p': packed = true; break;
			case 'C': color = true; break;
			case 'l':
				lineRate = parseNumber<int>(3, "-l/--lineRate arg must be at least 3");
				break;
			case 'i':
				linesPerSide = parseNumber<int>(1, "-i/--linesPerSide arg must be at least 1");
				break;
			case 'o':
				offRate = parseNumber<int>(0, "-o/--offRate arg must be at least 0");
				break;
			case '3':
				justRef = true;
				break;
			case 't':
				ftabChars = parseNumber<int>(1, "-t/--ftabChars arg must be at least 1");
				break;
			case 'n':
				// all f-s is used to mean "not set", so put 'e' on end
				bmax = 0xfffffffe;
				break;
			case 'h':
			case ARG_USAGE:
				printUsage(cout);
				throw 0;
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
			case ARG_NTOA: nsToAs = true; break;
			case ARG_NEW_REVERSE: reverseType = REF_READ_REVERSE; break;
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
				printUsage(cerr);
				throw 1;
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
	bool bisulfite = false;
	RefReadInParams refparams(color, reverse ? reverseType : REF_READ_FORWARD, nsToAs, bisulfite);
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
			if (f == NULL) {
				cerr << "Error: could not open "<< infiles[i] << endl;
				throw 1;
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
	vector<uint32_t> plens;
	std::pair<size_t, size_t> sztot;
	{
		if(verbose) cout << "Reading reference sizes" << endl;
		Timer _t(cout, "  Time reading reference sizes: ", verbose);
		if(!reverse && (writeRef || justRef)) {
			// For forward reference, dump it to .3.ebwt and .4.ebwt
			// files
			string file3 = outfile + ".3.ebwt";
			string file4 = outfile + ".4.ebwt";
			// Open output stream for the '.3.ebwt' file which will
			// hold the size records.
			ofstream fout3(file3.c_str(), ios::binary);
			if(!fout3.good()) {
				cerr << "Could not open index file for writing: \"" << file3 << "\"" << endl
					 << "Please make sure the directory exists and that permissions allow writing by" << endl
					 << "Bowtie." << endl;
				throw 1;
			}
			BitpairOutFileBuf bpout(file4.c_str());
			// Read in the sizes of all the unambiguous stretches of
			// the genome into a vector of RefRecords.  The input
			// streams are reset once it's done.
			writeU32(fout3, 1, bigEndian); // endianness sentinel
			if(color) {
				refparams.color = false;
				// Make sure the .3.ebwt and .4.ebwt files contain
				// nucleotides; not colors
				int numSeqs = 0;
				std::pair<size_t, size_t> sztot2 =
					fastaRefReadSizes(is, szs, plens, refparams, &bpout, numSeqs);
				refparams.color = true;
				writeU32(fout3, szs.size(), bigEndian); // write # records
				for(size_t i = 0; i < szs.size(); i++) {
					szs[i].write(fout3, bigEndian);
				}
				szs.clear();
				plens.clear();
				// Now read in the colorspace size records; these are
				// the ones that were indexed
				int numSeqs2 = 0;
				sztot = fastaRefReadSizes(is, szs, plens, refparams, NULL, numSeqs2);
				assert_geq(numSeqs, numSeqs2);
				//assert_eq(sztot2.second, sztot.second + numSeqs);
			} else {
				int numSeqs = 0;
				sztot = fastaRefReadSizes(is, szs, plens, refparams, &bpout, numSeqs);
				writeU32(fout3, szs.size(), bigEndian); // write # records
				for(size_t i = 0; i < szs.size(); i++) szs[i].write(fout3, bigEndian);
			}
			if(sztot.first == 0) {
				cerr << "Error: No unambiguous stretches of characters in the input.  Aborting..." << endl;
				throw 1;
			}
			assert_gt(sztot.first, 0);
			assert_gt(sztot.second, 0);
			bpout.close();
			fout3.close();
#ifndef NDEBUG
			if(sanityCheck) {
				BitPairReference bpr(
					outfile, // ebwt basename
					color,   // expect color?
					true,    // sanity check?
					&infiles,// files to check against
					NULL,    // sequences to check against
					format == CMDLINE, // whether infiles contains strings
					true,    // load sequence?
					false,   // use memory-mapped files
					false,   // use shared memory
					false,   // sweep through memory-mapped memory
					false,   // be talkative
					false);  // be talkative
			}
#endif
		} else {
			// Read in the sizes of all the unambiguous stretches of the
			// genome into a vector of RefRecords
			int numSeqs = 0;
			sztot = fastaRefReadSizes(is, szs, plens, refparams, NULL, numSeqs);
#ifndef NDEBUG
			if(refparams.color) {
				refparams.color = false;
				vector<RefRecord> szs2;
				vector<uint32_t> plens2;
				int numSeqs2 = 0;
				std::pair<size_t, size_t> sztot2 =
					fastaRefReadSizes(is, szs2, plens2, refparams, NULL, numSeqs2);
				assert_leq(numSeqs, numSeqs2);
				// One less color than base
				//assert_geq(sztot2.second, sztot.second + numSeqs);
				refparams.color = true;
			}
#endif
		}
	}
	if(justRef) return;
	assert_gt(sztot.first, 0);
	assert_gt(sztot.second, 0);
	assert_gt(szs.size(), 0);
	// Construct Ebwt from input strings and parameters
	Ebwt<TStr> ebwt(refparams.color ? 1 : 0,
	                lineRate,
	                linesPerSide,
	                offRate,      // suffix-array sampling rate
	                -1,           // ISA sampling rate
	                ftabChars,    // number of chars in initial arrow-pair calc
	                outfile,      // basename for .?.ebwt files
	                !reverse,     // fw
	                !entireSA,    // useBlockwise
	                bmax,         // block size for blockwise SA builder
	                bmaxMultSqrt, // block size as multiplier of sqrt(len)
	                bmaxDivN,     // block size as divisor of len
	                noDc? 0 : dcv,// difference-cover period
	                is,           // list of input streams
	                szs,          // list of reference sizes
	                plens,        // list of not-all-gap reference sequence lengths
	                sztot.first,  // total size of all unambiguous ref chars
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
		ebwt.loadIntoMemory(
			refparams.color ? 1 : 0,
			-1,
			false,
			false);
		TStr s2; ebwt.restore(s2);
		ebwt.evictFromMemory();
		{
			TStr joinedss = Ebwt<TStr>::join(
				is,          // list of input streams
				szs,         // list of reference sizes
				sztot.first, // total size of all unambiguous ref chars
				refparams,   // reference read-in parameters
				seed);       // pseudo-random number generator seed
			if(refparams.reverse == REF_READ_REVERSE) {
				reverseInPlace(joinedss);
			}
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

static const char *argv0 = NULL;

extern "C" {
/**
 * main function.  Parses command-line arguments.
 */
int bowtie_build(int argc, const char **argv) {
	try {
		// Reset all global state, including getopt state
		opterr = optind = 1;
		resetOptions();

		string infile;
		vector<string> infiles;
		string outfile;

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
}
