#ifdef EBWT_SEARCH_MAIN

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
#include "endian.h"
#include "packed_io.h"
#include "ebwt.h"
#include "params.h"
#include "sequence_io.h"
#include "tokenize.h"
#include "hit.h"
#include "pat.h"
#include "inexact_extend.h"

using namespace std;
using namespace seqan;

static int verbose				= 0; // be talkative
static int allowed_diffs		= -1; // 3' difference budget for extension mode
static int kmer					= -1; // 5' kmer length for extension mode
static int sanityCheck			= 0;  // enable expensive sanity checks
static int format				= FASTQ; // default read format is FASTQ
static string origString		= ""; // reference text, or filename(s)
static int revcomp				= 1; // search for reverse complements?
static int seed					= 0; // srandom() seed
static int timing				= 0; // whether to report basic timing data
static bool oneHit				= true;  // for multihits, report just one
static bool concise				= false; // report hits in id+/-:<x,y,z> format
static bool arrowMode			= false; // report SA arrows instead of locs
static int showVersion			= 0; // just print version and quit?
static int ipause				= 0; // pause before maching?
static int binOut				= 0; // write hits in binary
static int qUpto				= -1; // max # of queries to read
static int skipSearch			= 0; // abort before searching
static int qSameLen				= 0; // abort before searching
static int trim5				= 0; // amount to trim from 5' end
static int trim3				= 0; // amount to trim from 3' end
static int printStats			= 0; // whether to print statistics
static int reportOpps			= 0; // whether to report # of other mappings
static int offRate				= -1; // keep default offRate
static int mismatches			= 0; // allow 0 mismatches by default
static char *patDumpfile		= NULL; // filename to dump patterns to
static bool solexa_quals		= false; //quality strings are solexa qualities, instead of phred
static int maqLike				= 1; // do maq-like searching
static int seedLen              = 24; // seed length (not configurable in maq)
static int seedMms              = 2;  // # mismatches allowed in seed (maq's -n)
static int qualThresh           = 70; // max qual-weighted hamming dist (maq's -e)

static const char *short_options = "fqbmlcu:rvsat013:5:o:k:d:e:n:";

#define ARG_ORIG 256
#define ARG_SEED 257
#define ARG_DUMP_PATS 258
#define ARG_ARROW 259
#define ARG_CONCISE 260
#define ARG_SOLEXA_QUALS 261

static struct option long_options[] = {
	{"verbose",      no_argument,       0,            'v'},
	{"sanity",       no_argument,       0,            's'},
	{"exact",        no_argument,       0,            '0'},
	{"1mm",          no_argument,       0,            '1'},
	{"pause",        no_argument,       &ipause,      1},
	{"orig",         required_argument, 0,            ARG_ORIG},
	{"allhits",      no_argument,       0,            'a'},
	{"binout",       no_argument,       0,            'b'},
	{"concise",      no_argument,       0,            ARG_CONCISE},
	{"solexa-quals", no_argument,       0,            ARG_SOLEXA_QUALS},
	{"time",         no_argument,       0,            't'},
	{"trim3",        required_argument, 0,            '3'},
	{"trim5",        required_argument, 0,            '5'},
	{"seed",         required_argument, 0,            ARG_SEED},
	{"qupto",        required_argument, 0,            'u'},
	{"offrate",      required_argument, 0,            'o'},
	{"skipsearch",   no_argument,       &skipSearch,  1},
	{"qsamelen",     no_argument,       &qSameLen,    1},
	{"stats",        no_argument,       &printStats,  1},
	{"reportopps",   no_argument,       &reportOpps,  1},
	{"version",      no_argument,       &showVersion, 1},
	{"maq",          no_argument,       &maqLike,     1},
	{"dumppats",     required_argument, 0,            ARG_DUMP_PATS},
	{"revcomp",      no_argument,       0,            'r'},
	{"err",          required_argument, 0,            'e'},
	{"seedmms",      required_argument, 0,            'n'},
	{"kmer",         required_argument, 0,            'k'},
	{"3prime-diffs", required_argument, 0,            'd'},
	{"arrows",       no_argument,       0,            ARG_ARROW},
	{0, 0, 0, 0} // terminator
};

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
	out << "Usage: ebwt_search [options]* <ebwt_base> <query_in> [<hit_outfile>]" << endl
	    << "  <ebwt_base>        ebwt filename minus trailing .1.ebwt/.2.ebwt" << endl
	    << "  <query_in>         comma-separated list of files containing query reads" << endl
	    << "                     (or the sequences themselves, if -c is specified)" << endl
	    << "  <hit_outfile>      file to write hits to (default: stdout)" << endl
	    << "Options:" << endl
	    << "  -q                 query input files are FASTQ .fq/.fastq (default)" << endl
	    << "  -f                 query input files are (multi-)FASTA .fa/.mfa" << endl
	    //<< "  -m                 query input files are Maq .bfq" << endl
	    //<< "  -l                 query input files are Solexa _seq.txt" << endl
	    << "  -c                 query sequences given on command line (as <query_in>)" << endl
	    << "  -e/--err <int>     max sum of mismatch qualities (default: 70)" << endl
	    << "  -n/--seedmms <int> max mismatches in 5' 24 bps (0, 1 or 2, default: 2)" << endl
	    << "  -0/--exact         report end-to-end exact hits; ignore quals, -e, -n" << endl
	    << "  -1/--1mm           report end-to-end 1-mismatch hits; ignore quals, -e, -n" << endl
	    << "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads" << endl
	    << "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads" << endl
	    << "  -u/--qupto <int>   stop after the first <int> reads" << endl
	    //<< "  --maq              maq-like matching (forces -r, -k 24)" << endl
	    //<< "  -r/--revcomp       also search for rev. comp. of each query (default: off)" << endl
		//<< "  -k/--kmer [int]    match on the 5' #-mer and then extend hits with a more sensitive alignment (default: 22bp)" << endl
		//<< "  -d/--3prime-diffs  # of differences in the 3' end, when used with -k above (default: 4)" << endl
	    //<< "  -b/--binout        write hits in binary format (must specify <hit_outfile>)" << endl
	    << "  -t/--time          print wall-clock time taken by search phases" << endl
		<< "  --solexa-quals     convert FASTQ qualities from solexa-scaled to phred" << endl
	    //<< "  -s/--sanity        enable sanity checks (increases runtime and mem usage!)" << endl
	    //<< "  --orig <str>       specify original string (for sanity-checking)" << endl
	    //<< "  --qsamelen         die with error if queries don't all have the same length" << endl
	    //<< "  --stats            write statistics after hits" << endl
	    //<< "  --reportopps       report # of other potential mapping targets for each hit" << endl
	    //<< "  -a/--allhits       if query has >1 hit, give all hits (default: 1 random hit)" << endl
	    //<< "  --arrows           report hits as top/bottom offsets into SA" << endl
	    << "  --concise          write hits in a concise format" << endl
	    //<< "  --dumppats <file>  dump all patterns read to a file" << endl
	    << "  -o/--offrate <int> override offRate of Ebwt; must be <= value in index" << endl
	    << "  --seed <int>       seed for random number generator" << endl
	    << "  -v/--verbose       verbose output (for debugging)" << endl
	    << "  --version          print version information and quit" << endl
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
	   		case 'm': format = BFQ; break;
	   		case 'l': format = SOLEXA; break;
	   		case 'c': format = CMDLINE; break;
	   		case '0': maqLike = 0; mismatches = 0; break;
	   		case '1': maqLike = 0; mismatches = 1; break;
	   		case 'r': revcomp = 1; break;
	   		case ARG_ARROW: arrowMode = true; break;
	   		case ARG_CONCISE: concise = true; break;
			case ARG_SOLEXA_QUALS: solexa_quals = true; break;
	   		case ARG_SEED:
	   			seed = parseInt(0, "--seed arg must be at least 0");
	   			break;
	   		case 'u':
	   			qUpto = parseInt(1, "-u/--qUpto arg must be at least 1");
	   			break;
	   		case '3':
	   			trim3 = parseInt(0, "-3/--trim3 arg must be at least 0");
	   			break;
	   		case '5':
	   			trim5 = parseInt(0, "-5/--trim5 arg must be at least 0");
	   			break;
	   		case 'o':
	   			offRate = parseInt(1, "-o/--offRate arg must be at least 1");
	   			break;
	   		case 'e':
	   			qualThresh = parseInt(1, "-e/--err arg must be at least 1");
	   			break;
	   		case 'n':
	   			seedMms = parseInt(0, "-n/--seedMms arg must be at least 0");
	   			break;
	   		case 'a': oneHit = false; break;
	   		case 'v': verbose = true; break;
	   		case 's': sanityCheck = true; break;
	   		case 't': timing = true; break;
	   		case 'b': binOut = true; break;
			case 'k': 
				if (optarg != NULL)
					kmer = parseInt(1, "-k/--kmer must be at least 1");
				else
					kmer = 22;
				break;
			case 'd': 
				if (optarg)
					allowed_diffs = parseInt(0, "-d/--3prime-diffs must be at least 0");
				else
					allowed_diffs = 4;
				break;
	   		case ARG_DUMP_PATS:
	   			patDumpfile = optarg;
	   			break;
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
	if(revcomp && qUpto != -1) {
		// If revomp is enabled, double qUpto to reflect the fact that
		// it's a cap on reads, not pattermns
		qUpto <<= 1;
	}
	if(maqLike) {
		revcomp = true;
	}
	if(maqLike && !oneHit) {
		// No support for -a in Maq mode (yet)
		cerr << "Cannot combine -a/--allhits with Maq-like (default) mode
		     << endl
		     << "Either omit -a/--allhits or also specify -0 or -1 for end-to-end mode"
		     << endl;
		exit(1);
	}
}

static char *argv0 = NULL;

/**
 * Search through a single (forward) Ebwt index for exact end-to-end
 * hits.  Assumes that index is already loaded into memory.
 */
template<typename TStr>
static void exactSearch(PatternSource<TStr>& patsrc,
                        HitSink& sink,
                        EbwtSearchStats<TStr>& stats,
                        EbwtSearchParams<TStr>& params,
                        Ebwt<TStr>& ebwt,
                        vector<TStr>& os)
{
	uint32_t patid = 0;
	uint64_t lastHits = 0llu;
	uint32_t lastLen = 0;
	assert(patsrc.hasMorePatterns());
	EbwtSearchState<TStr> s(ebwt, params, seed);
    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
    	params.setFw(!revcomp || !patsrc.nextIsReverseComplement());
    	params.setPatId(patid++);
    	assert(!revcomp || (params.patId() & 1) == 0 || !params.fw());
    	assert(!revcomp || (params.patId() & 1) == 1 ||  params.fw());
    	TStr* pat = NULL;
		String<char>* qual = NULL;
		String<char>* name = NULL;
		patsrc.nextPattern(&pat, &qual, &name);
    	assert(!empty(*pat));
    	if(lastLen == 0) lastLen = length(*pat);
    	if(qSameLen && length(*pat) != lastLen) {
    		throw runtime_error("All reads must be the same length");
    	}
    	s.newQuery(pat, name, qual);
    	params.stats().incRead(s, *pat);
	    ebwt.search(s, params);
	    // If the forward direction matched exactly, ignore the
	    // reverse complement
	    if(oneHit && revcomp && sink.numHits() > lastHits) {
	    	lastHits = sink.numHits();
	    	if(params.fw()) {
	    		assert(patsrc.nextIsReverseComplement());
	    		assert(patsrc.hasMorePatterns());
	    		// Ignore this pattern (the reverse complement of
	    		// the one we just matched)
	        	TStr* pat2;
	    		String<char>* qual2;
	    		String<char>* name2;
				patsrc.nextPattern(&pat2, &qual2, &name2);
		    	assert(!empty(*pat2));
		    	patid++;
		    	if(qSameLen && length(*pat2) != lastLen) {
		    		throw runtime_error("All reads must be the same length");
		    	}
		    	params.setFw(false);
		    	params.stats().incRead(s, *pat2);
	    		assert(!patsrc.nextIsReverseComplement());
	    	}
	    }
    	// Optionally sanity-check results by confirming with a
    	// different matcher that the pattern occurs in exactly
    	// those locations reported.  
    	if(sanityCheck && !oneHit && !arrowMode && !os.empty()) {
    	    vector<Hit>& results = sink.retainedHits();
		    vector<U32Pair> results2;
		    results2.reserve(256);
		    for(unsigned int i = 0; i < os.size(); i++) {
	    		// Forward
	    		Finder<TStr> finder(os[i]);
	    		Pattern<TStr, Horspool> pattern(*pat);
	    		while (find(finder, pattern)) {
	    			results2.push_back(make_pair(i, position(finder)));
	    		}
		    }
    		sort(results.begin(), results.end());
    		if(oneHit) {
	    		assert_leq(results.size(), results2.size());
	    		for(int i = 0; i < (int)results.size(); i++) {
	    			bool foundMatch = false;
		    		for(int j = i; j < (int)results2.size(); j++) {
		    			if(results[i].h.first == results2[j].first &&
		    			   results[i].h.second == results2[j].second)
		    			{
		    				foundMatch = true;
		    				break;
		    			}
		    		}
		    		assert(foundMatch);
	    		}
    		} else {
	    		assert_eq(results.size(), results2.size());
	    		for(int i = 0; i < (int)results.size(); i++) {
	    			assert_eq(results[i].h.first, results2[i].first);
	    			assert_eq(results[i].h.second, results2[i].second);
	    		}
    		}
    		if(verbose) {
    			cout << "Passed orig/result sanity-check ("
    			     << results2.size() << " results checked) for pattern "
    			     << patid << endl;
    		}
    		sink.clearRetainedHits();
    	}
    }
}



/**
 * Search through a single (forward) Ebwt index for exact query hits in the 
 * 5' end of each read, and then extend that hit by shift-and to allow for 3'
 * mismatches.  Assumes that index is already loaded into memory.
 */
template<typename TStr>
static void exactSearchWithExtension(vector<String<Dna, Packed<> > >& packed_texts,
									 PatternSource<TStr>& patsrc,
									 HitSink& sink,
									 EbwtSearchStats<TStr>& stats,
									 EbwtSearchParams<TStr>& params,
									 Ebwt<TStr>& ebwt,
									 vector<TStr>& os)
{
	uint32_t patid = 0;
	uint64_t lastHits= 0llu;
	uint32_t lastLen = 0;
	assert(patsrc.hasMorePatterns());
	
	if (allowed_diffs == -1)
		allowed_diffs = default_allowed_diffs;
	
	ExactSearchWithLowQualityThreePrime<TStr> extend_policy(packed_texts, 
															false, 
															kmer,/* global, override */ 
															allowed_diffs /*global, override*/);
	
    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
    	params.setFw(!revcomp || !patsrc.nextIsReverseComplement());
    	params.setPatId(patid++);
    	assert(!revcomp || (params.patId() & 1) == 0 || !params.fw());
    	assert(!revcomp || (params.patId() & 1) == 1 ||  params.fw());
    	TStr* pat;
		String<char>* qual;
		String<char>* name;
		patsrc.nextPattern(&pat, &qual, &name);
    	assert(!empty(*pat));
    	if(lastLen == 0) lastLen = length(*pat);
    	if(qSameLen && length(*pat) != lastLen) {
    		throw runtime_error("All reads must be the same length");
    	}
		
    	// FIXME: not accumulating stats for search with extension
    	//params.stats().incRead(s, pat);
	    extend_policy.search(ebwt, stats, params, *pat, *name, *qual, sink);
		
	    // If the forward direction matched exactly, ignore the
	    // reverse complement
	    if(oneHit && revcomp && sink.numHits() > lastHits) {
	    	lastHits = sink.numHits();
	    	if(params.fw()) {
	    		assert(patsrc.nextIsReverseComplement());
	    		assert(patsrc.hasMorePatterns());
	    		// Ignore this pattern (the reverse complement of
	    		// the one we just matched)
				
		    	TStr* pat2;
				String<char>* qual2;
				String<char>* name2;
				patsrc.nextPattern(&pat2, &qual2, &name2);
		    	assert(!empty(*pat2));
		    	patid++;
		    	if(qSameLen && length(*pat2) != lastLen) {
		    		throw runtime_error("All reads must be the same length");
		    	}
		    	params.setFw(false);
				
				//FIXME:
		    	//params.stats().incRead(s, pat2);
	    		assert(!patsrc.nextIsReverseComplement());
	    	}
	    }
    }
}


/**
 * Given a pattern, a list of reference texts, and some other state,
 * find all hits for that pattern in all texts using a naive seed- 
 * and-extend algorithm where seeds are found using Horspool.
 */
template<typename TStr1, typename TStr2>
static bool findSanityHits(const TStr1& pat,
                           uint32_t patid,
                           bool fw,
                           vector<TStr2>& os,
                           vector<Hit>& sanityHits,
                           bool allowExact,
                           bool transpose)
{
	typedef TStr1 TStr;
	typedef typename Value<TStr>::Type TVal;
	bool ebwtFw = !transpose;
	bool fivePrimeOnLeft = (ebwtFw == fw);
    uint32_t plen = length(pat);
	TStr half;
	reserve(half, plen);
	uint32_t bump = 0;
	if(!transpose) bump = 1;
	// Grab the unrevisitable region of pat
	for(size_t i = ((plen+bump)>>1); i < plen; i++) {
		appendValue(half, (TVal)pat[i]);
	}
    uint32_t hlen = length(half); // length of seed (right) half
    assert_leq(hlen, plen);
    uint32_t ohlen = plen - hlen; // length of other (left) half
    assert_leq(ohlen, plen);
	Pattern<TStr, Horspool> pattern(half);
	for(size_t i = 0; i < os.size(); i++) {
		TStr2 o = os[i];
		if(transpose) {
			for(size_t j = 0; j < length(o)>>1; j++) {
				TVal tmp = o[j];
				o[j] = o[length(o)-j-1];
				o[length(o)-j-1] = tmp;
			}
		}
		Finder<TStr> finder(o);
		while (find(finder, pattern)) {
			uint32_t pos = position(finder);
			bitset<max_read_bp> diffs = 0;
			if(pos >= ohlen) {
				// Extend toward the left end of the pattern, counting
				// mismatches
				for(uint32_t j = 0; j < ohlen && diffs.count() <= 1; j++) {
					if(o[pos-j-1] != pat[ohlen-j-1]) {
						uint32_t off = ohlen-j-1;
						if(fivePrimeOnLeft) {
							diffs.set(off);
						} else {
							// The 3' end is on on the left end of the
							// pattern, but the diffs vector should
							// encode mismatches w/r/t the 5' end, so
							// we flip
							diffs.set(plen-off-1);
						}
					}
				}
			}
			// If the extend yielded 1 or fewer mismatches, keep it
			if((diffs == 0 && allowExact) || diffs.count() == 1) {
				uint32_t off = pos - ohlen;
				if(transpose) {
					off = length(o) - off;
					off -= length(pat);
				}
				// A hit followed by a transpose can sometimes fall
				// off the beginning of the text
				if(off < (0xffffffff - length(pat))) {
					Hit h(make_pair(i, off), 
						  patid, 
						  "",
						  pat, 
						  "" /*no need for qualities*/, 
						  fw, 
						  diffs);
					sanityHits.push_back(h);
				}
			}
		}
	}
	return true;
}

/**
 * Assert that the sanityHits array has been exhausted, presumably
 * after having been reconciled against actual hits with
 * reconcileHits().  Only used in allHits mode.
 */
template<typename TStr>
static bool checkSanityExhausted(const TStr& pat,
                                 uint32_t patid,
                                 bool fw,
                                 vector<Hit>& sanityHits,
                                 bool transpose)
{
    // If caller specified mustExhaust, then we additionally check
    // whether every sanityHit has now been matched up with some Ebwt
    // hit.  If not, that means that Ebwt may have missed a hit, so
    // we assert.
    size_t unfoundHits = 0;
	for(size_t j = 0; j < sanityHits.size(); j++) {
		uint32_t patid = sanityHits[j].patId;
		bool fw = sanityHits[j].fw;
		cout << "Did not find sanity hit: "
		     << (patid>>revcomp) << (fw? "+":"-")
		     << ":<" << sanityHits[j].h.first << ","
		     << sanityHits[j].h.second << ","
		     << sanityHits[j].mms << ">" << endl;
		cout << "  transpose: " << transpose << endl;
		unfoundHits++;
	}
	assert_eq(0, unfoundHits); // Ebwt missed a true hit?
	return true;
}

/**
 * Assert that every hit in the hits array also occurs in the
 * sanityHits array.
 */
template<typename TStr1, typename TStr2>
static bool reconcileHits(const TStr1& pat,
                          uint32_t patid,
                          bool fw,
                          vector<TStr2>& os,
                          vector<Hit>& hits,
                          vector<Hit>& sanityHits,
                          bool allowExact,
                          bool transpose)
{
	typedef TStr1 TStr;
	typedef typename Value<TStr>::Type TVal;
    // Sanity-check each result by checking whether it occurs
	// in the sanityHits array-of-vectors
    for(size_t i = 0; i < hits.size(); i++) {
    	const Hit& h = hits[i];
    	vector<Hit>::iterator itr;
    	bool found = false;
    	// Scan through the sanityHits vector corresponding to
    	// this hit text
    	for(itr = sanityHits.begin(); itr != sanityHits.end(); itr++) {
    		// If offset into hit text matches
			assert_gt(sanityHits.size(), 0);
    		if(h.h.first == itr->h.first && h.h.second == itr->h.second) {
    			// Assert that number of mismatches matches
    			if(h.fw != itr->fw || h.mms != itr->mms) {
    				cout << endl;
    				cout << "actual hit: fw=" << h.fw << endl;
    				cout << "sanity hit: fw=" << itr->fw << endl;
    			}
    			assert_eq(h.fw, itr->fw);
    			assert_eq(h.mms, itr->mms);
    			found = true;
    			sanityHits.erase(itr); // Retire this sanity hit
    			break;
    		}
    	}
    	// Assert that the Ebwt hit was covered by a sanity-check hit
    	if(!found) {
    		cout << "Ebwt hit not found in sanity-check hits:" << endl
    		     << "  " << pat << endl;
    		cout << "  ";
    		cout << endl;
    		cout << (patid>>revcomp) << (fw? "+":"-") << ":<"
    		     << h.h.first << "," << h.h.second << "," << h.mms << ">" << endl;
    		cout << "transpose: " << transpose << endl;
    		cout << "Candidates:" << endl;
        	for(itr = sanityHits.begin(); itr != sanityHits.end(); itr++) {
        		cout << "  " << itr->h.first << " (" << itr->h.second << ")" << endl;
        	}
    	}
    	assert(found); 
    }
    return true;
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
template<typename TStr>
static void mismatchSearch(PatternSource<TStr>& patsrc,
                           HitSink& sink,
                           EbwtSearchStats<TStr>& stats,
                           EbwtSearchParams<TStr>& params,
                           Ebwt<TStr>& ebwtFw,
                           Ebwt<TStr>& ebwtBw,
                           vector<TStr>& os)
{
	typedef typename Value<TStr>::Type TVal;
	assert(ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	assert(patsrc.hasMorePatterns());
    patsrc.setReverse(false); // reverse patterns
    params.setEbwtFw(true); // let search parameters reflect the forward index
	vector<Hit> sanityHits;
	uint32_t patid = 0;
	uint64_t lastHits = 0llu;
	uint32_t lastLen = 0; // for checking if all reads have same length
	String<uint8_t> doneMask;
    params.setEbwtFw(true);
	uint32_t numQs = ((qUpto == -1) ? 4 * 1024 * 1024 : qUpto);
	fill(doneMask, numQs, 0); // 4 MB, masks 32 million reads
	{
	Timer _t(cout, "Time for 1-mismatch forward search: ", timing);
	EbwtSearchState<TStr> s(ebwtFw, params, seed);
    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
    	bool exactOnly = false;
    	bool sfw = !revcomp || !patsrc.nextIsReverseComplement();
    	params.setFw(sfw);
    	uint32_t spatid = patid;
    	params.setPatId(spatid);
    	assert(!revcomp || (params.patId() & 1) == 0 || !params.fw());
    	assert(!revcomp || (params.patId() & 1) == 1 ||  params.fw());
    	TStr* pat = NULL;
		String<char>* qual = NULL;
		String<char>* name = NULL;
		patsrc.nextPattern(&pat, &qual, &name);
    	assert(!empty(*pat));
    	if(lastLen == 0) lastLen = length(*pat);
    	if(qSameLen && length(*pat) != lastLen) {
    		throw runtime_error("All reads must be the same length");
    	}
    	// Create state for a search on in the forward index
    	s.newQuery(pat, name, qual);
    	params.stats().incRead(s, *pat);
    	// Are there provisional hits?
    	if(sink.numProvisionalHits() > 0) {
    		// Shouldn't be any provisional hits unless we're doing
    		// pick-one and this is a reverse complement
    		assert(oneHit);
    		assert(!params.fw());
    		exactOnly = true;
    		// There is a provisional inexact match for the forward
    		// orientation of this pattern, so just try exact
    		ebwtFw.search(s, params);
    	    if(sink.numHits() > lastHits) {
    	    	// Got one or more exact hits from the reverse
    	    	// complement; reject provisional hits
    	    	sink.rejectProvisionalHits();
    	    } else {
    	    	// No exact hits from reverse complement; accept
    	    	// provisional hits, thus avoiding doing an inexact
    	    	// match on the reverse complement.
        	    ASSERT_ONLY(size_t numRetained = sink.retainedHits().size());
    	    	sink.acceptProvisionalHits();
    	    	assert_eq(sink.retainedHits().size(), numRetained);
    	    	assert_gt(sink.numHits(), lastHits);
    	    }
    	    assert_eq(0, sink.numProvisionalHits());
    	} else {
    		ebwtFw.search1MismatchOrBetter(s, params);
    	}
    	bool gotHits = sink.numHits() > lastHits;
	    // Set a bit indicating this pattern is done and needn't be
	    // considered by the 1-mismatch loop
	    if(oneHit && gotHits) {
	    	assert_eq(0, sink.numProvisionalHits());
	    	uint32_t mElt = patid >> 3;
	    	if(mElt > length(doneMask)) {
	    		// Add 50% more elements, initialized to 0
	    		fill(doneMask, mElt + patid>>4, 0);
	    	}
			
			// Set a bit indicating this pattern is done and needn't be
			// considered by the 1-mismatch loop
	    	doneMask[mElt] |= (1 << (patid & 7));
	    	if(revcomp && params.fw()) {
	    		assert(patsrc.hasMorePatterns());
	    		assert(patsrc.nextIsReverseComplement());
	    		// Ignore this pattern (the reverse complement of
	    		// the one we just matched)
		    	TStr* pat2;
				String<char>* qual2;
				String<char>* name2;
				patsrc.nextPattern(&pat2, &qual2, &name2);
		    	assert(!empty(*pat2));
		    	patid++;
		    	// Set a bit indicating this pattern is done
		    	doneMask[patid >> 3] |= (1 << (patid & 7));
		    	if(qSameLen && length(*pat2) != lastLen) {
		    		throw runtime_error("All reads must be the same length");
		    	}
		    	params.setFw(false);
		    	params.stats().incRead(s, *pat2);
	    		assert(!patsrc.nextIsReverseComplement());
	    	} else if(revcomp) {
    	    	// The reverse-complement version hit, so retroactively
	    		// declare the forward version done
    	    	uint32_t mElt = (patid-1) >> 3;
    	    	if(mElt > length(doneMask)) {
    	    		// Add 50% more elements, initialized to 0
    	    		fill(doneMask, mElt + patid>>4, 0);
    	    	}
    	    	doneMask[mElt] |= (1 << ((patid-1) & 7));
	    	}
	    }
	    // Check all hits against a naive oracle
    	if(sanityCheck && !os.empty() && !arrowMode) {
    	    vector<Hit>& hits = sink.retainedHits();
    	    // Accumulate hits found using a naive seed-and-extend into
    	    // sanityHits
    		findSanityHits(*pat, spatid, sfw, os, sanityHits, true, false);
    		if(hits.size() > 0) {
    			// We hit, check that oracle also got our hits
        	    assert(!oneHit || hits.size() == 1);
    			if(oneHit && hits[0].mms.count() > 0) {
					// If our oneHit hit is inexact, then there had
    				// better be no exact sanity hits
    				for(size_t i = 0; i < sanityHits.size(); i++) {
    					assert_gt(sanityHits[i].mms.count(), 0);
    				}
    			}
    			reconcileHits(*pat, spatid, sfw, os, hits, sanityHits, true, false);
    		} else if(!exactOnly) {
    			// If we tried exact and inexact and didn't hit, then
    			// oracle shouldn't have hit
        		assert_eq(0, sanityHits.size());
    		} else {
    			// If we tried exact only and didn't hit, then oracle
    			// shouldn't have any exact
				for(size_t i = 0; i < sanityHits.size(); i++) {
					assert_gt(sanityHits[i].mms.count(), 0);
				}
    		}
    		if(oneHit) {
    			// Ignore the rest of the oracle hits
    			sanityHits.clear(); 
    		} else {
    			// If in allHit mode, check that we covered *all* the
    			// hits produced by the oracle
    			checkSanityExhausted(*pat, spatid, sfw, sanityHits, false);
    		}
    		assert_eq(0, sanityHits.size());
    	    // Check that orientation of hits squares with orientation
    	    // of the pattern searched
    	    for(size_t i = 0; i < hits.size(); i++) {
    	    	assert_eq(sfw, hits[i].fw);
    	    }
    	    sink.clearRetainedHits();
    	}
	    patid++;
    	lastHits = sink.numHits();
    }
	}
	// Release most of the memory associated with the forward Ebwt
    ebwtFw.evictFromMemory();
	{
		// Load the rest of (vast majority of) the backward Ebwt into
		// memory
		Timer _t(cout, "Time loading Backward Ebwt: ", timing);
		ebwtBw.loadIntoMemory();
	}
    patsrc.reset();          // reset pattern source to 1st pattern
    patsrc.setReverse(true); // reverse patterns
    params.setEbwtFw(false); // let search parameters reflect the reverse index
	// Sanity-check the restored version of the Ebwt
	if(sanityCheck && !os.empty()) {
		TStr rest; ebwtBw.restore(rest);
		uint32_t restOff = 0;
		for(size_t i = 0; i < os.size(); i++) {
			uint32_t olen = length(os[i]);
			for(size_t j = 0; j < olen; j++) {
				assert_eq(os[i][olen-j-1], rest[restOff]);
				restOff++;
			}
			uint32_t leftover = (restOff & ~ebwtBw.eh().chunkMask());
			uint32_t diff = ebwtBw.eh().chunkLen() - leftover;
			if(leftover != 0) restOff += diff;
			assert_eq(0, restOff & ~ebwtBw.eh().chunkMask());
		}
	}
	assert(patsrc.hasMorePatterns());
	assert(!patsrc.nextIsReverseComplement());
	patid = 0;       // start again from id 0
	//lastHits = 0llu; // start again from 0 hits
	{
	Timer _t(cout, "Time for 1-mismatch backward search: ", timing);
	EbwtSearchState<TStr> s(ebwtBw, params, seed);
    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
    	bool sfw = !revcomp || !patsrc.nextIsReverseComplement();
    	params.setFw(sfw);
    	uint32_t spatid = patid;
    	params.setPatId(spatid);
    	assert(!revcomp || (params.patId() & 1) == 0 || !params.fw());
    	assert(!revcomp || (params.patId() & 1) == 1 ||  params.fw());
    	TStr* pat = NULL;
		String<char>* qual = NULL;
		String<char>* name = NULL;
		patsrc.nextPattern(&pat, &qual, &name);
    	assert(!empty(*pat));
    	s.newQuery(pat, name, qual);
    	params.stats().incRead(s, *pat);
		// Skip if previous phase determined this read is "done"; this
		// should only happen in oneHit mode
    	if((doneMask[patid >> 3] & (1 << (patid & 7))) != 0) {
    		assert(oneHit);
    		patid++;
    		continue;
    	}
		patid++;
    	// Try to match with one mismatch while suppressing exact hits
		//cerr << "searching for "<< pat<<endl;
    	ebwtBw.search1MismatchOrBetter(s, params, false /* suppress exact */);
    	sink.acceptProvisionalHits(); // automatically approve provisional hits
	    // If the forward direction matched with one mismatch, ignore
	    // the reverse complement
	    if(oneHit && revcomp && sink.numHits() > lastHits && params.fw()) {
    		assert(patsrc.nextIsReverseComplement());
    		assert(patsrc.hasMorePatterns());
    		// Ignore this pattern (the reverse complement of
    		// the one we just matched)
	    	TStr* pat2;
			String<char>* qual2;
			String<char>* name2;
			patsrc.nextPattern(&pat2, &qual2, &name2);

	    	assert(!empty(*pat2));
	    	patid++;
	    	params.setFw(false);
	    	params.stats().incRead(s, *pat2);
    		assert(!patsrc.nextIsReverseComplement());
	    }
	    // Check all hits against a naive oracle
    	if(sanityCheck && !os.empty() && !arrowMode) {
    	    vector<Hit>& hits = sink.retainedHits();
    	    // Accumulate hits found using a naive seed-and-extend into
    	    // sanityHits
    		findSanityHits(*pat, spatid, sfw, os, sanityHits, false, true);
    		if(hits.size() > 0) {
    			// We hit, check that oracle also got our hits
    			reconcileHits(*pat, spatid, sfw, os, hits, sanityHits, false, true);
    		} else {
    			// If we didn't hit, then oracle shouldn't have hit
        		assert_eq(0, sanityHits.size());
    		}
    		if(oneHit) {
    			// Ignore the rest of the oracle hits
    			sanityHits.clear(); 
    		} else {
    			// If in allHit mode, check that we covered *all* the
    			// hits produced by the oracle
    			checkSanityExhausted(*pat, spatid, sfw, sanityHits, true);
    		}
    		assert_eq(0, sanityHits.size());
    	    // Check that orientation of hits squares with orientation
    	    // of the pattern searched
    	    for(size_t i = 0; i < hits.size(); i++) {
    	    	assert_eq(sfw, hits[i].fw);
    	    }
    	    sink.clearRetainedHits();
    	}
    	lastHits = sink.numHits();
    }
	}
}

/**
 * Search through a pair of Ebwt indexes, one for the forward direction
 * and one for the backward direction, for exact query hits and hits
 * with at most one mismatch.
 * 
 * Forward Ebwt (ebwtFw) is already loaded into memory and backward
 * Ebwt (ebwtBw) is not loaded into memory.
 */

template<typename TStr>
static void mismatchSearchWithExtension(vector<String<Dna, Packed<> > >& packed_texts,
										PatternSource<TStr>& patsrc,
										HitSink& sink,
										EbwtSearchStats<TStr>& stats,
										EbwtSearchParams<TStr>& params,
										Ebwt<TStr>& ebwtFw,
										Ebwt<TStr>& ebwtBw,
										vector<TStr>& os)
{
	typedef typename Value<TStr>::Type TVal;
	assert(ebwtFw.isInMemory());
	assert(!ebwtBw.isInMemory());
	assert(patsrc.hasMorePatterns());
    patsrc.setReverse(false); // reverse patterns
    params.setEbwtFw(true); // let search parameters reflect the forward index
	
	uint32_t patid = 0;
	uint64_t lastHits = 0llu;
	uint32_t lastLen = 0; // for checking if all reads have same length
	String<uint8_t> doneMask;
    params.setEbwtFw(true);
	uint32_t numQs = ((qUpto == -1) ? 4 * 1024 * 1024 : qUpto);
	
	if (allowed_diffs == -1)
		allowed_diffs = default_allowed_diffs;
	
	ExactSearchWithLowQualityThreePrime<TStr> extend_exact(packed_texts, 
														   false, 
														   kmer,/* global, override */ 
														   allowed_diffs /*global, override*/,
														   revcomp /*global, override*/);
	
	OneMismatchSearchWithLowQualityThreePrime<TStr> extend_one_mismatch(packed_texts, 
																		false, 
																		kmer,/* global, override */ 
																		allowed_diffs /*global, override*/,
																		revcomp /*global, override*/);
	
	fill(doneMask, numQs, 0); // 4 MB, masks 32 million reads
	{
		Timer _t(cout, "Time for 1-mismatch forward search: ", timing);
		// Create state for a search on in the forward index
		EbwtSearchState<TStr> s(ebwtFw, params, seed);
		while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
			bool exactOnly = false;
			bool sfw = !revcomp || !patsrc.nextIsReverseComplement();
			params.setFw(sfw);
			uint32_t spatid = patid;
			params.setPatId(spatid);
			assert(!revcomp || (params.patId() & 1) == 0 || !params.fw());
			assert(!revcomp || (params.patId() & 1) == 1 ||  params.fw());
	    	TStr* pat = NULL;
			String<char>* qual = NULL;
			String<char>* name = NULL;
			patsrc.nextPattern(&pat, &qual, &name);
			assert(!empty(*pat));
			if(lastLen == 0) lastLen = length(*pat);
			if(qSameLen && length(*pat) != lastLen) {
				throw runtime_error("All reads must be the same length");
			}
			s.newQuery(pat, qual, name);
			params.stats().incRead(s, *pat);
			// Are there provisional hits?
			if(sink.numProvisionalHits() > 0) {
				// Shouldn't be any provisional hits unless we're doing
				// pick-one and this is a reverse complement
				assert(oneHit);
				assert(!params.fw());
				exactOnly = true;
				// There is a provisional inexact match for the forward
				// orientation of this pattern, so just try exact

				extend_exact.search(ebwtFw, stats, params, *pat, *name, *qual, sink);
				
				if(sink.numHits() > lastHits) {
					// Got one or more exact hits from the reverse
					// complement; reject provisional hits
					sink.rejectProvisionalHits();
				} else {
					// No exact hits from reverse complement; accept
					// provisional hits, thus avoiding doing an inexact
					// match on the reverse complement.
					ASSERT_ONLY(size_t numRetained = sink.retainedHits().size());
					sink.acceptProvisionalHits();
					assert_eq(sink.retainedHits().size(), numRetained);
					assert_gt(sink.numHits(), lastHits);
				}
				assert_eq(0, sink.numProvisionalHits());
			} else {
				//ebwtFw.search1MismatchOrBetter(s, params);
				extend_one_mismatch.search(ebwtFw, stats, params, *pat, *name, *qual, sink);
			}
			bool gotHits = sink.numHits() > lastHits;
			// Set a bit indicating this pattern is done and needn't be
			// considered by the 1-mismatch loop
			if(oneHit && gotHits) {
				assert_eq(0, sink.numProvisionalHits());
				uint32_t mElt = patid >> 3;
				if(mElt > length(doneMask)) {
					// Add 50% more elements, initialized to 0
					fill(doneMask, mElt + patid>>4, 0);
				}
				
				// Set a bit indicating this pattern is done and needn't be
				// considered by the 1-mismatch loop
				doneMask[mElt] |= (1 << (patid & 7));
				if(revcomp && params.fw()) {
					assert(patsrc.hasMorePatterns());
					assert(patsrc.nextIsReverseComplement());
					// Ignore this pattern (the reverse complement of
					// the one we just matched)
					TStr* pat2;
					String<char>* qual2;
					String<char>* name2;
					patsrc.nextPattern(&pat2, &qual2, &name2);
					assert(!empty(*pat2));
					patid++;
					// Set a bit indicating this pattern is done
					doneMask[patid >> 3] |= (1 << (patid & 7));
					if(qSameLen && length(*pat2) != lastLen) {
						throw runtime_error("All reads must be the same length");
					}
					params.setFw(false);
					params.stats().incRead(s, *pat2);
					assert(!patsrc.nextIsReverseComplement());
				} else if(revcomp) {
					// The reverse-complement version hit, so retroactively
					// declare the forward version done
					uint32_t mElt = (patid-1) >> 3;
					if(mElt > length(doneMask)) {
						// Add 50% more elements, initialized to 0
						fill(doneMask, mElt + patid>>4, 0);
					}
					doneMask[mElt] |= (1 << ((patid-1) & 7));
				}
			}
			patid++;
			lastHits = sink.numHits();
		}
	}
	// Release most of the memory associated with the forward Ebwt
    ebwtFw.evictFromMemory();
	{
		// Load the rest of (vast majority of) the backward Ebwt into
		// memory
		Timer _t(cout, "Time loading Backward Ebwt: ", timing);
		ebwtBw.loadIntoMemory();
	}
    patsrc.reset();          // reset pattern source to 1st pattern
    patsrc.setReverse(true); // reverse patterns
    params.setEbwtFw(false); // let search parameters reflect the reverse index
	
	assert(patsrc.hasMorePatterns());
	assert(!patsrc.nextIsReverseComplement());
	patid = 0;       // start again from id 0
	//lastHits = 0llu; // start again from 0 hits
	{
		Timer _t(cout, "Time for 1-mismatch backward search: ", timing);
		EbwtSearchState<TStr> s(ebwtBw, params, seed);
		while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
			bool sfw = !revcomp || !patsrc.nextIsReverseComplement();
			params.setFw(sfw);
			uint32_t spatid = patid;
			params.setPatId(spatid);
			assert(!revcomp || (params.patId() & 1) == 0 || !params.fw());
			assert(!revcomp || (params.patId() & 1) == 1 ||  params.fw());
	    	TStr* pat = NULL;
			String<char>* qual = NULL;
			String<char>* name = NULL;
			patsrc.nextPattern(&pat, &qual, &name);
			assert(!empty(*pat));
			s.newQuery(pat, name, qual);
			params.stats().incRead(s, *pat);
			// Skip if previous phase determined this read is "done"; this
			// should only happen in oneHit mode
			if((doneMask[patid >> 3] & (1 << (patid & 7))) != 0) {
				assert(oneHit);
				patid++;
				continue;
			}
			patid++;
			// Try to match with one mismatch while suppressing exact hits
			//ebwtBw.search1MismatchOrBetter(s, params, false /* suppress exact */);
			extend_one_mismatch.suppressExact(true);
			extend_one_mismatch.search(ebwtBw, stats, params, *pat, *name, *qual, sink);
						
			sink.acceptProvisionalHits(); // automatically approve provisional hits
			// If the forward direction matched with one mismatch, ignore
			// the reverse complement
			if(oneHit && revcomp && sink.numHits() > lastHits && params.fw()) {
				assert(patsrc.nextIsReverseComplement());
				assert(patsrc.hasMorePatterns());
				// Ignore this pattern (the reverse complement of
				// the one we just matched)
				TStr* pat2;
				String<char>* qual2;
				String<char>* name2;
				patsrc.nextPattern(&pat2, &qual2, &name2);
				assert(!empty(*pat2));
				patid++;
				params.setFw(false);
				params.stats().incRead(s, *pat2);
				assert(!patsrc.nextIsReverseComplement());
			}
			lastHits = sink.numHits();
		}
	}
}

#define SWITCH_TO_FW_INDEX() { \
	/* Evict the transposed index from memory if necessary */ \
	if(ebwtBw.isInMemory()) ebwtBw.evictFromMemory(); \
	assert(!ebwtBw.isInMemory()); \
	/* Load the forward index into memory if necessary */ \
	if(!ebwtFw.isInMemory()) { \
		Timer _t(cout, "Time loading forward Bowtie Index: ", timing); \
		ebwtFw.loadIntoMemory(); \
	} \
	assert(ebwtFw.isInMemory()); \
	patsrc.reset(); /* rewind pattern source to first pattern */ \
	assert(patsrc.hasMorePatterns()); \
	patsrc.setReverse(false); /* tell pattern source not to reverse patterns */ \
	params.setEbwtFw(true); /* tell search params that we're in the forward domain */ \
}

#define SWITCH_TO_BW_INDEX() { \
	/* Evict the forward index from memory if necessary */ \
	if(ebwtFw.isInMemory()) ebwtFw.evictFromMemory(); \
	assert(!ebwtFw.isInMemory()); \
	/* Load the forward index into memory if necessary */ \
	if(!ebwtBw.isInMemory()) { \
		Timer _t(cout, "Time loading transposed Bowtie Index: ", timing); \
		ebwtBw.loadIntoMemory(); \
	} \
	assert(ebwtBw.isInMemory()); \
	patsrc.reset(); /* rewind pattern source to first pattern */ \
	assert(patsrc.hasMorePatterns()); \
	patsrc.setReverse(true); /* tell pattern source to reverse patterns */ \
	params.setEbwtFw(false); /* tell search params that we're in the transposed domain */ \
}

#define GET_BOTH_PATTERNS(pf, qf, nf, pr, qr, nr) { \
	patsrc.nextPattern(&pf, &qf, &nf); \
	assert(pf != NULL); \
	assert(patsrc.nextIsReverseComplement()); \
	patsrc.nextPattern(&pr, &qr, &nr); \
	assert(pr != NULL); \
	assert(pr != patFw); \
	assert(!patsrc.nextIsReverseComplement()); \
	assert(isReverseComplement(*pf, *pr)); \
}

#define GET_FW_PATTERN(pf, qf, nf) { \
	patsrc.nextPattern(&pf, &qf, &nf); \
	assert(pf != NULL); \
	assert(patsrc.nextIsReverseComplement()); \
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
        int qualWobble,                 /// currently unused
        int seedMms,                    /// max # mismatches allowed in seed
                                        /// (like maq map's -n option)
                                        /// Can only be 1 or 2, default: 1
        PatternSource<TStr>& patsrc,    /// pattern source
        HitSink& sink,                  /// hit sink
        EbwtSearchStats<TStr>& stats,   /// statistics (mostly unused)
        EbwtSearchParams<TStr>& params, /// search parameters
        Ebwt<TStr>& ebwtFw,             /// index of original text
        Ebwt<TStr>& ebwtBw,             /// index of transposed text
        vector<TStr>& os)    /// text strings, if available (empty otherwise)
{
	typedef typename Value<TStr>::Type TVal;
	vector<Hit> sanityHits;
	uint32_t numPats;
	assert(revcomp);
	// Load forward index
	SWITCH_TO_FW_INDEX();
	uint32_t numQs = ((qUpto == -1) ? 10 * 1024 * 1024 : qUpto);
	vector<bool> doneMask(numQs, false);
	uint32_t s = seedLen;
	uint32_t s3 = seedLen >> 1; // length of 3' half of seed
	uint32_t s5 = (seedLen >> 1) + (s & 1); // length of 5' half of seed
	{
		// Phase 1: Consider cases 1R and 2R
		Timer _t(cout, "Time for seeded quality search Phase 1 of 4: ", timing);
		BacktrackManager<TStr> bt(ebwtFw, params,
		                          (seedMms > 0)?  s5 : s,// unrevOff, 
		                          (seedMms == 2)? s5 : s,// 1revOff
		                          s,                     // 2revOff
		                          0, 0,                  // itop, ibot
		                          qualCutoff,            // qualThresh
		                          0,                     // qualWobble
		                          0,                     // reportSeedlings (don't)
		                          NULL,                  // seedlings
		                          NULL,                  // mutations
		                          verbose,               // verbose
		                          true,                  // oneHit
		                          seed,                  // seed
		                          &os);
		uint32_t patid = 0;
		uint32_t lastLen = 0; // for checking if all reads have same length
		params.setFw(false);  // we'll only look at reverse complements this round
		TStr* patFw = NULL; String<char>* qualFw = NULL; String<char>* nameFw = NULL;
		TStr* patRc = NULL; String<char>* qualRc = NULL; String<char>* nameRc = NULL;
	    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
	    	if(patid>>1 > doneMask.capacity()) {
	    		// Expand doneMask
	    		doneMask.resize(doneMask.capacity()*2, 0);
	    	}
	    	GET_BOTH_PATTERNS(patFw, qualFw, nameFw, patRc, qualRc, nameRc);
			size_t plen = length(*patFw);
			if(qSameLen) {
				if(lastLen == 0) lastLen = plen;
				else assert_eq(lastLen, plen);
			}
			// Set up backtracker with reverse complement
			bt.setQuery(patRc, qualRc, nameRc);
			params.setPatId(patid+1);
			ASSERT_ONLY(uint64_t numHits = sink.numHits());
			bool hit = bt.backtrack();
			assert(hit  || numHits == sink.numHits());
			assert(!hit || numHits <  sink.numHits());
			if(hit) {
				// If we reach here, then we obtained a hit for case
				// 1R, 2R or 3R and can stop considering this read
				doneMask[patid>>1] = true;
			} else {
				// If we reach here, then cases 1R, 2R, and 3R have
				// been eliminated and the read needs further
				// examination
			}
			patid += 2;
	    }
	    numPats = patid;
	    assert_leq(numPats>>1, doneMask.size());
	}
	// Unload forward index and load transposed index
	SWITCH_TO_BW_INDEX();
	String<uint8_t> seedlingsRc;
	reserve(seedlingsRc, 10 * 1024 * 1024, Exact());
	{
		// Phase 2: Consider cases 1F, 2F and 3F and generate seedlings
		// for case 4R
		Timer _t(cout, "Time for seeded quality search Phase 2 of 4: ", timing);
		// BacktrackManager to search for hits for cases 1F, 2F, 3F
		BacktrackManager<TStr> btf(ebwtBw, params,
                                   (seedMms > 0)?  s5 : s,// unrevOff
                                   (seedMms == 2)? s5 : s,// 1revOff
                                   s,                     // 2revOff
		                           0, 0,                  // itop, ibot
		                           qualCutoff,            // qualThresh
		                           0,                     // qualWobble
		                           0,                     // reportSeedlings (no)
		                           NULL,                  // seedlings
			                       NULL,                  // mutations
		                           verbose,               // verbose
		                           true,                  // oneHit
			                       seed+1,                // seed
			                       &os);
		// BacktrackManager to search for seedlings for case 4R
		BacktrackManager<TStr> btr(ebwtBw, params,
		                           s3,                    // unrevOff
		                           (seedMms == 2)? s3 : s,// 1revOff
		                           s,                     // 2revOff
		                           0, 0,                  // itop, ibot
		                           qualCutoff,            // qualThresh (none)
		                           0,                     // qualWobble
		                           seedMms,               // reportSeedlings (yes)
		                           &seedlingsRc,          // seedlings
			                       NULL,                  // mutations
		                           verbose,               // verbose
		                           true,                  // oneHit
			                       seed+2,                // seed
			                       &os);
		uint32_t patid = 0;
		TStr* patFw = NULL; String<char>* qualFw = NULL; String<char>* nameFw = NULL;
		TStr* patRc = NULL; String<char>* qualRc = NULL; String<char>* nameRc = NULL;
	    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
	    	assert_lt((patid>>1), doneMask.capacity());
	    	assert_lt((patid>>1), doneMask.size());
	    	if(doneMask[patid>>1]) {
				patsrc.skipPattern();
				patsrc.skipPattern();
	    		patid += 2;
	    		continue;
	    	}
	    	
			// If we reach here, then cases 1R, 2R, and 3R have been
	    	// eliminated.  The next most likely cases are 1F, 2F and
	    	// 3F...
			params.setFw(true);  // looking at forward strand
	    	GET_BOTH_PATTERNS(patFw, qualFw, nameFw, patRc, qualRc, nameRc);
			btf.setQuery(patFw, qualFw, nameFw);
			params.setPatId(patid);
			ASSERT_ONLY(uint64_t numHits = sink.numHits());
			
			// Do a 12/24 backtrack on the forward-strand read using
			// the transposed index.  This will find all case 1F, 2F
			// and 3F hits.
			bool hit = btf.backtrack();
			assert(hit  || numHits == sink.numHits());
			assert(!hit || numHits <  sink.numHits());
			if(hit) {
				// The reverse complement hit, so we're done with this
				// read
				doneMask[patid>>1] = true;
				patid += 2;
				continue;
			}
			// No need to collect seedlings if we're not allowing
			// mismatches in the 5' seed half
			if(seedMms == 0) {
				patid += 2;
				continue;
			}
			
			// If we reach here, then cases 1F, 2F, 3F, 1R, 2R, and 3R
			// have been eliminated, leaving us with cases 4F and 4R
			// (the cases with 1 mismatch in the 5' half of the seed)
			params.setFw(false);  // looking at reverse-comp strand
			params.setPatId(patid+1);
			btr.setQuery(patRc, qualRc, nameRc);
			btr.setQlen(s); // just look at the seed
			uint32_t numSeedlings = length(seedlingsRc);
			// Do a 12/24 seedling backtrack on the reverse-comp read
			// using the transposed index.  This will find seedlings
			// for case 4R
			btr.backtrack();
			hit = length(seedlingsRc) > numSeedlings;
			append(seedlingsRc, 0xff);
#ifndef NDEBUG
			// Sanity-check the generated seedling hits
			append(seedlingsRc, (patid>>1));
			{
				uint32_t id2 = numSeedlings;
				char lastQual = 0;
				while(seedlingsRc[id2++] != 0xff) {
					uint8_t pos = seedlingsRc[id2-1];
					assert_lt(pos, s5);
					// Get the character to mutate it to
					uint8_t chr = seedlingsRc[id2++];
					assert_lt(chr, 4);
					uint8_t oldChar = (uint8_t)(*patRc)[pos];
					uint8_t oldQual = (uint8_t)(*qualRc)[pos]-33;
					assert_leq(oldQual, 40);
					if(seedMms < 2) {
						assert_geq((*qualRc)[pos], lastQual);
					}
					lastQual = (*qualRc)[pos];
					assert_neq(oldChar, chr);
					if(seedlingsRc[id2] == 0xfe) {
						id2++;
					}
				}
			}
#endif
			patid += 2;
	    }
	    assert_eq(numPats, patid);
	}
	if(seedMms == 0) {
		// If we're not allowing any mismatches in the seed, then there
		// is no need to continue to phases 3 and 4
		return;
	}
	// Unload transposed index and load forward index
	SWITCH_TO_FW_INDEX();
	String<uint8_t> seedlingsFw;
	reserve(seedlingsFw, 10 * 1024 * 1024, Exact());
	{
		// Phase 3: Consider cases 3R and 4R and generate seedlings for
		// case 4F
		Timer _t(cout, "Time for seeded quality search Phase 3 of 4: ", timing);
		// BacktrackManager to search for seedlings for case 4F
		BacktrackManager<TStr> btf(ebwtFw, params,
                                   s3,                    // unrevOff
                                   (seedMms == 2)? s3 : s,// 1revOff
                                   s,                     // 2revOff
		                           0, 0,                  // itop, ibot
		                           qualCutoff,            // qualThresh (none)
		                           0,                     // qualWobble
		                           seedMms,               // reportSeedlings (do)
		                           &seedlingsFw,          // seedlings
			                       NULL,                  // mutations
		                           verbose,               // verbose
		                           true,                  // oneHit
			                       seed+3,                // seed
			                       &os);
		// BacktrackManager to search for hits for case 4R; default
		// parameters are appropriate when searching with 1 mismatch
		// allowed in the entire seed.  If 2 mismatches are allowed,
		// then, depending on whether the seedling has 1 or 2
		// mismatches, we may need to move unrevOff up to s2.
		BacktrackManager<TStr> btr(ebwtFw, params,
		                           s,       // unrevOff
		                           s,       // 1revOff
		                           s,       // 2revOff
		                           0, 0,    // itop, ibot
		                           qualCutoff, // qualThresh
		                           0,       // qualWobble
		                           0,       // reportSeedlings (don't)
		                           NULL,    // seedlings
			                       NULL,    // mutations
		                           verbose, // verbose
		                           true,    // oneHit
			                       seed+4,  // seed
			                       &os);
		BacktrackManager<TStr> btr2(ebwtFw, params,
		                           0,       // unrevOff
		                           s5,      // 1revOff
		                           s,       // 2revOff
		                           0, 0,    // itop, ibot
		                           qualCutoff, // qualThresh
		                           0,       // qualWobble
		                           0,       // reportSeedlings (don't)
		                           NULL,    // seedlings
			                       NULL,    // mutations
		                           verbose, // verbose
		                           true,    // oneHit
			                       seed+5,  // seed
			                       &os,
			                       true);   // halfAndHalf
		uint32_t patid = 0;
		uint32_t seedlingId = 0;
		TStr* patFw = NULL; String<char>* qualFw = NULL; String<char>* nameFw = NULL;
		TStr* patRc = NULL; String<char>* qualRc = NULL; String<char>* nameRc = NULL;
		ASSERT_ONLY(uint32_t seedlingsRcLen = length(seedlingsRc));
	    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
	    	assert_lt((patid>>1), doneMask.capacity());
	    	assert_lt((patid>>1), doneMask.size());
	    	if(doneMask[patid>>1]) {
				patsrc.skipPattern();
				patsrc.skipPattern();
	    		patid += 2;
	    		continue;
	    	}
	    	GET_BOTH_PATTERNS(patFw, qualFw, nameFw, patRc, qualRc, nameRc);
			params.setFw(false);  // looking at reverse-comp strand
			params.setPatId(patid+1);
			btr.setQuery(patRc, qualRc, nameRc);

			// Given the seedlines generated in phase 2, check for hits
			// for case 4R
			uint32_t plen = length(*patRc);
			if(seedlingsRc[seedlingId++] != 0xff) {
				assert_lt(seedlingId, seedlingsRcLen);
				seedlingId--;
				bool hit = false;
#ifndef NDEBUG
				{
					// Sanity-check the content of the seedlings array
					// for this pattern; don't search yet
					uint32_t id2 = seedlingId;
					char lastQual = 0;
					while(seedlingsRc[id2++] != 0xff) {
						uint8_t pos = seedlingsRc[id2-1];
						assert_lt(pos, s5);
						assert_lt(pos, plen);
						// Get the character to mutate it to
						uint8_t chr = seedlingsRc[id2++];
						assert_lt(chr, 4);
						uint8_t oldChar = (uint8_t)(*patRc)[plen-1-pos];
						uint8_t oldQual = (uint8_t)(*qualRc)[plen-1-pos]-33;
						assert_leq(oldQual, 40);
						if(seedMms < 2) {
							assert_geq((*qualRc)[plen-1-pos], lastQual);
						}
						lastQual = (*qualRc)[plen-1-pos];
						assert_neq(oldChar, chr);
						if(seedlingsRc[id2] == 0xfe) {
							// Skip minor separator
							id2++;
						}
					}
					// Make sure that we and the seedlings array are
					// talking about the same patid
					assert_eq((patid>>1) & 0xff, seedlingsRc[id2]);
				}
#endif
				while(seedlingsRc[seedlingId++] != 0xff) {
					seedlingId--; // point to that non-0xff character
					String<QueryMutation> muts;
					reserve(muts, 4);
					assert_gt(plen, 0);
					uint8_t oldQuals = 0;
					do {
						// Get the query position (offset from the 5' end)
						// to be mutated
						uint8_t pos = seedlingsRc[seedlingId++];
						uint8_t tpos = plen-1-pos;
						// Get the character to mutate it to
						uint8_t chr = seedlingsRc[seedlingId++];
						uint8_t oldChar = (uint8_t)(*patRc)[tpos];
						oldQuals += ((uint8_t)(*qualRc)[tpos]-33);
						assert_leq(oldQuals, qualCutoff);
						append(muts, QueryMutation(tpos, oldChar, chr));
					} while(seedlingsRc[seedlingId++] == 0xfe);
					seedlingId--; // point to that non-0xfe character
					assert_gt(length(muts), 0);
					assert_leq(length(muts), 2);
					// Set the backtracking thresholds appropriately
					// Now begin the backtracking, treating the first
					// 24 bases as unrevisitable
					ASSERT_ONLY(uint64_t numHits = sink.numHits());
					ASSERT_ONLY(TStr tmp = (*patRc));
					btr.setMuts(&muts);
					hit = btr.backtrack(oldQuals);
					btr.setMuts(NULL);
					assert_eq(tmp, (*patRc));
					assert(hit  || numHits == sink.numHits());
					assert(!hit || numHits <  sink.numHits());
					if(hit) {
						// The reverse complement hit, so we're done with this
						// read
						doneMask[patid>>1] = true;
						while(seedlingsRc[seedlingId++] != 0xff);
						break;
					}
				}
				assert_leq(seedlingId, seedlingsRcLen);
				assert_eq(0xff, seedlingsRc[seedlingId-1]);
			}
			assert_eq(0xff, seedlingsRc[seedlingId-1]);
			assert_eq((patid>>1) & 0xff, seedlingsRc[seedlingId]);
			ASSERT_ONLY(seedlingId++); // skip patid sanity marker
			
			// Case 4R yielded a hit; mark this pattern as done and
			// continue to next pattern
	    	if(doneMask[patid>>1]) {
	    		patid += 2;
	    		continue;
	    	}
	    	
	    	// If we're in two-mismatch mode, then now is the time to
	    	// try the final case that might apply to the reverse
	    	// complement pattern: 1 mismatch in each of the 3' and 5'
	    	// halves of the seed.
	    	if(seedMms == 2) {
				btr2.setQuery(patRc, qualRc, nameRc);
				ASSERT_ONLY(uint64_t numHits = sink.numHits());
				bool hit = btr2.backtrack();
				assert(hit  || numHits == sink.numHits());
				assert(!hit || numHits <  sink.numHits());
				if(hit) {
					doneMask[patid>>1] = true;
		    		patid += 2;
					continue;
				}
	    	}
	    	
#ifndef NDEBUG
			// The reverse-complement version of the read doesn't hit
	    	// at all!  Check with the oracle to make sure it agrees.
	    	if(sanityCheck && os.size() > 0) {
				vector<Hit> hits;
				uint32_t twoRevOff = s;
				uint32_t oneRevOff = (seedMms <= 1) ? s : 0;
				uint32_t unrevOff   = (seedMms == 0) ? s : 0;
				bool newName = false;
				if(nameRc == NULL) {
					nameRc = new String<char>("no_name");
					newName = true;
				}
				BacktrackManager<TStr>::naiveOracle(
				        os,
						*patRc,
						plen,
				        *qualRc,
				        *nameRc,
				        patid+1,    // patid
				        hits,
				        qualCutoff, 
				        unrevOff,
				        oneRevOff,
				        twoRevOff,
				        false,      // fw
				        true);      // ebwtFw
				if(hits.size() > 0) {
					// Print offending hit obtained by oracle
					BacktrackManager<TStr>::printHit(
						os,
						hits[0],
						*patRc,
						plen,
					    unrevOff,
					    oneRevOff,
					    twoRevOff,
					    true);      // ebwtFw
				}
				if(newName) {
					delete nameRc;
				}
				assert_eq(0, hits.size());
	    	}
#endif

			// If we reach here, then cases 1F, 2F, 3F, 1R, 2R, 3R and
			// 4R have been eliminated leaving only 4F.
			params.setFw(true);  // looking at forward strand
			params.setPatId(patid);
			btf.setQuery(patFw, qualFw, nameFw);
			btf.setQlen(24); // just look at the seed
			ASSERT_ONLY(uint32_t numSeedlings = length(seedlingsFw));
			// Do a 12/24 seedling backtrack on the forward read
			// using the normal index.  This will find seedlings
			// for case 4F
			btf.backtrack();
			append(seedlingsFw, 0xff);
#ifndef NDEBUG
			append(seedlingsFw, (patid>>1));
			{
				uint32_t id2 = numSeedlings;
				char lastQual = 0;
				while(seedlingsFw[id2++] != 0xff) {
					uint8_t pos = seedlingsFw[id2-1];
					assert_lt(pos, s5);
					// Get the character to mutate it to
					uint8_t chr = seedlingsFw[id2++];
					assert_lt(chr, 4);
					uint8_t oldChar = (uint8_t)(*patFw)[pos];
					uint8_t oldQual = (uint8_t)(*qualFw)[pos]-33;
					assert_leq(oldQual, 40);
					if(seedMms < 2) {
						assert_geq((*qualFw)[pos], lastQual);
					}
					lastQual = (*qualFw)[pos];
					assert_neq(oldChar, chr);
					if(seedlingsFw[id2] == 0xfe) {
						id2++;
					}
				}
			}
#endif
			patid += 2;
	    }
	    assert_eq(seedlingId, seedlingsRcLen);
	    assert(numPats == patid || numPats+2 == patid);
	}
	// Unload forward index and load transposed index
	SWITCH_TO_BW_INDEX();
	{
		// Phase 4: Consider case 4F
		Timer _t(cout, "Time for seeded quality search Phase 4 of 4: ", timing);
		// BacktrackManager to search for hits for case 4F; default
		// parameters are appropriate when searching with 1 mismatch
		// allowed in the entire seed.  If 2 mismatches are allowed,
		// then, depending on whether the seedling has 1 or 2
		// mismatches, we may need to move unrevOff up to s2.
		BacktrackManager<TStr> btf(ebwtBw, params,
                                   s,       // unrevOff
                                   s,       // 1revOff
                                   s,       // 2revOff
		                           0, 0,    // itop, ibot
		                           qualCutoff, // qualThresh
		                           0,       // qualWobble
		                           0,       // reportSeedlings (don't)
		                           NULL,    // seedlings
			                       NULL,    // mutations
		                           verbose, // verbose
		                           true,    // oneHit
		                           seed+6,  // seed
		                           &os);
		BacktrackManager<TStr> btf2(ebwtBw, params,
                                   0,       // unrevOff
                                   s5,      // 1revOff
                                   s,       // 2revOff
		                           0, 0,    // itop, ibot
		                           qualCutoff, // qualThresh
		                           0,       // qualWobble
		                           0,       // reportSeedlings (don't)
		                           NULL,    // seedlings
			                       NULL,    // mutations
		                           verbose, // verbose
		                           true,    // oneHit
		                           seed+7,  // seed
		                           &os,
		                           true);   // halfAndHalf
		uint32_t patid = 0;
		uint32_t seedlingId = 0;
		uint32_t seedlingsFwLen = length(seedlingsFw);
		params.setFw(true);  // looking only at forward strand
		TStr* patFw = NULL; String<char>* qualFw = NULL; String<char>* nameFw = NULL;
	    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) {
	    	assert_lt((patid>>1), doneMask.capacity());
	    	assert_lt((patid>>1), doneMask.size());
	    	if(doneMask[patid>>1]) {
				patsrc.skipPattern();
				patsrc.skipPattern();
	    		patid += 2;
	    		continue;
	    	}
	    	GET_FW_PATTERN(patFw, qualFw, nameFw);
	    	patsrc.skipPattern();
			assert(!patsrc.nextIsReverseComplement());
			params.setPatId(patid);
			params.setFw(true);
			btf.setQuery(patFw, qualFw, nameFw);

			// Given the seedlines generated in phase 3, check for hits
			// for case 4F
			uint32_t plen = length(*patFw);
			if(seedlingsFw[seedlingId++] != 0xff) {
				assert_lt(seedlingId, seedlingsFwLen);
				seedlingId--;
				bool hit = false;
#ifndef NDEBUG
				{
					// Sanity-check the content of the seedlings array
					// for this pattern; don't search yet
					uint32_t id2 = seedlingId;
					char lastQual = 0;
					while(seedlingsFw[id2++] != 0xff) {
						uint8_t pos = seedlingsFw[id2-1];
						assert_lt(pos, s5);
						assert_lt(pos, plen);
						uint8_t tpos = plen-1-pos;
						// Get the character to mutate it to
						uint8_t chr = seedlingsFw[id2++];
						assert_lt(chr, 4);
						uint8_t oldChar = (uint8_t)(*patFw)[tpos];
						uint8_t oldQual = (uint8_t)(*qualFw)[tpos]-33;
						assert_leq(oldQual, 40);
						if(seedMms < 2) {
							assert_geq((*qualFw)[pos], lastQual);
						}
						lastQual = (*qualFw)[tpos];
						assert_neq(oldChar, chr);
						if(seedlingsFw[id2] == 0xfe) {
							// Skip minor separator
							id2++;
						}
					}
					// Make sure that we and the seedlings array are
					// talking about the same patid
					assert_eq((patid>>1) & 0xff, seedlingsFw[id2]);
				}
#endif
				while(seedlingsFw[seedlingId++] != 0xff) {
					seedlingId--; // point to that non-0xff character
					String<QueryMutation> muts;
					reserve(muts, 4);
					assert_gt(plen, 0);
					uint8_t oldQuals = 0;
					do {
						// Get the query position (offset from the 5' end)
						// to be mutated
						uint8_t pos = seedlingsFw[seedlingId++];
						uint8_t tpos = plen-1-pos;
						// Get the character to mutate it to
						uint8_t chr = seedlingsFw[seedlingId++];
						uint8_t oldChar = (uint8_t)(*patFw)[tpos];
						oldQuals += ((uint8_t)(*qualFw)[tpos]-33);
						assert_leq(oldQuals, qualCutoff);
						append(muts, QueryMutation(tpos, oldChar, chr));
					} while(seedlingsFw[seedlingId++] == 0xfe);
					seedlingId--; // point to that non-0xfe character
					assert_gt(length(muts), 0);
					assert_leq(length(muts), 2);
					// Now begin the backtracking, treating the first
					// 24 bases as unrevisitable
					ASSERT_ONLY(uint64_t numHits = sink.numHits());
					btf.setMuts(&muts);
					hit = btf.backtrack(oldQuals);
					btf.setMuts(NULL);
					assert(hit  || numHits == sink.numHits());
					assert(!hit || numHits <  sink.numHits());
					if(hit) {
						// The reverse complement hit, so we're done with this
						// read
						doneMask[patid>>1] = true;
						while(seedlingsFw[seedlingId++] != 0xff);
						break;
					}
				}
				assert_leq(seedlingId, seedlingsFwLen);
				if(seedlingId == seedlingsFwLen) {
					break; // break out of the pattern loop
				}
				assert_eq(0xff, seedlingsFw[seedlingId-1]);
			} // if(seedlingsFw[seedlingId++] != 0xff)
			assert_eq(0xff, seedlingsFw[seedlingId-1]);
			assert_eq((patid>>1) & 0xff, seedlingsFw[seedlingId]);
			ASSERT_ONLY(seedlingId++); // skip patid sanity marker

			// Case 4F yielded a hit; mark this pattern as done and
			// continue to next pattern
	    	if(doneMask[patid>>1]) {
	    		patid += 2;
	    		continue;
	    	}

	    	// If we're in two-mismatch mode, then now is the time to
	    	// try the final case that might apply to the reverse
	    	// complement pattern: 1 mismatch in each of the 3' and 5'
	    	// halves of the seed.
	    	if(seedMms == 2) {
				ASSERT_ONLY(uint64_t numHits = sink.numHits());
				btf2.setQuery(patFw, qualFw, nameFw);
				bool hit = btf2.backtrack();
				assert(hit  || numHits == sink.numHits());
				assert(!hit || numHits <  sink.numHits());
				if(hit) {
					doneMask[patid>>1] = true;
		    		patid += 2;
					continue;
				}
	    	}

#ifndef NDEBUG
	    	
			// The forward version of the read doesn't hit at all!
			// Check with the oracle to make sure it agrees.
	    	if(sanityCheck && os.size() > 0) {
				vector<Hit> hits;
				uint32_t twoRevOff = s;
				uint32_t oneRevOff = (seedMms <= 1) ? s : 0;
				uint32_t unrevOff   = (seedMms == 0) ? s : 0;
				bool newName = false;
				if(nameFw == NULL) {
					nameFw = new String<char>("no_name");
					newName = true;
				}
				BacktrackManager<TStr>::naiveOracle(
				        os,
						*patFw,
						plen,
				        *qualFw,
				        *nameFw,
				        patid,      // patid
				        hits,
				        qualCutoff, 
				        unrevOff,
				        oneRevOff,
				        twoRevOff,
				        true,       // fw
				        false);     // ebwtFw
				if(hits.size() > 0) {
					// Print offending hit obtained by oracle
					BacktrackManager<TStr>::printHit(
						os,
						hits[0],
						*patFw,
						plen,
					    unrevOff,
					    oneRevOff,
					    twoRevOff,
					    false);     // ebwtFw
				}
				if(newName) {
					delete nameFw;
				}
				assert_eq(0, hits.size());
			}
#endif
			patid += 2;
	    } // while(patsrc.hasMorePatterns()...
	    assert_eq(seedlingId, seedlingsFwLen);
	} // end of Phase 4
	return;
}

template<typename TStr>
static void driver(const char * type,
                   const string& infile,
                   const string& query,
                   const vector<string>& queries,
                   const string& outfile)
{
	vector<TStr> ps;
	vector<TStr> os;
	// Read original string(s) from command-line if given (for sanity checking)
	if(sanityCheck && !origString.empty()) {
		if(origString.substr(origString.length()-4) == ".mfa" ||
		   origString.substr(origString.length()-3) == ".fa")
		{
			vector<string> origFiles;
			tokenize(origString, ",", origFiles);
			readSequenceFiles<TStr, Fasta>(origFiles, os);
		} else {
			readSequenceString(origString, os);
		}
	}
	// Seed random number generator
	srandom(seed);
	// Create a pattern source for the queries
	PatternSource<TStr> *patsrc = NULL;
	switch(format) {
		case FASTA:   patsrc = new FastaPatternSource<TStr> (queries, revcomp, false, patDumpfile, trim3, trim5); break;
		case FASTQ:   patsrc = new FastqPatternSource<TStr> (queries, revcomp, false, patDumpfile, trim3, trim5, solexa_quals); break;
	    //case BFQ:     patsrc = new BfqPatternSource<TStr>   (queries, revcomp, false, patDumpfile, trim3, trim5); break;
	    //case SOLEXA:  patsrc = new SolexaPatternSource<TStr>(queries, revcomp, false, patDumpfile, trim3, trim5); break;
		case CMDLINE: patsrc = new VectorPatternSource<TStr>(queries, revcomp, false, patDumpfile, trim3, trim5); break;
		default: assert(false);
	}
	// Check that input is non-empty
	if(!patsrc->hasMorePatterns()) {
		cerr << "Error: Empty input!  Check that file format is correct." << endl;
		exit(1);
	}
	if(skipSearch) return;
	// Open hit output file
	ostream *fout;
	if(!outfile.empty()) {
		fout = new ofstream(outfile.c_str(), ios::binary);
	} else {
		fout = &cout;
	}
	// Initialize Ebwt object and read in header
    Ebwt<TStr> ebwt(infile, /* overriding: */ offRate, verbose, sanityCheck);
    assert_geq(ebwt.eh().offRate(), offRate);
    Ebwt<TStr>* ebwtBw = NULL;
    // We need the transposed index if mismatches are allowed
    if(mismatches > 0 || maqLike) {
    	ebwtBw = new Ebwt<TStr>(infile + ".rev", /* overriding: */ offRate, verbose, sanityCheck);
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
		Timer _t(cout, "Time loading Ebwt: ", timing);
	    ebwt.loadIntoMemory();
	}
	// Sanity-check the restored version of the Ebwt
	if(!maqLike && sanityCheck && !os.empty()) {
		TStr rest; ebwt.restore(rest);
		uint32_t restOff = 0;
		for(size_t i = 0; i < os.size(); i++) {
			for(size_t j = 0; j < length(os[i]); j++) {
				assert_eq(os[i][j], rest[restOff]);
				restOff++;
			}
			uint32_t leftover = (restOff & ~ebwt.eh().chunkMask());
			uint32_t diff = ebwt.eh().chunkLen() - leftover;
			if(leftover != 0) restOff += diff;
			assert_eq(0, restOff & ~ebwt.eh().chunkMask());
		}
	}
    // If sanity-check is enabled and an original text string
    // was specified, sanity-check the Ebwt by confirming that
    // its detransformation equals the original.
	if(!maqLike && sanityCheck && !os.empty()) {
		TStr rs; ebwt.restore(rs);
		TStr joinedo = Ebwt<TStr>::join(os, ebwt.eh().chunkRate(), seed);
		assert_eq(joinedo, rs);
	}
	{
		Timer _t(cout, "Time searching: ", timing);
		// Set up hit sink; if sanityCheck && !os.empty() is true,
		// then instruct the sink to "retain" hits in a vector in
		// memory so that we can easily sanity check them later on
		HitSink *sink;
		if(binOut)
		{
			sink = new BufferedBinaryHitSink(
					*fout,
					revcomp,
					reportOpps,
					sanityCheck && !os.empty());
		}
		else if(!concise)
		{
			sink = new VerboseHitSink(*fout,
									 revcomp,
									 sanityCheck && !os.empty());
		}
		else
		{
			sink = new PrettyHitSink(
					*fout,
					revcomp,
					reportOpps,
					sanityCheck && !os.empty());
		}
		EbwtSearchStats<TStr> stats;
		EbwtSearchParams<TStr> params(*sink,   // HitSink
		                              stats,   // EbwtSearchStats
		                              // Policy for how to resolve multiple hits
		                              (oneHit? MHP_PICK_1_RANDOM : MHP_CHASE_ALL),
		                              os,      //
		                              revcomp,
		                              true,
		                              true,
		                              arrowMode);
		if(maqLike) {
			seededQualCutoffSearch(seedLen,
			                       qualThresh,
			                       15,      // qualWobble
			                       seedMms,
			                       *patsrc,
			                       *sink,
			                       stats,
			                       params,
			                       ebwt,    // forward index
			                       *ebwtBw, // transposed index (not optional)
			                       os);     // references, if available
		}
		else if(mismatches > 0) 
		{
			// Search with mismatches
			if (kmer != -1 || allowed_diffs != -1)
			{
				vector<String<Dna, Packed<> > > ss;
				unpack(infile + ".3.ebwt", ss, NULL);
				mismatchSearchWithExtension(ss, *patsrc, *sink, stats, params, ebwt, *ebwtBw, os);
			}
			else
			{
				mismatchSearch(*patsrc, *sink, stats, params, ebwt, *ebwtBw, os);
			}
			
		} else {
			if (kmer != -1)
			{
				vector<String<Dna, Packed<> > > ss;
				unpack(infile + ".3.ebwt", ss, NULL);
				// Search for hits on the 5' end, and then try to extend them
				// with a dynamic programming algorithm
				exactSearchWithExtension(ss, *patsrc, *sink, stats, params, ebwt, os);
			}
			else
			{
				// Search without mismatches
				exactSearch(*patsrc, *sink, stats, params, ebwt, os);
			}
		}
	    sink->finish(); // end the hits section of the hit file
	    if(printStats) {
		    // Write some high-level searching parameters and inputs
	    	// to the hit file
		    sink->out() << "Binary name: " << argv0 << endl;
		    sink->out() << "  Checksum: " << (uint64_t)(EBWT_SEARCH_HASH) << endl;
		    sink->out() << "Ebwt file base: " << infile << endl;
			sink->out() << "Sanity checking: " << (sanityCheck? "on":"off") << endl;
			sink->out() << "Verbose: " << (verbose? "on":"off") << endl;
		    sink->out() << "Queries: " << endl;
		    for(size_t i = 0; i < queries.size(); i++) {
		    	sink->out() << "  " << queries[i] << endl;
		    }
		    params.write(sink->out()); // write searching parameters
		    stats.write(sink->out());  // write searching statistics
		    _t.write(sink->out());     // write timing info
	    }
	    sink->flush();
		if(!outfile.empty()) {
			((ofstream*)fout)->close();
		}
		delete sink;
	}
}

/**
 * main function.  Parses command-line arguments.
 */
int main(int argc, char **argv) {
	string infile;  // read serialized Ebwt from this file
	string query;   // read query string(s) from this file
	vector<string> queries;
	string outfile; // write query results to this file
	parseOptions(argc, argv);
	argv0 = argv[0];
	if(showVersion) {
		// TODO: handle versioning better
		cout << argv0 << " version 0.1 (beta)" << endl;
		cout << "Hash: " << EBWT_SEARCH_HASH << endl;
		return 0;
	}
	Timer _t(cout, "Overall time: ", timing);

	// Get input filename
	if(optind >= argc) {
		cerr << "No input sequence, query, or output file specified!" << endl;
		printUsage(cerr);
		return 1;
	}
	infile = argv[optind++];
	
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
	if(outfile.empty() && binOut) {
		cerr << "When --binOut is specified, an output file must also be specified" << endl;
		printUsage(cerr);
		return 1;
	}

	// Optionally summarize
	if(verbose) {
		cout << "Input ebwt file: \"" << infile << "\"" << endl;
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
	driver<String<Dna, Alloc<> > >("DNA", infile, query, queries, outfile);
    return 0;
}

#if 0
template<typename TStr>
static void prioritySearch(PatternSource<TStr>& patsrc,
						   Ebwt<TStr>& ebwt,
						   vector<String<Dna, Packed<> > >* ss,
						   bool revcomp)

{
	uint32_t patid = 0;
	assert(patsrc.hasMorePatterns());
	
	list<SearchPolicy<TStr>*> search_heirarchy;
	
	
	// ExactSearchWithLowQualityThreePrime uses Landau Vishkin extension,
	// which assumes that a 5 prime hit has a difference on the 3 prime end.
	// Exact end to end matches will screw it up.  This is easy to change, but
	// I haven't done so yet for performance reasons.  For now, prioritySearch
	// screens out all exact matches first.
	
	search_heirarchy.push_back(new ExactSearch<TStr>());
	search_heirarchy.push_back(new ExactSearchWithLowQualityThreePrime<TStr>(*ss));
	
	HitSink hit_sink(cout, true);
	
	// Since I'm manually accumulating hits against both strands for a given
	// policy, I don't want the hit_sink to >> my patId for me.  I'll manage
	// the ids and corresponding hits myself.
	PrettyHitSink report_sink(cout, false, false);
	
	EbwtSearchStats<TStr> stats;
	EbwtSearchParams<TStr> params(hit_sink,
								  stats,
								  MHP_PICK_1_RANDOM,
								  false);
	
    while(patsrc.hasMorePatterns() && patid < (uint32_t)qUpto) 
	{
		
		// Grab a pattern...
		const TStr& pat = patsrc.nextPattern();
		assert(!empty(pat));
		
		TStr pat_rc;
		
		// ...and if we are doing RC, it's RC pattern from the pattern source.
		if (revcomp)
		{
			assert (patsrc.hasMorePatterns() && 
					patsrc.nextIsReverseComplement());
			pat_rc = patsrc.nextPattern();
			assert(!empty(pat_rc));
			//cerr << "\tfwd: " << pat << endl;
			//cerr << "\trev: " << pat_rc << endl;
		}
		
		params.setPatId(patid++);
		
		// Search for hits, stopping at the first (i.e. highest priority) 
		// matching policy for which there is at least one hit.
		for (typename list<SearchPolicy<TStr>*>::iterator pol_itr = search_heirarchy.begin();
			 pol_itr != search_heirarchy.end(); 
			 ++pol_itr)
		{
			SearchPolicy<TStr>* policy = *pol_itr;
			params.setFw(true);
			policy->search(ebwt, stats, params, pat, hit_sink);
			if (revcomp)
			{
				params.setFw(false);
				policy->search(ebwt, stats, params, pat_rc, hit_sink);
			}
			
			vector<Hit>& hits = hit_sink.retainedHits();
			if (hits.size())
			{
				// FIXME: this is a bullshit reporting policy.  We need a real 
				// one.
				Hit& hit = hits[0];
				report_sink.reportHit(hit.h, hit.patId, hit.patSeq, hit.fw, hit.mms, hit.oms);
				break;
			}
		}
		
		hit_sink.clearRetainedHits();
	}
	
	report_sink.finish();
	
	for (typename list<SearchPolicy<TStr>*>::iterator pol_itr = search_heirarchy.begin();
		 pol_itr != search_heirarchy.end(); 
		 ++pol_itr)
	{
		delete *pol_itr;
	}
}
#endif

#endif
