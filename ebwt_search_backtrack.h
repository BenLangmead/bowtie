#ifndef EBWT_SEARCH_BACKTRACK_H_
#define EBWT_SEARCH_BACKTRACK_H_

#include <stdexcept>
#include <seqan/sequence.h>
#include "pat.h"
#include "qual.h"
#include "ebwt_search_util.h"
#include "range.h"
#include "range_source.h"
#include "aligner_metrics.h"
#include "search_globals.h"

/**
 * Class that coordinates quality- and quantity-aware backtracking over
 * some range of a read sequence.
 *
 * The creator can configure the BacktrackManager to treat different
 * stretches of the read differently.
 */
class GreedyDFSRangeSource {

	typedef std::pair<int, int> TIntPair;
	typedef seqan::String<seqan::Dna> DnaString;

public:
	GreedyDFSRangeSource(
			const Ebwt<DnaString>* ebwt,
			const EbwtSearchParams<DnaString>& params,
			const BitPairReference* refs,
			uint32_t qualThresh,  /// max acceptable q-distance
			const int maxBts, /// maximum # backtracks allowed
			uint32_t reportPartials = 0,
			bool reportExacts = true,
			bool reportRanges = false,
			PartialAlignmentManager* partials = NULL,
			String<QueryMutation>* muts = NULL,
			bool verbose = true,
			vector<String<Dna5> >* os = NULL,
			bool considerQuals = true,  // whether to consider quality values when making backtracking decisions
			bool halfAndHalf = false, // hacky way of supporting separate revisitable regions
			bool maqPenalty = true) :
		_refs(refs),
		_qry(NULL),
		_qlen(0),
		_qual(NULL),
		_name(NULL),
		_ebwt(ebwt),
		_params(params),
		_unrevOff(0),
		_1revOff(0),
		_2revOff(0),
		_3revOff(0),
		_maqPenalty(maqPenalty),
		_qualThresh(qualThresh),
		_pairs(NULL),
		_elims(NULL),
		_mms(),
		_refcs(),
		_chars(NULL),
		_reportPartials(reportPartials),
		_reportExacts(reportExacts),
		_reportRanges(reportRanges),
		_partials(partials),
		_muts(muts),
		_os(os),
		_sanity(_os != NULL && _os->size() > 0),
		_considerQuals(considerQuals),
		_halfAndHalf(halfAndHalf),
		_5depth(0),
		_3depth(0),
		_numBts(0),
		_totNumBts(0),
		_maxBts(maxBts),
		_precalcedSideLocus(false),
		_preLtop(),
		_preLbot(),
		_verbose(verbose),
		_ihits(0llu)
	{ }

	~GreedyDFSRangeSource() {
		if(_pairs != NULL) delete[] _pairs;
		if(_elims != NULL) delete[] _elims;
		if(_chars != NULL) delete[] _chars;
	}

	/**
	 * Set a new query read.
	 */
	void setQuery(ReadBuf& r) {
		const bool fw = _params.fw();
		const bool ebwtFw = _ebwt->fw();
		if(ebwtFw) {
			_qry  = fw ? &r.patFw : &r.patRc;
			_qual = fw ? &r.qual  : &r.qualRev;
		} else {
			_qry  = fw ? &r.patFwRev : &r.patRcRev;
			_qual = fw ? &r.qualRev  : &r.qual;
		}
		_name = &r.name;
		// Reset _qlen
		if(length(*_qry) > _qlen) {
			try {
				_qlen = length(*_qry);
				// Resize _pairs
				if(_pairs != NULL) { delete[] _pairs; }
				_pairs = new TIndexOffU[_qlen*_qlen*8];
				// Resize _elims
				if(_elims != NULL) { delete[] _elims; }
				_elims = new uint8_t[_qlen*_qlen];
				memset(_elims, 0, _qlen*_qlen);
				// Resize _chars
				if(_chars != NULL) { delete[] _chars; }
				_chars = new char[_qlen];
				assert(_pairs != NULL && _elims != NULL && _chars != NULL);
			} catch(std::bad_alloc& e) {
				ThreadSafe _ts(&gLock);
				cerr << "Unable to allocate memory for depth-first "
				     << "backtracking search; new length = " << length(*_qry)
				     << endl;
				throw 1;
			}
		} else {
			// New length is less than old length, so there's no need
			// to resize any data structures.
			assert(_pairs != NULL && _elims != NULL && _chars != NULL);
			_qlen = length(*_qry);
		}
		_mms.clear();
		_refcs.clear();
		assert_geq(length(*_qual), _qlen);
		if(_verbose) {
			cout << "setQuery(_qry=" << (*_qry) << ", _qual=" << (*_qual) << ")" << endl;
		}
		// Initialize the random source using new read as part of the
		// seed.
		_color = r.color;
		_seed = r.seed;
		_patid = r.patid;
		_primer = r.primer;
		_trimc = r.trimc;
		_rand.init(r.seed);
	}

	/**
	 * Apply a batch of mutations to this read, possibly displacing a
	 * previous batch of mutations.
	 */
	void setMuts(String<QueryMutation>* muts) {
		if(_muts != NULL) {
			// Undo previous mutations
			assert_gt(length(*_muts), 0);
			undoPartialMutations();
		}
		_muts = muts;
		if(_muts != NULL) {
			assert_gt(length(*_muts), 0);
			applyPartialMutations();
		}
	}

	/**
	 * Set backtracking constraints.
	 */
	void setOffs(uint32_t depth5,   // depth of far edge of hi-half
	             uint32_t depth3,   // depth of far edge of lo-half
	             uint32_t unrevOff, // depth above which we cannot backtrack
	             uint32_t revOff1,  // depth above which we may backtrack just once
	             uint32_t revOff2,  // depth above which we may backtrack just twice
	             uint32_t revOff3)  // depth above which we may backtrack just three times
	{
		_5depth   = depth5;
		_3depth   = depth3;
		assert_geq(depth3, depth5);
		_unrevOff = unrevOff;
		_1revOff  = revOff1;
		_2revOff  = revOff2;
		_3revOff  = revOff3;
	}

	/**
	 * Reset number of backtracks to 0.
	 */
	void resetNumBacktracks() {
		_totNumBts = 0;
	}

	/**
	 * Return number of backtracks since the last time the count was
	 * reset.
	 */
	uint32_t numBacktracks() {
		return _totNumBts;
	}

	/**
	 * Set whether to report exact hits.
	 */
	void setReportExacts(int stratum) {
		_reportExacts = stratum;
	}

	/**
	 * Set the Bowtie index to search against.
	 */
	void setEbwt(const Ebwt<String<Dna> >* ebwt) {
		_ebwt = ebwt;
	}

	/**
	 * Return the current range
	 */
	Range& range() {
		return _curRange;
	}

	/**
	 * Set _qlen.  Don't let it exceed length of query.
	 */
	void setQlen(uint32_t qlen) {
		assert(_qry != NULL);
		_qlen = min<uint32_t>((uint32_t)length(*_qry), qlen);
	}

	/// Return the maximum number of allowed backtracks in a given call
	/// to backtrack()
	uint32_t maxBacktracks() {
		return _maxBts;
	}

	/**
	 * Initiate the recursive backtracking routine starting at the
	 * extreme right-hand side of the pattern.  Use the ftab to match
	 * the first several characters in one chomp, as long as doing so
	 * does not "jump over" any legal backtracking targets.
	 *
	 * Return true iff the HitSink has indicated that we're done with
	 * this read.
	 */
	bool backtrack(uint32_t ham = 0) {
		assert_gt(length(*_qry), 0);
		assert_leq(_qlen, length(*_qry));
		assert_geq(length(*_qual), length(*_qry));
		const Ebwt<String<Dna> >& ebwt = *_ebwt;
		int ftabChars = ebwt._eh._ftabChars;
		int nsInSeed = 0; int nsInFtab = 0;
		if(!tallyNs(nsInSeed, nsInFtab)) {
			// No alignments are possible because of the distribution
			// of Ns in the read in combination with the backtracking
			// constraints.
			return false;
		}
		bool ret;
		// m = depth beyond which ftab must not extend or else we might
		// miss some legitimate paths
		uint32_t m = min<uint32_t>(_unrevOff, (uint32_t)_qlen);
		if(nsInFtab == 0 && m >= (uint32_t)ftabChars) {
			uint32_t ftabOff = calcFtabOff();
			TIndexOffU top = ebwt.ftabHi(ftabOff);
			TIndexOffU bot = ebwt.ftabLo(ftabOff+1);
			if(_qlen == (TIndexOffU)ftabChars && bot > top) {
				// We have a match!
				if(_reportPartials > 0) {
					// Oops - we're trying to find seedlings, so we've
					// gone too far; start again
					ret = backtrack(0,   // depth
					                0,   // top
					                0,   // bot
					                ham,
					                nsInFtab > 0);
				} else {
					// We have a match!
					ret = reportAlignment(0, top, bot, ham);
				}
			} else if (bot > top) {
				// We have an arrow pair from which we can backtrack
				ret = backtrack(ftabChars, // depth
				                top,       // top
				                bot,       // bot
				                ham,
				                nsInFtab > 0);
			} else {
				// The arrows are already closed; give up
				ret = false;
			}
		} else {
			// The ftab *does* extend past the unrevisitable portion;
			// we can't use it in this case, because we might jump past
			// a legitimate mismatch
			ret = backtrack(0,   // depth
			                0,   // top
			                0,   // bot
			                ham,
			                // disable ftab jumping if there is more
			                // than 1 N in it
			                nsInFtab > 0);
		}
		if(finalize()) ret = true;
		return ret;
	}

	/**
	 * If there are any buffered results that have yet to be committed,
	 * commit them.  This happens when looking for partial alignments.
	 */
	bool finalize() {
		bool ret = false;
		if(_reportPartials > 0) {
			// We're in partial alignment mode; take elements of the
			// _partialBuf and install them in the _partials database
			assert(_partials != NULL);
			if(_partialsBuf.size() > 0) {
#ifndef NDEBUG
				for(size_t i = 0; i < _partialsBuf.size(); i++) {
					assert(_partialsBuf[i].repOk(_qualThresh, (uint32_t)_qlen, (*_qual), _maqPenalty));
				}
#endif
				_partials->addPartials(_params.patId(), _partialsBuf);
				_partialsBuf.clear();
				ret = true;
			} else {
				assert(!ret);
			}
		}
		assert_eq(0, _partialsBuf.size());
		return ret;
	}

	/**
	 * Starting at the given "depth" relative to the 5' end, and the
	 * given top and bot indexes (where top=0 and bot=0 means it's up
	 * to us to calculate the initial range), and initial weighted
	 * hamming distance iham, find a hit using randomized, quality-
	 * aware backtracking.
	 */
	bool backtrack(uint32_t depth,
	               TIndexOffU top,
	               TIndexOffU bot,
	               uint32_t iham = 0,
	               bool disableFtab = false)
	{
		HitSinkPerThread& sink = _params.sink();
		_ihits = sink.retainedHits().size();

		// Initiate the recursive, randomized quality-aware backtracker
		// with a stack depth of 0 (no backtracks so far)
		_bailedOnBacktracks = false;
		bool done = backtrack(0, depth, _unrevOff, _1revOff, _2revOff, _3revOff,
		                      top, bot, iham, iham, _pairs, _elims, disableFtab);

		_totNumBts += _numBts;
		_numBts = 0;
		_precalcedSideLocus = false;
		_bailedOnBacktracks = false;
		return done;
	}

	/**
	 * Recursive routine for progressing to the next backtracking
	 * decision given some initial conditions.  If a hit is found, it
	 * is recorded and true is returned.  Otherwise, if there are more
	 * backtracking opportunities, the function will call itself
	 * recursively and return the result.  As soon as there is a
	 * mismatch and no backtracking opportunities, false is returned.
	 */
	bool backtrack(uint32_t  stackDepth, // depth of the recursion stack; = # mismatches so far
	               uint32_t  depth,    // next depth where a post-pair needs to be calculated
	               uint32_t  unrevOff, // depths < unrevOff are unrevisitable
	               uint32_t  oneRevOff,// depths < oneRevOff are 1-revisitable
	               uint32_t  twoRevOff,// depths < twoRevOff are 2-revisitable
	               uint32_t  threeRevOff,// depths < threeRevOff are 3-revisitable
	               TIndexOffU  top,      // top arrow in pair prior to 'depth'
	               TIndexOffU  bot,      // bottom arrow in pair prior to 'depth'
	               uint32_t  ham,      // weighted hamming distance so far
	               uint32_t  iham,     // initial weighted hamming distance
	               TIndexOffU* pairs,    // portion of pairs array to be used for this backtrack frame
	               uint8_t*  elims,    // portion of elims array to be used for this backtrack frame
	               bool disableFtab = false)
	{
		// Can't have already exceeded weighted hamming distance threshold
		assert_leq(stackDepth, depth);
		assert_gt(length(*_qry), 0);
		assert_leq(_qlen, length(*_qry));
		assert_geq(length(*_qual), length(*_qry));
		assert(_qry != NULL);
		assert(_qual != NULL);
		assert(_name != NULL);
		assert(_qlen != 0);
		assert_leq(ham, _qualThresh);
		assert_lt(depth, _qlen); // can't have run off the end of qry
		assert_geq(bot, top);    // could be that both are 0
		assert(pairs != NULL);
		assert(elims != NULL);
		assert_leq(stackDepth, _qlen);
		const Ebwt<String<Dna> >& ebwt = *_ebwt;
		HitSinkPerThread& sink = _params.sink();
		uint64_t prehits = sink.numValidHits();
		if(_halfAndHalf) {
			assert_eq(0, _reportPartials);
			assert_gt(_3depth, _5depth);
		}
		if(_reportPartials) {
			assert(!_halfAndHalf);
		}
		if(_verbose) {
			cout << "  backtrack(stackDepth=" << stackDepth << ", "
			     << "depth=" << depth << ", "
			     << "top=" << top << ", "
			     << "bot=" << bot << ", "
			     << "ham=" << ham << ", "
			     << "iham=" << iham << ", "
			     << "pairs=" << pairs << ", "
			     << "elims=" << (void*)elims << "): \"";
			for(int i = (int)depth - 1; i >= 0; i--) {
				cout << _chars[i];
			}
			cout << "\"" << endl;
		}
		// Do this early on so that we can clear _precalcedSideLocus
		// before we have too many opportunities to bail and leave it
		// 'true'
		SideLocus ltop, lbot;
		if(_precalcedSideLocus) {
			ltop = _preLtop;
			lbot = _preLbot;
			_precalcedSideLocus = false;
		} else if(top != 0 || bot != 0) {
			SideLocus::initFromTopBot(top, bot, ebwt._eh, ebwt._ebwt, ltop, lbot);
		}
		// Check whether we've exceeded any backtracking limit
		if(_halfAndHalf) {
			if(_maxBts > 0 && _numBts == _maxBts) {
				_bailedOnBacktracks = true;
				return false;
			}
			_numBts++;
		}
		// # positions with at least one legal outgoing path
		uint32_t altNum = 0;
		// # positions tied for "best" outgoing qual
		uint32_t eligibleNum = 0;
		// total range-size for all eligibles
		TIndexOffU eligibleSz = 0;
		// If there is just one eligible slot at the moment (a common
		// case), these are its parameters
		uint32_t eli = 0;
		bool     elignore = true; // ignore the el values because they didn't come from a recent override
		TIndexOffU eltop = 0;
		TIndexOffU elbot = 0;
		uint32_t elham = ham;
		char     elchar = 0;
		int      elcint = 0;
		// The lowest quality value associated with any alternative
		// ranges; all alternative ranges with this quality are
		// eligible
		uint8_t lowAltQual = 0xff;
		uint32_t d = depth;
		uint32_t cur = (uint32_t)_qlen - d - 1; // current offset into _qry
		while(cur < _qlen) {
			// Try to advance further given that
			if(_verbose) {
				cout << "    cur=" << cur << " \"";
				for(int i = (int)d - 1; i >= 0; i--) {
					cout << _chars[i];
				}
				cout << "\"";
			}

			// If we're searching for a half-and-half solution, then
			// enforce the boundary-crossing constraints here.
			if(_halfAndHalf && !hhCheckTop(stackDepth, d, iham, _mms, prehits)) {
				return false;
			}

			bool curIsEligible = false;
			// Reset eligibleNum and eligibleSz if there are any
			// eligible pairs discovered at this spot
			bool curOverridesEligible = false;
			// Determine whether ranges at this location are
			// candidates for backtracking
			int c = (int)(*_qry)[cur];
			assert_leq(c, 4);
			uint8_t q = qualAt(cur);
			// The current query position is a legit alternative if it a) is
			// not in the unrevisitable region, and b) the quality ceiling (if
			// one exists) is not exceeded
			bool curIsAlternative =
				(d >= unrevOff) &&
			    (!_considerQuals ||
			     (ham + mmPenalty(_maqPenalty, q) <= _qualThresh));
			if(curIsAlternative) {
				if(_considerQuals) {
					// Is it the best alternative?
					if(q < lowAltQual) {
						// Ranges at this depth in this backtracking frame are
						// eligible, unless we learn otherwise.  Ranges previously
						// thought to be eligible are not any longer.
						curIsEligible = true;
						curOverridesEligible = true;
					} else if(q == lowAltQual) {
						// Ranges at this depth in this backtracking frame
						// are eligible, unless we learn otherwise
						curIsEligible = true;
					}
				} else {
					// When quality values are not considered, all positions
					// are eligible
					curIsEligible = true;
				}
			}
			if(curIsEligible) assert(curIsAlternative);
			if(curOverridesEligible) assert(curIsEligible);
			if(curIsAlternative && !curIsEligible) {
				assert_gt(eligibleSz, 0);
				assert_gt(eligibleNum, 0);
			}
			if(_verbose) {
				cout << " alternative: " << curIsAlternative;
				cout << ", eligible: " << curIsEligible;
				if(curOverridesEligible) cout << "(overrides)";
				cout << endl;
			}

			// If c is 'N', then it's guaranteed to be a mismatch
			if(c == 4 && d > 0) {
				// Force the 'else if(curIsAlternative)' branch below
				top = bot = 1;
			} else if(c == 4) {
				// We'll take the 'if(top == 0 && bot == 0)' branch below
				assert_eq(0, top);
				assert_eq(0, bot);
			}
			// Calculate the ranges for this position
			if(top == 0 && bot == 0) {
				// Calculate first quartet of ranges using the _fchr[]
				// array
				               pairs[0 + 0] = ebwt._fchr[0];
				pairs[0 + 4] = pairs[1 + 0] = ebwt._fchr[1];
				pairs[1 + 4] = pairs[2 + 0] = ebwt._fchr[2];
				pairs[2 + 4] = pairs[3 + 0] = ebwt._fchr[3];
				pairs[3 + 4]                = ebwt._fchr[4];
				// Update top and bot
				if(c < 4) {
					top = pairTop(pairs, d, c); bot = pairBot(pairs, d, c);
					assert_geq(bot, top);
				}
			} else if(curIsAlternative) {
				// Clear pairs
				memset(&pairs[d*8], 0, 8 * OFF_SIZE);
				// Calculate next quartet of ranges
				ebwt.mapLFEx(ltop, lbot, &pairs[d*8], &pairs[(d*8)+4]);
				// Update top and bot
				if(c < 4) {
					top = pairTop(pairs, d, c); bot = pairBot(pairs, d, c);
					assert_geq(bot, top);
				}
			} else {
				// This query character is not even a legitimate
				// alternative (because backtracking here would blow
				// our mismatch quality budget), so no need to do the
				// bookkeeping for the entire quartet, just do c
				if(c < 4) {
					if(top+1 == bot) {
						bot = top = ebwt.mapLF1(top, ltop, c);
						if(bot != OFF_MASK) bot++;
					} else {
						top = ebwt.mapLF(ltop, c); bot = ebwt.mapLF(lbot, c);
						assert_geq(bot, top);
					}
				}
			}
			if(top != bot) {
				// Calculate loci from row indices; do it now so that
				// those prefetches are fired off as soon as possible.
				// This eventually calls SideLocus.initfromRow().
				SideLocus::initFromTopBot(top, bot, ebwt._eh, ebwt._ebwt, ltop, lbot);
			}
			// Update the elim array
			eliminate(elims, d, c);

			if(curIsAlternative) {
				// Given the just-calculated range quartet, update
				// elims, altNum, eligibleNum, eligibleSz
				for(int i = 0; i < 4; i++) {
					if(i == c) continue;
					assert_leq(pairTop(pairs, d, i), pairBot(pairs, d, i));
					TIndexOffU spread = pairSpread(pairs, d, i);
					if(spread == 0) {
						// Indicate this char at this position is
						// eliminated as far as this backtracking frame is
						// concerned, since its range is empty
						elims[d] |= (1 << i);
						assert_lt(elims[d], 16);
					}
					if(spread > 0 && ((elims[d] & (1 << i)) == 0)) {
						// This char at this position is an alternative
						if(curIsEligible) {
							if(curOverridesEligible) {
								// Only now that we know there is at least
								// one potential backtrack target at this
								// most-eligible position should we reset
								// these eligibility parameters
								lowAltQual = q;
								eligibleNum = 0;
								eligibleSz = 0;
								curOverridesEligible = false;
								// Remember these parameters in case
								// this turns out to be the only
								// eligible target
								eli = d;
								eltop = pairTop(pairs, d, i);
								elbot = pairBot(pairs, d, i);
								assert_eq(elbot-eltop, spread);
								elham = mmPenalty(_maqPenalty, q);
								elchar = "acgt"[i];
								elcint = i;
								elignore = false;
							}
							eligibleSz += spread;
							eligibleNum++;
						}
						assert_gt(eligibleSz, 0);
						assert_gt(eligibleNum, 0);
						altNum++;
					}
				}
			}
			if(altNum > 0) {
				assert_gt(eligibleSz, 0);
				assert_gt(eligibleNum, 0);
			}
			assert_leq(eligibleNum, eligibleSz);
			assert_leq(eligibleNum, altNum);
			assert_lt(elims[d], 16);
			assert(sanityCheckEligibility(depth, d, unrevOff, lowAltQual, eligibleSz, eligibleNum, pairs, elims));

			// Achieved a match, but need to keep going
			bool backtrackDespiteMatch = false;
			bool reportedPartial = false;
			if(cur == 0 &&  // we've consumed the entire pattern
			   top < bot && // there's a hit to report
			   stackDepth < _reportPartials && // not yet used up our mismatches
			   _reportPartials > 0)  // there are still legel backtracking targets
			{
				assert(!_halfAndHalf);
				if(altNum > 0) backtrackDespiteMatch = true;
				if(stackDepth > 0) {
					// This is a legit seedling; report it
					reportPartial(stackDepth);
					reportedPartial = true;
				}
				// Now continue on to find legitimate seedlings with
				// more mismatches than this one
			}
			// Check whether we've obtained an exact alignment when
			// we've been instructed not to report exact alignments
			bool invalidExact = false;
			if(cur == 0 && stackDepth == 0 && bot > top && !_reportExacts) {
				invalidExact = true;
				backtrackDespiteMatch = true;
			}
			// Set this to true if the only way to make legal progress
			// is via one or more additional backtracks.  This is
			// helpful in half-and-half mode.
			bool mustBacktrack = false;
			bool invalidHalfAndHalf = false;
			if(_halfAndHalf) {
				ASSERT_ONLY(uint32_t lim = (_3revOff == _2revOff)? 2 : 3);
				if((d == (_5depth-1)) && top < bot) {
					// We're crossing the boundary separating the hi-half
					// from the non-seed portion of the read.
					// We should induce a mismatch if we haven't mismatched
					// yet, so that we don't waste time pursuing a match
					// that was covered by a previous phase
					assert_eq(0, _reportPartials);
					assert_leq(stackDepth, lim-1);
					invalidHalfAndHalf = (stackDepth == 0);
					if(stackDepth == 0 && altNum > 0) {
						backtrackDespiteMatch = true;
						mustBacktrack = true;
					} else if(stackDepth == 0) {
						// We're returning from the bottommost frame
						// without having found any hits; let's
						// sanity-check that there really aren't any
						return false;
					}
				}
				else if((d == (_3depth-1)) && top < bot) {
					// We're crossing the boundary separating the lo-half
					// from the non-seed portion of the read
					assert_eq(0, _reportPartials);
					assert_leq(stackDepth, lim);
					assert_gt(stackDepth, 0);
					// Count the mismatches in the lo and hi halves
					uint32_t loHalfMms = 0, hiHalfMms = 0;
					assert_geq(_mms.size(), stackDepth);
					for(size_t i = 0; i < stackDepth; i++) {
						uint32_t d = (uint32_t)_qlen - _mms[i] - 1;
						if     (d < _5depth) hiHalfMms++;
						else if(d < _3depth) loHalfMms++;
						else assert(false);
					}
					assert_leq(loHalfMms + hiHalfMms, lim);
					invalidHalfAndHalf = (loHalfMms == 0 || hiHalfMms == 0);
					if((stackDepth < 2 || invalidHalfAndHalf) && altNum > 0) {
						// We backtracked fewer times than necessary;
						// force a backtrack
						mustBacktrack = true;
						backtrackDespiteMatch = true;
					} else if(stackDepth < 2) {
						return false;
					}
				}
				if(d < _5depth-1) {
					assert_leq(stackDepth, lim-1);
				}
				else if(d >= _5depth && d < _3depth-1) {
					assert_gt(stackDepth, 0);
					assert_leq(stackDepth, lim);
				}
			}
			// This is necessary for the rare case where we're about
			// to declare success because bot > top and we've consumed
			// the final character, but all hits between top and bot
			// are spurious.  This check ensures that we keep looking
			// for non-spurious hits in that case.
			if(cur == 0 &&            // we made it to the left-hand-side of the read
			   bot > top &&           // there are alignments to report
			   !invalidHalfAndHalf && // alignment isn't disqualified by half-and-half requirement
			   !invalidExact &&       // alignment isn't disqualified by no-exact-hits setting
			   !reportedPartial)      // for when it's a partial alignment we've already reported
			{
				bool ret = reportAlignment(stackDepth, top, bot, ham);
				if(!ret) {
					// reportAlignment returned false, so enter the
					// backtrack loop and keep going
					top = bot;
				} else {
					// reportAlignment returned true, so stop
					return true;
				}
			}
			//
			// Mismatch with alternatives
			//
			while((top == bot || backtrackDespiteMatch) && altNum > 0) {
				if(_verbose) cout << "    top (" << top << "), bot ("
				                 << bot << ") with " << altNum
				                 << " alternatives, eligible: "
				                 << eligibleNum << ", " << eligibleSz
				                 << endl;
				assert_gt(eligibleSz, 0);
				assert_gt(eligibleNum, 0);
				// Mismatch!  Must now choose where we are going to
				// take our quality penalty.  We can only look as far
				// back as our last decision point.
				assert(sanityCheckEligibility(depth, d, unrevOff, lowAltQual, eligibleSz, eligibleNum, pairs, elims));
				// Pick out the arrow pair we selected and target it
				// for backtracking
				ASSERT_ONLY(uint32_t eligiblesVisited = 0);
				size_t i = d, j = 0;
				assert_geq(i, depth);
				TIndexOffU bttop = 0;
				TIndexOffU btbot = 0;
				uint32_t btham = ham;
				char     btchar = 0;
				int      btcint = 0;
				uint32_t icur = 0;
				// The common case is that eligibleSz == 1
				if(eligibleNum > 1 || elignore) {
					ASSERT_ONLY(bool foundTarget = false);
					// Walk from left to right
					for(; i >= depth; i--) {
						assert_geq(i, unrevOff);
						icur = (uint32_t)(_qlen - i - 1); // current offset into _qry
						uint8_t qi = qualAt(icur);
						assert_lt(elims[i], 16);
						if((qi == lowAltQual || !_considerQuals) && elims[i] != 15) {
							// This is the leftmost eligible position with at
							// least one remaining backtrack target
							TIndexOffU posSz = 0;
							// Add up the spreads for A, C, G, T
							for(j = 0; j < 4; j++) {
								if((elims[i] & (1 << j)) == 0) {
									assert_gt(pairSpread(pairs, i, j), 0);
									posSz += pairSpread(pairs, i, j);
								}
							}
							// Generate a random number
							assert_gt(posSz, 0);
							uint32_t r = _rand.nextU32() % posSz;
							for(j = 0; j < 4; j++) {
								if((elims[i] & (1 << j)) == 0) {
									// This range has not been eliminated
									ASSERT_ONLY(eligiblesVisited++);
									uint32_t spread = pairSpread(pairs, i, j);
									if(r < spread) {
										// This is our randomly-selected
										// backtrack target
										ASSERT_ONLY(foundTarget = true);
										bttop = pairTop(pairs, i, j);
										btbot = pairBot(pairs, i, j);
										btham += mmPenalty(_maqPenalty, qi);
										btcint = (uint32_t)j;
										btchar = "acgt"[j];
										assert_leq(btham, _qualThresh);
										break; // found our target; we can stop
									}
									r -= spread;
								}
							}
							assert(foundTarget);
							break; // escape left-to-right walk
						}
					}
					assert_leq(i, d);
					assert_lt(j, 4);
					assert_leq(eligiblesVisited, eligibleNum);
					assert(foundTarget);
					assert_neq(0, btchar);
					assert_gt(btbot, bttop);
					assert_leq(btbot-bttop, eligibleSz);
				} else {
					// There was only one eligible target; we can just
					// copy its parameters
					assert_eq(1, eligibleNum);
					assert(!elignore);
					i = eli;
					bttop = eltop;
					btbot = elbot;
					btham += elham;
					j = btcint = elcint;
					btchar = elchar;
					assert_neq(0, btchar);
					assert_gt(btbot, bttop);
					assert_leq(btbot-bttop, eligibleSz);
				}
				// This is the earliest that we know what the next top/
				// bot combo is going to be
				SideLocus::initFromTopBot(bttop, btbot,
				                          ebwt._eh, ebwt._ebwt,
				                          _preLtop, _preLbot);
				icur = (uint32_t)(_qlen - i - 1); // current offset into _qry
				// Slide over to the next backtacking frame within
				// pairs and elims; won't interfere with our frame or
				// any of our parents' frames
				TIndexOffU *newPairs = pairs + (_qlen*8);
				uint8_t  *newElims = elims + (_qlen);
				// If we've selected a backtracking target that's in
				// the 1-revisitable region, then we ask the recursive
				// callee to consider the 1-revisitable region as also
				// being unrevisitable (since we just "used up" all of
				// our visits)
				uint32_t btUnrevOff  = unrevOff;
				uint32_t btOneRevOff = oneRevOff;
				uint32_t btTwoRevOff = twoRevOff;
				uint32_t btThreeRevOff = threeRevOff;
				assert_geq(i, unrevOff);
				assert_geq(oneRevOff, unrevOff);
				assert_geq(twoRevOff, oneRevOff);
				assert_geq(threeRevOff, twoRevOff);
				if(i < oneRevOff) {
					// Extend unrevisitable region to include former 1-
					// revisitable region
					btUnrevOff = oneRevOff;
					// Extend 1-revisitable region to include former 2-
					// revisitable region
					btOneRevOff = twoRevOff;
					// Extend 2-revisitable region to include former 3-
					// revisitable region
					btTwoRevOff = threeRevOff;
				}
				else if(i < twoRevOff) {
					// Extend 1-revisitable region to include former 2-
					// revisitable region
					btOneRevOff = twoRevOff;
					// Extend 2-revisitable region to include former 3-
					// revisitable region
					btTwoRevOff = threeRevOff;
				}
				else if(i < threeRevOff) {
					// Extend 2-revisitable region to include former 3-
					// revisitable region
					btTwoRevOff = threeRevOff;
				}
				// Note the character that we're backtracking on in the
				// mm array:
				if(_mms.size() <= stackDepth) {
					assert_eq(_mms.size(), stackDepth);
					_mms.push_back(icur);
				} else {
					_mms[stackDepth] = icur;
				}
				assert_eq(1, dna4Cat[(int)btchar]);
				if(_refcs.size() <= stackDepth) {
					assert_eq(_refcs.size(), stackDepth);
					_refcs.push_back(btchar);
				} else {
					_refcs[stackDepth] = btchar;
				}
#ifndef NDEBUG
				for(uint32_t j = 0; j < stackDepth; j++) {
					assert_neq(_mms[j], icur);
				}
#endif
				_chars[i] = btchar;
				assert_leq(i+1, _qlen);
				bool ret;
				if(i+1 == _qlen) {
					ret = reportAlignment(stackDepth+1, bttop, btbot, btham);
				} else if(_halfAndHalf &&
				          !disableFtab &&
				          _2revOff == _3revOff &&
				          i+1 < (uint32_t)ebwt._eh._ftabChars &&
				          (uint32_t)ebwt._eh._ftabChars <= _5depth)
				{
					// The ftab doesn't extend past the unrevisitable portion,
					// so we can go ahead and use it
					// Rightmost char gets least significant bit-pairs
					int ftabChars = ebwt._eh._ftabChars;
					TIndexOffU ftabOff = (TIndexOffU)(int)(*_qry)[_qlen - ftabChars];
					assert_lt(ftabOff, 4);
					assert_lt(ftabOff, ebwt._eh._ftabLen-1);
					for(int j = ftabChars - 1; j > 0; j--) {
						ftabOff <<= 2;
						if(_qlen-j == icur) {
							ftabOff |= btcint;
						} else {
							assert_lt((int)(*_qry)[_qlen-j], 4);
							ftabOff |= (int)(*_qry)[_qlen-j];
						}
						assert_lt(ftabOff, ebwt._eh._ftabLen-1);
					}
					assert_lt(ftabOff, ebwt._eh._ftabLen-1);
					TIndexOffU ftabTop = ebwt.ftabHi(ftabOff);
					TIndexOffU ftabBot = ebwt.ftabLo(ftabOff+1);
					assert_geq(ftabBot, ftabTop);
					if(ftabTop == ftabBot) {
						ret = false;
					} else {
						assert(!_precalcedSideLocus);
						assert_leq(iham, _qualThresh);
						ret = backtrack(stackDepth+1,
						                ebwt._eh._ftabChars,
						                btUnrevOff,  // new unrevisitable boundary
						                btOneRevOff, // new 1-revisitable boundary
						                btTwoRevOff, // new 2-revisitable boundary
						                btThreeRevOff, // new 3-revisitable boundary
						                ftabTop,  // top arrow in range prior to 'depth'
						                ftabBot,  // bottom arrow in range prior to 'depth'
						                btham,  // weighted hamming distance so far
						                iham,   // initial weighted hamming distance
						                newPairs,
						                newElims);
					}
				} else {
					// We already called initFromTopBot for the range
					// we're going to continue from
					_precalcedSideLocus = true;
					assert_leq(iham, _qualThresh);
					// Continue from selected alternative range
					ret = backtrack(stackDepth+1,// added 1 mismatch to alignment
					                (uint32_t)i+1, // start from next position after
					                btUnrevOff,  // new unrevisitable boundary
					                btOneRevOff, // new 1-revisitable boundary
					                btTwoRevOff, // new 2-revisitable boundary
					                btThreeRevOff, // new 3-revisitable boundary
					                bttop,  // top arrow in range prior to 'depth'
					                btbot,  // bottom arrow in range prior to 'depth'
					                btham,  // weighted hamming distance so far
					                iham,   // initial weighted hamming distance
					                newPairs,
					                newElims);
				}
				if(ret) {
					assert_gt(sink.numValidHits(), prehits);
					return true; // return, signaling that we're done
				}
				if(_bailedOnBacktracks ||
				  (_halfAndHalf && (_maxBts > 0) && (_numBts >= _maxBts)))
				{
					_bailedOnBacktracks = true;
					return false;
				}
				// No hit was reported; update elims[], eligibleSz,
				// eligibleNum, altNum
				_chars[i] = (*_qry)[icur];
				assert_neq(15, elims[i]);
				ASSERT_ONLY(uint8_t oldElim = elims[i]);
				elims[i] |= (1 << j);
				assert_lt(elims[i], 16);
				assert_gt(elims[i], oldElim);
				eligibleSz -= (btbot-bttop);
				eligibleNum--;
				elignore = true;
				assert_geq(eligibleNum, 0);
				altNum--;
				assert_geq(altNum, 0);
				if(altNum == 0) {
					// No alternative backtracking points; all legal
					// backtracking targets have been exhausted
					assert_eq(0, altNum);
					assert_eq(0, eligibleSz);
					assert_eq(0, eligibleNum);
					return false;
				}
				else if(eligibleNum == 0 && _considerQuals) {
					// Find the next set of eligible backtrack points
					// by re-scanning this backtracking frame (from
					// 'depth' up to 'd')
					lowAltQual = 0xff;
					for(size_t k = d; k >= depth && k <= _qlen; k--) {
						size_t kcur = _qlen - k - 1; // current offset into _qry
						uint8_t kq = qualAt(kcur);
						if(k < unrevOff) break; // already visited all revisitable positions
						bool kCurIsAlternative = (ham + mmPenalty(_maqPenalty, kq) <= _qualThresh);
						bool kCurOverridesEligible = false;
						if(kCurIsAlternative) {
							if(kq < lowAltQual) {
								// This target is more eligible than
								// any targets that came before, so we
								// set it to supplant/override them
								kCurOverridesEligible = true;
							}
							if(kq <= lowAltQual) {
								// Position is eligible
								for(int l = 0; l < 4; l++) {
									if((elims[k] & (1 << l)) == 0) {
										// Not yet eliminated
										TIndexOffU spread = pairSpread(pairs, k, l);
										if(kCurOverridesEligible) {
											// Clear previous eligible results;
											// this one's better
											lowAltQual = kq;
											kCurOverridesEligible = false;
											// Keep these parameters in
											// case this target turns
											// out to be the only
											// eligible target and we
											// can avoid having to
											// recalculate them
											eligibleNum = 0;
											eligibleSz = 0;
											eli = (uint32_t)k;
											eltop = pairTop(pairs, k, l);
											elbot = pairBot(pairs, k, l);
											assert_eq(elbot-eltop, spread);
											elham = mmPenalty(_maqPenalty, kq);
											elchar = "acgt"[l];
											elcint = l;
											elignore = false;
										}
										eligibleNum++;
										assert_gt(spread, 0);
										eligibleSz += spread;
									}
								}
							}
						}
					}
				}
				assert_gt(eligibleNum, 0);
				assert_leq(eligibleNum, altNum);
				assert_gt(eligibleSz, 0);
				assert_geq(eligibleSz, eligibleNum);
				assert(sanityCheckEligibility(depth, d, unrevOff, lowAltQual, eligibleSz, eligibleNum, pairs, elims));
				// Try again
			} // while(top == bot && altNum > 0)
			if(mustBacktrack || invalidHalfAndHalf || invalidExact) {
				return false;
			}
			// Mismatch with no alternatives
			if(top == bot && altNum == 0) {
				assert_eq(0, altNum);
				assert_eq(0, eligibleSz);
				assert_eq(0, eligibleNum);
				return false;
			}
			// Match!
			_chars[d] = (*_qry)[cur];
			d++; cur--;
		} // while(cur < _qlen)
		assert_eq(0xffffffff, cur);
		assert_gt(bot, top);
		if(_reportPartials > 0) {
			// Stack depth should not exceed given hamming distance
			assert_leq(stackDepth, _reportPartials);
		}
		bool ret = false;
		if(stackDepth >= _reportPartials) {
			ret = reportAlignment(stackDepth, top, bot, ham);
		}
		return ret;
	}

	/**
	 * Pretty print a hit along with the backtracking constraints.
	 */
	void printHit(const Hit& h) {
		::printHit(*_os, h, *_qry, _qlen, _unrevOff, _1revOff, _2revOff, _3revOff, _ebwt->fw());
	}

	/**
	 * Return true iff we're enforcing a half-and-half constraint
	 * (forced edits in both seed halves).
	 */
	bool halfAndHalf() const {
		return _halfAndHalf;
	}

protected:

	/**
	 * Return true iff we're OK to continue after considering which
	 * half-seed boundary we're passing through, together with the
	 * number of mismatches accumulated so far.  Return false if we
	 * should stop because a half-and-half constraint is violated.  If
	 * we're not currently passing a half-seed boundary, just return
	 * true.
	 */
	bool hhCheck(uint32_t stackDepth, uint32_t depth,
	             const std::vector<uint32_t>& mms, bool empty)
	{
		ASSERT_ONLY(uint32_t lim = (_3revOff == _2revOff)? 2 : 3);
		if((depth == (_5depth-1)) && !empty) {
			// We're crossing the boundary separating the hi-half
			// from the non-seed portion of the read.
			// We should induce a mismatch if we haven't mismatched
			// yet, so that we don't waste time pursuing a match
			// that was covered by a previous phase
			assert_eq(0, _reportPartials);
			assert_leq(stackDepth, lim-1);
			return stackDepth > 0;
		} else if((depth == (_3depth-1)) && !empty) {
			// We're crossing the boundary separating the lo-half
			// from the non-seed portion of the read
			assert_eq(0, _reportPartials);
			assert_leq(stackDepth, lim);
			assert_gt(stackDepth, 0);
			// Count the mismatches in the lo and hi halves
			uint32_t loHalfMms = 0, hiHalfMms = 0;
			for(size_t i = 0; i < stackDepth; i++) {
				uint32_t depth = (uint32_t)(_qlen - mms[i] - 1);
				if     (depth < _5depth) hiHalfMms++;
				else if(depth < _3depth) loHalfMms++;
				else assert(false);
			}
			assert_leq(loHalfMms + hiHalfMms, lim);
			bool invalidHalfAndHalf = (loHalfMms == 0 || hiHalfMms == 0);
			return (stackDepth >= 2 && !invalidHalfAndHalf);
		}
		if(depth < _5depth-1) {
			assert_leq(stackDepth, lim-1);
		}
		else if(depth >= _5depth && depth < _3depth-1) {
			assert_gt(stackDepth, 0);
			assert_leq(stackDepth, lim);
		}
		return true;
	}

	/**
	 * Calculate the stratum of the partial (or full) alignment
	 * currently under consideration.  Stratum is equal to the number
	 * of mismatches in the seed portion of the alignment.
	 */
	int calcStratum(const std::vector<TIndexOffU>& mms, uint32_t stackDepth) {
		int stratum = 0;
		for(size_t i = 0; i < stackDepth; i++) {
			if(mms[i] >= (_qlen - _3revOff)) {
				// This mismatch falls within the seed; count it
				// toward the stratum to report
				stratum++;
				// Don't currently support more than 3
				// mismatches in the seed
				assert_leq(stratum, 3);
			}
		}
		return stratum;
	}

	/**
	 * Mark character c at depth d as being eliminated with respect to
	 * future backtracks.
	 */
	void eliminate(uint8_t *elims, uint32_t d, int c) {
		if(c < 4) {
			elims[d] = (1 << c);
			assert_gt(elims[d], 0);
			assert_lt(elims[d], 16);
		} else {
			elims[d] = 0;
		}
		assert_lt(elims[d], 16);
	}

	/**
	 * Return true iff the state of the backtracker as encoded by
	 * stackDepth, d and iham is compatible with the current half-and-
	 * half alignment mode.  prehits is for sanity checking when
	 * bailing.
	 */
	bool hhCheckTop(uint32_t stackDepth,
	                uint32_t d,
	                uint32_t iham,
	                const std::vector<TIndexOffU>& mms,
	                uint64_t prehits = 0xffffffffffffffffllu)
	{
		assert_eq(0, _reportPartials);
		// Crossing from the hi-half into the lo-half
		if(d == _5depth) {
			if(_3revOff == _2revOff) {
				// Total of 2 mismatches allowed: 1 hi, 1 lo
				// The backtracking logic should have prevented us from
				// backtracking more than once into this region
				assert_leq(stackDepth, 1);
				// Reject if we haven't encountered mismatch by this point
				if(stackDepth == 0) {
					return false;
				}
			} else { // if(_3revOff != _2revOff)
				// Total of 3 mismatches allowed: 1 hi, 1 or 2 lo
				// The backtracking logic should have prevented us from
				// backtracking more than twice into this region
				assert_leq(stackDepth, 2);
				// Reject if we haven't encountered mismatch by this point
				if(stackDepth < 1) {
					return false;
				}
			}
		} else if(d == _3depth) {
			// Crossing from lo-half to outside of the seed
			if(_3revOff == _2revOff) {
				// Total of 2 mismatches allowed: 1 hi, 1 lo
				// The backtracking logic should have prevented us from
				// backtracking more than twice within this region
				assert_leq(stackDepth, 2);
				// Must have encountered two mismatches by this point
				if(stackDepth < 2) {
					// We're returning from the bottommost frame
					// without having found any hits; let's
					// sanity-check that there really aren't any
					return false;
				}
			} else { // if(_3revOff != _2revOff)
				// Total of 3 mismatches allowed: 1 hi, 1 or 2 lo
				// Count the mismatches in the lo and hi halves
				int loHalfMms = 0, hiHalfMms = 0;
				assert_geq(mms.size(), stackDepth);
				for(size_t i = 0; i < stackDepth; i++) {
					TIndexOffU d = (TIndexOffU)(_qlen - mms[i] - 1);
					if     (d < _5depth) hiHalfMms++;
					else if(d < _3depth) loHalfMms++;
					else assert(false);
				}
				assert_leq(loHalfMms + hiHalfMms, 3);
				assert_gt(hiHalfMms, 0);
				if(loHalfMms == 0) {
					// We're returning from the bottommost frame
					// without having found any hits; let's
					// sanity-check that there really aren't any
					return false;
				}
				assert_geq(stackDepth, 2);
				// The backtracking logic should have prevented us from
				// backtracking more than twice within this region
				assert_leq(stackDepth, 3);
			}
		} else {
			// We didn't just cross a boundary, so do an in-between check
			if(d >= _5depth) {
				assert_geq(stackDepth, 1);
			} else if(d >= _3depth) {
				assert_geq(stackDepth, 2);
			}
		}
		return true;
	}

	/**
	 * Return the Phred quality value for the most likely base at
	 * offset 'off' in the read.
	 */
	inline uint8_t qualAt(size_t off) {
		return phredCharToPhredQual((*_qual)[off]);
	}

	/// Get the top offset for character c at depth d
	inline TIndexOffU pairTop(TIndexOffU* pairs, size_t d, size_t c) {
		return pairs[d*8 + c + 0];
	}

	/// Get the bot offset for character c at depth d
	inline TIndexOffU pairBot(TIndexOffU* pairs, size_t d, size_t c) {
		return pairs[d*8 + c + 4];
	}

	/// Get the spread between the bot and top offsets for character c
	/// at depth d
	inline TIndexOffU pairSpread(TIndexOffU* pairs, size_t d, size_t c) {
		assert_geq(pairBot(pairs, d, c), pairTop(pairs, d, c));
		return pairBot(pairs, d, c) - pairTop(pairs, d, c);
	}

	/**
	 * Tally how many Ns occur in the seed region and in the ftab-
	 * jumpable region of the read.  Check whether the mismatches
	 * induced by the Ns already violates the current policy.  Return
	 * false if the policy is already violated, true otherwise.
	 */
	bool tallyNs(int& nsInSeed, int& nsInFtab) {
		const Ebwt<String<Dna> >& ebwt = *_ebwt;
		int ftabChars = ebwt._eh._ftabChars;
		// Count Ns in the seed region of the read and short-circuit if
		// the configuration of Ns guarantees that there will be no
		// valid alignments given the backtracking constraints.
		for(size_t i = 0; i < _3revOff; i++) {
			if((int)(*_qry)[_qlen-i-1] == 4) {
				nsInSeed++;
				if(nsInSeed == 1) {
					if(i < _unrevOff) {
						return false; // Exceeded mm budget on Ns alone
					}
				} else if(nsInSeed == 2) {
					if(i < _1revOff) {
						return false; // Exceeded mm budget on Ns alone
					}
				} else if(nsInSeed == 3) {
					if(i < _2revOff) {
						return false; // Exceeded mm budget on Ns alone
					}
				} else {
					assert_gt(nsInSeed, 3);
					return false;     // Exceeded mm budget on Ns alone
				}
			}
		}
		// Calculate the number of Ns there are in the region that
		// would get jumped over if the ftab were used.
		for(size_t i = 0; i < (size_t)ftabChars && i < _qlen; i++) {
			if((int)(*_qry)[_qlen-i-1] == 4) nsInFtab++;
		}
		return true;
	}

	/**
	 * Calculate the offset into the ftab for the rightmost 'ftabChars'
	 * characters of the current query. Rightmost char gets least
	 * significant bit-pair.
	 */
	uint32_t calcFtabOff() {
		const Ebwt<String<Dna> >& ebwt = *_ebwt;
		int ftabChars = ebwt._eh._ftabChars;
		uint32_t ftabOff = (*_qry)[_qlen - ftabChars];
		assert_lt(ftabOff, 4);
		assert_lt(ftabOff, ebwt._eh._ftabLen-1);
		for(int i = ftabChars - 1; i > 0; i--) {
			ftabOff <<= 2;
			assert_lt((uint32_t)(*_qry)[_qlen-i], 4);
			ftabOff |= (uint32_t)(*_qry)[_qlen-i];
			assert_lt(ftabOff, ebwt._eh._ftabLen-1);
		}
		assert_lt(ftabOff, ebwt._eh._ftabLen-1);
		return ftabOff;
	}

	/**
	 * Mutate the _qry string according to the contents of the _muts
	 * array, which represents a partial alignment.
	 */
	void applyPartialMutations() {
		if(_muts == NULL) {
			// No mutations to apply
			return;
		}
		for(size_t i = 0; i < length(*_muts); i++) {
			const QueryMutation& m = (*_muts)[i];
			assert_lt(m.pos, _qlen);
			assert_leq(m.oldBase, 4);
			assert_lt(m.newBase, 4);
			assert_neq(m.oldBase, m.newBase);
			assert_eq((uint32_t)((*_qry)[m.pos]), (uint32_t)m.oldBase);
			(*_qry)[m.pos] = (Dna5)(int)m.newBase; // apply it
		}
	}

	/**
	 * Take partial-alignment mutations present in the _muts list and
	 * place them on the _mm list so that they become part of the
	 * reported alignment.
	 */
	void promotePartialMutations(int stackDepth) {
		if(_muts == NULL) {
			// No mutations to undo
			return;
		}
		size_t numMuts = length(*_muts);
		assert_leq(numMuts, _qlen);
		for(size_t i = 0; i < numMuts; i++) {
			// Entries in _mms[] are in terms of offset into
			// _qry - not in terms of offset from 3' or 5' end
			assert_lt(stackDepth + i, _qlen);
			// All partial-alignment mutations should fall
			// within bounds
			assert_lt((*_muts)[i].pos, _qlen);
			// All partial-alignment mutations should fall
			// within unrevisitable region
			assert_lt(_qlen - (*_muts)[i].pos - 1, _unrevOff);
#ifndef NDEBUG
			// Shouldn't be any overlap between mismatched positions
			// and positions that mismatched in the partial alignment.
			for(size_t j = 0; j < stackDepth + i; j++) {
				assert_neq(_mms[j], (uint32_t)(*_muts)[i].pos);
			}
#endif
			if(_mms.size() <= stackDepth + i) {
				assert_eq(_mms.size(), stackDepth + i);
				_mms.push_back((*_muts)[i].pos);
			} else {
				_mms[stackDepth + i] = (*_muts)[i].pos;
			}
			if(_refcs.size() <= stackDepth + i) {
				assert_eq(_refcs.size(), stackDepth + i);
				_refcs.push_back("ACGT"[(*_muts)[i].newBase]);
			} else {
				_refcs[stackDepth + i] = "ACGT"[(*_muts)[i].newBase];
			}
		}
	}

	/**
	 * Undo mutations to the _qry string, returning it to the original
	 * read.
	 */
	void undoPartialMutations() {
		if(_muts == NULL) {
			// No mutations to undo
			return;
		}
		for(size_t i = 0; i < length(*_muts); i++) {
			const QueryMutation& m = (*_muts)[i];
			assert_lt(m.pos, _qlen);
			assert_leq(m.oldBase, 4);
			assert_lt(m.newBase, 4);
			assert_neq(m.oldBase, m.newBase);
			assert_eq((uint32_t)((*_qry)[m.pos]), (uint32_t)m.newBase);
			(*_qry)[m.pos] = (Dna5)(int)m.oldBase; // undo it
		}
	}

	/**
	 * Report a range of alignments with # mismatches = stackDepth and
	 * with the mutations (also mismatches) contained in _muts.  The
	 * range is delimited by top and bot.  Returns true iff one or more
	 * full alignments were successfully reported and the caller can
	 * stop searching.
	 */
	bool reportAlignment(uint32_t stackDepth, TIndexOffU top,
			TIndexOffU bot, uint16_t cost)
	{
#ifndef NDEBUG
		// No two elements of _mms[] should be the same
		assert_geq(_mms.size(), stackDepth);
		for(size_t i = 0; i < stackDepth; i++) {
			for(size_t j = i+1; j < stackDepth; j++) {
				assert_neq(_mms[j], _mms[i]);
			}
			// All elements of _mms[] should fall within bounds
			assert_lt(_mms[i], _qlen);
		}
#endif
		if(_reportPartials) {
			assert_leq(stackDepth, _reportPartials);
			if(stackDepth > 0) {
				// Report this partial alignment.  A partial alignment
				// is defined purely by its mismatches; top and bot are
				// ignored.
				reportPartial(stackDepth);
			}
			return false; // keep going - we want to find all partial alignments
		}
		int stratum = 0;
		if(stackDepth > 0) {
			stratum = calcStratum(_mms, stackDepth);
		}
		assert_lt(stratum, 4);
		assert_geq(stratum, 0);
		bool hit;
		// If _muts != NULL then this alignment extends a partial
		// alignment, so we have to account for the differences present
		// in the partial.
		if(_muts != NULL) {
			// Undo partial-alignment mutations to get original _qry
			ASSERT_ONLY(String<Dna5> tmp = (*_qry));
			undoPartialMutations();
			assert_neq(tmp, (*_qry));
			// Add the partial-alignment mutations to the _mms[] array
			promotePartialMutations(stackDepth);
			// All muts are in the seed, so they count toward the stratum
			size_t numMuts = length(*_muts);
			stratum += numMuts;
			cost |= (stratum << 14);
			assert_geq(cost, (uint32_t)(stratum << 14));
			// Report the range of full alignments
			hit = reportFullAlignment((uint32_t)(stackDepth + numMuts), top, bot, stratum, cost);
			// Re-apply partial-alignment mutations
			applyPartialMutations();
			assert_eq(tmp, (*_qry));
		} else {
			// Report the range of full alignments
			cost |= (stratum << 14);
			assert_geq(cost, (uint32_t)(stratum << 14));
			hit = reportFullAlignment(stackDepth, top, bot, stratum, cost);
		}
		return hit;
	}

	/**
	 * Report a range of full alignments with # mismatches = stackDepth.
	 * The range is delimited by top and bot.  Returns true if one or
	 * more alignments were successfully reported.  Returns true iff
	 * one or more full alignments were successfully reported and the
	 * caller can stop searching.
	 */
	bool reportFullAlignment(uint32_t stackDepth,
			TIndexOffU top,
			TIndexOffU bot,
	                         int stratum,
	                         uint16_t cost)
	{
		assert_gt(bot, top);
		if(stackDepth == 0 && !_reportExacts) {
			// We are not reporting exact hits (usually because we've
			// already reported them as part of a previous invocation
			// of the backtracker)
			return false;
		}
		assert(!_reportRanges);
		TIndexOffU spread = bot - top;
		// Pick a random spot in the range to begin report
		TIndexOffU r = top + (_rand.nextU<TIndexOffU>() % spread);
		for(TIndexOffU i = 0; i < spread; i++) {
			TIndexOffU ri = r + i;
			if(ri >= bot) ri -= spread;
			// reportChaseOne takes the _mms[] list in terms of
			// their indices into the query string; not in terms
			// of their offset from the 3' or 5' end.
			assert_geq(cost, (uint32_t)(stratum << 14));
			if(_ebwt->reportChaseOne((*_qry), _qual, _name,
			                         _color, _primer, _trimc, colorExEnds,
			                         snpPhred, _refs, _mms, _refcs,
			                         stackDepth, ri, top, bot,
			                         (uint32_t)_qlen, stratum, cost, _patid,
			                         _seed, _params))
			{
				// Return value of true means that we can stop
				return true;
			}
			// Return value of false means that we should continue
			// searching.  This could happen if we the call to
			// reportChaseOne() reported a hit, but the user asked for
			// multiple hits and we haven't reached the ceiling yet.
			// This might also happen if the call to reportChaseOne()
			// didn't report a hit because the alignment was spurious
			// (i.e. overlapped some padding).
		}
		// All range elements were examined and we should keep going
		return false;
	}

	/**
	 * Report the partial alignment represented by the current stack
	 * state (_mms[] and stackDepth).
	 */
	bool reportPartial(uint32_t stackDepth) {
		// Sanity-check stack depth
		if(_3revOff != _2revOff) {
			assert_leq(stackDepth, 3);
		} else if(_2revOff != _1revOff) {
			assert_leq(stackDepth, 2);
		} else {
			assert_leq(stackDepth, 1);
		}

		// Possibly report
		assert_gt(_reportPartials, 0);
		assert(_partials != NULL);
		ASSERT_ONLY(uint32_t qualTot = 0);
		PartialAlignment al;
		al.u64.u64 = 0xffffffffffffffffllu;
		assert_leq(stackDepth, 3);
		assert_gt(stackDepth, 0);

		// First mismatch
		assert_gt(_mms.size(), 0);
		assert_lt(_mms[0], _qlen);
		// First, append the mismatch position in the read
		al.entry.pos0 = (uint16_t)_mms[0]; // pos
		ASSERT_ONLY(uint8_t qual0 = mmPenalty(_maqPenalty, phredCharToPhredQual((*_qual)[_mms[0]])));
		ASSERT_ONLY(qualTot += qual0);
		uint32_t ci = (uint32_t)(_qlen - _mms[0] - 1);
		// _chars[] is index in terms of RHS-relative depth
		int c = (int)(Dna5)_chars[ci];
		assert_lt(c, 4);
		assert_neq(c, (int)(*_qry)[_mms[0]]);
		// Second, append the substituted character for the position
		al.entry.char0 = c;

		if(stackDepth > 1) {
			assert_gt(_mms.size(), 1);
			// Second mismatch
			assert_lt(_mms[1], _qlen);
			// First, append the mismatch position in the read
			al.entry.pos1 = (uint16_t)_mms[1]; // pos
			ASSERT_ONLY(uint8_t qual1 = mmPenalty(_maqPenalty, phredCharToPhredQual((*_qual)[_mms[1]])));
			ASSERT_ONLY(qualTot += qual1);
			ci = (uint32_t)(_qlen - _mms[1] - 1);
			// _chars[] is index in terms of RHS-relative depth
			c = (int)(Dna5)_chars[ci];
			assert_lt(c, 4);
			assert_neq(c, (int)(*_qry)[_mms[1]]);
			// Second, append the substituted character for the position
			al.entry.char1 = c;
			if(stackDepth > 2) {
				assert_gt(_mms.size(), 2);
				// Second mismatch
				assert_lt(_mms[2], _qlen);
				// First, append the mismatch position in the read
				al.entry.pos2 = (uint16_t)_mms[2]; // pos
				ASSERT_ONLY(uint8_t qual2 = mmPenalty(_maqPenalty, phredCharToPhredQual((*_qual)[_mms[2]])));
				ASSERT_ONLY(qualTot += qual2);
				ci = (uint32_t)(_qlen - _mms[2] - 1);
				// _chars[] is index in terms of RHS-relative depth
				c = (int)(Dna5)_chars[ci];
				assert_lt(c, 4);
				assert_neq(c, (int)(*_qry)[_mms[2]]);
				// Second, append the substituted character for the position
				al.entry.char2 = c;
			} else {
				// Signal that the '2' slot is empty
				al.entry.pos2 = 0xffff;
			}
		} else {
			// Signal that the '1' slot is empty
			al.entry.pos1 = 0xffff;
		}

		assert_leq(qualTot, _qualThresh);
		assert(validPartialAlignment(al));
#ifndef NDEBUG
		assert(al.repOk(_qualThresh, (uint32_t)_qlen, (*_qual), _maqPenalty));
		for(size_t i = 0; i < _partialsBuf.size(); i++) {
			assert(validPartialAlignment(_partialsBuf[i]));
			assert(!samePartialAlignment(_partialsBuf[i], al));
		}
#endif
		_partialsBuf.push_back(al);
		return true;
	}

	/**
	 * Check that the given eligibility parameters (lowAltQual,
	 * eligibleSz, eligibleNum) are correct, given the appropriate
	 * inputs (pairs, elims, depth, d, unrevOff)
	 */
	bool sanityCheckEligibility(uint32_t  depth,
	                            uint32_t  d,
	                            uint32_t  unrevOff,
	                            uint32_t  lowAltQual,
	                            uint32_t  eligibleSz,
	                            uint32_t  eligibleNum,
	                            TIndexOffU* pairs,
	                            uint8_t*  elims)
	{
		// Sanity check that the lay of the land is as we
		// expect given eligibleNum and eligibleSz
		size_t i = max(depth, unrevOff), j = 0;
		uint32_t cumSz = 0;
		uint32_t eligiblesVisited = 0;
		for(; i <= d; i++) {
			uint32_t icur = (uint32_t)(_qlen - i - 1); // current offset into _qry
			uint8_t qi = qualAt(icur);
			assert_lt(elims[i], 16);
			if((qi == lowAltQual || !_considerQuals) && elims[i] != 15) {
				// This is an eligible position with at least
				// one remaining backtrack target
				for(j = 0; j < 4; j++) {
					if((elims[i] & (1 << j)) == 0) {
						// This pair has not been eliminated
						assert_gt(pairBot(pairs, i, j), pairTop(pairs, i, j));
						cumSz += pairSpread(pairs, i, j);
						eligiblesVisited++;
					}
				}
			}
		}
		assert_eq(cumSz, eligibleSz);
		assert_eq(eligiblesVisited, eligibleNum);
		return true;
	}

	const BitPairReference* _refs; // reference sequences (or NULL if not colorspace)
	String<Dna5>*       _qry;    // query (read) sequence
	size_t              _qlen;   // length of _qry
	String<char>*       _qual;   // quality values for _qry
	String<char>*       _name;   // name of _qry
	bool                _color;  // whether read is colorspace
	const Ebwt<String<Dna> >* _ebwt;   // Ebwt to search in
	const EbwtSearchParams<String<Dna> >& _params;   // Ebwt to search in
	uint32_t            _unrevOff; // unrevisitable chunk
	uint32_t            _1revOff;  // 1-revisitable chunk
	uint32_t            _2revOff;  // 2-revisitable chunk
	uint32_t            _3revOff;  // 3-revisitable chunk
	/// Whether to round qualities off Maq-style when calculating penalties
	bool                _maqPenalty;
	uint32_t            _qualThresh; // only accept hits with weighted
	                             // hamming distance <= _qualThresh
	TIndexOffU           *_pairs;  // ranges, leveled in parallel
	                             // with decision stack
	uint8_t            *_elims;  // which ranges have been
	                             // eliminated, leveled in parallel
	                             // with decision stack
	std::vector<TIndexOffU> _mms;  // array for holding mismatches
	std::vector<uint8_t> _refcs;  // array for holding mismatches
	// Entries in _mms[] are in terms of offset into
	// _qry - not in terms of offset from 3' or 5' end
	char               *_chars;  // characters selected so far
	// If > 0, report partial alignments up to this many mismatches
	uint32_t            _reportPartials;
	/// Do not report alignments with stratum < this limit
	bool                _reportExacts;
	/// When reporting a full alignment, report top/bot; don't chase
	/// any of the results
	bool                _reportRanges;
	/// Append partial alignments here
	PartialAlignmentManager *_partials;
	/// Set of mutations that apply for a partial alignment
	String<QueryMutation> *_muts;
	/// Reference texts (NULL if they are unavailable
	vector<String<Dna5> >* _os;
	/// Whether to use the _os array together with a naive matching
	/// algorithm to double-check reported alignments (or the lack
	/// thereof)
	bool                _sanity;
	/// Whether to consider quality values when deciding where to
	/// backtrack
	bool                _considerQuals;
	bool                _halfAndHalf;
	/// Depth of 5'-seed-half border
	uint32_t            _5depth;
	/// Depth of 3'-seed-half border
	uint32_t            _3depth;
	/// Default quals
	String<char>        _qualDefault;
	/// Number of backtracks in last call to backtrack()
	uint32_t            _numBts;
	/// Number of backtracks since last reset
	uint32_t            _totNumBts;
	/// Max # of backtracks to allow before giving up
	uint32_t            _maxBts;
	/// Whether we precalcualted the Ebwt locus information for the
	/// next top/bot pair
	bool    _precalcedSideLocus;
	/// Precalculated top locus
	SideLocus           _preLtop;
	/// Precalculated bot locus
	SideLocus           _preLbot;
	/// Flag to record whether a 'false' return from backtracker is due
	/// to having exceeded one or more backrtacking limits
	bool                _bailedOnBacktracks;
	/// Source of pseudo-random numbers
	RandomSource        _rand;
	/// Be talkative
	bool                _verbose;
	uint64_t            _ihits;
	// Holding area for partial alignments
	vector<PartialAlignment> _partialsBuf;
	// Current range to expose to consumers
	Range               _curRange;
	uint32_t            _patid;
	char                _primer;
	char                _trimc;
	uint32_t            _seed;
#ifndef NDEBUG
	std::set<TIndexOff> allTops_;
#endif
};

/**
 * Class that coordinates quality- and quantity-aware backtracking over
 * some range of a read sequence.
 *
 * The creator can configure the BacktrackManager to treat different
 * stretches of the read differently.
 */
class EbwtRangeSource : public RangeSource {
	typedef Ebwt<String<Dna> > TEbwt;
	typedef std::pair<int, int> TIntPair;
public:
	EbwtRangeSource(
			const TEbwt* ebwt,
			bool         fw,
			TIndexOffU   qualLim,
			bool         reportExacts,
			bool         verbose,
			bool         quiet,
			int          halfAndHalf,
			bool         partial,
			bool         maqPenalty,
			bool         qualOrder,
			AlignerMetrics *metrics = NULL) :
		RangeSource(),
		qry_(NULL),
		qlen_(0),
		qual_(NULL),
		name_(NULL),
		altQry_(NULL),
		altQual_(NULL),
		alts_(0),
		fuzzy_(false),
		ebwt_(ebwt),
		fw_(fw),
		offRev0_(0),
		offRev1_(0),
		offRev2_(0),
		offRev3_(0),
		maqPenalty_(maqPenalty),
		qualOrder_(qualOrder),
		qualLim_(qualLim),
		reportExacts_(reportExacts),
		halfAndHalf_(halfAndHalf),
		partial_(partial),
		depth5_(0),
		depth3_(0),
		verbose_(verbose),
		quiet_(quiet),
		skippingThisRead_(false),
		metrics_(metrics)
	{ curEbwt_ = ebwt_; }

	/**
	 * Set a new query read.
	 */
	virtual void setQuery(ReadBuf& r, Range *seedRange) {
		const bool ebwtFw = ebwt_->fw();
		if(ebwtFw) {
			qry_  = fw_ ? &r.patFw : &r.patRc;
			qual_ = fw_ ? &r.qual  : &r.qualRev;
			altQry_  = (String<Dna5>*)(fw_ ? r.altPatFw : r.altPatRc);
			altQual_ = (String<char>*)(fw_ ? r.altQual  : r.altQualRev);
		} else {
			qry_  = fw_ ? &r.patFwRev : &r.patRcRev;
			qual_ = fw_ ? &r.qualRev  : &r.qual;
			altQry_  = (String<Dna5>*)(fw_ ? r.altPatFwRev : r.altPatRcRev);
			altQual_ = (String<char>*)(fw_ ? r.altQualRev  : r.altQual);
		}
		alts_ = r.alts;
		name_ = &r.name;
		fuzzy_ = r.fuzzy;
		if(seedRange != NULL) seedRange_ = *seedRange;
		else                  seedRange_.invalidate();
		qlen_ = length(*qry_);
		skippingThisRead_ = false;
		// Apply edits from the partial alignment to the query pattern
		if(seedRange_.valid()) {
			qryBuf_ = *qry_;
			const size_t srSz = seedRange_.mms.size();
			assert_gt(srSz, 0);
			assert_eq(srSz, seedRange_.refcs.size());
			for(size_t i = 0; i < srSz; i++) {
				assert_lt(seedRange_.mms[i], qlen_);
				char rc = (char)seedRange_.refcs[i];
				assert(rc == 'A' || rc == 'C' || rc == 'G' || rc == 'T');
				ASSERT_ONLY(char oc = (char)qryBuf_[qlen_ - seedRange_.mms[i] - 1]);
				assert_neq(rc, oc);
				qryBuf_[qlen_ - seedRange_.mms[i] - 1] = (Dna5)rc;
				assert_neq((Dna5)rc, (*qry_)[qlen_ - seedRange_.mms[i] - 1]);
			}
			qry_ = &qryBuf_;
		}
		// Make sure every qual is a valid qual ASCII character (>= 33)
		for(size_t i = 0; i < length(*qual_); i++) {
			assert_geq((*qual_)[i], 33);
			for(int j = 0; j < alts_; j++) {
				assert_geq(altQual_[j][i], 33);
			}
		}
		assert_geq(length(*qual_), qlen_);
		this->done = false;
		this->foundRange = false;
		color_ = r.color;
		rand_.init(r.seed);
	}

	/**
	 * Set backtracking constraints.
	 */
	void setOffs(uint32_t depth5,   // depth of far edge of hi-half
	             uint32_t depth3,   // depth of far edge of lo-half
	             uint32_t unrevOff, // depth above which we cannot backtrack
	             uint32_t revOff1,  // depth above which we may backtrack just once
	             uint32_t revOff2,  // depth above which we may backtrack just twice
	             uint32_t revOff3)  // depth above which we may backtrack just three times
	{
		depth5_   = depth5;
		depth3_   = depth3;
		assert_geq(depth3_, depth5_);
		offRev0_ = unrevOff;
		offRev1_  = revOff1;
		offRev2_  = revOff2;
		offRev3_  = revOff3;
	}

	/**
	 * Return true iff this RangeSource is allowed to report exact
	 * alignments (exact = no edits).
	 */
	bool reportExacts() const {
		return reportExacts_;
	}

	/// Return the current range
	virtual Range& range() {
		return curRange_;
	}

	/**
	 * Set qlen_ according to parameter, except don't let it fall below
	 * the length of the query.
	 */
	void setQlen(uint32_t qlen) {
		assert(qry_ != NULL);
		qlen_ = min<uint32_t>((uint32_t)length(*qry_), qlen);
	}

	/**
	 * Initiate continuations so that the next call to advance() begins
	 * a new search.  Note that contMan is empty upon return if there
	 * are no valid continuations to begin with.  Also note that
	 * calling initConts() may result in finding a range (i.e., if we
	 * immediately jump to a valid range using the ftab).
	 */
	virtual void
	initBranch(PathManager& pm) {
		assert(curEbwt_ != NULL);
		assert_gt(length(*qry_), 0);
		assert_leq(qlen_, length(*qry_));
		assert_geq(length(*qual_), length(*qry_));
		const Ebwt<String<Dna> >& ebwt = *ebwt_;
		int ftabChars = ebwt._eh._ftabChars;
		this->foundRange = false;
		int nsInSeed = 0; int nsInFtab = 0;
		ASSERT_ONLY(allTops_.clear());
		if(skippingThisRead_) {
			this->done = true;
			return;
		}
		if(qlen_ < 4) {
			uint32_t maxmms = 0;
			if(offRev0_ != offRev1_) maxmms = 1;
			if(offRev1_ != offRev2_) maxmms = 2;
			if(offRev2_ != offRev3_) maxmms = 3;
			if(qlen_ <= maxmms) {
				if(!quiet_) {
					ThreadSafe _ts(&gLock);
					cerr << "Warning: Read (" << (*name_) << ") is less than " << (maxmms+1) << " characters long; skipping..." << endl;
				}
				this->done = true;
				skippingThisRead_ = true;
				return;
			}
		}
		if(!tallyNs(nsInSeed, nsInFtab)) {
			// No alignments are possible because of the distribution
			// of Ns in the read in combination with the backtracking
			// constraints.
			return;
		}
		// icost = total cost penalty (major bits = stratum, minor bits =
		// quality penalty) incurred so far by partial alignment
		uint16_t icost = (seedRange_.valid()) ? seedRange_.cost : 0;
		// iham = total quality penalty incurred so far by partial alignment
		uint16_t iham = (seedRange_.valid() && qualOrder_) ? (seedRange_.cost & ~0xc000): 0;
		assert_leq(iham, qualLim_);
		// m = depth beyond which ftab must not extend or else we might
		// miss some legitimate paths
		uint32_t m = min<uint32_t>(offRev0_, (uint32_t)qlen_);
		// Let skipInvalidExact = true if using the ftab would be a
		// waste because it would jump directly to an alignment we
		// couldn't use.
		bool ftabSkipsToEnd = (qlen_ == (uint32_t)ftabChars);
		bool skipInvalidExact = (!reportExacts_ && ftabSkipsToEnd);

		// If it's OK to use the ftab...
		if(nsInFtab == 0 && m >= (uint32_t)ftabChars && !skipInvalidExact) {
			// Use the ftab to jump 'ftabChars' chars into the read
			// from the right
			uint32_t ftabOff = calcFtabOff();
			TIndexOffU top = ebwt.ftabHi(ftabOff);
			TIndexOffU bot = ebwt.ftabLo(ftabOff+1);
			if(qlen_ == (uint32_t)ftabChars && bot > top) {
				// We found a range with 0 mismatches immediately.  Set
				// fields to indicate we found a range.
				assert(reportExacts_);
				curRange_.top     = top;
				curRange_.bot     = bot;
				curRange_.stratum = (icost >> 14);
				curRange_.cost    = icost;
				curRange_.numMms  = 0;
				curRange_.ebwt    = ebwt_;
				curRange_.fw      = fw_;
				curRange_.mms.clear(); // no mismatches
				curRange_.refcs.clear(); // no mismatches
				// Lump in the edits from the partial alignment
				addPartialEdits();
				assert(curRange_.repOk());
				// no need to do anything with curRange_.refcs
				this->foundRange  = true;
				//this->done = true;
				return;
			} else if (bot > top) {
				// We have a range to extend
				assert_leq(top, ebwt._eh._len);
				assert_leq(bot, ebwt._eh._len);
				Branch *b = pm.bpool.alloc();
				if(b == NULL) {
					assert(pm.empty());
					return;
				}
				if(!b->init(
				        pm.rpool, pm.epool, pm.bpool.lastId(), (uint32_t)qlen_,
				        offRev0_, offRev1_, offRev2_, offRev3_,
				        0, ftabChars, icost, iham, top, bot,
				        ebwt._eh, ebwt._ebwt))
				{
					// Negative result from b->init() indicates we ran
					// out of best-first chunk memory
					assert(pm.empty());
					return;
				}
				assert(!b->curtailed_);
				assert(!b->exhausted_);
				assert_gt(b->depth3_, 0);
				pm.push(b); // insert into priority queue
				assert(!pm.empty());
			} else {
				// The arrows are already closed within the
				// unrevisitable region; give up
			}
		} else {
			// We can't use the ftab, so we start from the rightmost
			// position and use _fchr
			Branch *b = pm.bpool.alloc();
			if(b == NULL) {
				assert(pm.empty());
				return;
			}
			if(!b->init(pm.rpool, pm.epool, pm.bpool.lastId(), (uint32_t)qlen_,
			        offRev0_, offRev1_, offRev2_, offRev3_,
			        0, 0, icost, iham, 0, 0, ebwt._eh, ebwt._ebwt))
			{
				// Negative result from b->init() indicates we ran
				// out of best-first chunk memory
				assert(pm.empty());
				return;
			}
			assert(!b->curtailed_);
			assert(!b->exhausted_);
			assert_gt(b->depth3_, 0);
			pm.push(b); // insert into priority queue
			assert(!pm.empty());
		}
		return;
	}

	/**
	 * Advance along the lowest-cost branch managed by the given
	 * PathManager.  Keep advancing until condition 'until' is
	 * satisfied.  Typically, the stopping condition 'until' is
	 * set to stop whenever pm's minCost changes.
	 */
	virtual void
	advanceBranch(int until, uint16_t minCost, PathManager& pm) {
		assert(curEbwt_ != NULL);

		// Let this->foundRange = false; we'll set it to true iff this call
		// to advance yielded a new valid-alignment range.
		this->foundRange = false;

		// Can't have already exceeded weighted hamming distance threshold
		assert_gt(length(*qry_), 0);
		assert_leq(qlen_, length(*qry_));
		assert_geq(length(*qual_), length(*qry_));
		assert(!pm.empty());

		do {
			assert(pm.repOk());
			// Get the highest-priority branch according to the priority
			// queue in 'pm'
			Branch* br = pm.front();
			// Shouldn't be curtailed or exhausted
			assert(!br->exhausted_);
			assert(!br->curtailed_);
			assert_gt(br->depth3_, 0);
			assert_leq(br->ham_, qualLim_);
			if(verbose_) {
				br->print((*qry_), (*qual_), minCost, cout, (halfAndHalf_>0), partial_, fw_, ebwt_->fw());
				if(!br->edits_.empty()) {
					cout << "Edit: ";
					for(size_t i = 0; i < br->edits_.size(); i++) {
						Edit e = br->edits_.get(i);
						cout << (curEbwt_->fw() ? (qlen_ - e.pos - 1) : e.pos)
							 << (char)e.chr;
						if(i < br->edits_.size()-1) cout << " ";
					}
					cout << endl;
				}

			}
			assert(br->repOk((uint32_t)qlen_));

			ASSERT_ONLY(int stratum = br->cost_ >> 14); // shift the stratum over
			assert_lt(stratum, 4);
			// Not necessarily true with rounding
			uint32_t depth = br->tipDepth();

			const Ebwt<String<Dna> >& ebwt = *ebwt_;

			if(halfAndHalf_ > 0) assert_gt(depth3_, depth5_);

			bool reportedPartial = false;
			bool invalidExact = false;
			bool empty = false;
			bool hit = false;
			uint16_t cost = br->cost_;
			uint32_t cur = 0;
			uint32_t nedits = 0;

			if(halfAndHalf_ && !hhCheckTop(br, depth, 0)) {
				// Stop extending this branch because it violates a half-
				// and-half constraint
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, (uint32_t)qlen_, depth3_, qualOrder_);
				goto bail;
			}

			cur = (uint32_t)(qlen_ - depth - 1); // current offset into qry_
			if(depth < qlen_) {
				// Determine whether ranges at this location are candidates
				// for backtracking
				int c = (int)(*qry_)[cur]; // get char at this position
				int nextc = -1;
				if(cur < qlen_-1) nextc = (int)(*qry_)[cur+1];
				assert_leq(c, 4);
				// If any uncalled base's penalty is still under
				// the ceiling, then this position is an alternative
				uint8_t q[4] = {'!', '!', '!', '!'};
				uint8_t bestq;
				// get unrounded penalties at this position
				if(fuzzy_) {
					bestq = penaltiesAt(cur, q, alts_, *qual_, altQry_, altQual_);
				} else {
					bestq = q[0] = q[1] = q[2] = q[3] =
						mmPenalty(maqPenalty_, qualAt(cur));
				}

				// The current query position is a legit alternative if it a) is
				// not in the unrevisitable region, and b) its selection would
				// not necessarily cause the quality ceiling (if one exists) to
				// be exceeded
				bool curIsAlternative = (depth >= br->depth0_) &&
				                        (br->ham_ + bestq <= qualLim_);
				ASSERT_ONLY(TIndexOffU obot = br->bot_);
				TIndexOffU otop = br->top_;

				// If c is 'N', then it's a mismatch
				if(c == 4 && depth > 0) {
					// Force the 'else if(curIsAlternative)' or 'else'
					// branches below
					br->top_ = br->bot_ = 1;
				} else if(c == 4) {
					// We'll take the 'if(br->top == 0 && br->bot == 0)'
					// branch below
					assert_eq(0, br->top_);
					assert_eq(0, br->bot_);
				}

				// Get the range state for the current position
				RangeState *rs = br->rangeState();
				assert(rs != NULL);
				// Calculate the ranges for this position
				if(br->top_ == 0 && br->bot_ == 0) {
					// Calculate first quartet of ranges using the _fchr[]
					// array
								  rs->tops[0] = ebwt._fchr[0];
					rs->bots[0] = rs->tops[1] = ebwt._fchr[1];
					rs->bots[1] = rs->tops[2] = ebwt._fchr[2];
					rs->bots[2] = rs->tops[3] = ebwt._fchr[3];
					rs->bots[3]               = ebwt._fchr[4];
					ASSERT_ONLY(int r =)
					br->installRanges(c, nextc, fuzzy_, qualLim_ - br->ham_, q);
					assert(r < 4 || c == 4);
					// Update top and bot
					if(c < 4) {
						br->top_ = rs->tops[c];
						br->bot_ = rs->bots[c];
					}
				} else if(curIsAlternative && (br->bot_ > br->top_ || c == 4)) {
					// Calculate next quartet of ranges.  We hope that the
					// appropriate cache lines are prefetched.
					assert(br->ltop_.valid());
								  rs->tops[0] =
					rs->bots[0] = rs->tops[1] =
					rs->bots[1] = rs->tops[2] =
					rs->bots[2] = rs->tops[3] =
					rs->bots[3]               = 0;
					if(br->lbot_.valid()) {
						if(metrics_ != NULL) metrics_->curBwtOps_++;
						ebwt.mapLFEx(br->ltop_, br->lbot_, (TIndexOffU*)rs->tops, (TIndexOffU*)rs->bots);
					} else {
#ifndef NDEBUG
						TIndexOffU tmptops[] = {0, 0, 0, 0};
						TIndexOffU tmpbots[] = {0, 0, 0, 0};
						SideLocus ltop, lbot;
						ltop.initFromRow(otop, ebwt_->_eh, ebwt_->_ebwt);
						lbot.initFromRow(obot, ebwt_->_eh, ebwt_->_ebwt);
						ebwt.mapLFEx(ltop, lbot, tmptops, tmpbots);
#endif
						if(metrics_ != NULL) metrics_->curBwtOps_++;
						int cc = ebwt.mapLF1((TIndexOffU&)otop, br->ltop_);
						br->top_ = otop;
						assert(cc == -1 || (cc >= 0 && cc < 4));
						if(cc >= 0) {
							assert_lt(cc, 4);
							rs->tops[cc] = br->top_;
							rs->bots[cc] = (br->top_ + 1);
						}
#ifndef NDEBUG
						for(int i = 0; i < 4; i++) {
							assert_eq(tmpbots[i] - tmptops[i],
							          rs->bots[i] - rs->tops[i]);
						}
#endif
					}
					ASSERT_ONLY(int r =)
					br->installRanges(c, nextc, fuzzy_, qualLim_ - br->ham_, q);
					assert(r < 4 || c == 4);
					// Update top and bot
					if(c < 4) {
						br->top_ = rs->tops[c];
						br->bot_ = rs->bots[c];
					} else {
						br->top_ = br->bot_ = 1;
					}
				} else if(br->bot_ > br->top_) {
					// This read position is not a legitimate backtracking
					// alternative.  No need to do the bookkeeping for the
					// entire quartet, just do c.  We hope that the
					// appropriate cache lines are prefetched before now;
					// otherwise, we're about to take an expensive cache
					// miss.
					assert(br->ltop_.valid());
					rs->eliminated_ = true; // eliminate all alternatives leaving this node
					assert(br->eliminated(br->len_));
					if(c < 4) {
						if(br->top_ + 1 == br->bot_) {
							if(metrics_ != NULL) metrics_->curBwtOps_++;
							br->bot_ = br->top_ = ebwt.mapLF1(br->top_, br->ltop_, c);
							if(br->bot_ != OFF_MASK) br->bot_++;
						} else {
							if(metrics_ != NULL) metrics_->curBwtOps_++;
							br->top_ = ebwt.mapLF(br->ltop_, c);
							assert(br->lbot_.valid());
							if(metrics_ != NULL) metrics_->curBwtOps_++;
							br->bot_ = ebwt.mapLF(br->lbot_, c);
						}
					}
				} else {
					rs->eliminated_ = true;
				}
				assert(rs->repOk());
				// br->top_ and br->bot_ now contain the next top and bot
			} else {
				// The continuation had already processed the whole read
				assert_eq(qlen_, depth);
				cur = 0;
			}
			empty = (br->top_ == br->bot_);
			hit = (cur == 0 && !empty);

			// Check whether we've obtained an exact alignment when
			// we've been instructed not to report exact alignments
			nedits = (uint32_t)br->edits_.size();
			invalidExact = (hit && nedits == 0 && !reportExacts_);
			assert_leq(br->ham_, qualLim_);

			// Set this to true if the only way to make legal progress
			// is via one or more additional backtracks.
			if(halfAndHalf_ && !hhCheck(br, depth, empty)) {
				// This alignment doesn't satisfy the half-and-half
				// requirements; reject it
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, (uint32_t)qlen_, depth3_, qualOrder_);
				goto bail;
			}

			if(hit &&            // there is a range to report
			   !invalidExact &&  // not disqualified by no-exact-hits setting
			   !reportedPartial) // not an already-reported partial alignment
			{
				if(verbose_) {
					if(partial_) {
						cout << " Partial alignment:" << endl;
					} else {
						cout << " Final alignment:" << endl;
					}
					br->len_++;
					br->print((*qry_), (*qual_), minCost, cout, halfAndHalf_ > 0, partial_, fw_, ebwt_->fw());
					br->len_--;
					cout << endl;
				}
				assert_gt(br->bot_, br->top_);
				assert_leq(br->ham_, qualLim_);
				assert_leq((uint32_t)(br->cost_ & ~0xc000), qualLim_);
				if(metrics_ != NULL) metrics_->setReadHasRange();
				curRange_.top     = br->top_;
				curRange_.bot     = br->bot_;
				curRange_.cost    = br->cost_;
				curRange_.stratum = (br->cost_ >> 14);
				curRange_.numMms  = nedits;
				curRange_.fw      = fw_;
				curRange_.mms.clear();
				curRange_.refcs.clear();
				for(size_t i = 0; i < nedits; i++) {
					curRange_.mms.push_back((uint32_t)(qlen_ - br->edits_.get(i).pos - 1));
					curRange_.refcs.push_back((char)br->edits_.get(i).chr);
				}
				addPartialEdits();
				curRange_.ebwt    = ebwt_;
				this->foundRange  = true;
	#ifndef NDEBUG
				TIndexOff top2 = (TIndexOff)br->top_;
				top2++; // ensure it's not 0
				if(ebwt_->fw()) top2 = -top2;
				assert(allTops_.find(top2) == allTops_.end());
				allTops_.insert(top2);
	#endif
				assert(curRange_.repOk());
				// Must curtail because we've consumed the whole pattern
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, (uint32_t)qlen_, depth3_, qualOrder_);
			} else if(empty || cur == 0) {
				// The branch couldn't be extended further
				if(metrics_ != NULL) metrics_->curBacktracks_++;
				pm.curtail(br, (uint32_t)qlen_, depth3_, qualOrder_);
			} else {
				// Extend the branch by one position; no change to its cost
				// so there's no need to reconsider where it lies in the
				// priority queue
				assert_neq(0, cur);
				br->extend();
			}
		bail:
			// Make sure the front element of the priority queue is
			// extendable (i.e. not curtailed) and then prep it.
			if(!pm.splitAndPrep(rand_, (uint32_t)qlen_, qualLim_, depth3_,
			                    qualOrder_, fuzzy_,
			                    ebwt_->_eh, ebwt_->_ebwt, ebwt_->_fw))
			{
				pm.reset(0);
				assert(pm.empty());
			}
			if(pm.empty()) {
				// No more branches
				break;
			}
			assert(!pm.front()->curtailed_);
			assert(!pm.front()->exhausted_);

			if(until == ADV_COST_CHANGES && pm.front()->cost_ != cost) break;
			else if(until == ADV_STEP) break;

		} while(!this->foundRange);
		if(!pm.empty()) {
			assert(!pm.front()->curtailed_);
			assert(!pm.front()->exhausted_);
		}
	}

	/**
	 * Return true iff we're enforcing a half-and-half constraint
	 * (forced edits in both seed halves).
	 */
	int halfAndHalf() const {
		return halfAndHalf_;
	}

protected:

	/**
	 * Lump all the seed-alignment edits from the seedRange_ range
	 * found previously to the curRange_ range just found.
	 */
	void addPartialEdits() {
		// Lump in the edits from the partial alignment
		if(seedRange_.valid()) {
			const size_t srSz = seedRange_.mms.size();
			for(size_t i = 0; i < srSz; i++) {
				curRange_.mms.push_back((uint32_t)(qlen_ - seedRange_.mms[i] - 1));
				curRange_.refcs.push_back(seedRange_.refcs[i]);
			}
			curRange_.numMms += srSz;
		}
	}

	/**
	 * Return true iff we're OK to continue after considering which
	 * half-seed boundary we're passing through, together with the
	 * number of mismatches accumulated so far.  Return false if we
	 * should stop because a half-and-half constraint is violated.  If
	 * we're not currently passing a half-seed boundary, just return
	 * true.
	 */
	bool hhCheck(Branch *b, uint32_t depth, bool empty) {
		const uint32_t nedits = (uint32_t)b->edits_.size();
		ASSERT_ONLY(uint32_t lim3 = (offRev3_ == offRev2_)? 2 : 3);
		ASSERT_ONLY(uint32_t lim5 = (offRev1_ == offRev0_)? 2 : 1);
		if((depth == (depth5_-1)) && !empty) {
			// We're crossing the boundary separating the hi-half
			// from the non-seed portion of the read.
			// We should induce a mismatch if we haven't mismatched
			// yet, so that we don't waste time pursuing a match
			// that was covered by a previous phase
			assert_leq(nedits, lim5);
			return nedits > 0;
		} else if((depth == (depth3_-1)) && !empty) {
			// We're crossing the boundary separating the lo-half
			// from the non-seed portion of the read
			assert_leq(nedits, lim3);
			assert_gt(nedits, 0);
			// Count the mismatches in the lo and hi halves
			uint32_t loHalfMms = 0, hiHalfMms = 0;
			for(size_t i = 0; i < nedits; i++) {
				uint32_t depth = b->edits_.get(i).pos;
				if     (depth < depth5_) hiHalfMms++;
				else if(depth < depth3_) loHalfMms++;
				else assert(false);
			}
			assert_leq(loHalfMms + hiHalfMms, lim3);
			bool invalidHalfAndHalf = (loHalfMms == 0 || hiHalfMms == 0);
			return (nedits >= (uint32_t)halfAndHalf_ && !invalidHalfAndHalf);
		}
#ifndef NDEBUG
		if(depth < depth5_-1) {
			assert_leq(nedits, lim5);
		}
		else if(depth >= depth5_ && depth < depth3_-1) {
			assert_gt(nedits, 0);
			assert_leq(nedits, lim3);
		}
#endif
		return true;
	}

	/**
	 * Return true iff the state of the backtracker as encoded by
	 * stackDepth, d and iham is compatible with the current half-and-
	 * half alignment mode.  prehits is for sanity checking when
	 * bailing.
	 */
	bool hhCheckTop(Branch* b,
	                uint32_t d,
	                uint32_t iham,
	                uint64_t prehits = 0xffffffffffffffffllu)
	{
		// Crossing from the hi-half into the lo-half
		ASSERT_ONLY(uint32_t lim3 = (offRev3_ == offRev2_)? 2 : 3);
		ASSERT_ONLY(uint32_t lim5 = (offRev1_ == offRev0_)? 2 : 1);
		const uint32_t nedits = (uint32_t)b->edits_.size();
		if(d == depth5_) {
			assert_leq(nedits, lim5);
			if(nedits == 0) {
				return false;
			}
		} else if(d == depth3_) {
			assert_leq(nedits, lim3);
			if(nedits < (uint32_t)halfAndHalf_) {
				return false;
			}
		}
#ifndef NDEBUG
		else {
			// We didn't just cross a boundary, so do an in-between check
			if(d >= depth5_) {
				assert_geq(nedits, 1);
			} else if(d >= depth3_) {
				assert_geq(nedits, lim3);
			}
		}
#endif
		return true;
	}

	/**
	 * Return the Phred-scale quality value at position 'off'
	 */
	inline uint8_t qualAt(size_t off) {
		return phredCharToPhredQual((*qual_)[off]);
	}

	/**
	 * Tally how many Ns occur in the seed region and in the ftab-
	 * jumpable region of the read.  Check whether the mismatches
	 * induced by the Ns already violates the current policy.  Return
	 * false if the policy is already violated, true otherwise.
	 */
	bool tallyNs(int& nsInSeed, int& nsInFtab) {
		const Ebwt<String<Dna> >& ebwt = *ebwt_;
		int ftabChars = ebwt._eh._ftabChars;
		// Count Ns in the seed region of the read and short-circuit if
		// the configuration of Ns guarantees that there will be no
		// valid alignments given the backtracking constraints.
		for(size_t i = 0; i < offRev3_; i++) {
			if((int)(*qry_)[qlen_-i-1] == 4) {
				nsInSeed++;
				if(nsInSeed == 1) {
					if(i < offRev0_) {
						return false; // Exceeded mm budget on Ns alone
					}
				} else if(nsInSeed == 2) {
					if(i < offRev1_) {
						return false; // Exceeded mm budget on Ns alone
					}
				} else if(nsInSeed == 3) {
					if(i < offRev2_) {
						return false; // Exceeded mm budget on Ns alone
					}
				} else {
					assert_gt(nsInSeed, 3);
					return false;     // Exceeded mm budget on Ns alone
				}
			}
		}
		// Calculate the number of Ns there are in the region that
		// would get jumped over if the ftab were used.
		for(size_t i = 0; i < (size_t)ftabChars && i < qlen_; i++) {
			if((int)(*qry_)[qlen_-i-1] == 4) nsInFtab++;
		}
		return true;
	}

	/**
	 * Calculate the offset into the ftab for the rightmost 'ftabChars'
	 * characters of the current query. Rightmost char gets least
	 * significant bit-pair.
	 */
	uint32_t calcFtabOff() {
		const Ebwt<String<Dna> >& ebwt = *ebwt_;
		int ftabChars = ebwt._eh._ftabChars;
		uint32_t ftabOff = (*qry_)[qlen_ - ftabChars];
		assert_lt(ftabOff, 4);
		assert_lt(ftabOff, ebwt._eh._ftabLen-1);
		for(int i = ftabChars - 1; i > 0; i--) {
			ftabOff <<= 2;
			assert_lt((uint32_t)(*qry_)[qlen_-i], 4);
			ftabOff |= (uint32_t)(*qry_)[qlen_-i];
			assert_lt(ftabOff, ebwt._eh._ftabLen-1);
		}
		assert_lt(ftabOff, ebwt._eh._ftabLen-1);
		return ftabOff;
	}

	String<Dna5>*       qry_;    // query (read) sequence
	String<Dna5>        qryBuf_; // for composing modified qry_ strings
	size_t              qlen_;   // length of _qry
	String<char>*       qual_;   // quality values for _qry
	String<char>*       name_;   // name of _qry
	bool                color_;  // true -> read is colorspace
	String<Dna5>*       altQry_; // alternate basecalls
	String<char>*       altQual_; // quality values for alternate basecalls
	int                 alts_;   // max # alternatives
	bool                fuzzy_;  // alternate scoring scheme?
	const Ebwt<String<Dna> >*   ebwt_;   // Ebwt to search in
	bool                fw_;
	uint32_t            offRev0_; // unrevisitable chunk
	uint32_t            offRev1_;  // 1-revisitable chunk
	uint32_t            offRev2_;  // 2-revisitable chunk
	uint32_t            offRev3_;  // 3-revisitable chunk
	/// Whether to round qualities off Maq-style when calculating penalties
	bool                maqPenalty_;
	/// Whether to order paths on our search in a way that takes
	/// qualities into account.  If this is false, the effect is that
	/// the first path reported is guaranteed to be in the best
	/// stratum, but it's not guaranteed to have the best quals.
	bool                qualOrder_;
	/// Reject alignments where sum of qualities at mismatched
	/// positions is greater than qualLim_
	uint32_t            qualLim_;
	/// Report exact alignments iff this is true
	bool                reportExacts_;
	/// Whether to use the _os array together with a naive matching
	/// algorithm to double-check reported alignments (or the lack
	/// thereof)
	int                 halfAndHalf_;
	/// Whether we're generating partial alignments for a longer
	/// alignment in the opposite index.
	bool                partial_;
	/// Depth of 5'-seed-half border
	uint32_t            depth5_;
	/// Depth of 3'-seed-half border
	uint32_t            depth3_;
	/// Source of pseudo-random numbers
	RandomSource        rand_;
	/// Be talkative
	bool                verbose_;
	/// Suppress unnecessary output
	bool                quiet_;
	// Current range to expose to consumers
	Range               curRange_;
	// Range for the partial alignment we're extending (NULL if we
	// aren't extending a partial)
	Range               seedRange_;
	// Starts as false; set to true as soon as we know we want to skip
	// all further processing of this read
	bool                skippingThisRead_;
	// Object encapsulating metrics
	AlignerMetrics*     metrics_;
#ifndef NDEBUG
	std::set<TIndexOff> allTops_;
#endif
};

/**
 * Concrete factory for EbwtRangeSource objects.
 */
class EbwtRangeSourceFactory {
	typedef Ebwt<String<Dna> > TEbwt;
public:
	EbwtRangeSourceFactory(
			const TEbwt* ebwt,
			bool         fw,
			uint32_t     qualThresh,
			bool         reportExacts,
			bool         verbose,
			bool         quiet,
			bool         halfAndHalf,
			bool         seeded,
			bool         maqPenalty,
			bool         qualOrder,
			AlignerMetrics *metrics = NULL) :
			ebwt_(ebwt),
			fw_(fw),
			qualThresh_(qualThresh),
			reportExacts_(reportExacts),
			verbose_(verbose),
			quiet_(quiet),
			halfAndHalf_(halfAndHalf),
			seeded_(seeded),
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			metrics_(metrics) { }

	/**
	 * Return new EbwtRangeSource with predefined params.s
	 */
	EbwtRangeSource *create() {
		return new EbwtRangeSource(ebwt_, fw_, qualThresh_,
		                           reportExacts_, verbose_, quiet_,
		                           halfAndHalf_, seeded_, maqPenalty_,
		                           qualOrder_, metrics_);
	}

protected:
	const TEbwt* ebwt_;
	bool         fw_;
	uint32_t     qualThresh_;
	bool         reportExacts_;
	bool         verbose_;
	bool         quiet_;
	bool         halfAndHalf_;
	bool         seeded_;
	bool         maqPenalty_;
	bool         qualOrder_;
	AlignerMetrics *metrics_;
};

/**
 * What boundary within the alignment to "pin" a particular
 * backtracking constraint to.
 */
enum SearchConstraintExtent {
	PIN_TO_BEGINNING = 1, // depth 0; i.e., constraint is inactive
	PIN_TO_LEN,           // constraint applies to while alignment
	PIN_TO_HI_HALF_EDGE,  // constraint applies to hi-half of seed region
	PIN_TO_SEED_EDGE      // constraint applies to entire seed region
};

/**
 * Concrete RangeSourceDriver that deals properly with
 * GreedyDFSRangeSource by calling setOffs() with the appropriate
 * parameters when initializing it;
 */
class EbwtRangeSourceDriver :
	public SingleRangeSourceDriver<EbwtRangeSource>
{
public:
	EbwtRangeSourceDriver(
			EbwtSearchParams<String<Dna> >& params,
			EbwtRangeSource* rs,
			bool fw,
			bool seed,
			bool maqPenalty,
			bool qualOrder,
			HitSink& sink,
			HitSinkPerThread* sinkPt,
			uint32_t seedLen,
			bool nudgeLeft,
			SearchConstraintExtent rev0Off,
			SearchConstraintExtent rev1Off,
			SearchConstraintExtent rev2Off,
			SearchConstraintExtent rev3Off,
			vector<String<Dna5> >& os,
			bool verbose,
			bool quiet,
			bool mate1,
			ChunkPool* pool,
			int *btCnt) :
			SingleRangeSourceDriver<EbwtRangeSource>(
					params, rs, fw, sink, sinkPt, os, verbose,
					quiet, mate1, 0, pool, btCnt),
			seed_(seed),
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			rs_(rs), seedLen_(seedLen),
			nudgeLeft_(nudgeLeft),
			rev0Off_(rev0Off), rev1Off_(rev1Off),
			rev2Off_(rev2Off), rev3Off_(rev3Off),
			verbose_(verbose), quiet_(quiet)
	{
		if(seed_) assert_gt(seedLen, 0);
	}

	virtual ~EbwtRangeSourceDriver() { }

	bool seed() const { return seed_; }

	bool ebwtFw() const { return rs_->curEbwt()->fw(); }

	/**
	 * Called every time setQuery() is called in the parent class,
	 * after setQuery() has been called on the RangeSource but before
	 * initConts() has been called.
	 */
	virtual void initRangeSource(const String<char>& qual, bool fuzzy,
	                             int alts, const String<char>* altQuals)
	{
		// If seedLen_ is huge, then it will always cover the whole
		// alignment
		assert_eq(len_, seqan::length(qual));
		uint32_t s = (seedLen_ > 0 ? min(seedLen_, len_) : len_);
		uint32_t sLeft  = s >> 1;
		uint32_t sRight = s >> 1;
		// If seed has odd length, then nudge appropriate half up by 1
		if((s & 1) != 0) { if(nudgeLeft_) sLeft++; else sRight++; }
		uint32_t rev0Off = cextToDepth(rev0Off_, sRight, s, len_);
		uint32_t rev1Off = cextToDepth(rev1Off_, sRight, s, len_);
		uint32_t rev2Off = cextToDepth(rev2Off_, sRight, s, len_);
		uint32_t rev3Off = cextToDepth(rev3Off_, sRight, s, len_);
		// Truncate the pattern if necessary
		uint32_t qlen = (uint32_t)seqan::length(qual);
		if(seed_) {
			if(len_ > s) {
				rs_->setQlen(s);
				qlen = s;
			}
			assert(!rs_->reportExacts());
		}
		// If there are any Ns in the unrevisitable region, then this
		// driver is guaranteed to yield no fruit.
		uint16_t minCost = 0;
		if(rs_->reportExacts()) {
			// Keep minCost at 0
		} else if (!rs_->halfAndHalf() && rev0Off < s) {
			// Exacts not allowed, so there must be at least 1 mismatch
			// outside of the unrevisitable area
			minCost = 1 << 14;
			if(qualOrder_) {
				uint8_t lowQual = 0xff;
				for(uint32_t d = rev0Off; d < s; d++) {
					uint8_t lowAtPos;
					if(fuzzy) {
						lowAtPos = loPenaltyAt(qlen-d-1, alts, qual, altQuals);
					} else {
						lowAtPos = qual[qlen - d - 1];
					}
					if(lowAtPos < lowQual) lowQual = lowAtPos;
				}
				assert_lt(lowQual, 0xff);
				if(fuzzy) {
					minCost += lowQual;
				} else {
					minCost += mmPenalty(maqPenalty_, phredCharToPhredQual(lowQual));
				}
			}
		} else if(rs_->halfAndHalf() && sRight > 0 && sRight < (s-1)) {
			// Half-and-half constraints are active, so there must be
			// at least 1 mismatch in both halves of the seed
			assert(rs_->halfAndHalf());
			minCost = (seed_ ? 3 : 2) << 14;
			if(qualOrder_) {
				assert(rs_->halfAndHalf() == 2 || rs_->halfAndHalf() == 3);
				uint8_t lowQual1 = 0xff;
				for(uint32_t d = 0; d < sRight; d++) {
					uint8_t lowAtPos;
					if(fuzzy) {
						lowAtPos = loPenaltyAt(qlen-d-1, alts, qual, altQuals);
					} else {
						lowAtPos = qual[qlen - d - 1];
					}
					if(lowAtPos < lowQual1) lowQual1 = lowAtPos;
				}
				assert_lt(lowQual1, 0xff);
				if(fuzzy) {
					minCost += lowQual1;
				} else {
					minCost += mmPenalty(maqPenalty_, phredCharToPhredQual(lowQual1));
				}
				uint8_t lowQual2_1 = 0xff;
				uint8_t lowQual2_2 = 0xff;
				for(uint32_t d = sRight; d < s; d++) {
					uint8_t lowAtPos;
					if(fuzzy) {
						lowAtPos = loPenaltyAt(qlen-d-1, alts, qual, altQuals);
					} else {
						lowAtPos = qual[qlen - d - 1];
					}
					if(lowAtPos < lowQual2_1) {
						if(lowQual2_1 != 0xff) {
							lowQual2_2 = lowQual2_1;
						}
						lowQual2_1 = lowAtPos;
					} else if(lowAtPos < lowQual2_2) {
						lowQual2_2 = lowAtPos;
					}
				}
				assert_lt(lowQual2_1, 0xff);
				if(fuzzy) {
					minCost += lowQual2_1;
					if(rs_->halfAndHalf() > 2 && lowQual2_2 != 0xff) {
						minCost += lowQual2_2;
					}
				} else {
					minCost += mmPenalty(maqPenalty_, phredCharToPhredQual(lowQual2_1));
					if(rs_->halfAndHalf() > 2 && lowQual2_2 != 0xff) {
						minCost += mmPenalty(maqPenalty_, phredCharToPhredQual(lowQual2_2));
					}
				}
			}
		}
		if(verbose_) cout << "initRangeSource minCost: " << minCost << endl;
		this->minCostAdjustment_ = minCost;
		rs_->setOffs(sRight,   // depth of far edge of hi-half (only matters where half-and-half is possible)
		             s,        // depth of far edge of lo-half (only matters where half-and-half is possible)
		             rev0Off,  // depth above which we cannot backtrack
		             rev1Off,  // depth above which we may backtrack just once
		             rev2Off,  // depth above which we may backtrack just twice
		             rev3Off); // depth above which we may backtrack just three times
	}

protected:

	/**
	 * Convert a search constraint extent to an actual depth into the
	 * read.
	 */
	inline uint32_t cextToDepth(SearchConstraintExtent cext,
	                            uint32_t sRight,
	                            uint32_t s,
	                            uint32_t len)
	{
		if(cext == PIN_TO_SEED_EDGE)    return s;
		if(cext == PIN_TO_HI_HALF_EDGE) return sRight;
		if(cext == PIN_TO_BEGINNING)    return 0;
		if(cext == PIN_TO_LEN)          return len;
		cerr << "Bad SearchConstraintExtent: " << cext;
		throw 1;
	}

	bool seed_;
	bool maqPenalty_;
	bool qualOrder_;
	EbwtRangeSource* rs_;
	uint32_t seedLen_;
	bool nudgeLeft_;
	SearchConstraintExtent rev0Off_;
	SearchConstraintExtent rev1Off_;
	SearchConstraintExtent rev2Off_;
	SearchConstraintExtent rev3Off_;
	bool verbose_;
	bool quiet_;
};

/**
 * Concrete RangeSourceDriver that deals properly with
 * GreedyDFSRangeSource by calling setOffs() with the appropriate
 * parameters when initializing it;
 */
class EbwtRangeSourceDriverFactory {
public:
	EbwtRangeSourceDriverFactory(
			EbwtSearchParams<String<Dna> >& params,
			EbwtRangeSourceFactory* rs,
			bool fw,
			bool seed,
			bool maqPenalty,
			bool qualOrder,
			HitSink& sink,
			HitSinkPerThread* sinkPt,
			uint32_t seedLen,
			bool nudgeLeft,
			SearchConstraintExtent rev0Off,
			SearchConstraintExtent rev1Off,
			SearchConstraintExtent rev2Off,
			SearchConstraintExtent rev3Off,
			vector<String<Dna5> >& os,
			bool verbose,
			bool quiet,
			bool mate1,
			ChunkPool* pool,
			int *btCnt = NULL) :
			params_(params),
			rs_(rs),
			fw_(fw),
			seed_(seed),
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			sink_(sink),
			sinkPt_(sinkPt),
			seedLen_(seedLen),
			nudgeLeft_(nudgeLeft),
			rev0Off_(rev0Off),
			rev1Off_(rev1Off),
			rev2Off_(rev2Off),
			rev3Off_(rev3Off),
			os_(os),
			verbose_(verbose),
			quiet_(quiet),
			mate1_(mate1),
			pool_(pool),
			btCnt_(btCnt)
	{ }

	~EbwtRangeSourceDriverFactory() {
		delete rs_; rs_ = NULL;
	}

	/**
	 * Return a newly-allocated EbwtRangeSourceDriver with the given
	 * parameters.
	 */
	EbwtRangeSourceDriver *create() const {
		return new EbwtRangeSourceDriver(
				params_, rs_->create(), fw_, seed_, maqPenalty_,
				qualOrder_, sink_, sinkPt_, seedLen_, nudgeLeft_,
				rev0Off_, rev1Off_, rev2Off_, rev3Off_, os_, verbose_,
				quiet_, mate1_, pool_, btCnt_);
	}

protected:
	EbwtSearchParams<String<Dna> >& params_;
	EbwtRangeSourceFactory* rs_;
	bool fw_;
	bool seed_;
	bool maqPenalty_;
	bool qualOrder_;
	HitSink& sink_;
	HitSinkPerThread* sinkPt_;
	uint32_t seedLen_;
	bool nudgeLeft_;
	SearchConstraintExtent rev0Off_;
	SearchConstraintExtent rev1Off_;
	SearchConstraintExtent rev2Off_;
	SearchConstraintExtent rev3Off_;
	vector<String<Dna5> >& os_;
	bool verbose_;
	bool quiet_;
	bool mate1_;
	ChunkPool* pool_;
	int *btCnt_;
};

/**
 * A RangeSourceDriver that manages two child EbwtRangeSourceDrivers,
 * one for searching for seed strings with mismatches in the hi-half,
 * and one for extending those seed strings toward the 3' end.
 */
class EbwtSeededRangeSourceDriver : public RangeSourceDriver<EbwtRangeSource> {
	typedef RangeSourceDriver<EbwtRangeSourceDriver>* TRangeSrcDrPtr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
public:
	EbwtSeededRangeSourceDriver(
			EbwtRangeSourceDriverFactory* rsFact,
			EbwtRangeSourceDriver* rsSeed,
			bool fw,
			uint32_t seedLen,
			bool verbose,
			bool quiet,
			bool mate1) :
			RangeSourceDriver<EbwtRangeSource>(true, 0),
			rsFact_(rsFact), rsFull_(false, NULL, verbose, quiet, true),
			rsSeed_(rsSeed), patsrc_(NULL), seedLen_(seedLen), fw_(fw),
			mate1_(mate1), seedRange_(0)
	{
		assert(rsSeed_->seed());
	}

	virtual ~EbwtSeededRangeSourceDriver() {
		delete rsFact_; rsFact_ = NULL;
		delete rsSeed_; rsSeed_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQueryImpl(PatternSourcePerThread* patsrc, Range *partial) {
		this->done = false;
		rsSeed_->setQuery(patsrc, partial);
		this->minCostAdjustment_ = max(rsSeed_->minCostAdjustment_, rsSeed_->minCost);
		this->minCost = this->minCostAdjustment_;
		rsFull_.clearSources();
		rsFull_.setQuery(patsrc, partial);
		rsFull_.minCost = this->minCost;
		assert_gt(rsFull_.minCost, 0);
		patsrc_ = patsrc;
		// The minCostAdjustment comes from the seed range source
		// driver, based on Ns and quals in the hi-half
		this->foundRange = false;
		ASSERT_ONLY(allTops_.clear());
		assert_eq(this->minCost, min<uint16_t>(rsSeed_->minCost, rsFull_.minCost));
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advance(int until) {
		assert(!this->foundRange);
		until = max<int>(until, ADV_COST_CHANGES);
		ASSERT_ONLY(uint16_t preCost = this->minCost);
		advanceImpl(until);
		if(this->foundRange) {
			assert_eq(range().cost, preCost);
		}
#ifndef NDEBUG
		if(this->foundRange) {
			// Assert that we have not yet dished out a range with this
			// top offset
			assert_gt(range().bot, range().top);
			assert(range().ebwt != NULL);
			TIndexOff top = (TIndexOff)range().top;
			top++; // ensure it's not 0
			if(!range().ebwt->fw()) top = -top;
			assert(allTops_.find(top) == allTops_.end());
			allTops_.insert(top);
		}
#endif
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advanceImpl(int until) {
		assert(!this->done);
		assert(!this->foundRange);
		assert_gt(rsFull_.minCost, 0);
		// Advance the seed range source
		if(rsSeed_->done && rsFull_.done &&
		   !rsSeed_->foundRange && !rsFull_.foundRange)
		{
			this->done = true;
			return;
		}
		if(rsSeed_->done && !rsSeed_->foundRange) {
			rsSeed_->minCost = 0xffff;
			if(rsFull_.minCost > this->minCost) {
				this->minCost = rsFull_.minCost;
				// Cost changed, so return
				return;
			}
		}
		if(rsFull_.done && !rsFull_.foundRange) {
			rsFull_.minCost = 0xffff;
			if(rsSeed_->minCost > this->minCost) {
				this->minCost = rsSeed_->minCost;
				// Cost changed, so return
				return;
			}
		}
		assert(rsSeed_->minCost != 0xffff || rsFull_.minCost != 0xffff);
		// Extend a partial alignment
		ASSERT_ONLY(uint16_t oldMinCost = this->minCost);
		assert_eq(this->minCost, min<uint16_t>(rsSeed_->minCost, rsFull_.minCost));
		bool doFull = rsFull_.minCost <= rsSeed_->minCost;
		if(!doFull) {
			// Advance the partial-alignment generator
			assert_eq(rsSeed_->minCost, this->minCost);
			if(!rsSeed_->foundRange) {
				rsSeed_->advance(until);
			}
			if(rsSeed_->foundRange) {
				assert_eq(this->minCost, rsSeed_->range().cost);
				assert_eq(oldMinCost, rsSeed_->range().cost);
				seedRange_ = &rsSeed_->range();
				rsSeed_->foundRange = false;
				assert_geq(seedRange_->cost, this->minCostAdjustment_);
				this->minCostAdjustment_ = seedRange_->cost;
				assert_gt(seedRange_->numMms, 0);
				// Keep the range for the hi-half partial alignment so
				// that the driver can (a) modify the pattern string
				// and (b) modify results from the RangeSource to
				// include these edits.
				EbwtRangeSourceDriver *partial = rsFact_->create();
				partial->minCost = seedRange_->cost;
				rsFull_.minCost = seedRange_->cost;
				rsFull_.addSource(partial, seedRange_);
				if(rsFull_.foundRange) {
					this->foundRange = true;
					rsFull_.foundRange = false;
					assert(rsFull_.range().repOk());
					assert_eq(range().cost, oldMinCost);
				}
			}
			if(rsSeed_->minCost > this->minCost) {
				this->minCost = rsSeed_->minCost;
				if(!rsFull_.done) {
					this->minCost = min(this->minCost, rsFull_.minCost);
					assert_eq(this->minCost, min<uint16_t>(rsSeed_->minCost, rsFull_.minCost));
				}
			}
		} else {
			// Extend a full alignment
			assert(!rsFull_.done);
			assert(!rsFull_.foundRange);
			uint16_t oldFullCost = rsFull_.minCost;
			if(!rsFull_.foundRange) {
				rsFull_.advance(until);
			}
			// Found a minimum-cost range
			if(rsFull_.foundRange) {
				this->foundRange = true;
				rsFull_.foundRange = false;
				assert(rsFull_.range().repOk());
				assert_eq(range().cost, oldMinCost);
			}
			assert_geq(rsFull_.minCost, oldFullCost);
			// Did the min cost change?
			if(rsFull_.minCost > oldFullCost) {
				// If a range was found, hold on to it and save it for
				// later.  Update the minCost.
				assert(!rsSeed_->done || rsSeed_->minCost == 0xffff);
				this->minCost = min(rsFull_.minCost, rsSeed_->minCost);
			}
		}
	}

	/**
	 * Return the range found.
	 */
	virtual Range& range() {
		Range& r = rsFull_.range();
		r.fw = fw_;
		r.mate1 = mate1_;
		return r;
	}

	/**
	 * Return whether we're generating ranges for the first or the
	 * second mate.
	 */
	virtual bool mate1() const {
		return mate1_;
	}

	/**
	 * Return true iff current pattern is forward-oriented.
	 */
	virtual bool fw() const {
		return fw_;
	}

protected:

	EbwtRangeSourceDriverFactory* rsFact_;
	TCostAwareRangeSrcDr rsFull_;
	EbwtRangeSourceDriver* rsSeed_;
	PatternSourcePerThread* patsrc_;
	uint32_t seedLen_;
	bool fw_;
	bool mate1_;
	bool generating_;
	Range *seedRange_;
};

#endif /*EBWT_SEARCH_BACKTRACK_H_*/
