#ifndef EBWT_SEARCH_BACKTRACK_H_
#define EBWT_SEARCH_BACKTRACK_H_

#include "pat.h"
#include "qual.h"
#include "ebwt_search_util.h"
#include "range.h"
#include "range_source.h"

/**
 * Class that coordinates quality- and quantity-aware backtracking over
 * some range of a read sequence.
 *
 * The creator can configure the BacktrackManager to treat different
 * stretches of the read differently.
 */
class GreedyDFSRangeSource : public RangeSource {
	typedef std::pair<int, int> TIntPair;
public:
	GreedyDFSRangeSource(
			const Ebwt<String<Dna> >* __ebwt,
			const EbwtSearchParams<String<Dna> >& __params,
			uint32_t __qualThresh,  /// max acceptable q-distance
			const BacktrackLimits& __maxBts, /// maximum # backtracks allowed
			uint32_t __reportPartials = 0,
			bool __reportExacts = true,
			bool __reportRanges = false,
			PartialAlignmentManager* __partials = NULL,
			String<QueryMutation>* __muts = NULL,
			bool __verbose = true,
			uint32_t seed = 0,
			vector<String<Dna5> >* __os = NULL,
			bool __considerQuals = true,  // whether to consider quality values when making backtracking decisions
			bool __halfAndHalf = false, // hacky way of supporting separate revisitable regions
			bool __maqPenalty = true) :
	    RangeSource(),
	    _started(false),
		_qry(NULL),
		_qlen(0),
		_qual(NULL),
		_name(NULL),
		_ebwt(__ebwt),
		_params(__params),
		_unrevOff(0),
		_1revOff(0),
		_2revOff(0),
		_3revOff(0),
		_maqPenalty(__maqPenalty),
		_qualThresh(__qualThresh),
		_pairs(NULL),
		_elims(NULL),
		_mms(),
		_refcs(),
		_chars(NULL),
		_reportPartials(__reportPartials),
		_reportExacts(__reportExacts),
		_reportRanges(__reportRanges),
		_partials(__partials),
		_muts(__muts),
		_os(__os),
		_sanity(_os != NULL && _os->size() > 0),
		_considerQuals(__considerQuals),
		_halfAndHalf(__halfAndHalf),
		_5depth(0),
		_3depth(0),
		_numBts(0),
		_totNumBts(0),
		_maxBts(__maxBts.maxBts),
		_precalcedSideLocus(false),
		_preLtop(),
		_preLbot(),
		_maxBts0(__maxBts.maxBts0),
		_maxBts1(__maxBts.maxBts1),
		_maxBts2(__maxBts.maxBts2),
		_rand(seed),
		_randSeed(seed),
		_verbose(__verbose),
		_ihits(0llu)
	{
		curEbwt_ = __ebwt;
	}

	virtual ~GreedyDFSRangeSource() {
		if(_pairs != NULL)          delete[] _pairs;
		if(_elims != NULL)          delete[] _elims;
		if(_chars != NULL)          delete[] _chars;
	}

	/**
	 * Set a new query read.
	 */
	virtual void setQuery(String<Dna5>* __qry,
	                      String<char>* __qual,
	                      String<char>* __name)
	{
		_qry = __qry;
		_qual = __qual;
		_name = __name;
		assert(_qry != NULL);
		assert(_qual != NULL);
		assert(_name != NULL);
		// Make sure every qual is a valid qual ASCII character (>= 33)
		for(size_t i = 0; i < length(*_qual); i++) {
			assert_geq((*_qual)[i], 33);
		}
		// Reset _qlen
		if(length(*_qry) > _qlen) {
			_qlen = length(*_qry);
			// Resize _pairs
	 		if(_pairs != NULL) { delete[] _pairs; }
	 		_pairs = new uint32_t[_qlen*_qlen*8];
	 		// Resize _elims
	 		if(_elims != NULL) { delete[] _elims; }
	 		_elims = new uint8_t[_qlen*_qlen];
	 		memset(_elims, 0, _qlen*_qlen);
			// Resize _chars
	 		if(_chars != NULL) { delete[] _chars; }
			_chars = new char[_qlen];
			_refcs.resize(_qlen);
			assert(_pairs != NULL);
			assert(_elims != NULL);
			assert(_chars != NULL);
		} else {
			// New length is less than old length, so there's no need
			// to resize any data structures.
			assert(_pairs != NULL);
			assert(_elims != NULL);
			assert(_chars != NULL);
			_qlen = length(*_qry);
		}
		_mms.clear();
		assert_geq(length(*_qual), _qlen);
		if(_verbose) {
			cout << "setQuery(_qry=" << (*_qry) << ", _qual=" << (*_qual) << ")" << endl;
		}
		_started = false;
		this->done = false;
		this->foundRange = false;
		// Initialize the random source using new read as part of the
		// seed.
		_rand.init(genRandSeed(*_qry, *_qual, *_name) + _randSeed);
	}

	/**
	 * Apply a batch of mutations to this read, possibly displacing a
	 * previous batch of mutations.
	 */
	void setMuts(String<QueryMutation>* __muts) {
		if(_muts != NULL) {
			// Undo previous mutations
			assert_gt(length(*_muts), 0);
			undoPartialMutations();
		}
		_muts = __muts;
		if(_muts != NULL) {
			assert_gt(length(*_muts), 0);
			applyPartialMutations();
		}
	}

	/**
	 * Set backtracking constraints.
	 */
	void setOffs(uint32_t __5depth,   // depth of far edge of hi-half
	             uint32_t __3depth,   // depth of far edge of lo-half
	             uint32_t __unrevOff, // depth above which we cannot backtrack
	             uint32_t __1revOff,  // depth above which we may backtrack just once
	             uint32_t __2revOff,  // depth above which we may backtrack just twice
	             uint32_t __3revOff)  // depth above which we may backtrack just three times
	{
		_5depth   = __5depth;
		_3depth   = __3depth;
		assert_geq(_3depth, _5depth);
		_unrevOff = __unrevOff;
		_1revOff  = __1revOff;
		_2revOff  = __2revOff;
		_3revOff  = __3revOff;
	}

	/// Reset number of backtracks to 0
	void resetNumBacktracks() {
		_totNumBts = 0;
	}

	/// Return number of backtracks since last reset
	uint32_t numBacktracks() {
		return _totNumBts;
	}

	void setReportExacts(int stratum) {
		_reportExacts = stratum;
	}

	/**
	 * Return true iff this RangeSource is allowed to report exact
	 * alignments (exact = no edits).
	 */
	bool reportExats() const {
		return _reportExacts;
	}

	void setEbwt(const Ebwt<String<Dna> >* ebwt) {
		_ebwt = ebwt;
	}

	/// Return the current range
	virtual Range& range() {
		return _curRange;
	}

	/**
	 * Set _qlen according to parameter, except don't let it fall below
	 * the length of the query.
	 */
	void setQlen(uint32_t qlen) {
		assert(_qry != NULL);
		_qlen = min<uint32_t>(length(*_qry), qlen);
	}

	/// Return the maximum number of allowed backtracks in a given call
	/// to backtrack()
	uint32_t maxBacktracks() {
		return _maxBts;
	}

	/**
	 * Set whether this manager will report hits as matrix ranges
	 * (true) or as fully-resolved hits (false).
	 */
	void setreportRanges(bool a) {
		_reportRanges = a;
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
		uint32_t m = min<uint32_t>(_unrevOff, _qlen);
		if(nsInFtab == 0 && m >= (uint32_t)ftabChars) {
			uint32_t ftabOff = calcFtabOff();
			uint32_t top = ebwt.ftabHi(ftabOff);
			uint32_t bot = ebwt.ftabLo(ftabOff+1);
			if(_qlen == (uint32_t)ftabChars && bot > top) {
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
					ret = reportAlignment(0, top, bot);
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
		this->done = true;
		if(finalize()) ret = true;
		return ret;
	}

	/**
	 * Initiate continuations so that the next call to advance() begins
	 * a new search.  Note that contMan is empty upon return if there
	 * are no valid continuations to begin with.  Also note that
	 * calling initConts() may result in finding a range (i.e., if we
	 * immediately jump to a valid range using the ftab).
	 */
	virtual void
	initBranch(PathManager& pm, uint32_t ham) {
		assert(curEbwt_ != NULL);
		assert_gt(length(*_qry), 0);
		assert_leq(_qlen, length(*_qry));
		assert_geq(length(*_qual), length(*_qry));
		const Ebwt<String<Dna> >& ebwt = *_ebwt;
		int ftabChars = ebwt._eh._ftabChars;
		this->foundRange = false;
		int nsInSeed = 0; int nsInFtab = 0;
		ASSERT_ONLY(allTops_.clear());
		if(!tallyNs(nsInSeed, nsInFtab)) {
			// No alignments are possible because of the distribution
			// of Ns in the read in combination with the backtracking
			// constraints.
			return;
		}
		// m = depth beyond which ftab must not extend or else we might
		// miss some legitimate paths
		uint32_t m = min<uint32_t>(_unrevOff, _qlen);
		// Let skipInvalidExact = true if using the ftab would be a
		// waste because it would jump directly to an alignment we
		// couldn't use.
		bool ftabSkipsToEnd = (_qlen == (uint32_t)ftabChars);
		bool skipInvalidExact = (!_reportExacts && ftabSkipsToEnd);
		bool skipInvalidPartial = (_reportPartials > 0 && ftabSkipsToEnd);
		// If it's OK to use the ftab...
		if(nsInFtab == 0 && m >= (uint32_t)ftabChars &&
		   !skipInvalidExact && !skipInvalidPartial)
		{
			// Use the ftab to jump 'ftabChars' chars into the read
			// from the right
			uint32_t ftabOff = calcFtabOff();
			uint32_t top = ebwt.ftabHi(ftabOff);
			uint32_t bot = ebwt.ftabLo(ftabOff+1);
			if(_qlen == (uint32_t)ftabChars && bot > top) {
				// We found a range with 0 mismatches immediately.  Set
				// fields to indicate we found a range.
				assert_eq(0, _reportPartials);
				assert(_reportExacts);
				_curRange.top     = top;
				_curRange.bot     = bot;
				_curRange.stratum = calcStratum(_mms, 0);
				_curRange.numMms  = 0;
				_curRange.ebwt    = _ebwt;
				_curRange.fw      = _params.fw();
				_curRange.mms.clear(); // no mismatches
				// no need to do anything with _curRange.refcs
				this->foundRange  = true;
				return;
			} else if (bot > top) {
				// We have a range to extend
				assert_leq(top, ebwt._eh._len);
				assert_leq(bot, ebwt._eh._len);
				Branch *b = pm.bpool.alloc();
				b->init(pm.rpool, _qlen,
				        _unrevOff, _1revOff, _2revOff, _3revOff,
				        0, ftabChars, ham, top, bot,
				        ebwt._eh, ebwt._ebwt, 0);
				pm.push(b); // insert into priority queue
			} else {
				// The arrows are already closed within the
				// unrevisitable region; give up
			}
		} else {
			// We can't use the ftab, so we start from the rightmost
			// position and use _fchr
			Branch *b = pm.bpool.alloc();
			b->init(pm.rpool, _qlen,
			        _unrevOff, _1revOff, _2revOff, _3revOff,
			        0, 0, ham, 0, 0, ebwt._eh, ebwt._ebwt, 0);
			pm.push(b); // insert into priority queue
		}
		return;
	}

	virtual void
	advanceBranch(PathManager& pm) {
		// Restore alignment state from the frontmost continuation.
		// The frontmost continuation should in principle be the most
		// promising partial alignment found so far.  In the case of
		// the greedy DFS backtracker, which partial alignment is most
		// promising is determined in a greedy, depth-first, non-
		// optimal fashion.
		assert(curEbwt_ != NULL);

		// Let this->foundRange = false; we'll set it to true iff this call
		// to advance yielded a new valid-alignment range.
		this->foundRange = false;

		// Can't have already exceeded weighted hamming distance threshold
		assert_gt(length(*_qry), 0);
		assert_leq(_qlen, length(*_qry));
		assert_geq(length(*_qual), length(*_qry));

		do {
			// Get the highest-priority branch according to the priority
			// queue in 'pm'
			Branch* br = pm.front();
			// Shouldn't be curtailed or exhausted
			assert(!br->exhausted_);
			assert(!br->curtailed_);
			if(_verbose) {
				br->print((*_qry), (*_qual), cout, _halfAndHalf, _params.fw(), _ebwt->fw());
			}
			assert(br->repOk(_qlen));

			ASSERT_ONLY(int stratum = br->cost_ >> 14); // shift the stratum over
			assert_lt(stratum, 4);
			uint16_t ham = br->cost_ & ~0xc000; // AND off the stratum
			assert_leq(ham, 0x3000);
			assert_leq(ham, _qualThresh);
			uint32_t depth = br->tipDepth();

			const Ebwt<String<Dna> >& ebwt = *_ebwt;

			if(_halfAndHalf) {
				assert_eq(0, _reportPartials);
				assert_gt(_3depth, _5depth);
			}
			if(_reportPartials) {
				assert(!_halfAndHalf);
				assert_leq(br->numEdits_, _reportPartials);
			}

			// Compiler complains if our goto's jump over these
			// initializations, so stick them here
			bool backtrackDespiteMatch = false;
			bool reportedPartial = false;
			bool invalidExact = false;
			bool empty = false;
			bool hit = false;
			uint32_t cur = 0;

			if(_halfAndHalf && !hhCheckTop(br, depth, 0 /*iham*/)) {
				// Stop extending this branch because it violates a half-
				// and-half constraint
				pm.curtail(br, _3depth);
				goto bail;
			}

			cur = _qlen - depth - 1; // current offset into _qry
			if(depth < _qlen) {
				// Determine whether ranges at this location are candidates
				// for backtracking
				int c = (int)(*_qry)[cur]; // get char at this position
				int nextc = -1;
				if(cur < _qlen-1) nextc = (int)(*_qry)[cur+1];
				assert_leq(c, 4);
				uint8_t q = qualAt(cur);   // get qual at this position

				// The current query position is a legit alternative if it a) is
				// not in the unrevisitable region, and b) its selection would
				// not cause the quality ceiling (if one exists) to be exceeded
				bool curIsAlternative =
					 (depth >= br->depth0_) &&
					 (!_considerQuals ||
					  (ham + mmPenalty(_maqPenalty, q) <= _qualThresh));

#ifndef NDEBUG
				uint32_t obot = br->bot_;
#endif
				uint32_t otop = br->top_;

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
				// Should be all 0s
				assert_eq(0, rs->tops[0]); assert_eq(0, rs->bots[0]);
				assert_eq(0, rs->tops[1]); assert_eq(0, rs->bots[1]);
				assert_eq(0, rs->tops[2]); assert_eq(0, rs->bots[2]);
				assert_eq(0, rs->tops[3]); assert_eq(0, rs->bots[3]);

				// Calculate the ranges for this position
				if(br->top_ == 0 && br->bot_ == 0) {
					// Calculate first quartet of ranges using the _fchr[]
					// array
								  rs->tops[0] = ebwt._fchr[0];
					rs->bots[0] = rs->tops[1] = ebwt._fchr[1];
					rs->bots[1] = rs->tops[2] = ebwt._fchr[2];
					rs->bots[2] = rs->tops[3] = ebwt._fchr[3];
					rs->bots[3]               = ebwt._fchr[4];
					ASSERT_ONLY(int r =) br->installRanges(c, nextc, q);
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
						ebwt.mapLFEx(br->ltop_, br->lbot_, rs->tops, rs->bots);
					} else {
#ifndef NDEBUG
						uint32_t tmptops[] = {0, 0, 0, 0};
						uint32_t tmpbots[] = {0, 0, 0, 0};
						SideLocus ltop, lbot;
						ltop.initFromRow(otop, _ebwt->_eh, _ebwt->_ebwt);
						lbot.initFromRow(obot, _ebwt->_eh, _ebwt->_ebwt);
						ebwt.mapLFEx(ltop, lbot, tmptops, tmpbots);
#endif
						int cc = ebwt.mapLF1(otop, br->ltop_);
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
					ASSERT_ONLY(int r =) br->installRanges(c, nextc, q);
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
							br->bot_ = br->top_ = ebwt.mapLF1(br->top_, br->ltop_, c);
							if(br->bot_ != 0xffffffff) br->bot_++;
						} else {
							br->top_ = ebwt.mapLF(br->ltop_, c);
							assert(br->lbot_.valid());
							br->bot_ = ebwt.mapLF(br->lbot_, c);
						}
					}
				} else {
					rs->eliminated_ = true;
				}
				// br->top_ and br->bot_ now contain the next top and bot
			} else {
				// The continuation had already processed the whole read
				assert_eq(_qlen, depth);
				cur = 0;
			}
			empty = (br->top_ == br->bot_);
			hit = (cur == 0 && !empty);

			// Is this a potential partial-alignment range?
			if(hit && _reportPartials > 0) {
				// This is a potential alignment range; set
				// backtrackDespiteMatch to 'true' so that we keep looking
				backtrackDespiteMatch = true;
				// We don't care to report exact partial alignments, only
				// ones with mismatches.
				if(br->numEdits_ > 0) {
					reportPartial(br->numEdits_);
					reportedPartial = true;
				}
				// Now continue on to find legitimate partial
				// mismatches
			}

			// Check whether we've obtained an exact alignment when
			// we've been instructed not to report exact alignments
			invalidExact = (hit && br->numEdits_ == 0 && !_reportExacts);

			// Set this to true if the only way to make legal progress
			// is via one or more additional backtracks.
			if(_halfAndHalf && !hhCheck(br, depth, empty)) {
				// This alignment doesn't satisfy the half-and-half
				// requirements; reject it
				pm.curtail(br, _3depth);
				goto bail;
			}

			if(hit &&            // there is a range to report
			   !invalidExact &&  // not disqualified by no-exact-hits setting
			   !reportedPartial) // not an already-reported partial alignment
			{
				if(_verbose) {
					br->len_++;
					br->print((*_qry), (*_qual), cout, _halfAndHalf, _params.fw(), _ebwt->fw());
					br->len_--;
				}
				assert_gt(br->bot_, br->top_);
				_curRange.top     = br->top_;
				_curRange.bot     = br->bot_;
				_curRange.stratum = (br->cost_ >> 14);
				_curRange.numMms  = br->numEdits_;
				_curRange.fw      = _params.fw();
				_curRange.mms.clear();
				_curRange.refcs.clear();
				for(size_t i = 0; i < br->numEdits_; i++) {
					_curRange.mms.push_back(_qlen - br->edits_[i].pos - 1);
					_curRange.refcs.push_back("ACGTN"[br->edits_[i].chr]);
				}
				_curRange.ebwt    = _ebwt;
				this->foundRange  = true;
	#ifndef NDEBUG
				int64_t top2 = (int64_t)br->top_;
				top2++; // ensure it's not 0
				if(_ebwt->fw()) top2 = -top2;
				assert(allTops_.find(top2) == allTops_.end());
				allTops_.insert(top2);
	#endif
				// Must curtail because we've consumed the whole pattern
				pm.curtail(br, _3depth);
			} else if(empty || cur == 0) {
				// The branch couldn't be extended further
				pm.curtail(br, _3depth);
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
			pm.splitAndPrep(_rand, _qlen, _3depth, _ebwt->_eh, _ebwt->_ebwt);

		} while(!this->foundRange && !pm.empty());
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
	               uint32_t top,
	               uint32_t bot,
	               uint32_t iham = 0,
	               bool disableFtab = false)
	{
		HitSinkPerThread& sink = _params.sink();
		bool oldRetain = sink.retainHits();
		_ihits = sink.retainedHits().size();
		if(_sanity) {
			// Save some info about hits retained at this point
			sink.setRetainHits(true);
		}

		// Initiate the recursive, randomized quality-aware backtracker
		// with a stack depth of 0 (no backtracks so far)
		_bailedOnBacktracks = false;
		bool done = backtrack(0, depth, _unrevOff, _1revOff, _2revOff, _3revOff,
		                      top, bot, iham, iham, _pairs, _elims, disableFtab);

		// Remainder of this function is sanity checking
		sink.setRetainHits(oldRetain); // restore old value
		if(_sanity) {
			// Check that the alignments produced for the read and
			// backtracking constraints are consistent with results
			// obtained by a naive oracle
			sanityCheckHits(sink.retainedHits(), sink.retainedStrata(), _ihits, iham);
		}
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
	               uint32_t  top,      // top arrow in pair prior to 'depth'
	               uint32_t  bot,      // bottom arrow in pair prior to 'depth'
	               uint32_t  ham,      // weighted hamming distance so far
	               uint32_t  iham,     // initial weighted hamming distance
	               uint32_t* pairs,    // portion of pairs array to be used for this backtrack frame
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
		// The total number of ranges that are acceptable
		// backtracking targets ("alternative" ranges)
		uint32_t altNum = 0;
		// The total number of ranges that are candidates to be
		// the *next* backtracking target because they are low quality
		// ("eligible" ranges)
		uint32_t eligibleNum = 0;
		// Total distance between all lowest-quality "alternative"
		// ranges that haven't yet been eliminated
		uint32_t eligibleSz = 0;
		// If there is just one eligible slot at the moment (a common
		// case), these are its parameters
		uint32_t eli = 0;
		bool     elignore = true; // ignore the el values because they didn't come from a recent override
		uint32_t eltop = 0;
		uint32_t elbot = 0;
		uint32_t elham = ham;
		char     elchar = 0;
		int      elcint = 0;
		// The lowest quality value associated with any alternative
		// ranges; all alternative ranges with this quality are
		// eligible
		uint8_t lowAltQual = 0xff;
		uint32_t d = depth;
		uint32_t cur = _qlen - d - 1; // current offset into _qry
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
				memset(&pairs[d*8], 0, 8 * 4);
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
						if(bot != 0xffffffff) bot++;
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
					uint32_t spread = pairSpread(pairs, d, i);
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
						if(sink.numValidHits() == prehits) {
							confirmNoHit(iham);
						}
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
						uint32_t d = _qlen - _mms[i] - 1;
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
						if(stackDepth == 0 && sink.numValidHits() == prehits) {
							// We're returning from the bottommost frame
							// without having found any hits; let's
							// sanity-check that there really aren't any
							confirmNoHit(iham);
						}
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
				uint64_t rhits = sink.numReportedHits();
				bool ret = reportAlignment(stackDepth, top, bot);
				if(_sanity && sink.numReportedHits() > rhits) {
					confirmHit(iham);
				}
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
				uint32_t bttop = 0;
				uint32_t btbot = 0;
				uint32_t btham = ham;
				char     btchar = 0;
				int      btcint = 0;
				uint32_t icur = 0;
				// The common case is that eligibleSz == 1
				if(eligibleNum > 1 || elignore) {
					bool foundTarget = false;
					// Walk from left to right
					for(; i >= depth; i--) {
						assert_geq(i, unrevOff);
						icur = _qlen - i - 1; // current offset into _qry
						uint8_t qi = qualAt(icur);
						assert_lt(elims[i], 16); // 1.26% in profile (next or prev?)
						if((qi == lowAltQual || !_considerQuals) && elims[i] != 15) {
							// This is the leftmost eligible position with at
							// least one remaining backtrack target
							uint32_t posSz = 0;
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
										foundTarget = true;
										bttop = pairTop(pairs, i, j);
										btbot = pairBot(pairs, i, j);
										btham += mmPenalty(_maqPenalty, qi);
										btcint = j;
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
				icur = _qlen - i - 1; // current offset into _qry
				// Slide over to the next backtacking frame within
				// pairs and elims; won't interfere with our frame or
				// any of our parents' frames
				uint32_t *newPairs = pairs + (_qlen*8);
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
				_refcs[stackDepth] = btchar;
#ifndef NDEBUG
				for(uint32_t j = 0; j < stackDepth; j++) {
					assert_neq(_mms[j], icur);
				}
#endif
				_chars[i] = btchar;
				assert_leq(i+1, _qlen);
				bool ret;
				if(i+1 == _qlen) {
					ret = reportAlignment(stackDepth+1, bttop, btbot);
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
					uint32_t ftabOff = (*_qry)[_qlen - ftabChars];
					assert_lt(ftabOff, 4);
					assert_lt(ftabOff, ebwt._eh._ftabLen-1);
					for(int j = ftabChars - 1; j > 0; j--) {
						ftabOff <<= 2;
						if(_qlen-j == icur) {
							ftabOff |= btcint;
						} else {
							assert_lt((uint32_t)(*_qry)[_qlen-j], 4);
							ftabOff |= (uint32_t)(*_qry)[_qlen-j];
						}
						assert_lt(ftabOff, ebwt._eh._ftabLen-1);
					}
					assert_lt(ftabOff, ebwt._eh._ftabLen-1);
					uint32_t ftabTop = ebwt.ftabHi(ftabOff);
					uint32_t ftabBot = ebwt.ftabLo(ftabOff+1);
					assert_geq(ftabBot, ftabTop);
					if(ftabTop == ftabBot) {
						ret = false;
					} else {
						assert(!_precalcedSideLocus);
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
					// Continue from selected alternative range
					ret = backtrack(stackDepth+1,// added 1 mismatch to alignment
				                    i+1,         // start from next position after
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
					// Not necessarily true that there is a hit to
					// confirm; it could be that true was returned
					// because we just ran up against the max, where
					// max > k
					//if(_sanity) confirmHit(iham);
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
					if(stackDepth == 0 && sink.numValidHits() == prehits) {
						// We're returning from the bottommost frame
						// without having found any hits; let's
						// sanity-check that there really aren't any
						confirmNoHit(iham);
					}
					return false;
				}
				else if(eligibleNum == 0 && _considerQuals) {
					// Find the next set of eligible backtrack points
					// by re-scanning this backtracking frame (from
					// 'depth' up to 'd')
					lowAltQual = 0xff;
					for(size_t k = d; k >= depth && k <= _qlen; k--) {
						uint32_t kcur = _qlen - k - 1; // current offset into _qry
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
										uint32_t spread = pairSpread(pairs, k, l);
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
											eli = k;
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
				// This full alignment was rejected for one of the above
				// reasons - return false to indicate we should keep
				// searching
				if(stackDepth == 0 && sink.numValidHits() == prehits) {
					// We're returning from the bottommost frame
					// without having found any hits; let's
					// sanity-check that there really aren't any
					confirmNoHit(iham);
				}
				return false;
			}
			// Mismatch with no alternatives
			if(top == bot && altNum == 0) {
				assert_eq(0, altNum);
				assert_eq(0, eligibleSz);
				assert_eq(0, eligibleNum);
				if(stackDepth == 0 && sink.numValidHits() == prehits) {
					// We're returning from the bottommost frame
					// without having found any hits; let's
					// sanity-check that there really aren't any
					confirmNoHit(iham);
				}
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
			ret = reportAlignment(stackDepth, top, bot);
		}
		if(!ret && stackDepth == 0 && sink.numValidHits() == prehits) {
			confirmNoHit(iham);
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
				uint32_t depth = _qlen - mms[i] - 1;
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
	 * Return true iff we're OK to continue after considering which
	 * half-seed boundary we're passing through, together with the
	 * number of mismatches accumulated so far.  Return false if we
	 * should stop because a half-and-half constraint is violated.  If
	 * we're not currently passing a half-seed boundary, just return
	 * true.
	 */
	bool hhCheck(Branch *b, uint32_t depth, bool empty) {
		ASSERT_ONLY(uint32_t lim = (_3revOff == _2revOff)? 2 : 3);
		if((depth == (_5depth-1)) && !empty) {
			// We're crossing the boundary separating the hi-half
			// from the non-seed portion of the read.
			// We should induce a mismatch if we haven't mismatched
			// yet, so that we don't waste time pursuing a match
			// that was covered by a previous phase
			assert_eq(0, _reportPartials);
			assert_leq(b->numEdits_, lim-1);
			return b->numEdits_ > 0;
		} else if((depth == (_3depth-1)) && !empty) {
			// We're crossing the boundary separating the lo-half
			// from the non-seed portion of the read
			assert_eq(0, _reportPartials);
			assert_leq(b->numEdits_, lim);
			assert_gt(b->numEdits_, 0);
			// Count the mismatches in the lo and hi halves
			uint32_t loHalfMms = 0, hiHalfMms = 0;
			for(size_t i = 0; i < b->numEdits_; i++) {
				uint32_t depth = b->edits_[i].pos;
				if     (depth < _5depth) hiHalfMms++;
				else if(depth < _3depth) loHalfMms++;
				else assert(false);
			}
			assert_leq(loHalfMms + hiHalfMms, lim);
			bool invalidHalfAndHalf = (loHalfMms == 0 || hiHalfMms == 0);
			return (b->numEdits_ >= 2 && !invalidHalfAndHalf);
		}
#ifndef NDEBUG
		if(depth < _5depth-1) {
			assert_leq(b->numEdits_, lim-1);
		}
		else if(depth >= _5depth && depth < _3depth-1) {
			assert_gt(b->numEdits_, 0);
			assert_leq(b->numEdits_, lim);
		}
#endif
		return true;
	}

	/**
	 * Calculate the stratum of the partial (or full) alignment
	 * currently under consideration.  Stratum is equal to the number
	 * of mismatches in the seed portion of the alignment.
	 */
	int calcStratum(const std::vector<uint32_t>& mms, uint32_t stackDepth) {
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
	                const std::vector<uint32_t>& mms,
	                uint64_t prehits = 0xffffffffffffffffllu)
	{
		HitSinkPerThread& sink = _params.sink();
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
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits)
					{
						// We're returning from the bottommost frame
						// without having reported any hits; let's
						// sanity-check that there really aren't any
						confirmNoHit(iham);
					}
					return false;
				}
			} else { // if(_3revOff != _2revOff)
				// Total of 3 mismatches allowed: 1 hi, 1 or 2 lo
				// The backtracking logic should have prevented us from
				// backtracking more than twice into this region
				assert_leq(stackDepth, 2);
				// Reject if we haven't encountered mismatch by this point
				if(stackDepth < 1) {
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits)
					{
						// We're returning from the bottommost frame
						// without having found any hits; let's
						// sanity-check that there really aren't any
						confirmNoHit(iham);
					}
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
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits &&
					   stackDepth == 0)
					{
						confirmNoHit(iham);
					}
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
					uint32_t d = _qlen - mms[i] - 1;
					if     (d < _5depth) hiHalfMms++;
					else if(d < _3depth) loHalfMms++;
					else assert(false);
				}
				assert_leq(loHalfMms + hiHalfMms, 3);
				assert_gt(hiHalfMms, 0);
				if(loHalfMms == 0) {
					// Didn't encounter any mismatches in the lo-half
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits &&
					   stackDepth == 0)
					{
						confirmNoHit(iham);
					}
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
		HitSinkPerThread& sink = _params.sink();
		assert_eq(0, _reportPartials);
		// Crossing from the hi-half into the lo-half
		if(d == _5depth) {
			if(_3revOff == _2revOff) {
				// Total of 2 mismatches allowed: 1 hi, 1 lo
				// The backtracking logic should have prevented us from
				// backtracking more than once into this region
				assert_leq(b->numEdits_, 1);
				// Reject if we haven't encountered mismatch by this point
				if(b->numEdits_ == 0) {
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits)
					{
						// We're returning from the bottommost frame
						// without having reported any hits; let's
						// sanity-check that there really aren't any
						confirmNoHit(iham);
					}
					return false;
				}
			} else { // if(_3revOff != _2revOff)
				// Total of 3 mismatches allowed: 1 hi, 1 or 2 lo
				// The backtracking logic should have prevented us from
				// backtracking more than twice into this region
				assert_leq(b->numEdits_, 2);
				// Reject if we haven't encountered mismatch by this point
				if(b->numEdits_ < 1) {
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits)
					{
						// We're returning from the bottommost frame
						// without having found any hits; let's
						// sanity-check that there really aren't any
						confirmNoHit(iham);
					}
					return false;
				}
			}
		} else if(d == _3depth) {
			// Crossing from lo-half to outside of the seed
			if(_3revOff == _2revOff) {
				// Total of 2 mismatches allowed: 1 hi, 1 lo
				// The backtracking logic should have prevented us from
				// backtracking more than twice within this region
				assert_leq(b->numEdits_, 2);
				// Must have encountered two mismatches by this point
				if(b->numEdits_ < 2) {
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits &&
					   b->numEdits_ == 0)
					{
						confirmNoHit(iham);
					}
					// We're returning from the bottommost frame
					// without having found any hits; let's
					// sanity-check that there really aren't any
					return false;
				}
			} else { // if(_3revOff != _2revOff)
				// Total of 3 mismatches allowed: 1 hi, 1 or 2 lo
				// Count the mismatches in the lo and hi halves
				int loHalfMms = 0, hiHalfMms = 0;
				for(size_t i = 0; i < b->numEdits_; i++) {
					uint32_t d = b->edits_[i].pos;
					if     (d < _5depth) hiHalfMms++;
					else if(d < _3depth) loHalfMms++;
					else assert(false);
				}
				assert_leq(loHalfMms + hiHalfMms, 3);
				assert_gt(hiHalfMms, 0);
				if(loHalfMms == 0) {
					// Didn't encounter any mismatches in the lo-half
					if(prehits != 0xffffffffffffffffllu &&
					   sink.numValidHits() == prehits &&
					   b->numEdits_ == 0)
					{
						confirmNoHit(iham);
					}
					// We're returning from the bottommost frame
					// without having found any hits; let's
					// sanity-check that there really aren't any
					return false;
				}
				assert_geq(b->numEdits_, 2);
				// The backtracking logic should have prevented us from
				// backtracking more than twice within this region
				assert_leq(b->numEdits_, 3);
			}
		} else {
			// We didn't just cross a boundary, so do an in-between check
			if(d >= _5depth) {
				assert_geq(b->numEdits_, 1);
			} else if(d >= _3depth) {
				assert_geq(b->numEdits_, 2);
			}
		}
		return true;
	}

	inline uint8_t qualAt(size_t off) {
		return phredCharToPhredQual((*_qual)[off]);
	}

	/// Get the top offset for character c at depth d
	inline uint32_t pairTop(uint32_t* pairs, size_t d, size_t c) {
		return pairs[d*8 + c + 0];
	}

	/// Get the bot offset for character c at depth d
	inline uint32_t pairBot(uint32_t* pairs, size_t d, size_t c) {
		return pairs[d*8 + c + 4];
	}

	/// Get the spread between the bot and top offsets for character c
	/// at depth d
	inline uint32_t pairSpread(uint32_t* pairs, size_t d, size_t c) {
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
			_refcs[stackDepth + i] = "ACGT"[(*_muts)[i].newBase];
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
	bool reportAlignment(uint32_t stackDepth, uint32_t top, uint32_t bot) {
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
			// Report the range of full alignments
			hit = reportFullAlignment(stackDepth + numMuts, top, bot, stratum);
			// Re-apply partial-alignment mutations
			applyPartialMutations();
			assert_eq(tmp, (*_qry));
		} else {
			// Report the range of full alignments
			hit = reportFullAlignment(stackDepth, top, bot, stratum);
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
	                         uint32_t top,
	                         uint32_t bot,
	                         int stratum)
	{
		assert_gt(bot, top);
		if(stackDepth == 0 && !_reportExacts) {
			// We are not reporting exact hits (usually because we've
			// already reported them as part of a previous invocation
			// of the backtracker)
			return false;
		}
		if(_reportRanges) {
			assert(_params.arrowMode());
			return _ebwt->report((*_qry), _qual, _name,
                    _mms, _refcs, stackDepth, 0,
                    top, bot, _qlen, stratum, _params);
		}
		uint32_t spread = bot - top;
		// Pick a random spot in the range to begin report
		uint32_t r = top + (_rand.nextU32() % spread);
		for(uint32_t i = 0; i < spread; i++) {
			uint32_t ri = r + i;
			if(ri >= bot) ri -= spread;
			// reportChaseOne takes the _mms[] list in terms of
			// their indices into the query string; not in terms
			// of their offset from the 3' or 5' end.
			if(_ebwt->reportChaseOne((*_qry), _qual, _name,
			                         _mms, _refcs, stackDepth, ri,
			                         top, bot, _qlen, stratum, _params))
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
		ASSERT_ONLY(qualTot += mmPenalty(_maqPenalty, ((*_qual)[_mms[0]] - 33)));
		uint32_t ci = _qlen - _mms[0] - 1;
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
			ASSERT_ONLY(qualTot += mmPenalty(_maqPenalty, ((*_qual)[_mms[1]] - 33)));
			ci = _qlen - _mms[1] - 1;
			// _chars[] is index in terms of RHS-relative depth
			c = (int)(Dna5)_chars[ci];
			assert_lt(c, 4);
			assert_neq(c, (int)(*_qry)[_mms[1]]);
			// Second, append the substituted character for the position
			al.entry.char1 = c;
		} else {
			// Signal that the '1' slot is empty
			al.entry.pos1 = 0xff;
		}

		if(stackDepth > 2) {
			assert_gt(_mms.size(), 2);
			// Second mismatch
			assert_lt(_mms[2], _qlen);
			// First, append the mismatch position in the read
			al.entry.pos2 = (uint16_t)_mms[2]; // pos
			ASSERT_ONLY(qualTot += mmPenalty(_maqPenalty, ((*_qual)[_mms[2]] - 33)));
			ci = _qlen - _mms[2] - 1;
			// _chars[] is index in terms of RHS-relative depth
			c = (int)(Dna5)_chars[ci];
			assert_lt(c, 4);
			assert_neq(c, (int)(*_qry)[_mms[2]]);
			// Second, append the substituted character for the position
			al.entry.char2 = c;
		} else {
			// Signal that the '2' slot is empty
			al.entry.pos2 = 0xff;
		}
		assert(validPartialAlignment(al));
#ifndef NDEBUG
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
	                            uint32_t* pairs,
	                            uint8_t*  elims)
	{
		// Sanity check that the lay of the land is as we
		// expect given eligibleNum and eligibleSz
		size_t i = max(depth, unrevOff), j = 0;
		uint32_t cumSz = 0;
		uint32_t eligiblesVisited = 0;
		for(; i <= d; i++) {
			uint32_t icur = _qlen - i - 1; // current offset into _qry
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

	/**
	 * Confirm that no
	 */
	void confirmNoHit(uint32_t iham) {
		// Not smart enough to deal with seedling hits yet
		if(!_sanity || _reportPartials > 0) return;
		vector<Hit> oracleHits;
		vector<int> oracleStrata;
		// Invoke the naive oracle, which will place all qualifying
		// hits in the 'oracleHits' vector
		naiveOracle(oracleHits, oracleStrata, iham);
		if(oracleHits.size() > 0) {
			// Oops, the oracle found at least one hit; print
			// detailed info about the first oracle hit (for
			// debugging)
			cout << "Oracle hit " << oracleHits.size()
				 << " times, but backtracker did not hit" << endl;
			for(size_t i = 0; i < oracleHits.size() && i < 3; i++) {
				const Hit& h = oracleHits[i];
				cout << "  Oracle hit " << i << ": " << endl;
				if(_muts != NULL) {
					undoPartialMutations();
					cout << "  Unmutated Pat:  " << prefix(*_qry, _qlen) << endl;
					applyPartialMutations();
				}
				cout << "  Pat:            " << prefix(*_qry, _qlen) << endl;
				cout << "  Tseg:           ";
				bool ebwtFw = _ebwt->fw();
				if(ebwtFw) {
					for(size_t i = 0; i < _qlen; i++) {
						cout << (*_os)[h.h.first][h.h.second + i];
					}
				} else {
					for(int i = (int)_qlen-1; i >= 0; i--) {
						cout << (*_os)[h.h.first][h.h.second + i];
					}
				}
				cout << endl;
				cout << "  Quals:          " << prefix(*_qual, _qlen) << endl;
				cout << "  Bt:             ";
				for(int i = (int)_qlen-1; i >= 0; i--) {
					if     (i < (int)_unrevOff) cout << "0";
					else if(i < (int)_1revOff)  cout << "1";
					else if(i < (int)_2revOff)  cout << "2";
					else if(i < (int)_3revOff)  cout << "3";
					else cout << "X";
				}
				cout << endl;
			}
		}
		assert_eq(0, oracleHits.size());
	}

	/**
	 *
	 */
	void confirmHit(uint32_t iham) {
		// Not smart enough to deal with seedling hits yet
		if(!_sanity || _reportPartials > 0) return;
		vector<Hit> oracleHits;
		vector<int> oracleStrata;
		// Invoke the naive oracle, which will place all qualifying
		// hits in the 'oracleHits' vector
		naiveOracle(oracleHits, oracleStrata, iham);
		vector<Hit>& retainedHits = _params.sink().retainedHits();
		assert_gt(retainedHits.size(), 0);
		if(oracleHits.size() == 0) {
			::printHit(
					(*_os), retainedHits.back(), (*_qry), _qlen, _unrevOff,
					_1revOff, _2revOff, _3revOff, _ebwt->fw());
		}
		// If we found a hit, it had better be one of the ones
		// that the oracle found
		assert_gt(oracleHits.size(), 0);
		// Get the hit reported by the backtracker
		Hit& rhit = retainedHits.back();
		// Go through oracleHits and look for a match
		size_t i;
		for(i = 0; i < oracleHits.size(); i++) {
			const Hit& h = oracleHits[i];
    		if(h.h.first == rhit.h.first &&
    		   h.h.second == rhit.h.second &&
    		   h.fw == rhit.fw)
    		{
    			assert(h.mms == rhit.mms);
    			// It's a match - hit confirmed
    			break;
    		}
		}
		assert_lt(i, oracleHits.size()); // assert we found a matchup
	}

	/**
	 * Sanity-check the final set of alignments generated for this read
	 * with respect to the current backtracking constraints.
	 */
	void sanityCheckHits(const vector<Hit>& hits,
	                     const vector<int>& strata,
	                     size_t first, uint32_t iham)
	{
		if(!_sanity || _reportPartials > 0 || _bailedOnBacktracks || first == hits.size()) {
			return;
		}
		assert_eq(hits.size(), strata.size());
		HitSinkPerThread& sink = _params.sink();
		uint32_t maxHitsAllowed = sink.maxHits();
		vector<Hit> oracleHits; vector<int> oracleStrata;
		// Invoke the naive oracle, which will place all qualifying
		// hits and their strata in the 'oracleHits' and
		// 'oracleStrata' vectors
		naiveOracle(oracleHits, oracleStrata, iham);
		assert_eq(oracleHits.size(), oracleStrata.size());
		assert_gt(oracleHits.size(), 0);
		int lastStratum = -1;
		int bestStratum = -1;
		int worstStratum = -1;
		// For each hit generated for this read
		for(size_t i = first; i < hits.size(); i++) {
			const Hit& rhit = hits[i];
			int stratum = strata[i];
			// Keep a running measure of the "worst" stratum
			// observed in the reported hits
			if(worstStratum == -1) {
				worstStratum = stratum;
			} else if(stratum > worstStratum) {
				worstStratum = stratum;
			}
			// Keep a running measure of the "best" stratum
			// observed in the reported hits
			if(bestStratum == -1) {
				bestStratum = stratum;
			} else if(stratum < bestStratum) {
				bestStratum = stratum;
			}
			if(lastStratum == -1) {
				lastStratum = stratum;
			} else if(!sink.spanStrata()) {
				// Retained hits are not allowed to "span" strata,
				// so this hit's stratum had better be the same as
				// the last hit's
				assert_eq(lastStratum, stratum);
			}
			// Go through oracleHits and look for a match
			size_t i;
			bool found = false;
			for(i = 0; i < oracleHits.size(); i++) {
				const Hit& h = oracleHits[i];
				if(h.h.first  == rhit.h.first  &&
				   h.h.second == rhit.h.second &&
				   h.fw       == rhit.fw)
				{
					// Oracle hit i and reported hit j refer to the
					// same alignment
					assert(h.mms == rhit.mms);
					assert_eq(stratum, oracleStrata[i]);
					// Erase the elements in the oracle vectors
					oracleHits.erase(oracleHits.begin() + i);
					oracleStrata.erase(oracleStrata.begin() + i);
					found = true;
					break;
				}
			}
			// If the backtracker found a hit, it had better be one
			// of the ones that the oracle found
			assert(found);
		}
		// All of the oracle hits that corresponded to a reported
		// hit have now been eliminated
		assert_neq(-1, lastStratum);
		assert_neq(-1, bestStratum);
		assert_neq(-1, worstStratum);
		assert_eq(oracleHits.size(), oracleStrata.size());
		if(maxHitsAllowed == 0xffffffff && !sink.exceededOverThresh()) {
			if(sink.spanStrata()) {
				// All hits remaining must occur more times than the max
				map<uint32_t,uint32_t> readToCnt;
				for(size_t i = 0; i < oracleHits.size(); i++) {
					readToCnt[oracleHits[i].patId]++;
				}
				map<uint32_t,uint32_t>::iterator it;
				for(it = readToCnt.begin(); it != readToCnt.end(); it++) {
					assert_gt(it->second, sink.overThresh());
				}
			} else {
				// Must have matched all oracle hits at the best
				// stratum
				size_t numLeftovers = 0;
				for(size_t i = 0; i < oracleStrata.size(); i++) {
					if(oracleStrata[i] == bestStratum) {
						numLeftovers++;
					}
				}
				// Must have matched every oracle hit
				if(numLeftovers > 0 && numLeftovers <= sink.overThresh()) {
					for(size_t i = 0; i < oracleHits.size(); i++) {
						if(oracleStrata[i] == bestStratum) {
							printHit(oracleHits[i]);
						}
					}
					assert(0);
				}
			}
		}
		if(sink.best() && sink.overThresh() == 0xffffffff && !sink.exceededOverThresh()) {
			// Ensure that all oracle hits have a stratum at least
			// as bad as the worst stratum observed in the reported
			// hits
			for(size_t i = 0; i < oracleStrata.size(); i++) {
				assert_leq(bestStratum, oracleStrata[i]);
			}
		}
	}

	/**
	 * Naively search for hits for the current pattern under the
	 * current backtracking strategy and store hits in hits vector.
	 */
	void naiveOracle(vector<Hit>& hits,
	                 vector<int>& strata,
	                 uint32_t iham = 0, /// initial weighted hamming distance
	                 uint32_t unrevOff = 0xffffffff,
	                 uint32_t oneRevOff = 0xffffffff,
	                 uint32_t twoRevOff = 0xffffffff,
	                 uint32_t threeRevOff = 0xffffffff)
	{
		if(unrevOff  == 0xffffffff) unrevOff  = _unrevOff;
		if(oneRevOff == 0xffffffff) oneRevOff = _1revOff;
		if(twoRevOff == 0xffffffff) twoRevOff = _2revOff;
		if(threeRevOff == 0xffffffff) threeRevOff = _3revOff;
		bool ebwtFw = _ebwt->fw();
		bool fw = _params.fw();
		uint32_t patid = _params.patId();

		::naiveOracle(
		            (*_os),
		            (*_qry),
		            _qlen,
		            (*_qual),
		            (*_name),
		            patid,
		            hits,
		            strata,
		            _qualThresh,
		            unrevOff,
		            oneRevOff,
		            twoRevOff,
		            threeRevOff,
		            fw,    // transpose
		            ebwtFw,
		            iham,
		            _muts,
		            _maqPenalty,
		            _halfAndHalf,
		            _reportExacts,
		            false);
	}

	bool                _started;
	String<Dna5>*       _qry;    // query (read) sequence
	size_t              _qlen;   // length of _qry
	String<char>*       _qual;   // quality values for _qry
	String<char>*       _name;   // name of _qry
	const Ebwt<String<Dna> >*   _ebwt;   // Ebwt to search in
	const EbwtSearchParams<String<Dna> >& _params;   // Ebwt to search in
	uint32_t            _unrevOff; // unrevisitable chunk
	uint32_t            _1revOff;  // 1-revisitable chunk
	uint32_t            _2revOff;  // 2-revisitable chunk
	uint32_t            _3revOff;  // 3-revisitable chunk
	/// Whether to round qualities off Maq-style when calculating penalties
	bool                _maqPenalty;
	uint32_t            _qualThresh; // only accept hits with weighted
	                             // hamming distance <= _qualThresh
	uint32_t           *_pairs;  // ranges, leveled in parallel
	                             // with decision stack
	uint8_t            *_elims;  // which ranges have been
	                             // eliminated, leveled in parallel
	                             // with decision stack
	uint8_t            *_bts;    // how many backtracks remain in each
	                             // equivalence class
	std::vector<uint32_t> _mms;  // array for holding mismatches
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
	/// Number of backtracks permitted w/ stackDepth 0 before giving up
	/// in halfAndHalf mode
	uint32_t            _maxBts0;
	/// Number of backtracks permitted w/ stackDepth 1 before giving up
	/// in halfAndHalf mode
	uint32_t            _maxBts1;
	/// Number of backtracks permitted w/ stackDepth 2 before giving up
	/// in halfAndHalf mode
	uint32_t            _maxBts2;
	/// Flag to record whether a 'false' return from backtracker is due
	/// to having exceeded one or more backrtacking limits
	bool                _bailedOnBacktracks;
	/// Source of pseudo-random numbers
	RandomSource        _rand;
	/// Seed for random number generator
	uint32_t            _randSeed;
	/// Be talkative
	bool                _verbose;
	uint64_t            _ihits;
	// Holding area for partial alignments
	vector<PartialAlignment> _partialsBuf;
	// Current range to expose to consumers
	Range               _curRange;
#ifndef NDEBUG
	std::set<int64_t>   allTops_;
#endif
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
class GreedyDFSRangeSourceDriver :
	public SingleRangeSourceDriver<GreedyDFSRangeSource>
{
public:
	GreedyDFSRangeSourceDriver(
			EbwtSearchParams<String<Dna> >& params,
			GreedyDFSRangeSource* rs,
			bool fw,
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
			uint32_t seed) :
			SingleRangeSourceDriver<GreedyDFSRangeSource>(
					params, rs, fw, sink, sinkPt, os, verbose, seed),
			rs_(rs), seedLen_(seedLen),
			nudgeLeft_(nudgeLeft),
			rev0Off_(rev0Off), rev1Off_(rev1Off),
			rev2Off_(rev2Off), rev3Off_(rev3Off)
	{ }

	virtual ~GreedyDFSRangeSourceDriver() { }

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
		exit(1);
	}

	/**
	 * Called every time setQuery() is called in the parent class,
	 * after setQuery() has been called on the RangeSource but before
	 * initConts() has been called.
	 */
	virtual void initRangeSource() {
		// If seedLen_ is huge, then it will always cover the whole
		// alignment
		uint32_t s = (seedLen_ > 0 ? min(seedLen_, len_) : len_);
		uint32_t sLeft  = s >> 1;
		uint32_t sRight = s >> 1;
		// If seed has odd length, then nudge appropriate half up by 1
		if((s & 1) != 0) { if(nudgeLeft_) sLeft++; else sRight++; }
		uint32_t rev0Off = cextToDepth(rev0Off_, sRight, s, len_);
		uint32_t rev1Off = cextToDepth(rev1Off_, sRight, s, len_);
		uint32_t rev2Off = cextToDepth(rev2Off_, sRight, s, len_);
		uint32_t rev3Off = cextToDepth(rev3Off_, sRight, s, len_);
		// If there are any Ns in the unrevisitable region, then this
		// driver is guaranteed to yield no fruit.
		uint16_t minCost = 0;
		if(rs_->reportExats()) {
			// Keep minCost at 0
		} else if (!rs_->halfAndHalf()) {
			// Exacts not allowed, so there must be at least 1 mismatch
			// outside of the unrevisitable area
			minCost = 1 << 14;
			uint8_t lowQual = 0xff;
			for(uint32_t d = rev0Off; d < s; d++) {
				if((*this->qual_)[len_ - d - 1] < lowQual) {
					lowQual = (*this->qual_)[len_ - d - 1];
				}
			}
			assert_lt(lowQual, 0xff);
			assert_geq(lowQual, 33);
			minCost += (lowQual - 33);
		} else {
			// Half-and-half constraints are active, so there must be
			// at least 1 mismatch in both halves of the seed
			assert(rs_->halfAndHalf());
			minCost = 2 << 14;
			uint8_t lowQual1 = 0xff;
			for(uint32_t d = 0; d < sRight; d++) {
				if((*this->qual_)[len_ - d - 1] < lowQual1) {
					lowQual1 = (*this->qual_)[len_ - d - 1];
				}
			}
			assert_lt(lowQual1, 0xff);
			assert_geq(lowQual1, 33);
			uint8_t lowQual2 = 0xff;
			for(uint32_t d = sRight; d < s; d++) {
				if((*this->qual_)[len_ - d - 1] < lowQual2) {
					lowQual2 = (*this->qual_)[len_ - d - 1];
				}
			}
			assert_lt(lowQual2, 0xff);
			assert_geq(lowQual2, 33);
			minCost += (lowQual1 + lowQual2 - 66);
		}
		this->minCostAdjustment_ = minCost;
		rs_->setOffs(sRight,   // depth of far edge of hi-half (only matters where half-and-half is possible)
		             s,        // depth of far edge of lo-half (only matters where half-and-half is possible)
		             rev0Off,  // depth above which we cannot backtrack
		             rev1Off,  // depth above which we may backtrack just once
		             rev2Off,  // depth above which we may backtrack just twice
		             rev3Off); // depth above which we may backtrack just three times
	}

	GreedyDFSRangeSource* rs_;
	uint32_t seedLen_;
	bool nudgeLeft_;
	SearchConstraintExtent rev0Off_;
	SearchConstraintExtent rev1Off_;
	SearchConstraintExtent rev2Off_;
	SearchConstraintExtent rev3Off_;
};

#endif /*EBWT_SEARCH_BACKTRACK_H_*/
