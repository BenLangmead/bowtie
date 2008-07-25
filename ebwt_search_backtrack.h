#ifndef EBWT_SEARCH_BACKTRACK_H_
#define EBWT_SEARCH_BACKTRACK_H_

#define DEFAULT_SPREAD 64

/// Encapsulates a change made to a query base, i.e. "the 3rd base from
/// the 5' end was changed from an A to a T".  Useful when using
/// for matching seeded by "seedlings".
struct QueryMutation {
	QueryMutation(uint8_t _pos, uint8_t _oldBase, uint8_t _newBase) :
		pos(_pos), oldBase(_oldBase), newBase(_newBase)
	{
		assert_neq(oldBase, newBase);
		assert_lt(oldBase, 4);
		assert_lt(newBase, 4);
	}
	uint8_t pos;
	uint8_t oldBase;
	uint8_t newBase;
};

/**
 * Class that coordinates quality- and quantity-aware backtracking over
 * some range of a read sequence.
 * 
 * The creator can configure the BacktrackManager to treat different
 * stretches of the read differently.
 */
template<typename TStr>
class BacktrackManager {
public:
	BacktrackManager(const Ebwt<TStr>& __ebwt,
	                 const EbwtSearchParams<TStr>& __params,
	                 uint32_t __unrevOff, // size of unrevisitable chunk
	                 uint32_t __1revOff,  // size of 1-revisitable chunk
	                 uint32_t __itop,
	                 uint32_t __ibot,
	                 uint32_t __qualThresh,
	                 uint32_t __qualWobble,
	                 bool __reportSeedlings = false,
	                 String<uint8_t>* __seedlings = NULL,
	                 String<QueryMutation>* __muts = NULL,
	                 bool __verbose = true,
	                 bool __oneHit = true,
	                 uint32_t seed = 0,
	                 TStr* __qry = NULL,
	                 String<char>* __qual = NULL,
	                 String<char>* __name = NULL) :
		_qry(__qry),
		_qlen(0),
		_qual(__qual),
		_name(__name),
		_ebwt(__ebwt),
		_params(__params),
		_unrevOff(__unrevOff),
		_1revOff(__1revOff),
		_itop(__itop),
		_ibot(__ibot),
		_spread(DEFAULT_SPREAD),
		_maxStackDepth(DEFAULT_SPREAD),
		_qualThresh(__qualThresh),
		_qualWobble(__qualWobble),
		_oneHit(__oneHit),
		_reportOnHit(true),
		_pairs(NULL),
		_elims(NULL),
		_mms(NULL),
		_chars(NULL),
		_reportSeedlings(__reportSeedlings),
		_seedlings(__seedlings),
		_muts(__muts),
		_nameDefault("default"),
		_rand(RandomSource(seed)),
		_verbose(__verbose)
	{
	    // For a 40-bp query range, the _pairs array occupies
	    // 40 * 40 * 8 * 4 = 51,200 bytes, and _elims
	    // occupy 40 * 40 = 1,600 bytes
	    assert_geq(__1revOff, __unrevOff);
	    fill(_qualDefault, DEFAULT_SPREAD, (char)(40+33));
 		if(_qry != NULL) {
 			_qlen = length(*_qry);
 			_spread = length(*_qry);
 			if(_qual == NULL || length(*_qual) == 0) {
 				_qual = &_qualDefault;
 			}
 			assert_geq(length(*_qual), _qlen);
 			for(size_t i = 0; i < length(*_qual); i++) {
 				assert_geq((*_qual)[i], 33);
 				assert_leq((*_qual)[i], 73);
 			}
 			assert_geq(length(*_qual), _qlen);
 			if(_name == NULL || length(*_name) == 0) {
 				_name = &_nameDefault;
 			}
 			assert_leq(_spread, DEFAULT_SPREAD);
 	 		_maxStackDepth = length(*_qry) - min(_unrevOff, length(*_qry)) + 2;
 	 		_pairs  = new uint32_t[DEFAULT_SPREAD*_maxStackDepth*8];
 	 		_elims  = new uint8_t [DEFAULT_SPREAD*_maxStackDepth];
 			if(_muts != NULL) {
 				applyMutations();
 			}
 		}
		if(_itop != 0 || _ibot != 0) {
			assert_lt(_itop, _ibot);
		}
		_mms = new uint32_t[DEFAULT_SPREAD];
		_chars = new char[DEFAULT_SPREAD];
	}

	~BacktrackManager() {
		if(_pairs != NULL) {
			delete[] _pairs; _pairs = NULL;
		}
		if(_elims != NULL) {
			delete[] _elims; _elims = NULL;
		}
		if(_mms != NULL) {
			delete[] _mms; _mms = NULL;
		}
		if(_chars != NULL) {
			delete[] _chars; _chars = NULL;
		}
	}

	#define PAIR_TOP(d, c)    (pairs[d*8 + c + 0])
	#define PAIR_BOT(d, c)    (pairs[d*8 + c + 4])
	#define PAIR_SPREAD(d, c) (PAIR_BOT(d, c) - PAIR_TOP(d, c))
	#define QUAL(k)           ((uint8_t)(*_qual)[k] >= 33 ? ((uint8_t)(*_qual)[k] - 33) : 0)
	
	void setQuery(TStr* __qry,
	              String<char>* __qual,
	              String<char>* __name,
	              String<QueryMutation>* __muts = NULL)
	{
		_qry = __qry;
		_qual = __qual;
		_name = __name;
		if(_muts != NULL) {
			undoMutations();
		}
		_muts = __muts;
		if(_muts != NULL) {
			applyMutations();
		}
		assert(_qry != NULL);
		// Reset _qlen
		_qlen = length(*_qry);
		_spread = _qlen;
		assert_leq(_spread, DEFAULT_SPREAD);
		if(_qual == NULL || length(*_qual) == 0) {
			_qual = &_qualDefault;
		}
		assert_geq(length(*_qual), _qlen);
		for(size_t i = 0; i < length(*_qual); i++) {
			assert_geq((*_qual)[i], 33);
			assert_leq((*_qual)[i], 73);
		}
		if(_name == NULL || length(*_name) == 0) {
			_name = &_nameDefault;
		}
 		_maxStackDepth = length(*_qry) - min(_unrevOff, length(*_qry)) + 2;
 		if(_pairs == NULL) {
 			_pairs = new uint32_t[DEFAULT_SPREAD*_maxStackDepth*8];
 		}
 		if(_elims == NULL) {
 			_elims = new uint8_t[DEFAULT_SPREAD*_maxStackDepth];
 		}
		if(_verbose) {
			String<char> qual = (*_qual);
			if(length(qual) > length(*_qry)) {
				resize(qual, length(*_qry));
			}
			cout << "setQuery(_qry=" << (*_qry) << ", _qual=" << qual << ")" << endl;
		}
	}
	
	void setMuts(String<QueryMutation>* __muts) {
		if(_muts != NULL) {
			// Undo previous mutations
			assert_gt(length(*_muts), 0);
			undoMutations();
		}
		_muts = __muts;
		if(_muts != NULL) {
			assert_gt(length(*_muts), 0);
			applyMutations();
		}
	}

	/**
	 * Set _qlen according to parameter, except don't let it fall below
	 * the length of the query.
	 */
	void setQlen(uint32_t qlen) {
		assert(_qry != NULL);
		_qlen = min(length(*_qry), qlen);
	}

	/**
	 * Initiate the recursive backtracking routine starting at the
	 * extreme right-hand side of the pattern.  Use the ftab to match
	 * the first several characters in one chomp, as long as doing so
	 * does not "jump over" any legal backtracking targets.
	 */
	bool backtrack(vector<TStr>* os = NULL, uint32_t ham = 0) {
		assert_gt(length(*_qry), 0);
		assert_leq(_qlen, length(*_qry));
		assert_geq(length(*_qual), length(*_qry));
		int ftabChars = _ebwt._eh._ftabChars;
		uint32_t m = min<uint32_t>(_unrevOff, _qlen);
		if(m >= (uint32_t)ftabChars) {
			// The ftab doesn't extend past the unrevisitable portion,
			// so we can go ahead and use it
			// Rightmost char gets least significant bit-pair
			uint32_t ftabOff = (*_qry)[_qlen - ftabChars];
			assert_lt(ftabOff, _ebwt._eh._ftabLen-1);
			for(int i = ftabChars - 1; i > 0; i--) {
				ftabOff <<= 2;
				assert_lt((uint32_t)(*_qry)[_qlen-i], 4);
				ftabOff |= (uint32_t)(*_qry)[_qlen-i];
				assert_lt(ftabOff, _ebwt._eh._ftabLen-1);
			}
			assert_lt(ftabOff, _ebwt._eh._ftabLen-1);
			uint32_t top = _ebwt.ftabHi(ftabOff);
			uint32_t bot = _ebwt.ftabLo(ftabOff+1);
			if(_qlen == (uint32_t)ftabChars && bot > top) {
				// We have a match!
				if(_reportSeedlings) {
					// Oops - we're trying to find seedlings, so we've
					// gone too far; start again
					return backtrack(0,   // depth
					                 0,   // top
					                 0,   // bot
					                 os,
					                 ham);
				} else {
					// We have a match!
					return report(0, top, bot);
				}
			} else if (bot > top) {
				// We have an arrow pair from which we can backtrack
				return backtrack(ftabChars, // depth
				                 top,       // top
				                 bot,       // bot
				                 os,
				                 ham);
			}
			// The arrows are already closed; give up
			return false;
		} else {
			// The ftab *does* extend past the unrevisitable portion;
			// we can't use it in this case, because we might jump past
			// a legitimate mismatch
			return backtrack(0,   // depth
			                 0,   // top
			                 0,   // bot
			                 os,
			                 ham);
		}
	}

	/**
	 * Run the backtracker with initial conditions in place
	 */
	bool backtrack(uint32_t depth,
	               uint32_t top,
	               uint32_t bot,
	               vector<TStr>* os = NULL,
	               uint32_t iham = 0)
	{
		assert_gt(length(*_qry), 0);
		assert_leq(_qlen, length(*_qry));
		assert_geq(length(*_qual), length(*_qry));
		if(_verbose) cout << "backtrack(top=" << top << ", "
		                 << "bot=" << bot << ", "
		                 << "iham=" << iham << ", "
		                 << "_pairs" << _pairs << ", "
		                 << "_elims=" << (void*)_elims << ")" << endl;
		bool oldRetain = false;
		size_t oldRetainSz = 0;
		if(os != NULL && (*os).size() > 0) {
			oldRetain = _params.sink().retainHits();
			oldRetainSz = _params.sink().retainedHits().size();
			_params.sink().setRetainHits(true);
		}
		ASSERT_ONLY(uint64_t nhits = _params.sink().numHits());
		bool ret = backtrack(0, depth, _unrevOff, top, bot, iham, _pairs, _elims);
		if(ret) {
			assert_gt(_params.sink().numHits(), nhits);
		} else {
			assert_eq(_params.sink().numHits(), nhits);
		}
		// If these are end-to-end matches and we have the original
		// texts, then we can double-check them with a naive oracle.
		if(os != NULL && (*os).size() > 0 && !_reportSeedlings) {
			_params.sink().setRetainHits(oldRetain);
			vector<Hit> oracleHits;
			naiveOracle(*os, oracleHits);
			vector<Hit>& retainedHits = _params.sink().retainedHits();
			if(ret == false) {
				assert_eq(oldRetainSz, retainedHits.size());
				if(oracleHits.size() > 0) {
					const Hit& h = oracleHits[0];
					cout << "Oracle hit " << oracleHits.size()
					     << " times, but backtracker did not hit" << endl;
					cout << "First oracle hit: " << endl;
					cout << "  Pat:  " << (*_qry) << endl;
					cout << "  Tseg: ";
					bool ebwtFw = _params.ebwtFw();
					if(ebwtFw) {
						for(size_t i = 0; i < _qlen; i++) {
							cout << (*os)[h.h.first][h.h.second + i];
						}
					} else {
						for(int i = (int)_qlen-1; i >= 0; i--) {
							cout << (*os)[h.h.first][h.h.second + i];
						}
					}
					cout << endl;
					cout << "  Bt:   ";
					for(int i = (int)_qlen-1; i >= 0; i--) {
						if(i < (int)_unrevOff) cout << "0";
						else if(i < (int)_1revOff) cout << "1";
						else cout << "X";
					}
					cout << endl;
				}
				assert_eq(0, oracleHits.size());
			} else {
				assert_gt(oracleHits.size(), 0);
				assert_eq(oldRetainSz+1, retainedHits.size());
				// Get the most-recently-retained hit
				Hit& rhit = retainedHits.back();
				// Go through oracleHits and look for one that matches
				// up with rhit
				size_t i;
				for(i = 0; i < oracleHits.size(); i++) {
					const Hit& h = oracleHits[i];
		    		if(h.h.first == rhit.h.first && h.h.second == rhit.h.second) {
		    			// Assert that number of mismatches matches
		    			assert_eq(h.fw, rhit.fw);
		    			assert_eq(h.mms, rhit.mms);
		    			break;
		    		}
				}
				assert_lt(i, oracleHits.size()); // assert we found a matchup
			}
		}
		return ret;
	}

	/**
	 * Recursive routine for progressing to the next backtracking
	 * decision given some initial conditions.  If a hit is found, it
	 * is recorded and true is returned.  Otherwise, if there are more
	 * backtracking opportunities, the function will call itself
	 * recursively and return the result.  As soon as there is a
	 * mismatch and no backtracking opporunities, false is returned.
	 */
	bool backtrack(uint32_t  stackDepth,
	               uint32_t  depth,    // next depth where a post-pair needs to be calculated
	               uint32_t  unrevOff, // depths < unrevOff are unrevisitable 
	               uint32_t  top,      // top arrow in pair prior to 'depth'
	               uint32_t  bot,      // bottom arrow in pair prior to 'depth'
	               uint32_t  ham,      // weighted hamming distance so far
	               uint32_t* pairs,    // portion of pairs array to be used for this backtrack frame
	               uint8_t*  elims)    // portion of elims array to be used for this backtrack frame
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
		assert_leq(stackDepth, _maxStackDepth);
		if(_verbose) {
			cout << "  backtrack(stackDepth=" << stackDepth << ", "
                 << "depth=" << depth << ", "
                 << "top=" << top << ", "
                 << "bot=" << bot << ", "
                 << "ham=" << ham << ", "
                 << "pairs=" << pairs << ", "
                 << "elims=" << (void*)elims << "): \"";
			for(int i = (int)depth - 1; i >= 0; i--) {
				cout << _chars[i];
			}
			cout << "\"" << endl;
		}
		
		// The total number of arrow pairs that are acceptable
		// backtracking targets ("alternative" arrow pairs)
		uint32_t altNum = 0;
		// The total number of arrow pairs that are candidates to be
		// the *next* backtracking target because they are low quality
		// ("eligible" arrow pairs)
		uint32_t eligibleNum = 0;
		// Total distance between all lowest-quality "alternative"
		// arrow pairs that haven't yet been eliminated
		uint32_t eligibleSz = 0;
		// The lowest quality value associated with any alternative
		// arrow pairs; all alternative pairs with this quality are
		// eligible
		uint8_t lowAltQual = 0xff;
		uint32_t d = depth;
		uint32_t cur = _qlen - d - 1; // current offset into _qry
		SideLocus ltop, lbot;
		if(top != 0 || bot != 0) {
			SideLocus::initFromTopBot(top, bot, _ebwt._eh, _ebwt._ebwt, ltop, lbot);
		}
		// Advance along the read until we hit a mismatch
		while(cur < _qlen) {
			if(_verbose) {
				cout << "    cur=" << cur << " \"";
				for(int i = (int)d - 1; i >= 0; i--) {
					cout << _chars[i];
				}
				cout << "\"";
			}
			assert_lt((int)(*_qry)[cur], 4);
			bool curIsEligible = false;
			// Reset eligibleNum and eligibleSz if there are any
			// eligible pairs discovered at this spot
			bool curOverridesEligible = false;
			// Determine whether arrow pairs at this location are
			// candidates for backtracking
			int c = (int)(*_qry)[cur];
			uint8_t q = QUAL(cur);
			assert_lt((uint32_t)q, 100);
			bool curIsAlternative = (d >= unrevOff) && (ham + q <= _qualThresh);
			if(curIsAlternative) {
				// q is low enough to make this position an alternative,
				// but is it the best alternative?
				if(q < lowAltQual) {
					// Arrow pairs at this depth in this backtracking frame
					// are eligible, unless we learn otherwise.  Arrow
					// pairs previously thought to be eligible are not any
					// longer.
					curIsEligible = true;
					curOverridesEligible = true;
				} else if(q == lowAltQual) {
					// Arrow pairs at this depth in this backtracking frame
					// are eligible, unless we learn otherwise
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
			if(top == 0 && bot == 0) {
				// Calculate first quartet of pairs using the _fchr[]
				// array
				assert_eq(0, d);
				PAIR_TOP(0, 0)                  = _ebwt._fchr[0];
				PAIR_BOT(0, 0) = PAIR_TOP(0, 1) = _ebwt._fchr[1];
				PAIR_BOT(0, 1) = PAIR_TOP(0, 2) = _ebwt._fchr[2];
				PAIR_BOT(0, 2) = PAIR_TOP(0, 3) = _ebwt._fchr[3];
				PAIR_BOT(0, 3)                  = _ebwt._fchr[4];
				// Update top and bot
				top = PAIR_TOP(d, c); bot = PAIR_BOT(d, c);
			} else if(curIsAlternative) {
				// Clear pairs
				bzero(&pairs[d*8], 8 * 4);
				// Calculate next quartet of pairs
				_ebwt.mapLFEx(ltop, lbot, &pairs[d*8], &pairs[(d*8)+4]);
				// Update top and bot
				top = PAIR_TOP(d, c); bot = PAIR_BOT(d, c);
			} else {
				// This query character is not even a legitimate
				// alternative (because backtracking here would blow
				// our mismatch quality budget), so no need to do the
				// bookkeeping for the entire quartet, just do c
				top = _ebwt.mapLF(ltop, c); bot = _ebwt.mapLF(lbot, c);
			}
			if(top != bot) {
				// Calculate loci from row indices; do it now so that
				// those prefetches are fired off as soon as possible
				SideLocus::initFromTopBot(top, bot, _ebwt._eh, _ebwt._ebwt, ltop, lbot);
			}
			// Update the elim array
			elims[d] = (1 << c);
			assert_lt(elims[d], 16);
			assert_gt(elims[d], 0);
			if(curIsAlternative) {
				// Given the just-calculated quartet of arrow pairs, update
				// elims, altNum, eligibleNum, eligibleSz
				for(int i = 0; i < 4; i++) {
					assert_leq(PAIR_TOP(d, i), PAIR_BOT(d, i));
					uint32_t spread = PAIR_SPREAD(d, i);
					if(spread == 0) {
						// Indicate this char at this position is
						// eliminated as far as this backtracking frame is
						// concerned, since its arrow pair is closed
						elims[d] |= (1 << i);
					}
					if(i != c && spread > 0 && ((elims[d] & (1 << i)) == 0)) {
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
			// Mismatch with alternatives
			while((top == bot && altNum > 0) ||
			      (stackDepth == 0 && cur == 0 && _reportSeedlings && altNum > 0))
			{
				if(_verbose) cout << "    top (" << top << ") == bot ("
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
				uint32_t r = _rand.nextU32() % eligibleSz;
				bool foundTarget = false;
				uint32_t cumSz = 0;
				ASSERT_ONLY(uint32_t eligiblesVisited = 0);
				size_t i = depth, j = 0;
				uint32_t bttop = 0;
				uint32_t btbot = 0;
				uint32_t btham = ham;
				char btchar = 0;
				uint32_t icur = 0;
				for(; i <= d; i++) {
					if(i < unrevOff) {
						// i is an unrevisitable position, so don't consider it
						continue;
					}
					icur = _qlen - i - 1; // current offset into _qry
					uint8_t qi = QUAL(icur);
					assert_lt(elims[i], 16);
					assert_gt(elims[i], 0);
					if(qi == lowAltQual && elims[i] != 15) {
						// This is an eligible position with at least
						// one remaining backtrack target
						for(j = 0; j < 4; j++) {
							if((elims[i] & (1 << j)) == 0) {
								// This pair has not been eliminated
								assert_gt(PAIR_BOT(i, j), PAIR_TOP(i, j));
								cumSz += PAIR_SPREAD(i, j);
								ASSERT_ONLY(eligiblesVisited++);
								if(r < cumSz) {
									// This is our randomly-selected 
									// backtrack target
									foundTarget = true;
									bttop = PAIR_TOP(i, j);
									btbot = PAIR_BOT(i, j);
									btham += qi;
									btchar = "acgt"[j];
									assert_leq(btham, _qualThresh);
									break;
								}
							}
						}
						if(foundTarget) break;
					}
				}
				assert_leq(eligiblesVisited, eligibleNum);
				assert_leq(i, d);
				assert_lt(j, 4);
				assert_neq(0, btchar);
				assert_leq(cumSz, eligibleSz);
				assert(foundTarget);
				assert_gt(btbot, bttop);
				assert_leq(btbot-bttop, eligibleSz);
				// Slide over to the next backtacking frame within
				// pairs and elims; won't interfere with our frame or
				// any of our parents' frames
				uint32_t *newPairs = pairs + (_spread*8);
				uint8_t  *newElims = elims + (_spread);
				// If we've selected a backtracking target that's in
				// the 1-revisitable region, then we ask the recursive
				// callee to consider the 1-revisitable region as also
				// being unrevisitable (since we just "used up" all of
				// our visits)
				uint32_t btUnrevOff = unrevOff;
				if(i < _1revOff) {
					assert_geq(_1revOff, unrevOff);
					btUnrevOff = _1revOff;
				}
				// Note the character that we're backtracking on in the
				// mm array:
				_mms[stackDepth] = icur;
				_chars[i] = btchar;
				// Now backtrack to target
				ASSERT_ONLY(uint64_t numHits = _params.sink().numHits());
				assert_leq(i+1, _qlen);
				bool ret;
				if(i+1 == _qlen) {
					ret = report(stackDepth+1, bttop, btbot);
				} else {
					ret = backtrack(stackDepth+1,
				                    i+1,
				                    btUnrevOff,
				                    bttop,  // top arrow in pair prior to 'depth'
				                    btbot,  // bottom arrow in pair prior to 'depth'
				                    btham,  // weighted hamming distance so far
				                    newPairs,
				                    newElims);
				}
				if(ret) {
					assert_gt(_params.sink().numHits(), numHits);
					return true; // return, signaling that we've reported
				}
				assert_eq(_params.sink().numHits(), numHits);
				// No hit was reported; update elims[], eligibleSz,
				// eligibleNum, altNum
				_chars[i] = (*_qry)[icur];
				assert_neq(15, elims[i]);
				ASSERT_ONLY(uint8_t oldElim = elims[i]);
				elims[i] |= (1 << j);
				assert_gt(elims[i], oldElim);
				eligibleSz -= (btbot-bttop);
				eligibleNum--;
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
				else if(eligibleNum == 0) {
					// Find the next set of eligible backtrack points
					// by re-scanning this backtracking frame (from
					// 'depth' up to 'd')
					lowAltQual = 0xff;
					for(size_t k = depth; k <= d; k++) {
						uint32_t kcur = _qlen - k - 1; // current offset into _qry
						uint8_t kq = QUAL(kcur);
						bool kCurIsAlternative = (k >= unrevOff) && (ham + kq <= _qualThresh);
						bool kCurOverridesEligible = false;
						if(kCurIsAlternative) {
							if(kq < lowAltQual) {
								kCurOverridesEligible = true;
							}
							if(kq <= lowAltQual) {
								// Position is eligible
								for(int l = 0; l < 4; l++) {
									if((elims[k] & (1 << l)) == 0) {
										// Not yet eliminated
										if(kCurOverridesEligible) {
											// Clear previous eligible results;
											// this one's better
											lowAltQual = kq;
											kCurOverridesEligible = false;
											eligibleNum = 0;
											eligibleSz = 0;
										}
										eligibleNum++;
										uint32_t spread = PAIR_SPREAD(k, l);
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
			// Mismatch with no alternatives
			if(top == bot && altNum == 0) {
				assert_eq(0, altNum);
				assert_eq(0, eligibleSz);
				assert_eq(0, eligibleNum);
				// Mismatched with no backtracking opportunities;
				// return failure
				return false;
			}
			// Match!
			_chars[d] = (*_qry)[cur];
			d++; cur--;
		}
		assert_eq(0xffffffff, cur);
		assert_gt(bot, top);
		if(!_reportSeedlings || stackDepth > 0) {
			return report(stackDepth, top, bot);
		} else {
			return false;
		}
	}

protected:
	
	void applyMutations() {
		if(_muts == NULL) {
			// No mutations to apply
			return;
		}
		for(size_t i = 0; i < length(*_muts); i++) {
			const QueryMutation& m = (*_muts)[i];
			assert_lt(m.pos, _qlen);
			assert_lt(m.oldBase, 4);
			assert_lt(m.newBase, 4);
			assert_neq(m.oldBase, m.newBase);
			assert_eq((uint32_t)((*_qry)[m.pos]), (uint32_t)m.oldBase);
			(*_qry)[m.pos] = (Dna)(int)m.newBase; // apply it
		}
	}
	
	void undoMutations() {
		if(_muts == NULL) {
			// No mutations to undo
			return;
		}
		for(size_t i = 0; i < length(*_muts); i++) {
			const QueryMutation& m = (*_muts)[i];
			assert_lt(m.pos, _qlen);
			assert_lt(m.oldBase, 4);
			assert_lt(m.newBase, 4);
			assert_neq(m.oldBase, m.newBase);
			assert_eq((uint32_t)((*_qry)[m.pos]), (uint32_t)m.newBase);
			(*_qry)[m.pos] = (Dna)(int)m.oldBase; // undo it
		}
	}
	
	bool report(uint32_t stackDepth, uint32_t top, uint32_t bot) {
		if(_reportSeedlings) {
			reportSeedling(stackDepth);
			return false; // keep going
		} else {
			// Undo all the mutations
			ASSERT_ONLY(TStr tmp = (*_qry));
			undoMutations();
			bool hit;
			if(_muts != NULL) {
				assert_neq(tmp, (*_qry));
				for(size_t i = 0; i < length(*_muts); i++) {
					// Entries in _mms[] are in terms of offset into
					// _qry - not in terms of offset from 3' or 5' end
					_mms[stackDepth] = (*_muts)[i].pos;
				}
				hit = reportHit(stackDepth+1, top, bot);
			} else {
				hit = reportHit(stackDepth, top, bot);
			}
			applyMutations();
			assert_eq(tmp, (*_qry));
			return hit;
		}
	}

	/**
	 * Report a hit with # mismatches = stackDepth, at rows delimited
	 * by top and bot
	 */
	bool reportHit(uint32_t stackDepth, uint32_t top, uint32_t bot) {
		// Possibly report
		if(_oneHit) {
			uint32_t spread = bot - top;
			uint32_t r = top + (_rand.nextU32() % spread);
			for(uint32_t i = 0; i < spread; i++) {
				uint32_t ri = r + i;
				if(ri >= bot) ri -= spread;
				// reportChaseOne takes ths _mms[] list in terms of
				// their indices into the query string; not in terms
				// of their offset from the 3' or 5' end.
				if(_ebwt.reportChaseOne((*_qry), _qual, _name,
				                        _mms, stackDepth, ri,
				                        top, bot, _qlen, _params))
				{
					return true;
				}
			}
			return false;
		} else {
			// Not yet smart enough to report all hits
			assert(false);
			return false;
		}
	}

	/**
	 * Report a "seedling hit" - i.e. report the mismatch that got us
	 * here.  We only know how to deal with the mismatch in _mms[0].
	 */
	bool reportSeedling(uint32_t stackDepth) {
		// Possibly report
		assert(_reportSeedlings);
		assert(_seedlings != NULL);
		// Right now we only know how to report single-mismatch seedlings
		assert_eq(1, stackDepth);
		// Enrties in _mms[] hold the offset from the 5' end 
		assert_lt(_mms[0], _qlen);
		append((*_seedlings), (uint8_t)_mms[0]); // pos
		uint32_t i = _qlen - _mms[0] - 1;
		// _chars[] is index in terms of RHS-relative depth
		int c = (int)(Dna)_chars[i];
		assert_lt(c, 4);
		assert_neq(c, (int)(*_qry)[_mms[0]]);
		append((*_seedlings), (uint8_t)c); // chr
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
			uint8_t qi = QUAL(icur);
			assert_lt(elims[i], 16);
			assert_gt(elims[i], 0);
			if(qi == lowAltQual && elims[i] != 15) {
				// This is an eligible position with at least
				// one remaining backtrack target
				for(j = 0; j < 4; j++) {
					if((elims[i] & (1 << j)) == 0) {
						// This pair has not been eliminated
						assert_gt(PAIR_BOT(i, j), PAIR_TOP(i, j));
						cumSz += PAIR_SPREAD(i, j);
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
	 * Naively search for the same hits that should be found by 
	 */
	void naiveOracle(vector<TStr>& os, vector<Hit>& hits) {
		typedef typename Value<TStr>::Type TVal;
		bool ebwtFw = _params.ebwtFw();
		bool fw = _params.fw();
		bool fivePrimeOnLeft = (ebwtFw == fw);
		uint32_t patid = _params.patId();
	    uint32_t plen = length(*_qry);
		uint8_t *pstr = (uint8_t *)begin(*_qry, Standard());
	    // For each text...
		for(size_t i = 0; i < os.size(); i++) {
			// For each text position...
			if(length(os[i]) < plen) continue;
			TStr o = os[i];
			uint32_t olen = length(o);
			if(!ebwtFw) {
				for(size_t j = 0; j < olen>>1; j++) {
					TVal tmp = o[j];
					o[j] = o[olen-j-1];
					o[olen-j-1] = tmp;
				}
			}
			uint8_t *ostr = (uint8_t *)begin(o, Standard());
			// For each possible alignment of pattern against text
			for(size_t j = 0; j <= olen - plen; j++) {
				size_t rev1mm  = 0; // mismatches observed in the 1-revisitable region
				uint32_t ham = 0; // weighted hamming distance so far
				bitset<max_read_bp> diffs = 0; // mismatch bitvector
				// For each alignment column, from right to left
				bool success = true;
				for(int k = (int)plen-1; k >= 0; k--) {
					size_t kr = plen-1-k;
					if(pstr[k] != ostr[j+k]) {
						ham += QUAL(k);
						if(ham > _qualThresh) {
							// Alignment is invalid because it exceeds
							// our target weighted hamming distance
							// threshold
							success = false;
							break;
						}
						// What region does the mm fall into?
						if(kr < _unrevOff) {
							// Alignment is invalid because it contains
							// a mismatch in the unrevisitable region
							success = false;
							break;
						} else if(kr < _1revOff) {
							rev1mm++;
							if(rev1mm > 1) {
								// Alignment is invalid because it
								// contains more than 1 mismatch in the
								// 1-revisitable region
								success = false;
								break;
							}
						}
						if(fivePrimeOnLeft) {
							diffs.set(k);
						} else {
							// The 3' end is on on the left end of the
							// pattern, but the diffs vector should
							// encode mismatches w/r/t the 5' end, so
							// we flip
							diffs.set(plen-k-1);
						}
					}
				}
				if(success) {
					// It's a hit
					uint32_t off = j;
					if(!ebwtFw) {
						off = olen - off;
						off -= plen;
					}
					Hit h(make_pair(i, off), 
						  patid,  // read id
						  *_name, // read name
						  *_qry,  // read sequence
						  *_qual, // read qualities 
						  fw,     // forward/reverse-comp
						  diffs); // mismatch bitvector
					hits.push_back(h);
				} // For each pattern character
			} // For each alignment over current text
		} // For each text
	}
	
	bool sanityCheckPairs(uint32_t* pairs) {
		for(size_t i = 0; i < _spread; i++) {
			if(pairs[i*2] == 0 && pairs[i*2+1] == 0) {
				continue;
			}
			assert_lt(pairs[i*2], pairs[i*2+1]);
		}
		return true;
	}
	
	TStr*               _qry;    // query (read) sequence
	size_t              _qlen;   // length of _qry
	String<char>*       _qual;   // quality values for _qry
	String<char>*       _name;   // name of _qry
	const Ebwt<TStr>&   _ebwt;   // Ebwt to search in
	const EbwtSearchParams<TStr>& _params;   // Ebwt to search in
	uint32_t            _unrevOff; // unrevisitable chunk
	uint32_t            _1revOff;  // 1-revisitable chunk
	uint32_t            _itop;   // initial top arrow, or 0xffffffff if
	                             // we're starting from the beginning
	                             // of the query
	uint32_t            _ibot;   // initial bot arrow, or 0xffffffff if
	                             // we're starting from the beginning
	                             // of the query
	uint32_t            _spread; // size of window within which to
	                             // backtrack
	uint32_t            _maxStackDepth;
	uint32_t            _qualThresh; // only accept hits with weighted
	                             // hamming distance <= _qualThresh
	uint32_t            _qualWobble; // hits within _qualWobble weighted
	                             // hamming distance of each other are
	                             // considered equally good
	bool                _oneHit; // stop backtracking after finding 1
	                             // legit hit?  (doesn't really work -
	                             // we always operate in _oneHit mode)
	bool                _reportOnHit; // report as soon as we find a
	                             // hit? (as opposed to leaving it up
	                             // to the caller whether to report)
	uint32_t           *_pairs;  // arrow pairs, leveled in parallel
	                             // with decision stack
	uint8_t            *_elims;  // which arrow pairs have been
	                             // eliminated, leveled in parallel
	                             // with decision stack
	uint32_t           *_mms;    // array for holding mismatches
	// Entries in _mms[] are in terms of offset into
	// _qry - not in terms of offset from 3' or 5' end
	char               *_chars;  // characters selected so far
	bool                _reportSeedlings;
	String<uint8_t>    *_seedlings; // list in which to store seedlings
	String<QueryMutation> *_muts;
	String<char>        _nameDefault; // default name, for when it's
	                             // not specified by caller
	String<char>        _qualDefault; // default quals
	RandomSource        _rand;   // Source of pseudo-random numbers
	bool                _verbose;// be talkative
};

#endif /*EBWT_SEARCH_BACKTRACK_H_*/
