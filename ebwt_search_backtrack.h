#ifndef EBWT_SEARCH_BACKTRACK_H_
#define EBWT_SEARCH_BACKTRACK_H_

#define DEFAULT_SPREAD 64

/**
 * Class that coordinates quality- and quantity-aware backtracking over
 * some range of a read sequence.
 */
template<typename TStr>
class BacktrackManager {
public:
	BacktrackManager(const Ebwt<TStr>& __ebwt,
	                 const EbwtSearchParams<TStr>& __params,
	                 uint32_t __off,
	                 uint32_t __itop,
	                 uint32_t __ibot,
	                 uint32_t __iham,
	                 uint32_t __qualThresh,
	                 uint32_t __qualWobble,
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
		_off(__off),
		_itop(__itop),
		_ibot(__ibot),
		_iham(__iham),
		_spread(40),
		_qualThresh(__qualThresh),
		_qualWobble(__qualWobble),
		_oneHit(__oneHit),
		_reportOnHit(true),
		_pairs(NULL),
		_elims(NULL),
		_mms(NULL),
	    _nameDefault("default"),
		_rand(RandomSource(seed)),
		_verbose(__verbose)
	{
	    // For a 40-bp query range, the _pairs array occupies
	    // 40 * 40 * 8 * 4 = 51,200 bytes, and _elims
	    // occupy 40 * 40 = 1,600 bytes
	    fill(_qualDefault, DEFAULT_SPREAD, (char)(40+33));
 		if(_qry != NULL) {
 			_qlen = length(*_qry);
 			_spread = length(*_qry) - _off;
 			if(_qual == NULL || length(*_qual) == 0) {
 				_qual = &_qualDefault;
 			}
 			assert_geq(length(*_qual), _qlen);
 			if(_name == NULL || length(*_name) == 0) {
 				_name = &_nameDefault;
 			}
 			assert_leq(_spread, DEFAULT_SPREAD);
 		}
		_pairs  = new uint32_t[DEFAULT_SPREAD*DEFAULT_SPREAD*8];
		_elims  = new uint8_t [DEFAULT_SPREAD*DEFAULT_SPREAD];
		if(_itop != 0 || _ibot != 0) {
			assert_lt(_itop, _ibot);
		}
		_mms = new uint32_t[DEFAULT_SPREAD];
		bzero(_mms, sizeof(uint32_t) * DEFAULT_SPREAD);
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
	}

	#define PAIR_TOP(d, c)    (pairs[d*8 + c + 0])
	#define PAIR_BOT(d, c)    (pairs[d*8 + c + 4])
	#define PAIR_SPREAD(d, c) (PAIR_BOT(d, c) - PAIR_TOP(d, c))
	#define QUAL(k)           ((uint8_t)(*_qual)[k] - 33)
	
	void setQuery(TStr* __qry,
	              String<char>* __qual,
	              String<char>* __name)
	{
		_qry = __qry;
		_qual = __qual;
		_name = __name;
		assert(_qry != NULL);
		// Reset mismatch array
		bzero(_mms, sizeof(uint32_t) * DEFAULT_SPREAD);
		// Reset _qlen
		_qlen = length(*_qry);
		_spread = length(*_qry) - _off;
		assert_leq(_spread, DEFAULT_SPREAD);
		if(_qual == NULL || length(*_qual) == 0) {
			_qual = &_qualDefault;
		}
		assert_geq(length(*_qual), _qlen);
		if(_name == NULL || length(*_name) == 0) {
			_name = &_nameDefault;
		}
	}

	bool backtrack(vector<TStr>* os = NULL)	{
		return backtrack(_itop, _ibot, os);
	}

	/**
	 * Run the backtracker with initial conditions in place
	 */
	bool backtrack(uint32_t top,
	               uint32_t bot,
	               vector<TStr>* os = NULL)
	{
		if(_verbose) cout << "backtrack(_itop=" << _itop << ", "
		                 << "_ibot=" << _ibot << ", "
		                 << "_iham=" << _iham << ", "
		                 << "_pairs" << _pairs << ", "
		                 << "_elims=" << (void*)_elims << ")" << endl;
		bool oldRetain = false;
		size_t oldRetainSz = 0;
		if(os != NULL && (*os).size() > 0 && top == 0 && bot == 0) {
			oldRetain = _params.sink().retainHits();
			oldRetainSz = _params.sink().retainedHits().size();
			_params.sink().setRetainHits(true);
		}
		bool ret = backtrack(0, 0, top, bot, _iham, _pairs, _elims);
		// If these are end-to-end matches and we have the original
		// texts, then we can double-check them with a naive oracle.
		if(os != NULL && (*os).size() > 0 && top == 0 && bot == 0) {
			_params.sink().setRetainHits(oldRetain);
			vector<Hit> oracleHits;
			naiveOracle(*os, oracleHits);
			if(ret == false) {
				assert_eq(0, oracleHits.size());
			} else {
				assert_gt(oracleHits.size(), 0);
				vector<Hit>& retainedHits = _params.sink().retainedHits();
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
	 * Naively search for the same hits that should be found by 
	 */
	void naiveOracle(vector<TStr>& os, vector<Hit>& hits) {
		typedef typename Value<TStr>::Type TVal;
		bool ebwtFw = _params.ebwtFw();
		bool fw = _params.fw();
		bool fivePrimeOnLeft = (ebwtFw == fw);
		uint32_t patid = _params.patId();
	    uint32_t plen = length(_qry);
		uint8_t *pstr = (uint8_t *)begin(_qry, Standard());
	    // For each text...
		for(size_t i = 0; i < os.size(); i++) {
			// For each text position...
			TStr o = os[i];
			uint32_t olen = length(o);
			if(!ebwtFw) {
				for(size_t j = 0; j < olen>>1; j++) {
					TVal tmp = o[j];
					o[j] = o[olen-j-1];
					o[olen-j-1] = tmp;
				}
			}
			uint8_t *ostr = (uint8_t *)begin(os[i], Standard());
			for(size_t j = 0; j <= olen - plen; j++) {
				uint32_t ham = 0; // weighted hamming distance so far
				bitset<max_read_bp> diffs = 0; // mismatch bitvector
				for(size_t k = 0; k < plen; k++) {
					if(pstr[k] != ostr[j+k]) {
						ham += QUAL(k);
						if(ham > _qualThresh) break;
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
				if(ham <= _qualThresh) {
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
	
	/**
	 * Recursive routine for progressing to the next backtracking
	 * decision given some initial conditions.  If a hit is found, it
	 * is recorded and true is returned.  Otherwise, if there are more
	 * backtracking opportunities, the function will call itself
	 * recursively and return the result.  As soon as there is a
	 * mismatch and no backtracking opporunities, false is returned.
	 */
	bool backtrack(uint32_t  stackDepth,
	               uint32_t  depth,  // next depth where a post-pair needs to be calculated
	               uint32_t  top,    // top arrow in pair prior to 'depth'
	               uint32_t  bot,    // bottom arrow in pair prior to 'depth'
	               uint32_t  ham,    // weighted hamming distance so far
	               uint32_t* pairs,  // portion of pairs array to be used for this backtrack frame
	               uint8_t*  elims)
	{
		// Can't have already exceeded weighted hamming distance threshold
		assert_leq(stackDepth, depth);
		assert(_qry != NULL);
		assert(_qual != NULL);
		assert(_name != NULL);
		assert(_qlen != 0);
		assert_leq(ham, _qualThresh);
		assert_lt(depth, _qlen); // can't have run off the end of qry
		assert_geq(bot, top);    // could be that both are 0
		assert(pairs != NULL);
		assert(elims != NULL);
		assert_lt(_off, _qlen);
		if(_verbose) cout << "  backtrack(stackDepth=" << stackDepth << ", "
		                 << "depth=" << depth << ", "
		                 << "top=" << top << ", "
		                 << "bot=" << bot << ", "
		                 << "ham=" << ham << ", "
		                 << "pairs=" << pairs << ", "
		                 << "elims=" << (void*)elims << ")" << endl;
		
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
		uint32_t cur = _qlen - _off - d - 1; // current offset into _qry
		SideLocus ltop, lbot;
		if(top != 0 || bot != 0) {
			SideLocus::initFromTopBot(top, bot, _ebwt._eh, _ebwt._ebwt, ltop, lbot);
		}
		// Advance along the read until we hit a mismatch
		while(cur < _qlen) {
			if(_verbose) cout << "    cur=" << cur << endl;
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
			bool curIsAlternative = (ham + q <= _qualThresh);
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
					if(curIsAlternative) {
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
			assert(sanityCheckEligibility(depth, d, lowAltQual, eligibleSz, eligibleNum, pairs, elims));
			// Mismatch with alternatives
			while(top == bot && altNum > 0) {
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
				assert(sanityCheckEligibility(depth, d, lowAltQual, eligibleSz, eligibleNum, pairs, elims));
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
				for(; i <= d; i++) {
					uint32_t icur = _qlen - _off - i - 1; // current offset into _qry
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
				assert_leq(cumSz, eligibleSz);
				assert(foundTarget);
				assert_gt(btbot, bttop);
				assert_leq(btbot-bttop, eligibleSz);
				// Slide over to the next backtacking frame within
				// pairs and elims; won't interfere with our frame or
				// any of our parents' frames
				uint32_t *newPairs = pairs + (_spread*8);
				uint8_t  *newElims = elims + (_spread);
				// Note the character that we're backtracking on in the
				// mm array:
				_mms[stackDepth] = i;
				// Now backtrack to target
				ASSERT_ONLY(uint64_t numHits = _params.sink().numHits());
				assert_leq(i+1, _qlen);
				bool ret;
				if(i+1 == _qlen) {
					ret = report(stackDepth, bttop, btbot);
				} else {
					ret = backtrack(stackDepth+1,
				                     i+1,
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
					// No alternative backtracking points; either all
					// possible backtracking targets have been
					// exhausted or all backtracking targets that
					// remain are on spaces with quality values that
					// would exceed the threshold
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
						uint32_t kcur = _qlen - _off - k - 1; // current offset into _qry
						uint8_t kq = QUAL(kcur);
						bool kCurIsAlternative = (ham + kq <= _qualThresh);
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
				assert(sanityCheckEligibility(depth, d, lowAltQual, eligibleSz, eligibleNum, pairs, elims));
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
			d++; cur--;
		}
		assert_eq(0xffffffff, cur);
		assert_gt(bot, top);
		return report(stackDepth, top, bot);
	}
	
	bool report(uint32_t stackDepth, uint32_t top, uint32_t bot) {
		// Possibly report
		if(_oneHit) {
			uint32_t spread = bot - top;
			uint32_t r = top + (_rand.nextU32() % spread);
			for(uint32_t i = 0; i < spread; i++) {
				uint32_t ri = r + i;
				if(ri >= bot) ri -= spread;
				if(_reportOnHit) {
					if(_ebwt.reportChaseOne((*_qry), _qual, _name,
					                        _mms, stackDepth, ri,
					                        top, bot, _qlen, _params))
					{
						return true;
					}
				} else {
					// Not yet smart enough to not immediately report
					// hits 
					assert(false);
				}
			}
			return false;
		} else {
			// Not yet smart enough to report all hits
			assert(false);
			return false;
		}
	}

protected:

	bool sanityCheckEligibility(uint32_t  depth,
	                            uint32_t  d,
	                            uint32_t  lowAltQual,
	                            uint32_t  eligibleSz,
	                            uint32_t  eligibleNum,
	                            uint32_t* pairs,
	                            uint8_t*  elims)
	{
		// Sanity check that the lay of the land is as we
		// expect given eligibleNum and eligibleSz
		size_t i = depth, j = 0;
		uint32_t cumSz = 0;
		uint32_t eligiblesVisited = 0;
		for(; i <= d; i++) {
			uint32_t icur = _qlen - _off - i - 1; // current offset into _qry
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
	uint32_t            _off;    // number of charts from the end of
	                             // the query to start from 
	uint32_t            _itop;   // initial top arrow, or 0xffffffff if
	                             // we're starting from the beginning
	                             // of the query
	uint32_t            _ibot;   // initial bot arrow, or 0xffffffff if
	                             // we're starting from the beginning
	                             // of the query
	uint32_t            _iham;   // initial weighted hamming distance
	uint32_t            _spread; // size of window within which to
	                             // backtrack
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
	String<char>        _nameDefault; // default name, for when it's
	                             // not specified by caller
	String<char>        _qualDefault; // default quals
	RandomSource        _rand;   // Source of pseudo-random numbers
	bool                _verbose;// be talkative
};

#endif /*EBWT_SEARCH_BACKTRACK_H_*/
