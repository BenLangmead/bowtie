#ifndef EBWT_SEARCH_BACKTRACK_H_
#define EBWT_SEARCH_BACKTRACK_H_

#define DEFAULT_SPREAD 128
#include "pat.h"

/// An array that transforms Phred qualities into their maq-like
/// equivalents by dividing by ten and rounding to the nearest 1,
/// but saturating at 3.
static unsigned char qualRounds[] = {
	0, 0, 0, 0, 0,                //   0 -   4
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, //   5 -  14
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, //  15 -  24
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  25 -  34
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  35 -  44
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  45 -  54
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  55 -  64
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  65 -  74
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  75 -  84
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  85 -  94
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, //  95 - 104
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 105 - 114
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 115 - 124
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 125 - 134
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 135 - 144
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 145 - 154
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 155 - 164
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 165 - 174
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 175 - 184
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 185 - 194
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 195 - 204
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 205 - 214
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 215 - 224
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 225 - 234
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 235 - 244
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, // 245 - 254
	3                             // 255
};

/// Encapsulates a change made to a query base, i.e. "the 3rd base from
/// the 5' end was changed from an A to a T".  Useful when using
/// for matching seeded by "seedlings".
struct QueryMutation {
	QueryMutation() : pos(0), oldBase(0), newBase(0) { }
	QueryMutation(uint8_t _pos, uint8_t _oldBase, uint8_t _newBase) :
		pos(_pos), oldBase(_oldBase), newBase(_newBase)
	{
		assert_neq(oldBase, newBase);
		assert_leq(oldBase, 4);
		assert_lt(newBase, 4);
	}
	uint8_t pos;
	uint8_t oldBase; // original base from the read
	uint8_t newBase; // mutated to fit the reference in at least one place
};

/**
 * Encapsulates a partial alignment.  Supports up to 256 positions and
 * up to 3 substitutions.  The 'type' field of all the alternative
 * structs tells us whether this entry is a singleton entry, an offset
 * into the spillover list, a non-tail entry in the spillover list, or
 * a tail entry in the spillover list.
 */
typedef union {
	struct {
		uint32_t pos0  : 8; // mismatched pos 1
		uint32_t pos1  : 8; // mismatched pos 2
		uint32_t pos2  : 8; // mismatched pos 3
		uint32_t char0 : 2; // substituted char for pos 1
		uint32_t char1 : 2; // substituted char for pos 2
		uint32_t char2 : 2; // substituted char for pos 3
		uint32_t type  : 2; // type of entry; 0=singleton_entry,
		                    // 1=list_offset, 2=list_entry,
		                    // 3=list_tail
	} entry;
	struct {
		uint32_t off   : 30;// offset into list
		uint32_t type  : 2; // type of entry; 0=singleton,
                            // 1=list_offset, 2=list_entry,
                            // 3=list_tail
	} off; // offset into list
	struct {
		uint32_t unk   : 30;// padding
		uint32_t type  : 2; // type of entry; 0=singleton,
                            // 1=list_offset, 2=list_entry,
                            // 3=list_tail
	} unk;   // unknown
} PartialAlignment;

/**
 * A synchronized data structure for storing partial alignments
 * associated with patids, with particular attention to compactness.
 */
class PartialAlignmentManager {
public:
	PartialAlignmentManager(size_t listSz = 10 * 1024 * 1024) {
		MUTEX_INIT(_partialLock);
		// Reserve space for 10M partialsList entries = 40 MB
		_partialsList.reserve(listSz);
	}

	~PartialAlignmentManager() { }

	/**
	 * Add a set of partial alignments for a particular patid into the
	 * partial-alignment database.  This version locks the database,
	 * and so is safe to call if there are potential readers or
	 * writers currently running.
	 */
	void addPartials(uint32_t patid, const vector<PartialAlignment>& ps) {
		if(ps.size() == 0) return;
		MUTEX_LOCK(_partialLock);
		// Assert that the entry doesn't exist yet
		assert(_partialsMap.find(patid) == _partialsMap.end());
		if(ps.size() == 1) {
			_partialsMap[patid] = ps[0];
			_partialsMap[patid].entry.type = 0; // singleton
		} else {
			PartialAlignment al;
			al.off.off = _partialsList.size();
			al.off.type = 1; // list offset
			_partialsMap[patid] = al;
			assert_gt(ps.size(), 1);
			for(size_t i = 0; i < ps.size()-1; i++) {
				_partialsList.push_back(ps[i]);
				// list entry (non-tail)
				_partialsList.back().entry.type = 2;
			}
			_partialsList.push_back(ps.back());
			// list tail
			_partialsList.back().entry.type = 3;
		}
		// Assert that we added an entry
		assert(_partialsMap.find(patid) != _partialsMap.end());
		MUTEX_UNLOCK(_partialLock);
	}

	/**
	 * Get a set of partial alignments for a particular patid out of
	 * the partial-alignment database.
	 */
	void getPartials(uint32_t patid, vector<PartialAlignment>& ps) {
		assert_eq(0, ps.size());
		MUTEX_LOCK(_partialLock);
		getPartialsUnsync(patid, ps);
		MUTEX_UNLOCK(_partialLock);
	}

	/**
	 * Get a set of partial alignments for a particular patid out of
	 * the partial-alignment database.  This version does not attempt to
	 * lock the database.  This is more efficient than the synchronized
	 * version, but is unsafe if there are other threads that may be
	 * writing to the database.
	 */
	void getPartialsUnsync(uint32_t patid, vector<PartialAlignment>& ps) {
		assert_eq(0, ps.size());
		if(_partialsMap.find(patid) == _partialsMap.end()) {
			return;
		}
		PartialAlignment al = _partialsMap[patid];
		uint32_t type = al.unk.type;
		if(type == 0) {
			// singleton
			ps.push_back(al);
		} else {
			// list
			assert_eq(1, type);
			uint32_t off = al.off.off;
			do {
				assert_lt(off, _partialsList.size());
				ASSERT_ONLY(type = _partialsList[off].entry.type);
				assert(type == 2 || type == 3);
				ps.push_back(_partialsList[off]);
				ASSERT_ONLY(uint32_t pos0 = ps.back().entry.pos0);
				ASSERT_ONLY(uint32_t pos1 = ps.back().entry.pos1);
				ASSERT_ONLY(uint32_t pos2 = ps.back().entry.pos2);
				assert(pos1 == 0xff || pos0 != pos1);
				assert(pos2 == 0xff || pos0 != pos2);
				assert(pos2 == 0xff || pos1 != pos2);
			} while(_partialsList[off++].entry.type == 2);
			assert_eq(3, _partialsList[off-1].entry.type);
		}
		assert_gt(ps.size(), 0);
	}

	/// Call to clear the database when there is only one element in it
	void clear(uint32_t patid) {
		assert_eq(1, _partialsMap.count(patid));
		assert_eq(1, _partialsMap.size());
		_partialsMap.erase(patid);
		assert_eq(0, _partialsMap.size());
		_partialsList.clear();
	}

	size_t size() {
		return _partialsMap.size();
	}

	/**
	 * Convert a partial alignment into a QueryMutation string.
	 */
	static uint8_t toMutsString(const PartialAlignment& pal,
	                            const String<Dna5>& seq,
	                            const String<char>& quals,
	                            String<QueryMutation>& muts)
	{
		reserve(muts, 4);
		assert_eq(0, length(muts));
		uint32_t plen = length(seq);
		assert_gt(plen, 0);
		assert_neq(1, pal.unk.type);
		// Do first mutation
		uint8_t oldQuals = 0;
		uint32_t pos0 = pal.entry.pos0;
		uint8_t tpos0 = plen-1-pos0;
		uint32_t chr0 = pal.entry.char0;
		uint8_t oldChar = (uint8_t)seq[tpos0];
		oldQuals += qualRounds[quals[tpos0]-33]; // take quality hit
		appendValue(muts, QueryMutation(tpos0, oldChar, chr0)); // apply mutation
		if(pal.entry.pos1 != 0xff) {
			// Do second mutation
			uint32_t pos1 = pal.entry.pos1;
			uint8_t tpos1 = plen-1-pos1;
			uint32_t chr1 = pal.entry.char1;
			oldChar = (uint8_t)seq[tpos1];
			oldQuals += qualRounds[quals[tpos1]-33]; // take quality hit
			assert_neq(tpos1, tpos0);
			appendValue(muts, QueryMutation(tpos1, oldChar, chr1)); // apply mutation
			if(pal.entry.pos2 != 0xff) {
				// Do second mutation
				uint32_t pos2 = pal.entry.pos2;
				uint8_t tpos2 = plen-1-pos2;
				uint32_t chr2 = pal.entry.char2;
				oldChar = (uint8_t)seq[tpos2];
				oldQuals += qualRounds[quals[tpos2]-33]; // take quality hit
				assert_neq(tpos2, tpos0);
				assert_neq(tpos2, tpos1);
				append(muts, QueryMutation(tpos2, oldChar, chr2)); // apply mutation
			}
		}
		assert_gt(length(muts), 0);
		assert_leq(length(muts), 3);
		return oldQuals;
	}
private:
	/// Maps patids to partial alignments for that patid
	map<uint32_t, PartialAlignment> _partialsMap;
	/// Overflow for when a patid has more than 1 partial alignment
	vector<PartialAlignment> _partialsList;
	/// Lock for 'partialsMap' and 'partialsList'; necessary because
	/// search threads will be reading and writing them
	MUTEX_T _partialLock;
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
	BacktrackManager(const Ebwt<TStr>* __ebwt,
	                 const EbwtSearchParams<TStr>& __params,
	                 uint32_t __qualThresh,  /// max acceptable q-distance
	                 uint32_t __maxBts = 125, /// maximum # backtracks allowed
	                 uint32_t __reportPartials = 0,
	                 bool __reportExacts = true,
	                 PartialAlignmentManager* __partials = NULL,
	                 String<QueryMutation>* __muts = NULL,
	                 bool __verbose = true,
	                 uint32_t seed = 0,
	                 vector<String<Dna5> >* __os = NULL,
	                 bool __considerQuals = true, // whether to consider quality values when making backtracking decisions
	                 bool __halfAndHalf = false, // hacky way of supporting separate revisitable regions
	                 String<Dna5>* __qry = NULL,
	                 String<char>* __qual = NULL,
	                 String<char>* __name = NULL) :
		_qry(__qry),
		_qlen(0),
		_qual(__qual),
		_name(__name),
		_ebwt(__ebwt),
		_params(__params),
		_unrevOff(0),
		_1revOff(0),
		_2revOff(0),
		_3revOff(0),
		_spread(DEFAULT_SPREAD),
		_maxStackDepth(DEFAULT_SPREAD),
		_qualThresh(__qualThresh),
		_reportOnHit(true),
		_pairs(NULL),
		_elims(NULL),
		_mms(NULL),
		_refcs(NULL),
		_chars(NULL),
		_reportPartials(__reportPartials),
		_reportExacts(__reportExacts),
		_partials(__partials),
		_muts(__muts),
		_os(__os),
		_sanity(_os != NULL && _os->size() > 0),
		_considerQuals(__considerQuals),
		_halfAndHalf(__halfAndHalf),
		_5depth(0),
		_3depth(0),
		_nameDefault("default"),
		_hiDepth(0),
		_numBts(0),
		_totNumBts(0),
		_maxBts(__maxBts),
		_precalcedSideLocus(false),
		_preLtop(),
		_preLbot(),
		_btsAtDepths(NULL),
		_totBtsAtDepths(NULL),
		_maxBts0(28),
		_maxBts1(14),
		//_maxBts0(__maxBts0),
		//_maxBts1(__maxBts1),
		_rand(RandomSource(seed)),
		_verbose(__verbose)
	{
	    // For a 40-bp query range, the _pairs array occupies
	    // 40 * 40 * 8 * 4 = 51,200 bytes, and _elims
	    // occupy 40 * 40 = 1,600 bytes
		assert_geq(strlen(qualDefault), DEFAULT_SPREAD);
 		_qualDefault = qualDefault;
 		if(_qry != NULL) {
 			_qlen = length(*_qry);
 			_spread = length(*_qry);
 			if(_qual == NULL || length(*_qual) == 0) {
 				_qual = &_qualDefault;
 			}
 			assert_geq(length(*_qual), _qlen);
 			for(size_t i = 0; i < length(*_qual); i++) {
 				assert_geq((*_qual)[i], 33);
 			}
 			assert_geq(length(*_qual), _qlen);
 			if(_name == NULL || length(*_name) == 0) {
 				_name = &_nameDefault;
 			}
 			assert_leq(_spread, DEFAULT_SPREAD);
 			// TODO: try to make _maxStackDepth a tighter bound to save
 			// space, especially for _pairs
 	 		_maxStackDepth = length(*_qry);
 	 		_pairs  = new uint32_t[DEFAULT_SPREAD*_maxStackDepth*8];
 	 		_elims  = new uint8_t [DEFAULT_SPREAD*_maxStackDepth];
 	 		memset(_elims, 0, DEFAULT_SPREAD*_maxStackDepth);
 			if(_muts != NULL) {
 				applyMutations();
 			}
 		}
		_mms            = new uint32_t[DEFAULT_SPREAD];
		_refcs          = new char[DEFAULT_SPREAD];
		_btsAtDepths    = new uint32_t[DEFAULT_SPREAD];
		_totBtsAtDepths = new uint32_t[DEFAULT_SPREAD];
		_chars          = new char[DEFAULT_SPREAD];
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
		if(_refcs != NULL) {
			delete[] _refcs; _refcs = NULL;
		}
		if(_btsAtDepths != NULL) {
			delete[] _btsAtDepths; _btsAtDepths = NULL;
		}
		if(_totBtsAtDepths != NULL) {
			delete[] _totBtsAtDepths; _totBtsAtDepths = NULL;
		}
		if(_chars != NULL) {
			delete[] _chars; _chars = NULL;
		}
	}

	#define PAIR_TOP(d, c)    (pairs[d*8 + c + 0])
	#define PAIR_BOT(d, c)    (pairs[d*8 + c + 4])
	#define PAIR_SPREAD(d, c) (PAIR_BOT(d, c) - PAIR_TOP(d, c))
	#define PHRED_QUAL(k)     ((uint8_t)(*_qual)[k] >= 33 ? ((uint8_t)(*_qual)[k] - 33) : 0)
	#define PHRED_QUAL2(q, k) ((uint8_t)(q)[k] >= 33 ? ((uint8_t)(q)[k] - 33) : 0)
	#define QUAL(k)           PHRED_QUAL(k)
	#define QUAL2(q, k)       PHRED_QUAL2(q, k)

	void setQuery(String<Dna5>* __qry,
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

		if(_qual == NULL || empty(*_qual)) {
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
		// TODO: try to make _maxStackDepth a tighter bound to save
		// space, especially for _pairs
 		_maxStackDepth = length(*_qry);
 		if(_pairs == NULL) {
 			_pairs = new uint32_t[DEFAULT_SPREAD*_maxStackDepth*8];
 		}
 		if(_elims == NULL) {
 			_elims = new uint8_t[DEFAULT_SPREAD*_maxStackDepth];
 			memset(_elims, 0, DEFAULT_SPREAD*_maxStackDepth);
 		}
		if(_verbose) {
			String<char> qual = (*_qual);
			if(length(qual) > length(*_qry)) {
				resize(qual, length(*_qry));
			}
			cout << "setQuery(_qry=" << (*_qry) << ", _qual=" << qual << ")" << endl;
		}
	}

	/**
	 * Apply a batch of mutations to this read, possibly displacing a
	 * previous batch of mutations.
	 */
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

	void setOffs(uint32_t __5depth,
	             uint32_t __3depth,
	             uint32_t __unrevOff,
	             uint32_t __1revOff,
	             uint32_t __2revOff,
	             uint32_t __3revOff)
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

	void setEbwt(const Ebwt<TStr>* ebwt) {
		_ebwt = ebwt;
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
		const Ebwt<TStr>& ebwt = *_ebwt;
		int ftabChars = ebwt._eh._ftabChars;
		// m = depth beyond which ftab must not extend or else we might
		// miss some legitimate paths
		uint32_t m = min<uint32_t>(_unrevOff, _qlen);
		int nsInSeed = 0;
		int nsInFtab = 0;
		// Count Ns
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
		for(size_t i = 0; i < ebwt._eh._ftabLen && i < _qlen; i++) {
			if((int)(*_qry)[_qlen-i-1] == 4) nsInFtab++;
		}
		bool ret;
		if(nsInFtab == 0 && m >= (uint32_t)ftabChars) {
			// The ftab doesn't extend past the unrevisitable portion,
			// so we can go ahead and use it
			// Rightmost char gets least significant bit-pair
			uint32_t ftabOff = (*_qry)[_qlen - ftabChars];
			assert_lt(ftabOff, ebwt._eh._ftabLen-1);
			for(int i = ftabChars - 1; i > 0; i--) {
				ftabOff <<= 2;
				assert_lt((uint32_t)(*_qry)[_qlen-i], 4);
				ftabOff |= (uint32_t)(*_qry)[_qlen-i];
				assert_lt(ftabOff, ebwt._eh._ftabLen-1);
			}
			assert_lt(ftabOff, ebwt._eh._ftabLen-1);
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
		if(_reportPartials > 0) {
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
		// Initial sanity checking
		assert_gt(length(*_qry), 0);
		assert_leq(_qlen, length(*_qry));
		assert_geq(length(*_qual), length(*_qry));
		if(_verbose) cout << "backtrack(top=" << top << ", "
		                 << "bot=" << bot << ", "
		                 << "iham=" << iham << ", "
		                 << "_pairs" << _pairs << ", "
		                 << "_elims=" << (void*)_elims << ")" << endl;
		bool oldRetain = sink.retainHits();
		size_t oldRetainSz;
		if(_sanity) {
			// Save some info about hits retained at this point
			sink.setRetainHits(true);
			oldRetainSz = sink.retainedHits().size();
		}
		uint64_t nhits = sink.numHits();

		// Initiate the recursive, randomized quality-aware backtracker
		// with a stack depth of 0 (no backtracks so far)
		memset(_btsAtDepths, 0, DEFAULT_SPREAD * sizeof(uint32_t));
		memset(_totBtsAtDepths, 0, DEFAULT_SPREAD * sizeof(uint32_t));
		_hiDepth = 0; _bailedOnBacktracks = false;
		bool done = backtrack(0, depth, _unrevOff, _1revOff, _2revOff, _3revOff,
		                      top, bot, iham, iham, _pairs, _elims, disableFtab);
		bool hits = sink.numHits() > nhits;

		// Remainder of this function is sanity checking
		sink.setRetainHits(oldRetain); // restore old value
		// If we have the original texts, then we double-check the
		// backtracking result against the naive oracle
		// TODO: also check partial alignments
		if(_sanity &&
		   hits && // no-hit results are confirmed as part of the call to backtrack()
		   // ignore partial alignments; we'll check full alignments downstream
		   _reportPartials == 0 &&
		   !_bailedOnBacktracks)   // ignore excessive-backtracking copouts
		{
			uint32_t maxHitsAllowed = sink.maxHits();
			vector<Hit>& retainedHits = sink.retainedHits();
			vector<int>& retainedStrata = sink.retainedStrata();
			assert_leq(retainedHits.size() - oldRetainSz, maxHitsAllowed);
			assert_eq(retainedHits.size(), retainedStrata.size());
			vector<Hit> oracleHits;
			vector<int> oracleStrata;
			// Invoke the naive oracle, which will place all qualifying
			// hits in the 'oracleHits' vector
			naiveOracle(oracleHits, oracleStrata, iham);
			assert_eq(oracleHits.size(), oracleStrata.size());
			// If we found a hit, it had better be one of the ones
			// that the oracle found
			assert_gt(oracleHits.size(), 0);
			int lastStratum = -1;
			int worstStratum = -1;
			// Get the hit reported by the backtracker
			for(size_t j = oldRetainSz; j < retainedHits.size(); j++) {
				Hit& rhit = retainedHits[j];
				if(worstStratum == -1) {
					worstStratum = retainedStrata[j];
				} else if(retainedStrata[j] > worstStratum) {
					worstStratum = retainedStrata[j];
				}
				// Go through oracleHits and look for a match
				size_t i;
				bool found = false;
				for(i = 0; i < oracleHits.size(); i++) {
					const Hit& h = oracleHits[i];
					if(h.h.first == rhit.h.first &&
					   h.h.second == rhit.h.second)
					{
						assert_eq(h.fw, rhit.fw);
						assert(h.mms == rhit.mms);
						assert_eq(retainedStrata[j], oracleStrata[i]);
						// Erase the element in the oracle vector
						oracleHits.erase(oracleHits.begin() + i);
						if(lastStratum == -1) {
							lastStratum = oracleStrata[i];
						} else if(!sink.spanStrata()) {
							assert_eq(lastStratum, oracleStrata[i]);
						}
						oracleStrata.erase(oracleStrata.begin() + i);
						found = true;
						break;
					}
				}
				assert(found); // assert we found a matchup
			}
			assert_neq(-1, lastStratum);
			assert_neq(-1, worstStratum);
			assert_eq(oracleHits.size(), oracleStrata.size());
			if(maxHitsAllowed == 0xffffffff) {
				if(sink.spanStrata()) {
					// Must have matched every oracle hit
					assert_eq(0, oracleHits.size());
				} else {
					// Must have matched all oracle hits at the best
					// stratum
					for(size_t i = 0; i < oracleStrata.size(); i++) {
						assert_gt(oracleStrata[i], lastStratum);
					}
				}
			}
			if(sink.best()) {
				// Ensure that all oracle hits have a stratum at least
				// as bad as the worst stratum observed in the retained
				// hits
				for(size_t i = 0; i < oracleStrata.size(); i++) {
					assert_geq(oracleStrata[i], worstStratum);
				}
			}
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
	 * mismatch and no backtracking opporunities, false is returned.
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
		assert_leq(stackDepth, _maxStackDepth);
		const Ebwt<TStr>& ebwt = *_ebwt;
		HitSinkPerThread& sink = _params.sink();
		uint64_t prehits = sink.numHits();
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
			if(_numBts == _maxBts) {
				_bailedOnBacktracks = true;
				return false;
			}
			_numBts++;
		}
		// Update stack-depth high water mark
		if(stackDepth > _hiDepth) {
			_hiDepth = stackDepth;
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
		// arrow pairs; all alternative pairs with this quality are
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
			// enforce the boundary-crossing constraints here
			if(_halfAndHalf) {
				assert_eq(0, _reportPartials);
				// Crossing from the hi-half into the lo-half
				if(d == _5depth) {
					if(_3revOff == _2revOff) {
						// 1 and 1

						// The backtracking logic should have prevented us from
						// backtracking more than once into this region
						assert_leq(stackDepth, 1);
						// Reject if we haven't encountered mismatch by this point
						if(stackDepth == 0) {
							if(sink.numHits() == prehits) {
								// We're returning from the bottommost frame
								// without having reported any hits; let's
								// sanity-check that there really aren't any
								confirmNoHit(iham);
							}
							return false;
						}
					} else {
						// 1 and 1,2

						// The backtracking logic should have prevented us from
						// backtracking more than twice into this region
						assert_leq(stackDepth, 2);
						// Reject if we haven't encountered mismatch by this point
						if(stackDepth < 1) {
							if(sink.numHits() == prehits) {
								// We're returning from the bottommost frame
								// without having found any hits; let's
								// sanity-check that there really aren't any
								confirmNoHit(iham);
							}
							return false;
						}
					}
				}
				else if(d == _3depth) {
					if(_3revOff == _2revOff) {
						// 1 and 1
						// The backtracking logic should have prevented us from
						// backtracking more than twice within this region
						assert_leq(stackDepth, 2);
						// Must have encountered two mismatches by this point
						if(stackDepth < 2) {
							if(stackDepth == 0 && sink.numHits() == prehits) {
								confirmNoHit(iham);
							}
							// We're returning from the bottommost frame
							// without having found any hits; let's
							// sanity-check that there really aren't any
							return false;
						}
					} else {
						// 1 and 1,2
						// Count the mismatches in the lo and hi halves
						int loHalfMms = 0, hiHalfMms = 0;
						for(size_t i = 0; i < stackDepth; i++) {
							uint32_t d = _qlen - _mms[i] - 1;
							if     (d < _5depth) hiHalfMms++;
							else if(d < _3depth) loHalfMms++;
							else assert(false);
						}
						assert_leq(loHalfMms+hiHalfMms, 3);
						assert_gt(hiHalfMms, 0);
						if(loHalfMms == 0) {
							// Didn't encounter any mismatches in the lo-half
							if(stackDepth == 0 && sink.numHits() == prehits) {
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
				}
				// In-between sanity checks
				if(d >= _5depth) {
					assert_geq(stackDepth, 1);
				} else if(d >= _3depth) {
					assert_geq(stackDepth, 2);
				}
			}
			bool curIsEligible = false;
			// Reset eligibleNum and eligibleSz if there are any
			// eligible pairs discovered at this spot
			bool curOverridesEligible = false;
			// Determine whether arrow pairs at this location are
			// candidates for backtracking
			int c = (int)(*_qry)[cur];
			assert_leq(c, 4);
			uint8_t q = QUAL(cur);
			assert_lt((uint32_t)q, 100);
			// The current query position is a legit alternative if it a) is
			// not in the unrevisitable region, and b) there is a quality
			// ceiling and its selection would cause the ceiling to be exceeded
			bool curIsAlternative = (d >= unrevOff) &&
			                        (!_considerQuals || (ham + qualRounds[q] <= _qualThresh));
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
			// If c is 'N', then it's a mismatch
			if(c == 4 && d > 0) {
				top = bot = 1;
			} else if(c == 4) {
				assert_eq(0, top);
				assert_eq(0, bot);
			}
			// Calculate the ranges for this position
			if(top == 0 && bot == 0) {
				// Calculate first quartet of pairs using the _fchr[]
				// array
				PAIR_TOP(0, 0)                  = ebwt._fchr[0];
				PAIR_BOT(0, 0) = PAIR_TOP(0, 1) = ebwt._fchr[1];
				PAIR_BOT(0, 1) = PAIR_TOP(0, 2) = ebwt._fchr[2];
				PAIR_BOT(0, 2) = PAIR_TOP(0, 3) = ebwt._fchr[3];
				PAIR_BOT(0, 3)                  = ebwt._fchr[4];
				// Update top and bot
				if(c < 4) { top = PAIR_TOP(d, c); bot = PAIR_BOT(d, c); }
			} else if(curIsAlternative) {
				// Clear pairs
				memset(&pairs[d*8], 0, 8 * 4);
				// Calculate next quartet of pairs
				ebwt.mapLFEx(ltop, lbot, &pairs[d*8], &pairs[(d*8)+4]);
				// Update top and bot
				if(c < 4) { top = PAIR_TOP(d, c); bot = PAIR_BOT(d, c); }
			} else {
				// This query character is not even a legitimate
				// alternative (because backtracking here would blow
				// our mismatch quality budget), so no need to do the
				// bookkeeping for the entire quartet, just do c
				if(c < 4) { top = ebwt.mapLF(ltop, c); bot = ebwt.mapLF(lbot, c); }
			}
			if(top != bot) {
				// Calculate loci from row indices; do it now so that
				// those prefetches are fired off as soon as possible.
				// This eventually calls SideLocus.initfromRow().
				SideLocus::initFromTopBot(top, bot, ebwt._eh, ebwt._ebwt, ltop, lbot);
			}
			// Update the elim array
			if(c < 4) {
				elims[d] = (1 << c);
				assert_gt(elims[d], 0);
				assert_lt(elims[d], 16);
			} else {
				elims[d] = 0;
			}
			assert_lt(elims[d], 16);

			if(curIsAlternative) {
				// Given the just-calculated range quartet, update
				// elims, altNum, eligibleNum, eligibleSz
				for(int i = 0; i < 4; i++) {
					if(i == c) continue;
					assert_leq(PAIR_TOP(d, i), PAIR_BOT(d, i));
					uint32_t spread = PAIR_SPREAD(d, i);
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
								eltop = PAIR_TOP(d, i);
								elbot = PAIR_BOT(d, i);
								assert_eq(elbot-eltop, spread);
								elham = qualRounds[q];
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
						if(sink.numHits() == prehits) {
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
					for(size_t i = 0; i < stackDepth; i++) {
						uint32_t d = _qlen - _mms[i] - 1;
						if     (d < _5depth) hiHalfMms++;
						else if(d < _3depth) loHalfMms++;
						else assert(false);
					}
					assert_leq(loHalfMms + hiHalfMms, lim);
					invalidHalfAndHalf = (loHalfMms == 0 || hiHalfMms == 0);
					if(stackDepth < 2 && altNum > 0) {
						// We backtracked fewer times than necessary;
						// force a backtrack
						mustBacktrack = true;
						backtrackDespiteMatch = true;
					} else if(stackDepth < 2) {
						if(stackDepth == 0 && sink.numHits() == prehits) {
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
			   !invalidExact)         // alignment isn't disqualified by no-exact-hits setting
			{
				uint64_t hits = sink.numHits();
				bool ret = reportAlignment(stackDepth, top, bot);
				if(_sanity && sink.numHits() > hits) {
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
				//SideLocus toploci[4];
				//SideLocus botloci[4];
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
				ASSERT_ONLY(uint32_t eligiblesVisited = 0);
				size_t i = d, j = 0;
				assert_geq(i, depth);
				uint32_t bttop = 0;
				uint32_t btbot = 0;
				uint32_t btham = ham;
				char     btchar = 0;
				int      btcint = 0;
				//SideLocus *btltop = NULL;
				//SideLocus *btlbot = NULL;
				uint32_t icur = 0;
				// The common case is that eligibleSz == 1
				if(eligibleNum > 1 || elignore) {
					bool foundTarget = false;
					// Walk from left to right
					for(; i >= depth; i--) {
						assert_geq(i, unrevOff);
						icur = _qlen - i - 1; // current offset into _qry
						uint8_t qi = QUAL(icur);
						assert_lt(elims[i], 16); // 1.26% in profile (next or prev?)
						if((qi == lowAltQual || !_considerQuals) && elims[i] != 15) {
							// This is the leftmost eligible position with at
							// least one remaining backtrack target
							uint32_t posSz = 0;
							// Add up the spreads for A, C, G, T
							for(j = 0; j < 4; j++) {
								if((elims[i] & (1 << j)) == 0) {
									assert_gt(PAIR_BOT(i, j), PAIR_TOP(i, j));
									assert_gt(PAIR_SPREAD(i, j), 0);
									posSz += PAIR_SPREAD(i, j);
									// Calculate BWT locus and prefetch
									// appropriate cache lines
									//SideLocus::initFromTopBot(
									//	PAIR_TOP(i, j),
									//	PAIR_BOT(i, j),
									//	ebwt._eh,
									//	ebwt._ebwt,
									//	toploci[j],
									//	botloci[j]);
								}
							}
							// Generate a random number
							assert_gt(posSz, 0);
							uint32_t r = _rand.nextU32() % posSz;
							for(j = 0; j < 4; j++) {
								if((elims[i] & (1 << j)) == 0) {
									// This range has not been eliminated
									ASSERT_ONLY(eligiblesVisited++);
									uint32_t spread = PAIR_SPREAD(i, j);
									if(r < spread) {
										// This is our randomly-selected
										// backtrack target
										foundTarget = true;
										bttop = PAIR_TOP(i, j);
										btbot = PAIR_BOT(i, j);
										//SideLocus::initFromTopBot(
										//	bttop,
										//	btbot,
										//	ebwt._eh,
										//	ebwt._ebwt,
										//	toploci[j],
										//	botloci[j]);
										//btltop = &toploci[j];
										//btlbot = &botloci[j];
										btham += qualRounds[qi];
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
				uint32_t *newPairs = pairs + (_spread*8);
				uint8_t  *newElims = elims + (_spread);
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
				_mms[stackDepth] = icur;
				assert_eq(1, dna4Cat[(int)btchar]);
				_refcs[stackDepth] = btchar;
#ifndef NDEBUG
				for(int j = 0; j < (int)stackDepth; j++) {
					assert_neq(_mms[j], icur)
				}
#endif
				_chars[i] = btchar;
				// Now backtrack to target
				_btsAtDepths[stackDepth]++;
				_totBtsAtDepths[stackDepth]++;
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
					                    ftabTop,  // top arrow in pair prior to 'depth'
					                    ftabBot,  // bottom arrow in pair prior to 'depth'
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
				                    bttop,  // top arrow in pair prior to 'depth'
				                    btbot,  // bottom arrow in pair prior to 'depth'
				                    btham,  // weighted hamming distance so far
				                    iham,   // initial weighted hamming distance
				                    newPairs,
				                    newElims);
				}
				if(ret) {
					assert_gt(sink.numHits(), prehits);
					if(_sanity) confirmHit(iham);
					return true; // return, signaling that we're done
				}
				_btsAtDepths[stackDepth+1] = 0; // clear bts deeper in
				if(_numBts >= _maxBts) {
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
					if(stackDepth == 0 && sink.numHits() == prehits) {
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
						uint8_t kq = QUAL(kcur);
						if(k < unrevOff) break; // already visited all revisitable positions
						bool kCurIsAlternative = (ham + qualRounds[kq] <= _qualThresh);
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
										uint32_t spread = PAIR_SPREAD(k, l);
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
											eltop = PAIR_TOP(k, l);
											elbot = PAIR_BOT(k, l);
											assert_eq(elbot-eltop, spread);
											elham = qualRounds[kq];
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
				if(stackDepth == 0 && sink.numHits() == prehits) {
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
				if(stackDepth == 0 && sink.numHits() == prehits) {
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
		if(!ret && stackDepth == 0 && sink.numHits() == prehits) {
			confirmNoHit(iham);
		}
		return ret;
	}

	/**
	 * Print a hit along with information about the backtracking
	 * regions constraining the hit.
	 */
	static void printHit(const vector<String<Dna5> >& os,
	                     const Hit& h,
	                     const String<Dna5>& qry,
	                     size_t qlen,
	                     uint32_t unrevOff,
	                     uint32_t oneRevOff,
	                     uint32_t twoRevOff,
	                     uint32_t threeRevOff,
	                     bool ebwtFw)
	{
		// Print pattern sequence
		cout << "  Pat:  " << qry << endl;
		// Print text sequence
		cout << "  Tseg: ";
		if(ebwtFw) {
			for(size_t i = 0; i < qlen; i++) {
				cout << os[h.h.first][h.h.second + i];
			}
		} else {
			for(int i = (int)qlen-1; i >= 0; i--) {
				cout << os[h.h.first][h.h.second + i];
			}
		}
		cout << endl;
		cout << "  Bt:   ";
		for(int i = (int)qlen-1; i >= 0; i--) {
			if     (i < (int)unrevOff)    cout << "0";
			else if(i < (int)oneRevOff)   cout << "1";
			else if(i < (int)twoRevOff)   cout << "2";
			else if(i < (int)threeRevOff) cout << "3";
			else cout << "X";
		}
		cout << endl;
	}

	/**
	 * Naively search for the same hits that should be found by the
	 * backtracking search.  This is used only if --orig and -s have
	 * been specified for bowtie-debug.
	 */
	static void naiveOracle(const vector<String<Dna5> >& os,
	                        const String<Dna5>& qry,
	                        uint32_t qlen,
	                        const String<char>& qual,
	                        const String<char>& name,
	                        uint32_t patid,
	                        vector<Hit>& hits,
	                        vector<int>& strata,
	                        uint32_t qualThresh,
	                        uint32_t unrevOff,
	                        uint32_t oneRevOff,
	                        uint32_t twoRevOff,
	                        uint32_t threeRevOff,
	                        bool fw,
	                        bool ebwtFw,
	                        uint32_t iham = 0,
	                        String<QueryMutation>* muts = NULL,
	                        bool halfAndHalf = false,
	                        bool invert = false)
	{
		bool fivePrimeOnLeft = (ebwtFw == fw);
	    uint32_t plen = qlen;
		uint8_t *pstr = (uint8_t *)begin(qry, Standard());
	    // For each text...
		for(size_t i = 0; i < os.size(); i++) {
			// For each text position...
			if(length(os[i]) < plen) continue;
			uint32_t olen = length(os[i]);
			uint8_t *ostr = (uint8_t *)begin(os[i], Standard());
			// For each possible alignment of pattern against text
			for(size_t j = 0; j <= olen - plen; j++) {
				size_t rev1mm  = 0; // mismatches observed in the 1-revisitable region
				size_t rev2mm  = 0; // mismatches observed in the 2-revisitable region
				size_t rev3mm  = 0; // mismatches observed in the 3-revisitable region
				uint32_t ham = iham; // weighted hamming distance so far
				FixedBitset<max_read_bp> diffs; // mismatch bitvector
				vector<char> refcs; // reference characters for mms
				refcs.resize(qlen, 0);
				// For each alignment column, from right to left
				bool success = true;
				int ok, okInc;
				if(ebwtFw) {
					ok = j+(int)plen-1;
					okInc = -1;
				} else {
					ok = olen-(j+((int)plen-1))-1;
					okInc = 1;
				}
				bool rejectN = false;
				for(int k = (int)plen-1; k >= 0; k--) {
					size_t kr = plen-1-k;
					if((int)ostr[ok] == 4) {
						rejectN = true;
						break;
					}
					if(pstr[k] != ostr[ok]) {
						ham += qualRounds[QUAL2(qual, k)];
						if(ham > qualThresh) {
							// Alignment is invalid because it exceeds
							// our target weighted hamming distance
							// threshold
							success = false;
							break;
						}
						size_t koff = kr;
						if(invert) koff = (size_t)k;
						// What region does the mm fall into?
						if(koff < unrevOff) {
							// Alignment is invalid because it contains
							// a mismatch in the unrevisitable region
							success = false;
							break;
						} else if(koff < oneRevOff) {
							rev1mm++;
							if(rev1mm > 1) {
								// Alignment is invalid because it
								// contains more than 1 mismatch in the
								// 1-revisitable region
								success = false;
								break;
							}
						} else if(koff < twoRevOff) {
							rev2mm++;
							if(rev2mm > 2) {
								// Alignment is invalid because it
								// contains more than 2 mismatches in the
								// 2-revisitable region
								success = false;
								break;
							}
						} else if(koff < threeRevOff) {
							rev3mm++;
							if(rev3mm > 3) {
								// Alignment is invalid because it
								// contains more than 3 mismatches in the
								// 3-revisitable region
								success = false;
								break;
							}
						}
						if(halfAndHalf) {
							if(twoRevOff == threeRevOff) {
								// 1 and 1
								assert_eq(0, rev3mm);
								if(rev1mm > 1 || rev2mm > 1) {
									// Half-and-half alignment is invalid
									// because it contains more than 1 mismatch
									// in either one or the other half
									success = false;
									break;
								}
							} else {
								// 1 and 1,2
								assert_eq(unrevOff, oneRevOff);
								assert_eq(0, rev1mm);
								if(rev2mm > 2 || rev3mm > 2) {
									success = false;
									break;
								}
							}
						}
						if(fivePrimeOnLeft) {
							diffs.set(k);
							refcs[k] = (char)ostr[ok];
						} else {
							// The 3' end is on on the left end of the
							// pattern, but the diffs vector should
							// encode mismatches w/r/t the 5' end, so
							// we flip
							diffs.set(plen-k-1);
							refcs[plen-k-1] = (char)ostr[ok];
						}
					}
					ok += okInc;
				}
				if(rejectN) {
					// Rejected because the reference half of the
					// alignment contained one or more Ns
					continue;
				}
				if(halfAndHalf) {
					if(twoRevOff == threeRevOff) {
						if(rev1mm != 1 || rev2mm != 1) {
							success = false;
						}
					} else {
						if(rev2mm == 0 || rev3mm == 0 ||
						   rev2mm + rev3mm < 2 ||
						   rev2mm + rev3mm > 3)
						{
							success = false;
						}
					}
				}
				if(success) {
					// It's a hit
					uint32_t off = j;
					int stratum = (int)(rev1mm + rev2mm + rev3mm);
					if(!ebwtFw) {
						off = olen - off;
						off -= plen;
					}
					// Add in mismatches from _muts
					if(muts != NULL) {
						for(size_t i = 0; i < length(*muts); i++) {
							// Entries in _mms[] are in terms of offset into
							// _qry - not in terms of offset from 3' or 5' end
							if(fivePrimeOnLeft) {
								diffs.set((*muts)[i].pos);
								refcs[(*muts)[i].pos] = (*muts)[i].newBase;
							} else {
								diffs.set(plen - (*muts)[i].pos - 1);
								refcs[plen - (*muts)[i].pos - 1] = (*muts)[i].newBase;
							}
							stratum++;
						}
					}
					Hit h(make_pair(i, off),
						  patid,  // read id
						  name,   // read name
						  qry,    // read sequence
						  qual,   // read qualities
						  fw,     // forward/reverse-comp
						  diffs,  // mismatch bitvector
						  refcs);
					hits.push_back(h);
					strata.push_back(stratum);
				} // For each pattern character
			} // For each alignment over current text
		} // For each text
	}


protected:

	/**
	 * Mutate the _qry string according to the contents of the _muts
	 * array, which represents a partial alignment.
	 */
	void applyMutations() {
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
	 * Undo mutations to the _qry string, returning it to the original
	 * read.
	 */
	void undoMutations() {
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
		} else {
			int stratum = 0;
			if(stackDepth > 0) {
				for(size_t i = 0; i < stackDepth; i++) {
					if(_mms[i] >= (_qlen - _3revOff)) {
						// This mismatch falls within the seed; count it
						// toward the stratum to report
						stratum++;
						assert_leq(stratum, 3);
					}
				}
			}
			// Undo mutations
			ASSERT_ONLY(String<Dna5> tmp = (*_qry));
			undoMutations();
			bool hit;
			if(_muts != NULL) {
				// Add the mutations to the _mms[] array
				assert_neq(tmp, (*_qry));
				size_t numMuts = length(*_muts);
				assert_leq(numMuts, _qlen);
				for(size_t i = 0; i < numMuts; i++) {
					// Entries in _mms[] are in terms of offset into
					// _qry - not in terms of offset from 3' or 5' end
					assert_lt(stackDepth + i, DEFAULT_SPREAD);
					// All partial-alignment mutations should fall
					// within bounds
					assert_lt((*_muts)[i].pos, _qlen);
					// All partial-alignment mutations should fall
					// within unrevisitable region
					assert_lt(_qlen - (*_muts)[i].pos - 1, _unrevOff);
#ifndef NDEBUG
					for(size_t j = 0; j < stackDepth + i; j++) {
						assert_neq(_mms[j], (uint32_t)(*_muts)[i].pos);
					}
#endif
					_mms[stackDepth + i] = (*_muts)[i].pos;
					_refcs[stackDepth + i] = "ACGT"[(*_muts)[i].newBase];
					// all muts are in the seed, so they count toward the stratum
					stratum++;
				}
				// Report the range of full alignments
				hit = reportFullAlignment(stackDepth+numMuts, top, bot, stratum);
			} else {
				// Report the range of full alignments
				hit = reportFullAlignment(stackDepth, top, bot, stratum);
			}
			// Re-apply mutations
			applyMutations();
			assert_eq(tmp, (*_qry));
			return hit;
		}
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
		// Possibly report
		assert_gt(_reportPartials, 0);
		assert(_partials != NULL);
		ASSERT_ONLY(uint32_t qualTot = 0);
		PartialAlignment al;
		assert_leq(stackDepth, 3);
		assert_gt(stackDepth, 0);

		// First mismatch
		assert_lt(_mms[0], _qlen);
		// First, append the mismatch position in the read
		al.entry.pos0 = (uint8_t)_mms[0]; // pos
		ASSERT_ONLY(qualTot += qualRounds[((*_qual)[_mms[0]] - 33)]);
		uint32_t ci = _qlen - _mms[0] - 1;
		// _chars[] is index in terms of RHS-relative depth
		int c = (int)(Dna5)_chars[ci];
		assert_lt(c, 4);
		assert_neq(c, (int)(*_qry)[_mms[0]]);
		// Second, append the substituted character for the position
		al.entry.char0 = c;

		if(stackDepth > 1) {
			// Second mismatch
			assert_lt(_mms[1], _qlen);
			// First, append the mismatch position in the read
			al.entry.pos1 = (uint8_t)_mms[1]; // pos
			ASSERT_ONLY(qualTot += qualRounds[((*_qual)[_mms[1]] - 33)]);
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
			// Second mismatch
			assert_lt(_mms[2], _qlen);
			// First, append the mismatch position in the read
			al.entry.pos2 = (uint8_t)_mms[2]; // pos
			ASSERT_ONLY(qualTot += qualRounds[((*_qual)[_mms[2]] - 33)]);
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
			uint8_t qi = QUAL(icur);
			assert_lt(elims[i], 16);
			if((qi == lowAltQual || !_considerQuals) && elims[i] != 15) {
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
					undoMutations();
					cout << "  Unmutated Pat:  " << prefix(*_qry, _qlen) << endl;
					applyMutations();
				}
				cout << "  Pat:            " << prefix(*_qry, _qlen) << endl;
				cout << "  Tseg:           ";
				bool ebwtFw = _params.ebwtFw();
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
			BacktrackManager::printHit(
					(*_os), retainedHits.back(), (*_qry), _qlen, _unrevOff,
					_1revOff, _2revOff, _3revOff, _params.ebwtFw());
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
    		if(h.h.first == rhit.h.first && h.h.second == rhit.h.second) {
    			assert_eq(h.fw, rhit.fw);
    			assert(h.mms == rhit.mms);
    			// It's a match - hit confirmed
    			break;
    		}
		}
		assert_lt(i, oracleHits.size()); // assert we found a matchup
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
		bool ebwtFw = _params.ebwtFw();
		bool fw = _params.fw();
		uint32_t patid = _params.patId();

		naiveOracle((*_os),
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
		            _halfAndHalf);
	}

	String<Dna5>*       _qry;    // query (read) sequence
	size_t              _qlen;   // length of _qry
	String<char>*       _qual;   // quality values for _qry
	String<char>*       _name;   // name of _qry
	const Ebwt<TStr>*   _ebwt;   // Ebwt to search in
	const EbwtSearchParams<TStr>& _params;   // Ebwt to search in
	String<uint8_t>*    _btEquivs; // backtracking equivalence classes
	uint32_t            _unrevOff; // unrevisitable chunk
	uint32_t            _1revOff;  // 1-revisitable chunk
	uint32_t            _2revOff;  // 2-revisitable chunk
	uint32_t            _3revOff;  // 3-revisitable chunk
	uint32_t            _spread; // size of window within which to
	                             // backtrack
	uint32_t            _maxStackDepth;
	uint32_t            _qualThresh; // only accept hits with weighted
	                             // hamming distance <= _qualThresh
	bool                _reportOnHit; // report as soon as we find a
	                             // hit? (as opposed to leaving it up
	                             // to the caller whether to report)
	uint32_t           *_pairs;  // arrow pairs, leveled in parallel
	                             // with decision stack
	uint8_t            *_elims;  // which arrow pairs have been
	                             // eliminated, leveled in parallel
	                             // with decision stack
	uint8_t            *_bts;    // how many backtracks remain in each
	                             // equivalence class
	uint32_t           *_mms;    // array for holding mismatches
	char               *_refcs;  // array for holding mismatches
	// Entries in _mms[] are in terms of offset into
	// _qry - not in terms of offset from 3' or 5' end
	char               *_chars;  // characters selected so far
	// If > 0, report partial alignments up to this many mismatches
	uint32_t            _reportPartials;
	/// Do not report alignments with stratum < this limit
	bool                _reportExacts;
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
	/// Default name, for when it's not specified by caller
	String<char>        _nameDefault;
	/// Default quals
	String<char>        _qualDefault;
	/// Greatest stack depth seen since last reset
	uint32_t            _hiDepth;
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
	/// Keep track of # backtracks at various stack depths, but not
	/// cumulatively (just "in order to get where I am now")
	uint32_t*           _btsAtDepths;
	/// Keep track of # backtracks at various stack depths cumulatively
	uint32_t*           _totBtsAtDepths;
	/// Number of backtracks permitted w/ stackDepth 0 before giving up
	/// in halfAndHalf mode
	uint32_t            _maxBts0;
	/// Number of backtracks permitted w/ stackDepth 1 before giving up
	/// in halfAndHalf mode
	uint32_t            _maxBts1;
	/// Flag to record whether a 'false' return from backtracker is due
	/// to having exceeded one or more backrtacking limits
	bool                _bailedOnBacktracks;
	/// Source of pseudo-random numbers
	RandomSource        _rand;
	/// Be talkative
	bool                _verbose;
	// Holding area for partial alignments
	vector<PartialAlignment> _partialsBuf;
};

#endif /*EBWT_SEARCH_BACKTRACK_H_*/
