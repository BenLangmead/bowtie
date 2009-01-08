#ifndef EBWT_SEARCH_UTIL_H_
#define EBWT_SEARCH_UTIL_H_

#include <iostream>
#include <vector>
#include <map>
#include <stdint.h>
#include <seqan/sequence.h>
#include "hit.h"
#include "qual.h"

/// Encapsulates a change made to a query base, i.e. "the 3rd base from
/// the 5' end was changed from an A to a T".  Useful when using
/// for matching seeded by "seedlings".
struct QueryMutation {
	QueryMutation() : pos(0), oldBase(0), newBase(0) { }
	QueryMutation(uint16_t _pos, uint8_t _oldBase, uint8_t _newBase) :
		pos(_pos), oldBase(_oldBase), newBase(_newBase)
	{
		assert_neq(oldBase, newBase);
		assert_leq(oldBase, 4);
		assert_lt(newBase, 4);
	}
	uint16_t pos;
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
		uint64_t pos0  : 16;   // mismatched pos 1
		uint64_t pos1  : 16;   // mismatched pos 2
		uint64_t pos2  : 16;   // mismatched pos 3
		uint64_t char0 : 2;    // substituted char for pos 1
		uint64_t char1 : 2;    // substituted char for pos 2
		uint64_t char2 : 2;    // substituted char for pos 3
		uint64_t reserved : 8;
		uint64_t type  : 2;    // type of entry; 0=singleton_entry,
		                       // 1=list_offset, 2=list_entry,
		                       // 3=list_tail
	} entry;
	struct {
		uint64_t off   : 62;// offset into list
		uint64_t type  : 2; // type of entry; 0=singleton,
                            // 1=list_offset, 2=list_entry,
                            // 3=list_tail
	} off; // offset into list
	struct {
		uint64_t unk   : 62;// padding
		uint64_t type  : 2; // type of entry; 0=singleton,
                            // 1=list_offset, 2=list_entry,
                            // 3=list_tail
	} unk;   // unknown
	struct {
		uint64_t u64   : 64;
	} u64;
} PartialAlignment;

#ifndef NDEBUG
static bool sameHalfPartialAlignment(PartialAlignment pa1, PartialAlignment pa2) {
	if(pa1.unk.type == 1 || pa2.unk.type == 1) return false;
	assert_neq(0xff, pa1.entry.pos0);
	assert_neq(0xff, pa2.entry.pos0);

	// Make sure pa1's pos0 is represented in pa1
	if(pa1.entry.pos0 == pa2.entry.pos0) {
		if(pa1.entry.char0 != pa2.entry.char0) return false;
	} else if(pa1.entry.pos0 == pa2.entry.pos1) {
		if(pa1.entry.char0 != pa2.entry.char1) return false;
	} else if(pa1.entry.pos0 == pa2.entry.pos2) {
		if(pa1.entry.char0 != pa2.entry.char2) return false;
	} else {
		return false;
	}
	if(pa1.entry.pos1 != 0xff) {
		if       (pa1.entry.pos1 == pa2.entry.pos0) {
			if(pa1.entry.char1 != pa2.entry.char0) return false;
		} else if(pa1.entry.pos1 == pa2.entry.pos1) {
			if(pa1.entry.char1 != pa2.entry.char1) return false;
		} else if(pa1.entry.pos1 == pa2.entry.pos2) {
			if(pa1.entry.char1 != pa2.entry.char2) return false;
		} else {
			return false;
		}
	}
	if(pa1.entry.pos2 != 0xff) {
		if       (pa1.entry.pos2 == pa2.entry.pos0) {
			if(pa1.entry.char2 != pa2.entry.char0) return false;
		} else if(pa1.entry.pos2 == pa2.entry.pos1) {
			if(pa1.entry.char2 != pa2.entry.char1) return false;
		} else if(pa1.entry.pos2 == pa2.entry.pos2) {
			if(pa1.entry.char2 != pa2.entry.char2) return false;
		} else {
			return false;
		}
	}
	return true;
}

static bool samePartialAlignment(PartialAlignment pa1, PartialAlignment pa2) {
	return sameHalfPartialAlignment(pa1, pa2) && sameHalfPartialAlignment(pa2, pa1);
}

static bool validPartialAlignment(PartialAlignment pa) {
	if(pa.entry.pos0 != 0xff) {
		if(pa.entry.pos0 == pa.entry.pos1) return false;
		if(pa.entry.pos0 == pa.entry.pos2) return false;
	} else {
		if(pa.entry.pos1 != 0xff) return false;
		if(pa.entry.pos2 != 0xff) return false;
	}

	if(pa.entry.pos1 != 0xff) {
		if(pa.entry.pos1 == pa.entry.pos2) return false;
	} else {
		if(pa.entry.pos2 != 0xff) return false;
	}
	return true;
}
#endif

/**
 * Calculate a per-read random seed based on a combination of
 * the read data (incl. sequence, name, quals) and the global
 * seed in '_randSeed'.
 */
static inline uint32_t genRandSeed(const String<Dna5>& qry,
                                   const String<char>& qual,
                                   const String<char>& name)
{
	// Calculate a per-read random seed based on a combination of
	// the read data (incl. sequence, name, quals) and the global
	// seed
	uint32_t rseed = 0;
	size_t qlen = seqan::length(qry);
	// Throw all the characters of the read into the random seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qry[i];
		assert_leq(p, 4);
		size_t off = ((i & 15) << 1);
		rseed |= (p << off);
	}
	// Throw all the quality values for the read into the random
	// seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qual[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed |= (p << off);
	}
	// Throw all the characters in the read name into the random
	// seed
	size_t namelen = seqan::length(name);
	for(size_t i = 0; i < namelen; i++) {
		int p = (int)name[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed |= (p << off);
	}
	return rseed;
}

extern
void printHit(const vector<String<Dna5> >& os,
			  const Hit& h,
			  const String<Dna5>& qry,
			  size_t qlen,
			  uint32_t unrevOff,
			  uint32_t oneRevOff,
			  uint32_t twoRevOff,
			  uint32_t threeRevOff,
			  bool ebwtFw);

extern
void naiveOracle(const vector<String<Dna5> >& os,
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
				 uint32_t iham,
				 String<QueryMutation>* muts,
				 bool maqPenalty,
				 bool halfAndHalf,
				 bool reportExacts,
				 bool invert);

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
		size_t origPlSz = _partialsList.size();
		// Assert that the entry doesn't exist yet
		assert(_partialsMap.find(patid) == _partialsMap.end());
		if(ps.size() == 1) {
			_partialsMap[patid] = ps[0];
			_partialsMap[patid].entry.type = 0; // singleton
		} else {
#ifndef NDEBUG
			// Make sure there are not duplicate entries
			for(size_t i = 0; i < ps.size()-1; i++)
				for(size_t j = i+1; j < ps.size(); j++)
					assert(!samePartialAlignment(ps[i], ps[j]));
#endif
			// Insert a "pointer" record into the map that refers to
			// the stretch of the _partialsList vector that contains
			// the partial alignments.
			PartialAlignment al;
			al.u64.u64 = 0xffffffffffffffffllu;
			al.off.off = origPlSz;
			al.off.type = 1; // list offset
			_partialsMap[patid] = al; // install pointer
			assert_gt(ps.size(), 1);
			// Now add all the non-tail partial alignments (all but the
			// last) to the _partialsList
			for(size_t i = 0; i < ps.size()-1; i++) {
				assert(validPartialAlignment(ps[i]));
				_partialsList.push_back(ps[i]);
				// list entry (non-tail)
				_partialsList.back().entry.type = 2;
			}
			// Now add the tail (last) partial alignment and mark it as
			// such
			assert(validPartialAlignment(ps.back()));
			_partialsList.push_back(ps.back());
			// list tail
			_partialsList.back().entry.type = 3;
#ifndef NDEBUG
			// Make sure there are not duplicate entries
			assert_eq(_partialsList.size(), origPlSz + ps.size());
			for(size_t i = origPlSz; i < _partialsList.size()-1; i++) {
				for(size_t j = i+1; j < _partialsList.size(); j++) {
					assert(!samePartialAlignment(_partialsList[i], _partialsList[j]));
				}
			}
#endif
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
		PartialAlignment al;
		al.u64.u64 = _partialsMap[patid].u64.u64;
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
#ifndef NDEBUG
				// Make sure this entry isn't equal to any other entry
				for(size_t i = 0; i < ps.size(); i++) {
					assert(validPartialAlignment(ps[i]));
					assert(!samePartialAlignment(ps[i], _partialsList[off]));
				}
#endif
				assert(validPartialAlignment(_partialsList[off]));
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
		assert_eq(0, _partialsList.size());
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
	                            String<QueryMutation>& muts,
	                            bool maqPenalty = true)
	{
		reserve(muts, 4);
		assert_eq(0, length(muts));
		uint32_t plen = length(seq);
		assert_gt(plen, 0);
		assert_neq(1, pal.unk.type);
		// Do first mutation
		uint8_t oldQuals = 0;
		uint32_t pos0 = pal.entry.pos0;
		uint16_t tpos0 = plen-1-pos0;
		uint32_t chr0 = pal.entry.char0;
		uint8_t oldChar = (uint8_t)seq[tpos0];
		oldQuals += mmPenalty(maqPenalty, quals[tpos0]-33); // take quality hit
		appendValue(muts, QueryMutation(tpos0, oldChar, chr0)); // apply mutation
		if(pal.entry.pos1 != 0xff) {
			// Do second mutation
			uint32_t pos1 = pal.entry.pos1;
			uint16_t tpos1 = plen-1-pos1;
			uint32_t chr1 = pal.entry.char1;
			oldChar = (uint8_t)seq[tpos1];
			oldQuals += mmPenalty(maqPenalty, quals[tpos1]-33); // take quality hit
			assert_neq(tpos1, tpos0);
			appendValue(muts, QueryMutation(tpos1, oldChar, chr1)); // apply mutation
			if(pal.entry.pos2 != 0xff) {
				// Do second mutation
				uint32_t pos2 = pal.entry.pos2;
				uint16_t tpos2 = plen-1-pos2;
				uint32_t chr2 = pal.entry.char2;
				oldChar = (uint8_t)seq[tpos2];
				oldQuals += mmPenalty(maqPenalty, quals[tpos2]-33); // take quality hit
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

struct BacktrackLimits {

	BacktrackLimits() {
		maxBts  = 0;
		maxBts0 = 0;
		maxBts1 = 0;
		maxBts2 = 0;
	}

	BacktrackLimits(uint32_t _maxBts,
	                uint32_t _maxBts0,
	                uint32_t _maxBts1,
	                uint32_t _maxBts2)
	{
		maxBts  = _maxBts;
		maxBts0 = _maxBts0;
		maxBts1 = _maxBts1;
		maxBts2 = _maxBts2;
	}

	uint32_t maxBts;  // limit on overall number of backtracks per read
	uint32_t maxBts0; // limit on number of backtracks at depth 0
	uint32_t maxBts1; // limit on number of backtracks at depth 1
	uint32_t maxBts2; // limit on number of backtracks at depth 2
};

#endif /* EBWT_SEARCH_UTIL_H_ */
