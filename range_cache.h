/*
 * range_cache.h
 *
 * Classes that encapsulate the caching of
 */

#ifndef RANGE_CACHE_H_
#define RANGE_CACHE_H_

#include <stdint.h>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <map>
#include "ebwt.h"
#include "row_chaser.h"

#define RANGE_NOT_SET 0xffffffff
#define RANGE_CACHE_BAD_ALLOC 0xffffffff

/**
 * Manages a pool of memory used exclusively for range cache entries.
 * This manager is allocate-only; it exists mainly so that we can avoid
 * lots of new[]s and delete[]s.
 *
 * A given stretch of words may be one of two types: a cache entry, or
 * a cache entry wrapper.  A cache entry has a length and a list of
 * already-resolved reference positions.  A cache entry wrapper has a
 * pointer to a cache entry for a different range, along with an
 * integer indicating how many "jumps" to the left that range is from
 * the one that owns the wrapper.
 */
class RangeCacheMemPool {
public:
	RangeCacheMemPool(uint32_t lim) :
		lim_(lim), occ_(0), buf_(NULL), closed_(false)
	{
		if(lim > 0) {
			try {
				buf_ = new uint32_t[lim_];
				if(buf_ == NULL) throw std::bad_alloc();
			} catch(std::bad_alloc& e) {
				cerr << "Allocation error allocating " << lim
					 << " words of range-cache memory" << endl;
				exit(1);
			}
			assert(buf_ != NULL);
		}
	}

	~RangeCacheMemPool() {
		// Release all word memory!
		if(lim_ > 0) delete[] buf_;
	}

	/**
	 * Allocate numElts elements from the word pool.
	 */
	uint32_t alloc(uint32_t numElts) {
		assert_gt(numElts, 0);
		assert_leq(occ_, lim_);
		if(occ_ + numElts > lim_ || numElts >= 0x80000000) {
			return RANGE_CACHE_BAD_ALLOC;
		}
		assert_gt(lim_, 0);
		uint32_t ret = occ_;
		occ_ += numElts;
		assert_leq(occ_, lim_);
		if(lim_ - occ_ < 10) {
			// No more room - don't try anymore
			closed_ = true;
		}
		return ret;
	}

	/**
	 * Turn a pool-array index into a pointer; check that it doesn't
	 * fall outside the pool first.
	 */
	inline uint32_t *get(uint32_t off) {
		assert_gt(lim_, 0);
		assert_lt(off, lim_);
		uint32_t *ret = buf_ + off;
		assert_neq(0x80000000, ret[0]);
		return ret;
	}

	/**
	 * Return true iff there's no more room in the cache.
	 */
	inline bool closed() {
		return closed_;
	}

private:
	uint32_t lim_;  /// limit on number of 32-bit words to dish out in total
	uint32_t occ_;  /// number of occupied words
	uint32_t *buf_; /// buffer of 32-bit words
	bool closed_;   ///
};

/**
 * A view to a range of cached reference positions.
 */
class RangeCacheEntry {

	typedef Ebwt<String<Dna> > TEbwt;
	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef RowChaser<String<Dna> > TRowChaser;

public:
	/**
	 *
	 */
	RangeCacheEntry() :
		top_(0xffffffff), jumps_(0), len_(0), ents_(NULL), ebwt_(NULL)
	{ }

	/**
	 * Create a new RangeCacheEntry from the data in the pool at 'ents'.
	 */
	RangeCacheEntry(RangeCacheMemPool& pool, uint32_t top, uint32_t ent, TEbwt* ebwt) {
		init(pool, top, ent, ebwt);
	}

	/**
	 * Initialize a RangeCacheEntry from the data in the pool at 'ents'.
	 */
	void init(RangeCacheMemPool& pool, uint32_t top, uint32_t ent, TEbwt* ebwt) {
		assert(ebwt != NULL);
		top_ = top;
		ebwt_ = ebwt;
		uint32_t *ents = pool.get(ent);
		// Is hi bit set?
		if(ents[0] & 0x80000000) {
			// If so, the target is a wrapper and the non-hi bits
			// contain the # jumps
			jumps_ = (ents[0] & ~0x80000000);
			assert_gt(jumps_, 0);
			assert_leq(jumps_, ebwt_->_eh._len);
			// Get the target entry
			uint32_t *dest = pool.get(ents[1]);
			// Get the length from the target entry
			len_ = dest[0];
			assert_gt(len_, 0);
			assert_leq(len_, ebwt_->_eh._len);
			// Get the pointer to the entries themselves
			ents_ = dest + 1;
		} else {
			// Not a wrapper, so there are no jumps
			jumps_ = 0;
			// Get the length from the target entry
			len_  = ents[0];
			assert_gt(len_, 0);
			assert_leq(len_, ebwt_->_eh._len);
			// Get the pointer to the entries themselves
			ents_ = ents + 1;
		}
		assert(sanityCheckEnts());
	}

	/**
	 * Initialize a wrapper with given number of jumps and given target
	 * entry index.
	 */
	void init(RangeCacheMemPool& pool, uint32_t top, uint32_t jumps,
	          uint32_t ent, TEbwt* ebwt)
	{
		assert(ebwt != NULL);
		ebwt_ = ebwt;
		top_ = top;
		jumps_ = jumps;
		uint32_t *ents = pool.get(ent);
		// Must not be a wrapper
		assert_eq(0, ents[0] & 0x80000000);
		// Get the length from the target entry
		len_ = ents[0];
		assert_gt(len_, 0);
		assert_leq(len_, ebwt_->_eh._len);
		// Get the pointer to the entries themselves
		ents_ = ents + 1;
		assert(sanityCheckEnts());
	}

	uint32_t len() const   {
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		return len_;
	}

	uint32_t jumps() const {
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		return jumps_;
	}

	/**
	 *
	 */
	void reset() {
		ents_ = NULL;
	}

	/**
	 * Return true iff this object represents a valid cache entry.
	 */
	bool valid() const {
		return ents_ != NULL;
	}

	/**
	 * Install a result obtained by a client of this cache; be sure to
	 * adjust for how many jumps down the tunnel the cache entry is
	 * situated.
	 */
	void install(uint32_t elt, uint32_t val) {
		if(ents_ == NULL) {
			// This is not a valid cache entry; do nothing
			return;
		}
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		assert_leq(jumps_, val);
		if(elt < len_) {
			val -= jumps_;
			ASSERT_ONLY(uint32_t sanity = TRowChaser::toFlatRefOff(ebwt_, 1, top_ + elt));
			assert_eq(sanity, val);
			ents_[elt] = val;
		} else {
			// ignore install request
		}
	}

	/**
	 * Get an element from the cache, adjusted for tunnel jumps.
	 */
	inline uint32_t get(uint32_t elt) const {
		if(ents_ == NULL) {
			// This is not a valid cache entry; do nothing
			return RANGE_NOT_SET;
		}
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		if(elt < len_ && ents_[elt] != RANGE_NOT_SET) {
			uint32_t ret = ents_[elt] + jumps_;
			ASSERT_ONLY(uint32_t sanity = TRowChaser::toFlatRefOff(ebwt_, 1, top_ + elt));
			assert_eq(sanity, ret);
			return ret;
		} else {
			return RANGE_NOT_SET;
		}
	}

private:

	/**
	 * Check that len_ and the ents_ array both make sense.
	 */
	bool sanityCheckEnts() {
		assert_gt(len_, 0);
		assert_leq(len_, ebwt_->_eh._len);
		for(size_t i = 0; i < len_; i++) {
			if(ents_[i] == 0xffffffff) continue;
			assert_leq(ents_[i], ebwt_->_eh._len);
			for(size_t j = i+1; j < len_; j++) {
				if(ents_[j] == 0xffffffff) continue;
				assert_neq(ents_[i], ents_[j]);
			}
		}
		return true;
	}

	uint32_t top_;   /// top pointer for this range
	uint32_t jumps_; /// how many tunnel-jumps it is away from the requester
	uint32_t len_;   /// # of entries in cache entry
	uint32_t *ents_; /// ptr to entries, which are flat offs within joined ref
	//U32Pair *ents_;  /// pointer to entries, which are tidx,toff pairs
	TEbwt    *ebwt_; /// index that alignments are in
};

/**
 *
 */
class RangeCache {

	typedef Ebwt<String<Dna> > TEbwt;
	typedef std::vector<uint32_t> TU32Vec;
	typedef std::map<uint32_t, uint32_t> TMap;
	typedef std::map<uint32_t, uint32_t>::iterator TMapItr;

public:
	RangeCache(uint32_t lim, TEbwt* ebwt) :
		lim_(lim), map_(), pool_(lim), closed_(false), ebwt_(ebwt) { }

	/**
	 * Given top and bot offsets, retrieve the canonical cache entry
	 * that best covers that range.  The cache entry may not directly
	 * correspond to the top offset provided, rather, it might be an
	 * entry that lies "at the end of the tunnel" when top and bot are
	 * walked backward.
	 */
	bool lookup(uint32_t top, uint32_t bot, RangeCacheEntry& ent) {
		if(ebwt_ == NULL || lim_ == 0) return false;
		assert_gt(bot, top);
		TMapItr itr = map_.find(top);
		if(itr == map_.end()) {
			// No cache entry for the given 'top' offset
			if(closed_) {
				return false; // failed to get cache entry
			} else {
				if(pool_.closed()) {
					closed_ = true;
					return false; // failed to get cache entry
				}
			}
			// Use the tunnel
			return tunnel(top, bot, ent);
		} else {
			// There is a cache entry for the given 'top' offset
			uint32_t ret = itr->second;
			ent.init(pool_, top, ret, ebwt_);
			return true; // success
		}
	}

protected:
	/**
	 * Tunnel through to the first range that 1) includes all the same
	 * suffixes (though longer) as the given range, and 2) has a cache
	 * entry for it.
	 */
	bool tunnel(uint32_t top, uint32_t bot, RangeCacheEntry& ent) {
		assert_gt(bot, top);
		TU32Vec tops, bots;
		const uint32_t spread = bot - top;
		SideLocus tloc, bloc;
		SideLocus::initFromTopBot(top, bot, ebwt_->_eh, ebwt_->_ebwt, tloc, bloc);
		SideLocus::prefetchTopBot(tloc, bloc);
		uint32_t newtop = top, newbot = bot;
		uint32_t jumps = 0;
		// Walk left through the tunnel
		while(true) {
			jumps++;
			if(ebwt_->rowL(tloc) != ebwt_->rowL(bloc)) {
				// Different characters at top and bot positions of
				// BWT; this means that the calls to mapLF below are
				// guaranteed to yield rows in two different character-
				// sections of the BWT.
				break;
			}
			// Advance top and bot
			newtop = ebwt_->mapLF(tloc);
			newbot = ebwt_->mapLF(bloc);
			assert_geq(newbot, newtop);
			// If the new spread is the same as the old spread, we can
			// be confident that the new range includes all of the same
			// suffixes as the last range (though longer by 1 char)
			if((newbot - newtop) == spread) {
				// Check if newtop is already cached
				TMapItr itr = map_.find(newtop);
				if(itr != map_.end()) {
					// This range, which is further to the left in the
					// same tunnel as the query range, has a cache
					// entry already, so use that
					uint32_t idx = itr->second;
					uint32_t *ents = pool_.get(idx);
					// Allocate a new wrapper
					uint32_t newentIdx = pool_.alloc(2);
					if(newentIdx != RANGE_CACHE_BAD_ALLOC) {
						uint32_t *newent = pool_.get(newentIdx);
						if(ents[0] & 0x80000000) {
							// The cache entry we found was a wrapper; make
							// a new wrapper that points to that wrapper's
							// target, with the appropriate number of jumps
							jumps += (ents[0] & ~0x80000000);
							newent[0] = 0x80000000 | jumps;
							newent[1] = ents[1]; // same target
						} else {
							// The cache entry we found was a wrapper; make
							// a new wrapper that points to that wrapper's
							// target, with the appropriate number of jumps
							newent[0] = 0x80000000 | jumps;
							newent[1] = idx;
						}
						// Initialize the entry
						ent.init(pool_, top, newentIdx, ebwt_);
					} else {
						// Could not allocate a new entry
						ent.init(pool_, top, newentIdx, ebwt_);
					}
					// Cache the entry for 'top'
					map_.insert(make_pair(top, newentIdx));
					return true;
				}
				// Save this range
				tops.push_back(newtop);
				bots.push_back(newbot);
				SideLocus::initFromTopBot(newtop, newbot, ebwt_->_eh, ebwt_->_ebwt, tloc, bloc);
				SideLocus::prefetchTopBot(tloc, bloc);
			} else {
				// Not all the suffixes were preserved, so we can't
				// link the source range's cached result to this
				// range's cached results
				break;
			}
		}
		assert_eq(tops.size(), bots.size());
#ifndef NDEBUG
		for(size_t i = 0; i < tops.size(); i++) {
			assert_eq(spread, bots[i] - tops[i]);
		}
#endif
		// Try to create a new cache entry for the leftmost range in
		// the tunnel (which might be the query range)
		uint32_t newentIdx = pool_.alloc(spread + 1);
		if(newentIdx != RANGE_CACHE_BAD_ALLOC) {
			// Successfully allocated new range cache entry; install it
			uint32_t *newent = pool_.get(newentIdx);
			newent[0] = spread;
			memset(&newent[1], 0xff, spread << 2);
			uint32_t entTop = top;
			uint32_t jumps = 0;
			if(tops.size() > 0) {
				entTop = tops.back();
				jumps = tops.size();
			}
			// Cache the entry for the end of the tunnel
			map_.insert(make_pair(entTop, newentIdx));
			ent.init(pool_, entTop, jumps, newentIdx, ebwt_);
			// Cache a wrapper entry for the query range (if possible)
			uint32_t wrapentIdx = pool_.alloc(2);
			if(wrapentIdx != RANGE_CACHE_BAD_ALLOC) {
				uint32_t *wrapent = pool_.get(wrapentIdx);
				wrapent[0] = 0x80000000 | jumps;
				wrapent[1] = newentIdx;
			}
			return true;
		} else {
			// Could not allocate new range cache entry
			return false;
		}
	}

	uint32_t lim_;           /// Total number of key/val bytes to keep in cache
	TMap map_;               ///
	RangeCacheMemPool pool_; /// Memory pool
	bool closed_;            /// Out of space; no new entries
	TEbwt* ebwt_;            /// Index that alignments are in
};

#endif /* RANGE_CACHE_H_ */
