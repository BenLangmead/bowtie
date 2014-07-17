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

#define RANGE_NOT_SET OFF_MASK
#define RANGE_CACHE_BAD_ALLOC OFF_MASK

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
	RangeCacheMemPool(TIndexOffU lim /* max cache size in bytes */) :
		lim_(lim >> 2 /* convert to words */), occ_(0), buf_(NULL),
		closed_(false)
	{
		if(lim_ > 0) {
			try {
				buf_ = new TIndexOffU[lim_];
				if(buf_ == NULL) throw std::bad_alloc();
			} catch(std::bad_alloc& e) {
				cerr << "Allocation error allocating " << lim
					 << " words of range-cache memory" << endl;
				throw 1;
			}
			assert(buf_ != NULL);
			// Fill with 1s to signal that these elements are
			// uninitialized
			memset(buf_, 0xff, lim_ << 2 /* convert back to bytes */);
		}
	}

	~RangeCacheMemPool() {
		// Release all word memory!
		if(lim_ > 0) delete[] buf_;
	}

	/**
	 * Allocate numElts elements from the word pool.
	 */
	TIndexOffU alloc(TIndexOffU numElts) {
		assert_gt(numElts, 0);
		assert_leq(occ_, lim_);
		if(occ_ + numElts > lim_ || numElts >= CACHE_WRAPPER_BIT) {
			return RANGE_CACHE_BAD_ALLOC;
		}
		assert_gt(lim_, 0);
		TIndexOffU ret = occ_;
		assert(allocs_.find(ret) == allocs_.end());
		ASSERT_ONLY(allocs_.insert(ret));
		// Clear the first elt so that we don't think there's already
		// something there
#ifndef NDEBUG
		for(TIndexOffU i = 0; i < numElts; i++) {
			assert_eq(OFF_MASK, buf_[occ_ + i]);
		}
#endif
		buf_[occ_] = 0;
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
	inline TIndexOffU *get(TIndexOffU off) {
		assert_gt(lim_, 0);
		assert_lt(off, lim_);
		assert(allocs_.find(off) != allocs_.end());
		TIndexOffU *ret = buf_ + off;
		assert_neq(CACHE_WRAPPER_BIT, ret[0]);
		assert_neq(OFF_MASK, ret[0]);
		return ret;
	}

	/**
	 * Return true iff there's no more room in the cache.
	 */
	inline bool closed() {
		return closed_;
	}

private:
	TIndexOffU lim_;  /// limit on number of 32-bit words to dish out in total
	TIndexOffU occ_;  /// number of occupied words
	TIndexOffU *buf_; /// buffer of 32-bit words
	bool closed_;   ///
#ifndef NDEBUG
	std::set<TIndexOffU> allocs_; // elements allocated
#endif
};

/**
 * A view to a range of cached reference positions.
 */
class RangeCacheEntry {

	typedef Ebwt<String<Dna> > TEbwt;
	typedef RowChaser<String<Dna> > TRowChaser;

public:
	/**
	 *
	 */
	RangeCacheEntry(bool sanity = false) :
		top_(OFF_MASK), jumps_(0), len_(0), ents_(NULL), ebwt_(NULL),
		sanity_(sanity)
	{ }

	/**
	 * Create a new RangeCacheEntry from the data in the pool at 'ents'.
	 */
	RangeCacheEntry(RangeCacheMemPool& pool, TIndexOffU top,
			TIndexOffU ent, TEbwt* ebwt, bool sanity = false) :
	    sanity_(sanity)
	{
		init(pool, top, ent, ebwt);
	}

	/**
	 * Initialize a RangeCacheEntry from the data in the pool at 'ents'.
	 */
	void init(RangeCacheMemPool& pool, TIndexOffU top, TIndexOffU ent, TEbwt* ebwt) {
		assert(ebwt != NULL);
		top_ = top;
		ebwt_ = ebwt;
		TIndexOffU *ents = pool.get(ent);
		assert_neq(CACHE_WRAPPER_BIT, ents[0]);
		// Is hi bit set?
		if((ents[0] & CACHE_WRAPPER_BIT) != 0) {
			// If so, the target is a wrapper and the non-hi bits
			// contain the # jumps
			jumps_ = (ents[0] & ~CACHE_WRAPPER_BIT);
			assert_gt(jumps_, 0);
			assert_leq(jumps_, ebwt_->_eh._len);
			// Get the target entry
			TIndexOffU *dest = pool.get(ents[1]);
			// Get the length from the target entry
			len_ = dest[0];
			assert_leq(top_ + len_, ebwt_->_eh._len);
			assert_gt(len_, 0);
			assert_leq(len_, ebwt_->_eh._len);
			// Get the pointer to the entries themselves
			ents_ = dest + 1;
		} else {
			// Not a wrapper, so there are no jumps
			jumps_ = 0;
			// Get the length from the target entry
			len_  = ents[0];
			assert_leq(top_ + len_, ebwt_->_eh._len);
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
	void init(RangeCacheMemPool& pool, TIndexOffU top, TIndexOffU jumps,
			TIndexOffU ent, TEbwt* ebwt)
	{
		assert(ebwt != NULL);
		ebwt_ = ebwt;
		top_ = top;
		jumps_ = jumps;
		TIndexOffU *ents = pool.get(ent);
		// Must not be a wrapper
		assert_eq(0, ents[0] & CACHE_WRAPPER_BIT);
		// Get the length from the target entry
		len_ = ents[0];
		assert_gt(len_, 0);
		assert_leq(len_, ebwt_->_eh._len);
		// Get the pointer to the entries themselves
		ents_ = ents + 1;
		assert_leq(top_ + len_, ebwt_->_eh._len);
		assert(sanityCheckEnts());
	}

	TIndexOffU len() const   {
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		return len_;
	}

	TIndexOffU jumps() const {
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

	TEbwt *ebwt() {
		return ebwt_;
	}

	/**
	 * Install a result obtained by a client of this cache; be sure to
	 * adjust for how many jumps down the tunnel the cache entry is
	 * situated.
	 */
	void install(TIndexOffU elt, TIndexOffU val) {
		if(ents_ == NULL) {
			// This is not a valid cache entry; do nothing
			return;
		}
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		assert_leq(jumps_, val);
		assert_neq(OFF_MASK, val);
		assert_leq(top_ + len_, ebwt_->_eh._len);
		if(elt < len_) {
			val -= jumps_;
			if(verbose_) cout << "Installed reference offset: " << (top_ + elt) << endl;
			ASSERT_ONLY(TIndexOffU sanity = TRowChaser::toFlatRefOff(ebwt_, 1, top_ + elt));
			assert_eq(sanity, val);
#ifndef NDEBUG
			for(size_t i = 0; i < len_; i++) {
				if(i == elt) continue;
				assert_neq(val, ents_[i]);
			}
#endif
			ents_[elt] = val;
		} else {
			// ignore install request
			if(verbose_) cout << "Fell off end of cache entry for install: " << (top_ + elt) << endl;
		}
	}

	/**
	 * Get an element from the cache, adjusted for tunnel jumps.
	 */
	inline TIndexOffU get(TIndexOffU elt) const {
		if(ents_ == NULL) {
			// This is not a valid cache entry; do nothing
			return RANGE_NOT_SET;
		}
		assert(ents_ != NULL);
		assert(ebwt_ != NULL);
		assert_leq(top_ + len_, ebwt_->_eh._len);
		if(elt < len_ && ents_[elt] != RANGE_NOT_SET) {
			if(verbose_) cout << "Retrieved result from cache: " << (top_ + elt) << endl;
			TIndexOffU ret = ents_[elt] + jumps_;
			ASSERT_ONLY(TIndexOffU sanity = TRowChaser::toFlatRefOff(ebwt_, 1, top_ + elt));
			assert_eq(sanity, ret);
			return ret;
		} else {
			if(verbose_) cout << "Cache entry not set: " << (top_ + elt) << endl;
			return RANGE_NOT_SET;
		}
	}

	/**
	 * Check that len_ and the ents_ array both make sense.
	 */
	static bool sanityCheckEnts(TIndexOffU len, TIndexOffU *ents, TEbwt* ebwt) {
		assert_gt(len, 0);
		assert_leq(len, ebwt->_eh._len);
		if(len < 10) {
			for(size_t i = 0; i < len; i++) {
				if(ents[i] == OFF_MASK) continue;
				assert_leq(ents[i], ebwt->_eh._len);
				for(size_t j = i+1; j < len; j++) {
					if(ents[j] == OFF_MASK) continue;
					assert_neq(ents[i], ents[j]);
				}
			}
		} else {
			std::set<TIndexOffU> seen;
			for(size_t i = 0; i < len; i++) {
				if(ents[i] == OFF_MASK) continue;
				assert(seen.find(ents[i]) == seen.end());
				seen.insert(ents[i]);
			}
		}
		return true;
	}

	/**
	 * Check that len_ and the ents_ array both make sense.
	 */
	bool sanityCheckEnts() {
		return RangeCacheEntry::sanityCheckEnts(len_, ents_, ebwt_);
	}

private:

	TIndexOffU top_;   /// top pointer for this range
	TIndexOffU jumps_; /// how many tunnel-jumps it is away from the requester
	TIndexOffU len_;   /// # of entries in cache entry
	TIndexOffU *ents_; /// ptr to entries, which are flat offs within joined ref
	TEbwt    *ebwt_; /// index that alignments are in
	bool     verbose_; /// be talkative?
	bool     sanity_;  /// do consistency checks?
};

/**
 *
 */
class RangeCache {

	typedef Ebwt<String<Dna> > TEbwt;
	typedef std::vector<TIndexOffU> TUVec;
	typedef std::map<TIndexOffU, TIndexOffU> TMap;
	typedef std::map<TIndexOffU, TIndexOffU>::iterator TMapItr;

public:
	RangeCache(size_t lim, TEbwt* ebwt) :
		lim_(lim), map_(), pool_(lim), closed_(false), ebwt_(ebwt), sanity_(true) { }

	/**
	 * Given top and bot offsets, retrieve the canonical cache entry
	 * that best covers that range.  The cache entry may not directly
	 * correspond to the top offset provided, rather, it might be an
	 * entry that lies "at the end of the tunnel" when top and bot are
	 * walked backward.
	 */
	bool lookup(TIndexOffU top, TIndexOffU bot, RangeCacheEntry& ent) {
		if(ebwt_ == NULL || lim_ == 0) return false;
		assert_gt(bot, top);
		ent.reset();
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
			bool ret = tunnel(top, bot, ent);
			return ret;
		} else {
			// There is a cache entry for the given 'top' offset
			TIndexOffU ret = itr->second;
			ent.init(pool_, top, ret, ebwt_);
			return true; // success
		}
	}

	/**
	 * Exhaustively check all entries linked to from map_ to ensure
	 * they're well-formed.
	 */
	bool repOk() {
#ifndef NDEBUG
		for(TMapItr itr = map_.begin(); itr != map_.end(); itr++) {
			TIndexOffU top = itr->first;
			TIndexOffU idx = itr->second;
			TIndexOffU jumps = 0;
			assert_leq(top, ebwt_->_eh._len);
			TIndexOffU *ents = pool_.get(idx);
			if((ents[0] & CACHE_WRAPPER_BIT) != 0) {
				jumps = ents[0] & ~CACHE_WRAPPER_BIT;
				assert_leq(jumps, ebwt_->_eh._len);
				idx = ents[1];
				ents = pool_.get(idx);
			}
			TIndexOffU len = ents[0];
			assert_leq(top + len, ebwt_->_eh._len);
			RangeCacheEntry::sanityCheckEnts(len, ents + 1, ebwt_);
		}
#endif
		return true;
	}

protected:

	/**
	 * Tunnel through to the first range that 1) includes all the same
	 * suffixes (though longer) as the given range, and 2) has a cache
	 * entry for it.
	 */
	bool tunnel(TIndexOffU top, TIndexOffU bot, RangeCacheEntry& ent) {
		assert_gt(bot, top);
		TUVec tops;
		const TIndexOffU spread = bot - top;
		SideLocus tloc, bloc;
		SideLocus::initFromTopBot(top, bot, ebwt_->_eh, ebwt_->_ebwt, tloc, bloc);
		TIndexOffU newtop = top, newbot = bot;
		TIndexOffU jumps = 0;
		// Walk left through the tunnel
		while(true) {
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
			assert_leq(newbot - newtop, spread);
			// If the new spread is the same as the old spread, we can
			// be confident that the new range includes all of the same
			// suffixes as the last range (though longer by 1 char)
			if((newbot - newtop) == spread) {
				// Check if newtop is already cached
				TMapItr itr = map_.find(newtop);
				jumps++;
				if(itr != map_.end()) {
					// This range, which is further to the left in the
					// same tunnel as the query range, has a cache
					// entry already, so use that
					TIndexOffU idx = itr->second;
					TIndexOffU *ents = pool_.get(idx);
					if((ents[0] & CACHE_WRAPPER_BIT) != 0) {
						// The cache entry we found was a wrapper; make
						// a new wrapper that points to that wrapper's
						// target, with the appropriate number of jumps
						jumps += (ents[0] & ~CACHE_WRAPPER_BIT);
						idx = ents[1];
					}
					// Allocate a new wrapper
					TIndexOffU newentIdx = pool_.alloc(2);
					if(newentIdx != RANGE_CACHE_BAD_ALLOC) {
						// We successfully made a new wrapper entry;
						// now populate it and install it in map_
						TIndexOffU *newent = pool_.get(newentIdx); // get ptr to it
						assert_eq(0, newent[0]);
						newent[0] = CACHE_WRAPPER_BIT | jumps; // set jumps
						newent[1] = idx;                // set target
						assert(map_.find(top) == map_.end());
						map_[top] = newentIdx;
						if(sanity_) assert(repOk());
					}
					// Initialize the entry
					ent.init(pool_, top, jumps, idx, ebwt_);
					return true;
				}
				// Save this range
				tops.push_back(newtop);
				SideLocus::initFromTopBot(newtop, newbot, ebwt_->_eh, ebwt_->_ebwt, tloc, bloc);
			} else {
				// Not all the suffixes were preserved, so we can't
				// link the source range's cached result to this
				// range's cached results
				break;
			}
			assert_eq(jumps, tops.size());
		}
		assert_eq(jumps, tops.size());
		// Try to create a new cache entry for the leftmost range in
		// the tunnel (which might be the query range)
		TIndexOffU newentIdx = pool_.alloc(spread + 1);
		if(newentIdx != RANGE_CACHE_BAD_ALLOC) {
			// Successfully allocated new range cache entry; install it
			TIndexOffU *newent = pool_.get(newentIdx);
			assert_eq(0, newent[0]);
			// Store cache-range length in first word
			newent[0] = spread;
			assert_lt(newent[0], CACHE_WRAPPER_BIT);
			assert_eq(spread, newent[0]);
			TIndexOffU entTop = top;
			TIndexOffU jumps = 0;
			if(tops.size() > 0) {
				entTop = tops.back();
				jumps = tops.size();
			}
			// Cache the entry for the end of the tunnel
			assert(map_.find(entTop) == map_.end());
			map_[entTop] = newentIdx;
			if(sanity_) assert(repOk());
			ent.init(pool_, entTop, jumps, newentIdx, ebwt_);
			assert_eq(spread, newent[0]);
			if(jumps > 0) {
				assert_neq(entTop, top);
				// Cache a wrapper entry for the query range (if possible)
				TIndexOffU wrapentIdx = pool_.alloc(2);
				if(wrapentIdx != RANGE_CACHE_BAD_ALLOC) {
					TIndexOffU *wrapent = pool_.get(wrapentIdx);
					assert_eq(0, wrapent[0]);
					wrapent[0] = CACHE_WRAPPER_BIT | jumps;
					wrapent[1] = newentIdx;
					assert(map_.find(top) == map_.end());
					map_[top] = wrapentIdx;
					if(sanity_) assert(repOk());
				}
			}
			return true;
		} else {
			// Could not allocate new range cache entry
			return false;
		}
	}

	TIndexOffU lim_;           /// Total number of key/val bytes to keep in cache
	TMap map_;               ///
	RangeCacheMemPool pool_; /// Memory pool
	bool closed_;            /// Out of space; no new entries
	TEbwt* ebwt_;            /// Index that alignments are in
	bool sanity_;
};

#endif /* RANGE_CACHE_H_ */
