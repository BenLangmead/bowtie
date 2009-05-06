/*
 * pool.h
 */

#ifndef POOL_H_
#define POOL_H_

#include <iostream>
#include <vector>
#include <stdexcept>
#include <string.h>
#include <stdlib.h>
#include "bitset.h"

/**
 * Very simple allocator for fixed-size chunks of memory.  Chunk size
 * is set at construction time.  Heap memory is only allocated at
 * construction and deallocated at destruction.
 */
class ChunkPool {
public:
	/**
	 * Initialize a new pool with an initial size of about 'bytes'
	 * bytes.  Exit with an error message if we can't allocate it.
	 */
	ChunkPool(uint32_t chunkSz, uint32_t totSz) :
		pool_(NULL), cur_(0), chunkSz_(chunkSz), totSz_(totSz),
		lim_(totSz/chunkSz), bits_(lim_)
	{
		assert_gt(lim_, 0);
		try {
			if((pool_ = new int8_t[totSz_]) == NULL) {
				throw std::bad_alloc();
			}
		} catch(std::bad_alloc& e) {
			std::cerr << "Error: Could not allocate ChunkPool of "
			          << totSz << " bytes" << std::endl;
			exit(1);
		}
	}

	/**
	 * Delete all the pools.
	 */
	~ChunkPool() {
		if(pool_ != NULL) delete[] pool_;
	}

	/**
	 * Reset the pool, freeing all arrays that had been given out.
	 */
	void reset() {
		cur_ = 0;
		bits_.clear();
		assert_eq(0, bits_.test(0));
	}

	/**
	 * Return our current position.
	 */
	uint32_t pos() {
		return cur_;
	}

	/**
	 * Return our current position.
	 */
	uint32_t remaining() {
		assert_geq(lim_, cur_);
		return lim_ - cur_;
	}

	/**
	 * Allocate a single T from the pool.
	 */
	void* alloc() {
		assert_leq(cur_, lim_);
		if(cur_ == lim_) {
			return NULL;
		}
		bits_.set(cur_);
		return (void *)&pool_[cur_++ * chunkSz_];
	}

	/**
	 *
	 */
	void free(void *ptr) {
		uint32_t off = (uint32_t)((int8_t*)ptr - pool_);
		assert_eq(0, off % chunkSz_);
		off /= chunkSz_;
		bits_.clear(off);
	}

	/**
	 *
	 */
	uint32_t chunkSize() const {
		return chunkSz_;
	}

	/**
	 *
	 */
	uint32_t totalSize() const {
		return totSz_;
	}

protected:
	int8_t*  pool_; /// the memory pools
	uint32_t curPool_; /// pool we're current allocating from
	uint32_t cur_;  /// index of next free element of pool_
	const uint32_t chunkSz_;
	const uint32_t totSz_;
	uint32_t lim_;  /// # elements held in pool_
	FixedBitset2 bits_;
};

/**
 * Class for managing a pool of memory from which items of type T
 * (which must have a default constructor) are allocated.  Does not
 * support freeing or resizing individual items - just allocation and
 * then freeing all items at once.
 */
template<typename T>
class AllocOnlyPool {
	typedef std::pair<uint32_t, uint32_t> U32Pair;
public:
	/**
	 * Initialize a new pool with an initial size of about 'bytes'
	 * bytes.  Exit with an error message if we can't allocate it.
	 */
	AllocOnlyPool(ChunkPool* pool, const char *name) :
		pool_(pool), name_(name), curPool_(0), cur_(0)
	{
		assert(pool != NULL);
		lim_ = pool->chunkSize() / sizeof(T);
		assert_gt(lim_, 0);
	}

	/**
	 * Reset the pool, freeing all arrays that had been given out.
	 */
	void reset() {
#ifndef NDEBUG
		pools_.clear();
#endif
		cur_ = 0;
		curPool_ = 0;
		lastAlloc_ = NULL;
		lastAllocSz_ = 0;
	}

	/**
	 * Rewind to a an old position, essentially freeing everything past
	 * it.
	 */
	void rewind(U32Pair pos) {
		curPool_ = pos.first;
		cur_ = pos.second;
		ASSERT_ONLY(memset(&pools_[curPool_][cur_], 0, (lim_-cur_) * sizeof(T)));
	}

	/**
	 * Return our current position.
	 */
	U32Pair getPos() {
		return make_pair(curPool_, cur_);
	}

	/**
	 * Return the last RangeState allocated from the pool.
	 */
	T* lastAlloc() {
		return lastAlloc_;
	}

	/**
	 * Return the size of the last array of Ts allocated.
	 */
	uint32_t lastAllocSz() const {
		return lastAllocSz_;
	}

	/**
	 * Allocate a single T from the pool.
	 */
	T* alloc() {
		lazyInit();
		if(cur_ + 1 >= lim_) {
			allocNextPool();
		}
		lastAlloc_ = &pools_[curPool_][cur_];
		ASSERT_ONLY(lastAlloc_->allocPool = curPool_);
		ASSERT_ONLY(lastAlloc_->allocCur  = cur_);
		lastAllocSz_ = 1;
		cur_ ++;
		return lastAlloc_;
	}

	/**
	 * Allocate an array of Ts from the pool.
	 */
	T* alloc(uint32_t num) {
		lazyInit();
		if(cur_ + num >= lim_) {
			allocNextPool();
		}
		lastAlloc_ = &pools_[curPool_][cur_];
		ASSERT_ONLY(lastAlloc_->allocPool = curPool_);
		ASSERT_ONLY(lastAlloc_->allocCur  = cur_);
		lastAllocSz_ = num;
		cur_ += num;
		return lastAlloc_;
	}

	/**
	 * Return the current pool.
	 */
	uint32_t curPool() const {
		return curPool_;
	}

	/**
	 * Return the current position within the current pool.
	 */
	uint32_t cur() const {
		return cur_;
	}

	void free(T* t) {
		assert(t != NULL);
		if(t == lastAlloc_ && cur_ >= lastAllocSz_) {
			assert(lastAllocSz_ != 0);
			cur_ -= lastAllocSz_;
			ASSERT_ONLY(memset(&pools_[curPool_][cur_], 0, (lastAllocSz_) * sizeof(T)));
			lastAlloc_ = 0;
			lastAllocSz_ = 0;
		}
	}

#ifndef NDEBUG
	bool empty() const {
		assert(pools_.empty());
		assert_eq(0, cur_);
		assert_eq(0, curPool_);
		assert_eq(0, lastAllocSz_);
		assert(lastAlloc_ == NULL);
		return true;
	}
#endif

protected:

	void allocNextPool() {
		if(curPool_ >= pools_.size()-1) {
			T *pool;
			try {
				pool = new T[lim_];
				if((pool = (T*)pool_->alloc()) == NULL) throw std::bad_alloc();
			} catch(std::bad_alloc& e) {
				cerr << "Error: Could not allocate " << name_ << " pool #" << (curPool_+2) << " of " << (lim_ * sizeof(T)) << " bytes";
				exit(1);
			}
			ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
			pools_.push_back(pool);
		}
		curPool_++;
		cur_ = 0;
		ASSERT_ONLY(memset(pools_[curPool_], 0, lim_ * sizeof(T)));
	}

	void lazyInit() {
		if(cur_ == 0 && pools_.empty()) {
			T *tpool;
			try {
				if((tpool = (T*)pool_->alloc()) == NULL) throw std::bad_alloc();
			} catch(std::bad_alloc& e) {
				std::cerr << "Error: Could not allocate " << name_ << " pool #1" << std::endl;
				exit(1);
			}
			ASSERT_ONLY(memset(tpool, 0, lim_ * sizeof(T)));
			pools_.push_back(tpool);
			assert_eq(1, pools_.size());
		}
		assert(!pools_.empty());
	}

	ChunkPool*      pool_;
	const char     *name_;
	std::vector<T*> pools_; /// the memory pools
	uint32_t        curPool_; /// pool we're current allocating from
	uint32_t        lim_;  /// # elements held in pool_
	uint32_t        cur_;  /// index of next free element of pool_
	T *             lastAlloc_; /// last T array allocated
	uint32_t        lastAllocSz_; /// size of last T array allocated
};

#endif /* POOL_H_ */

