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
	ChunkPool(uint32_t chunkSz, uint32_t totSz, bool verbose_) :
		verbose(verbose_), pool_(NULL), cur_(0), chunkSz_(chunkSz),
		totSz_(totSz), lim_(totSz/chunkSz), bits_(lim_)
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
		assert_lt(cur_, lim_);
		uint32_t cur = cur_;
		while(bits_.test(cur)) {
			cur++;
			if(cur >= lim_) {
				cur = 0;
			}
			if(cur == cur_) {
				// Wrapped all the way around without finding a free
				// chunk
				return NULL;
			}
		}
		void * ptr = (void *)(&pool_[cur * chunkSz_]);
		bits_.set(cur);
		if(verbose) cout << "Freeing chunk with offset: " << cur << endl;
		cur_ = cur;
		return ptr;
	}

	/**
	 *
	 */
	void free(void *ptr) {
		uint32_t off = (uint32_t)((int8_t*)ptr - pool_);
		assert_eq(0, off % chunkSz_);
		off /= chunkSz_;
		if(verbose) cout << "Freeing chunk with offset: " << off << endl;
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

	bool verbose;

protected:
	int8_t*  pool_; /// the memory pools
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
		lastCurInPrevPool_ = lim_;
		assert_gt(lim_, 0);
	}

	/**
	 * Reset the pool, freeing all arrays that had been given out.
	 */
	void reset() {
		pools_.clear();
		cur_ = 0;
		curPool_ = 0;
	}

	/**
	 * Allocate a single T from the pool.
	 */
	T* alloc() {
		lazyInit();
		if(cur_ + 1 >= lim_) {
			allocNextPool();
		}
		cur_ ++;
		return &pools_[curPool_][cur_-1];
	}

	/**
	 * Allocate a single T from the pool and clear it.
	 */
	T* allocC() {
		T* t = alloc();
		if(t != NULL) {
			memset(t, 0, sizeof(T));
		}
		return t;
	}

	/**
	 * Allocate an array of Ts from the pool.
	 */
	T* alloc(uint32_t num) {
		lazyInit();
		if(cur_ + num >= lim_) {
			allocNextPool();
		}
		cur_ += num;
		return &pools_[curPool_][cur_-num];
	}

	/**
	 * Allocate an array of Ts and clear them.
	 */
	T* allocC(uint32_t num) {
		T* t = alloc(num);
		if(t != NULL) {
			memset(t, 0, sizeof(T) * num);
		}
		return t;
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

	/**
	 * Free a pointer allocated from this pool.  Fow, we only know how
	 * to free the topmost element.
	 */
	void free(T* t) {
		assert(t != NULL);
		if(pool_->verbose) {
			cout << "Freeing a " << name_ << endl;
		}
		if(cur_ > 0 && t == &pools_[curPool_][cur_-1]) {
			cur_--;
			ASSERT_ONLY(memset(&pools_[curPool_][cur_], 0, sizeof(T)));
			if(cur_ == 0 && curPool_ > 0) {
				assert_eq(curPool_+1, pools_.size());
				if(pool_->verbose) {
					cout << "Freeing a pool" << endl;
				}
				pool_->free(pools_.back());
				pools_.pop_back();
				curPool_--;
				cur_ = lastCurInPrevPool_;
				lastCurInPrevPool_ = lim_;
			}
		}
	}

	/**
	 * Free an array of pointers allocated from this pool.  For now, we
	 * only know how to free the topmost array.
	 */
	void free(T* t, uint32_t num) {
		assert(t != NULL);
		if(pool_->verbose) {
			cout << "Freeing a " << name_ << "; num " << num << endl;
		}
		if(num <= cur_ && t == &pools_[curPool_][cur_ - num]) {
			cur_ -= num;
			ASSERT_ONLY(memset(&pools_[curPool_][cur_], 0, num * sizeof(T)));
			if(cur_ == 0 && curPool_ > 0) {
				assert_eq(curPool_+1, pools_.size());
				if(pool_->verbose) {
					cout << "Freeing a pool" << endl;
				}
				pool_->free(pools_.back());
				pools_.pop_back();
				curPool_--;
				cur_ = lastCurInPrevPool_;
				lastCurInPrevPool_ = lim_;
			}
		}
	}

	/**
	 * Return a unique (with respect to every other object allocated
	 * from this pool) identifier for the last object that was just
	 * allocated.
	 */
	uint32_t lastId() const {
		return (curPool_ << 16) | cur_;
	}

#ifndef NDEBUG
	bool empty() const {
		assert(pools_.empty());
		assert_eq(0, cur_);
		assert_eq(0, curPool_);
		return true;
	}
#endif

protected:

	void allocNextPool() {
		lastCurInPrevPool_ = cur_;
		assert_eq(curPool_+1, pools_.size());
		T *pool;
		try {
			if((pool = (T*)pool_->alloc()) == NULL) {
				throw std::bad_alloc();
			}
		} catch(std::bad_alloc& e) {
			cerr << "Error: Could not allocate " << name_ << " pool #" << (curPool_+2) << " of " << (lim_ * sizeof(T)) << " bytes";
			exit(1);
		}
		ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
		pools_.push_back(pool);
		curPool_++;
		cur_ = 0;
	}

	void lazyInit() {
		if(cur_ == 0 && pools_.empty()) {
			T *pool;
			try {
				if((pool = (T*)pool_->alloc()) == NULL) {
					throw std::bad_alloc();
				}
			} catch(std::bad_alloc& e) {
				std::cerr << "Error: Could not allocate " << name_ << " pool #1" << std::endl;
				exit(1);
			}
			ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
			pools_.push_back(pool);
			assert_eq(1, pools_.size());
		}
		assert(!pools_.empty());
	}

	ChunkPool*      pool_;
	const char     *name_;
	std::vector<T*> pools_; /// the memory pools
	uint32_t        curPool_; /// pool we're current allocating from
	uint32_t        lastCurInPrevPool_;
	uint32_t        lim_;  /// # elements held in pool_
	uint32_t        cur_;  /// index of next free element of pool_
};

#endif /* POOL_H_ */

