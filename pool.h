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

/**
 * Class for managing a pool of memory from which items of type T
 * (which must have a default constructor) are allocated.  Does not
 * support freeing or resizing individual items - just allocation and
 * then freeing all items at once.
 */
template<typename T>
class AllocOnlyPool {
public:
	/**
	 * Initialize a new pool with an initial size of about 'bytes'
	 * bytes.  Exit with an error message if we can't allocate it.
	 */
	AllocOnlyPool(uint32_t bytes, const char *name) :
		name_(name), curPool_(0), cur_(0)
	{
		lim_ = bytes / sizeof(T);
		T *pool;
		try {
			pool = new T[lim_];
			if(pool == NULL) throw std::bad_alloc();
		} catch(std::bad_alloc& e) {
			std::cerr << "Error: Could not allocate " << name_ << " pool #1 of "
			          << bytes << " bytes" << std::endl;
			exit(1);
		}
		ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
		pools_.push_back(pool);
	}

	/**
	 * Delete all the pools.
	 */
	~AllocOnlyPool() {
		for(size_t i = 0; i < pools_.size(); i++) {
			delete[] pools_[i];
		}
	}

	/**
	 * Reset the pool, freeing all arrays that had been given out.
	 */
	void reset() {
#ifndef NDEBUG
		for(size_t i = 0; i < pools_.size(); i++) {
			memset(pools_[i], 0, (lim_) * sizeof(T));
		}
#endif
		cur_ = 0;
		curPool_ = 0;
		lastAlloc_ = NULL;
		lastAllocSz_ = 0;
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
		if(cur_ + 1 >= lim_) {
			if(curPool_ >= pools_.size()-1) {
				T *pool;
				try {
					pool = new T[lim_];
					if(pool == NULL) throw std::bad_alloc();
				} catch(std::bad_alloc& e) {
					cerr << "Error: Could not allocate " << name_ << " pool #" << (curPool_+2) << " of " << (lim_ * sizeof(T)) << " bytes";
					exit(1);
				}
				ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
				pools_.push_back(pool);
			}
			curPool_++;
			cur_ = 0;
		}
		lastAlloc_ = &pools_[curPool_][cur_];
		lastAllocSz_ = 1;
		cur_ ++;
		return lastAlloc_;
	}

	/**
	 * Allocate an array of Ts from the pool.
	 */
	T* alloc(uint32_t num) {
		if(cur_ + num >= lim_) {
			if(curPool_ >= pools_.size()-1) {
				T *pool;
				try {
					pool = new T[lim_];
					if(pool == NULL) throw std::bad_alloc();
				} catch(std::bad_alloc& e) {
					cerr << "Error: Could not allocate " << name_ << " pool #" << (curPool_+2) << " of " << (lim_ * sizeof(T)) << " bytes";
					exit(1);
				}
				ASSERT_ONLY(memset(pool, 0, lim_ * sizeof(T)));
				pools_.push_back(pool);
			}
			curPool_++;
			cur_ = 0;
		}
		lastAlloc_ = &pools_[curPool_][cur_];
		lastAllocSz_ = num;
		cur_ += num;
		return lastAlloc_;
	}

protected:
	const char     *name_;
	std::vector<T*> pools_; /// the memory pools
	uint32_t        curPool_; /// pool we're current allocating from
	uint32_t        lim_;  /// # elements held in pool_
	uint32_t        cur_;  /// index of next free element of pool_
	T *             lastAlloc_; /// last T array allocated
	uint32_t        lastAllocSz_; /// size of last T array allocated
};

#endif /* POOL_H_ */

