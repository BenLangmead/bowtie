#ifndef BITSET_H_
#define BITSET_H_

#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string.h>
#include <stdexcept>
#include "assert_helpers.h"
#include "threading.h"
#include "btypes.h"

/**
 * Given a words array and a size, allocate a new, larger array, moving
 * data from the old to the new array, and set all newly-allocated
 * words to 0.  Return the new, larger array, which can be substituted
 * for the old one.  The new array is larger than the old by about 50%.
 */
static inline TIndexOffU*
bitsetRealloc(TIndexOffU& sz, TIndexOffU* words, const char *errmsg = NULL) {
	TIndexOffU oldsz = sz;
	if(sz > 0) {
		sz += (sz >> 1) + BITSET_MASK; // Add 50% more elements, plus a bit
		sz &= ~BITSET_MASK;            // Make sure it's 32-aligned
	} else {
		sz = 1024; // Start off at 1024 bits to avoid many expansions
	}
	assert_gt(sz, oldsz);
	assert_eq(0, (sz & BITSET_MASK));
	TIndexOffU *newwords;
	try {
		newwords = new TIndexOffU[sz / WORD_SIZE /* convert to words */];
	} catch(std::bad_alloc& ba) {
		if(errmsg != NULL) {
			// Output given error message
			std::cerr << errmsg;
		}
		throw 1;
	}
	if(oldsz > 0) {
		// Move old values into new array
		memcpy(newwords, words, oldsz >> 3 /* convert to bytes */);
	}
	// Initialize all new words to 0
	memset(newwords + (oldsz / WORD_SIZE /*convert to words*/), 0,
	       (sz - oldsz) >> 3 /* convert to bytes */);
	return newwords; // return new array
}

/**
 * A simple synchronized bitset class.
 */
class SyncBitset {

public:
	/**
	 * Allocate enough words to accommodate 'sz' bits.  Output the given
	 * error message and quit if allocation fails.
	 */
	SyncBitset(TIndexOffU sz, const char *errmsg = NULL) : _errmsg(errmsg) {
		TIndexOffU nwords = (sz / WORD_SIZE)+1; // divide by 32 and add 1
		try {
			_words = new TIndexOffU[nwords];
		} catch(std::bad_alloc& ba) {
			if(_errmsg != NULL) {
				std::cerr << _errmsg;
			}
			throw 1;
		}
		assert(_words != NULL);
		memset(_words, 0, nwords * OFF_SIZE /* words to bytes */);
		_sz = nwords * WORD_SIZE /* words to bits */;
	}

	/**
	 * Free memory for words.
	 */
	~SyncBitset() {
		delete[] _words;
	}

	/**
	 * Test whether the given bit is set in an unsynchronized manner.
	 */
	bool testUnsync(TIndexOffU i) {
		if(i < _sz) {
			return ((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) != 0;
		}
		return false;
	}

	/**
	 * Test whether the given bit is set in a synchronized manner.
	 */
	bool test(TIndexOffU i) {
		bool ret;
		ThreadSafe _ts(&mutex_m);
		ret = testUnsync(i);
		return ret;
	}

	/**
	 * Set a bit in the vector that hasn't been set before.  Assert if
	 * it has been set.  Uses synchronization.
	 */
	void set(TIndexOffU i) {
		ThreadSafe _ts(&mutex_m);
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(TIndexOffU oldsz = _sz);
			expand();
			assert_gt(_sz, oldsz);
		}
		// Fast path
		assert_lt(i, _sz);
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 0);
		_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}

	/**
	 * Set a bit in the vector that might have already been set.  Uses
	 * synchronization.
	 */
	void setOver(TIndexOffU i) {
		ThreadSafe _ts(&mutex_m);
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(TIndexOffU oldsz = _sz);
			expand();
			assert_gt(_sz, oldsz);
		}
		// Fast path
		assert_lt(i, _sz);
		_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}


private:

	/**
	 * Expand the size of the _words array by 50% to accommodate more
	 * bits.
	 */
	void expand() {
		TIndexOffU *newwords = bitsetRealloc(_sz, _words, _errmsg);
		delete[] _words;   // delete old array
		_words = newwords; // install new array
	}

	const char *_errmsg; // error message if an allocation fails
	TIndexOffU _sz;        // size as # of bits
	MUTEX_T mutex_m;       // mutex
	TIndexOffU *_words;    // storage
};

/**
 * A simple unsynchronized bitset class.
 */
class Bitset {

public:
	Bitset(TIndexOffU sz, const char *errmsg = NULL) : _errmsg(errmsg) {
		TIndexOffU nwords = (sz / WORD_SIZE)+1;
		try {
			_words = new TIndexOffU[nwords];
		} catch(std::bad_alloc& ba) {
			if(_errmsg != NULL) {
				std::cerr << _errmsg;
			}
			throw 1;
		}
		assert(_words != NULL);
		memset(_words, 0, nwords * OFF_SIZE);
		_sz = nwords * WORD_SIZE;
		_cnt = 0;
	}

	Bitset(const Bitset& o) : _words(NULL) {
		this->operator=(o);
	}

	~Bitset() {
		delete[] _words;
	}

	/**
	 * Test whether the given bit is set.
	 */
	bool test(TIndexOffU i) const {
		bool ret = false;
		if(i < _sz) {
			ret = ((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) != 0;
		}
		return ret;
	}

	/**
	 * Set a bit in the vector that hasn't been set before.  Assert if
	 * it has been set.
	 */
	void set(TIndexOffU i) {
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(TIndexOffU oldsz = _sz);
			expand();
			assert_gt(_sz, oldsz);
		}
		// Fast path
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 0);
		_cnt++;
		_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}

	/**
	 * Set a bit in the vector that might have already been set.
	 */
	void setOver(TIndexOffU i) {
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(TIndexOffU oldsz = _sz);
			expand();
			assert_gt(_sz, oldsz);
		}
		// Fast path
		if(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 0) _cnt++;
		_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}

	/**
	 * Unset all entries.  Don't adjust size.
	 */
	void clear() {
		for(size_t i = 0; i < ((_sz+BITSET_MASK) / WORD_SIZE); i++) {
			_words[i] = 0;
		}
		_cnt = 0;
	}

	/**
	 * Return the number of set bits.
	 */
	TIndexOffU count() const {
		return _cnt;
	}

	/**
	 * Return true iff no bits are set.
	 */
	bool empty() const {
		return _cnt == 0;
	}

	/**
	 * Deep copy from given Bitset to this one.
	 */
	Bitset& operator=(const Bitset& o) {
		_errmsg = o._errmsg;
		_sz = o._sz;
		_cnt = o._cnt;
		if(_words != NULL) delete[] _words;
		_words = new TIndexOffU[(_sz+BITSET_MASK) / WORD_SIZE];
		for(size_t i = 0; i < (_sz+BITSET_MASK) / WORD_SIZE; i++) {
			_words[i] = o._words[i];
		}
		return *this;
	}

private:

	/**
	 * Expand the size of the _words array by 50% to accommodate more
	 * bits.
	 */
	void expand() {
		TIndexOffU *newwords = bitsetRealloc(_sz, _words, _errmsg);
		delete[] _words;   // delete old array
		_words = newwords; // install new array
	}

	TIndexOffU _cnt;       // number of set bits
	const char *_errmsg; // error message if an allocation fails
	TIndexOffU _sz;        // size as # of bits
	TIndexOffU *_words;    // storage
};

/**
 * A simple fixed-length unsynchronized bitset class.
 */
template<int LEN>
class FixedBitset {

public:
	FixedBitset() : _cnt(0), _size(0) {
		memset(_words, 0, ((LEN / WORD_SIZE)+1) * OFF_SIZE);
	}

	/**
	 * Unset all bits.
	 */
	void clear() {
		memset(_words, 0, ((LEN / WORD_SIZE)+1) * OFF_SIZE);
	}

	/**
	 * Return true iff the bit at offset i has been set.
	 */
	bool test(TIndexOffU i) const {
		bool ret = false;
		assert_lt(i, LEN);
		ret = ((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) != 0;
		return ret;
	}

	/**
	 * Set the bit at offset i.  Assert if the bit was already set.
	 */
	void set(TIndexOffU i) {
		// Fast path
		assert_lt(i, LEN);
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 0);
		_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
		_cnt++;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}

	/**
	 * Set the bit at offset i.  Do not assert if the bit was already
	 * set.
	 */
	void setOver(TIndexOffU i) {
		// Fast path
		assert_lt(i, LEN);
		_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
		_cnt++;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}

	TIndexOffU count() const { return _cnt; }
	TIndexOffU size() const  { return _size; }

	/**
	 * Return true iff this FixedBitset has the same bits set as
	 * FixedBitset 'that'.
	 */
	bool operator== (const FixedBitset<LEN>& that) const {
		for(TIndexOffU i = 0; i < (LEN / WORD_SIZE)+1; i++) {
			if(_words[i] != that._words[i]) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return true iff this FixedBitset does not have the same bits set
	 * as FixedBitset 'that'.
	 */
	bool operator!= (const FixedBitset<LEN>& that) const {
		for(TIndexOffU i = 0; i < (LEN / WORD_SIZE)+1; i++) {
			if(_words[i] != that._words[i]) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Return a string-ized version of this FixedBitset.
	 */
	std::string str() const {
		std::ostringstream oss;
		for(int i = (int)size()-1; i >= 0; i--) {
			oss << (test(i)? "1" : "0");
		}
		return oss.str();
	}

private:
	TIndexOffU _cnt;
	TIndexOffU _size;
	TIndexOffU _words[(LEN / WORD_SIZE)+1]; // storage
};

/**
 * A simple fixed-length unsynchronized bitset class.
 */
class FixedBitset2 {

public:
	FixedBitset2(TIndexOffU len) : len_(len), _cnt(0), _size(0) {
		_words = new TIndexOffU[((len_ / WORD_SIZE)+1)];
		memset(_words, 0, ((len_ / WORD_SIZE)+1) * OFF_SIZE);
	}

	~FixedBitset2() { delete[] _words; }

	/**
	 * Unset all bits.
	 */
	void clear() {
		memset(_words, 0, ((len_ / WORD_SIZE)+1) * OFF_SIZE);
		_cnt = 0;
		_size = 0;
	}

	/**
	 * Return true iff the bit at offset i has been set.
	 */
	bool test(TIndexOffU i) const {
		bool ret = false;
		assert_lt(i, len_);
		ret = ((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) != 0;
		return ret;
	}

	/**
	 * Set the bit at offset i.  Assert if the bit was already set.
	 */
	void set(TIndexOffU i) {
		// Fast path
		assert_lt(i, len_);
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 0);
		_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
		_cnt++;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}

	/**
	 * Clear the bit at offset i.  Assert if the bit was not already set.
	 */
	void clear(TIndexOffU i) {
		// Fast path
		assert_lt(i, len_);
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
		_words[i / WORD_SIZE] &= ~((TIndexOffU)1 << (i & BITSET_MASK));
		_cnt--;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 0);
	}

	/**
	 * Set the bit at offset i.  Do not assert if the bit was already
	 * set.
	 */
	void setOver(TIndexOffU i) {
		// Fast path
		assert_lt(i, len_);
		if(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 0) {
			_words[i / WORD_SIZE] |= ((TIndexOffU)1 << (i & BITSET_MASK));
			_cnt++;
		}
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i / WORD_SIZE] >> (i & BITSET_MASK)) & 1) == 1);
	}

	TIndexOffU count() const { return _cnt; }
	TIndexOffU size() const  { return _size; }

	/**
	 * Return true iff this FixedBitset has the same bits set as
	 * FixedBitset 'that'.
	 */
	bool operator== (const FixedBitset2& that) const {
		for(TIndexOffU i = 0; i < (len_ / WORD_SIZE)+1; i++) {
			if(_words[i] != that._words[i]) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return true iff this FixedBitset does not have the same bits set
	 * as FixedBitset 'that'.
	 */
	bool operator!= (const FixedBitset2& that) const {
		for(TIndexOffU i = 0; i < (len_ / WORD_SIZE)+1; i++) {
			if(_words[i] != that._words[i]) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Return a string-ized version of this FixedBitset.
	 */
	std::string str() const {
		std::ostringstream oss;
		for(int i = (int)size()-1; i >= 0; i--) {
			oss << (test(i)? "1" : "0");
		}
		return oss.str();
	}

private:
	const TIndexOffU len_;
	TIndexOffU _cnt;
	TIndexOffU _size;
	TIndexOffU *_words; // storage
};

#endif /* BITSET_H_ */
