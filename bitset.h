#ifndef BITSET_H_
#define BITSET_H_

#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string.h>
#include "assert_helpers.h"
#include "threading.h"

/**
 * A simple synchronized bitset class.
 */
class SyncBitset {

public:
	SyncBitset(size_t sz, const char *errmsg = NULL) : _errmsg(errmsg) {
		MUTEX_INIT(_lock);
		size_t nwords = (sz >> 5)+1;
		try {
			_words = new uint32_t[nwords];
		} catch(std::bad_alloc& ba) {
			if(_errmsg != NULL) {
				std::cerr << _errmsg;
			}
			exit(1);
		}
		memset(_words, 0, nwords * 4);
		_sz = nwords << 5;
	}

	~SyncBitset() {
		delete[] _words;
	}

	bool testUnsync(size_t i) {
		if(i < _sz) {
			return ((_words[i >> 5] >> (i & 0x1f)) & 1) != 0;
		}
		return false;
	}

	bool test(size_t i) {
		bool ret = false;
		MUTEX_LOCK(_lock);
		if(i < _sz) {
			ret = ((_words[i >> 5] >> (i & 0x1f)) & 1) != 0;
		}
		MUTEX_UNLOCK(_lock);
		return ret;
	}

	void set(size_t i) {
		MUTEX_LOCK(_lock);
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(size_t oldsz = _sz);
			expand();
			assert_eq(_sz, oldsz);
		}
		// Fast path
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0);
		_words[i >> 5] |= (1 << (i & 0x1f));
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
		MUTEX_UNLOCK(_lock);
	}

private:

	void expand() {
		size_t oldsz = _sz;
		_sz += (_sz>>1); // Add 50% more elements
		uint32_t *newwords;
		try {
			newwords = new uint32_t[_sz >> 5];
		} catch(std::bad_alloc& ba) {
			if(_errmsg != NULL) {
				std::cerr << _errmsg;
			}
			exit(1);
		}
		memcpy(newwords, _words, oldsz >> 3);
		delete[] _words;
		memset(newwords + (oldsz >> 5), 0, (oldsz >> 4));
		_words = newwords;
	}

	const char *_errmsg; // error message if an allocation fails
	size_t _sz;          // size as # of bits
	MUTEX_T _lock;       // mutex
	uint32_t *_words;    // storage
};

/**
 * A simple unsynchronized bitset class.
 */
class Bitset {

public:
	Bitset(size_t sz, const char *errmsg = NULL) : _errmsg(errmsg) {
		size_t nwords = (sz >> 5)+1;
		try {
			_words = new uint32_t[nwords];
		} catch(std::bad_alloc& ba) {
			if(_errmsg != NULL) {
				std::cerr << _errmsg;
			}
			exit(1);
		}
		memset(_words, 0, nwords * 4);
		_sz = nwords << 5;
	}

	~Bitset() {
		delete[] _words;
	}

	bool test(size_t i) {
		bool ret = false;
		if(i < _sz) {
			ret = ((_words[i >> 5] >> (i & 0x1f)) & 1) != 0;
		}
		return ret;
	}

	void set(size_t i) {
		while(i >= _sz) {
			// Slow path: bitset needs to be expanded before the
			// specified bit can be set
			ASSERT_ONLY(size_t oldsz = _sz);
			expand();
			assert_eq(_sz, oldsz);
		}
		// Fast path
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0);
		_words[i >> 5] |= (1 << (i & 0x1f));
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

private:

	void expand() {
		size_t oldsz = _sz;
		_sz += (_sz>>1); // Add 50% more elements
		uint32_t *newwords;
		try {
			newwords = new uint32_t[_sz >> 5];
		} catch(std::bad_alloc& ba) {
			if(_errmsg != NULL) {
				std::cerr << _errmsg;
			}
			exit(1);
		}
		memcpy(newwords, _words, oldsz >> 3);
		delete[] _words;
		memset(newwords + (oldsz >> 5), 0, (oldsz >> 4));
		_words = newwords;
	}

	const char *_errmsg; // error message if an allocation fails
	size_t _sz;          // size as # of bits
	uint32_t *_words;    // storage
};

/**
 * A simple fixed-length unsynchronized bitset class.
 */
template<int LEN>
class FixedBitset {

public:
	FixedBitset(const char *errmsg = NULL) : _errmsg(errmsg), _cnt(0), _size(0) {
		memset(_words, 0, ((LEN>>5)+1) * 4);
	}

	bool test(size_t i) const {
		bool ret = false;
		assert_lt(i, LEN);
		ret = ((_words[i >> 5] >> (i & 0x1f)) & 1) != 0;
		return ret;
	}

	void set(size_t i) {
		// Fast path
		assert_lt(i, LEN);
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 0);
		_words[i >> 5] |= (1 << (i & 0x1f));
		_cnt++;
		if(i >= _size) {
			_size = i+1;
		}
		assert(((_words[i >> 5] >> (i & 0x1f)) & 1) == 1);
	}

	size_t count() const { return _cnt; }
	size_t size() const  { return _size; }

	bool operator== (const FixedBitset<LEN>& that) const {
		for(size_t i = 0; i < (LEN>>5)+1; i++) {
			if(_words[i] != that._words[i]) {
				return false;
			}
		}
		return true;
	}

	bool operator!= (const FixedBitset<LEN>& that) const {
		for(size_t i = 0; i < (LEN>>5)+1; i++) {
			if(_words[i] != that._words[i]) {
				return true;
			}
		}
		return false;
	}

	std::string str() const {
		std::ostringstream oss;
		for(int i = (int)size()-1; i >= 0; i--) {
			oss << (test(i)? "1" : "0");
		}
		return oss.str();
	}

private:
	const char *_errmsg;         // error message if an allocation fails
	size_t _cnt;
	size_t _size;
	uint32_t _words[(LEN>>5)+1]; // storage
};

#endif /* BITSET_H_ */
