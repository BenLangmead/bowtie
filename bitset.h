#ifndef BITSET_H_
#define BITSET_H_

#include <stdint.h>
#include <string.h>
#include "assert_helpers.h"
#include "threading.h"

/**
 * A simple synchronized bitset class.
 */
class Bitset {

public:
	Bitset(size_t sz) {
		MUTEX_INIT(_lock);
		size_t nwords = (sz >> 5)+1;
		try {
			_words = new uint32_t[nwords];
		} catch(bad_alloc& ba) {
			cerr << "Could not allocate enough memory for the read mask; please subdivide reads and" << endl
			     << "run bowtie separately on each subset." << endl;
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
		} catch(bad_alloc& ba) {
			cerr << "Could not allocate enough memory for the read mask; please subdivide reads and" << endl
			     << "run bowtie separately on each subset." << endl;
			exit(1);
		}
		memcpy(newwords, _words, oldsz >> 3);
		delete[] _words;
		memset(newwords + (oldsz >> 5), 0, (oldsz >> 4));
		_words = newwords;
	}

	size_t _sz;       // size as # of bits
	MUTEX_T _lock;    // mutex
	uint32_t *_words; // storage
};

#endif /* BITSET_H_ */
