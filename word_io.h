#ifndef WORD_IO_H_
#define WORD_IO_H_

#include <stdint.h>
#include <iostream>
#include <fstream>
#include "assert_helpers.h"
#include "endian_swap.h"
#include "btypes.h"

/**
 * Write a 32/64 bit unsigned to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
template <typename T>
static inline void writeU(std::ostream& out, T x, bool toBigEndian) {
	T y = endianizeU<T>(x, toBigEndian);
	out.write((const char*)&y, sizeof(T));
}

/**
 * Write a 32/64 bit unsigned to an output stream using the native
 * endianness.
 */
template <typename T>
static inline void writeU(std::ostream& out, T x) {
	out.write((const char*)&x, sizeof(T));
}

/**
 * Write a 32/64 bit signed int to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
template <typename T>
static inline void writeI(std::ostream& out, T x, bool toBigEndian) {
	T y = endianizeI<T>(x, toBigEndian);
	out.write((const char*)&y, sizeof(T));
}

/**
 * Write a 32/64 bit unsigned to an output stream using the native
 * endianness.
 */
template <typename T>
static inline void writeI(std::ostream& out, T x) {
	out.write((const char*)&x, sizeof(T));
}

template <typename T>
static inline T readU(std::istream& in, bool swap) {
	T x;
	in.read((char *)&x, sizeof(T));
	assert_eq(sizeof(T), in.gcount());
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapU32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapU64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}


template <typename T>
static inline T readU(FILE* in, bool swap) {
	T x;
	if(fread((void *)&x, 1, sizeof(T), in) != sizeof(T)) {
		assert(false);
	}
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapU32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapU64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}

template <typename T>
static inline T readI(std::istream& in, bool swap) {
	T x;
	in.read((char *)&x, sizeof(T));
	assert_eq(sizeof(T), in.gcount());
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapI32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapI64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}

template <typename T>
static inline T readI(FILE* in, bool swap) {
	T x;
	if(fread((void *)&x, 1, sizeof(T), in) != sizeof(T)) {
		assert(false);
	}
	if(swap) {
		if(sizeof(T) == 4) {
			return endianSwapI32(x);
		} else if(sizeof(T) == 8) {
			return endianSwapI64(x);
		} else {
			assert(false);
		}
	} else {
		return x;
	}
}


#endif /*WORD_IO_H_*/
