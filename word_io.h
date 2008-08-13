#ifndef WORD_IO_H_
#define WORD_IO_H_

#include <stdint.h>
#include <iostream>
#include <fstream>
#include "assert_helpers.h"
#include "endian_swap.h"

/**
 * Write a 32-bit unsigned to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
static inline void writeU32(std::ostream& out, uint32_t x, bool toBigEndian) {
	uint32_t y = endianizeU32(x, toBigEndian);
	out.write((const char*)&y, 4);
}

/**
 * Write a 32-bit unsigned to an output stream using the native
 * endianness.
 */
static inline void writeU32(std::ostream& out, uint32_t x) {
	out.write((const char*)&x, 4);
}

/**
 * Write a 32-bit signed int to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
static inline void writeI32(std::ostream& out, int32_t x, bool toBigEndian) {
	int32_t y = endianizeI32(x, toBigEndian);
	out.write((const char*)&y, 4);
}

/**
 * Write a 32-bit unsigned to an output stream using the native
 * endianness.
 */
static inline void writeI32(std::ostream& out, int32_t x) {
	out.write((const char*)&x, 4);
}

/**
 * Read a 32-bit unsigned from an input stream, inverting endianness
 * if necessary.
 */
static inline uint32_t readU32(std::istream& in, bool toBigEndian) {
	uint32_t x;
	in.read((char *)&x, 4);
	assert_eq(4, in.gcount());
	return endianizeU32(x, toBigEndian);
}

/**
 * Read a 32-bit signed from an input stream, inverting endianness
 * if necessary.
 */
static inline int32_t readI32(std::istream& in, bool toBigEndian) {
	int32_t x;
	in.read((char *)&x, 4);
	assert_eq(4, in.gcount());
	return endianizeI32(x, toBigEndian);
}

#endif /*WORD_IO_H_*/
