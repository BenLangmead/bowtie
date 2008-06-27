#ifndef BITPACK_H_
#define BITPACK_H_

#include "assert_helpers.h"

//#define BITPACK_MACROS

#ifdef BITPACK_MACROS

/**
 * Use macros instead of inline functions; I find this useful for
 * debugging since the asserts report the line number 
 */

#define pack_2b_in_8b(t,e,o) { \
	assert_lt((t), 4); \
	assert_lt((o), 4); \
	(e) |= ((t) << ((o)*2)); \
}

#define unpack_2b_from_8b(e,o) { \
	assert_lt((o), 4); \
	(((e) >> ((o)*2)) & 0xff); \
}

#else

static inline void pack_2b_in_8b(int two, uint8_t& eight, int off) {
	assert_lt(two, 4);
	assert_lt(off, 4);
	eight |= (two << (off*2));
}

static inline int unpack_2b_from_8b(const uint8_t eight, const int off) {
	assert_lt(off, 4);
	return ((eight >> (off*2)) & 0x3);
}

static inline void pack_2b_in_32b(int two, uint32_t& thirty2, int off) {
	assert_lt(two, 4);
	assert_lt(off, 16);
	thirty2 |= (two << (off*2));
}

static inline int unpack_2b_from_32b(const uint32_t thirty2, const int off) {
	assert_lt(off, 16);
	return ((thirty2 >> (off*2)) & 0x3);
}

static inline void pack_3b_in_8b(int three, uint8_t& eight, int off) {
	assert_lt(three, 8);
	assert_lt(off, 2);
	eight |= (three << (off*3));
}

static inline void pack_3b_in_32b(int three, uint32_t& thirty2, int off) {
	assert_lt(three, 8);
	assert_lt(off, 10);
	thirty2 |= (three << (off*3));
}
#endif /*BITPACK_MACROS*/

#endif /*BITPACK_H_*/
