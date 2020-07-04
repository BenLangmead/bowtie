#ifndef ALPHABETS_H_
#define ALPHABETS_H_

#include <math.h>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <sstream>

#include "assert_helpers.h"

using namespace std;

/// Reverse a string in-place
template <typename TStr>
static inline void reverseInPlace(TStr& s) {
	size_t len = s.length();
	for(size_t i = 0; i < (len>>1); i++) {
		char tmp = s[i];
		s[i] = s[len-i-1];
		s[len-i-1] = tmp;
	}
}

/**
 * Return a new TStr containing the reverse-complement of s.  Ns go to
 * Ns.
 */
template<typename TStr>
static inline TStr reverseComplement(const TStr& s) {
	TStr s_rc;
	size_t slen = length(s);
	resize(s_rc, slen);
	for(size_t i = 0; i < slen; i++) {
		int sv = (int)s[slen-i-1];
		if(sv == 4) {
			s_rc[i] = (char)4;
		} else {
			s_rc[i] = (char)(sv ^ 3);
		}
	}
	return s_rc;
}

/**
 * Reverse-complement s in-place.  Ns go to Ns.
 */
template<typename TStr>
static inline void reverseComplementInPlace(TStr& s) {
	size_t len = s.length();
	size_t i;
	for(i = 0; i < (len>>1); i++) {
		int sv = (int)s[len-i-1];
		int sf = (int)s[i];
		if(sv == 4) {
			s[i] = (char)4;
		} else {
			s[i] = (char)(sv ^ 3);
		}
		if(sf == 4)  {
			s[len-i-1] = (char)4;
		} else {
			s[len-i-1] = (char)(sf ^ 3);
		}
	}
	if((len & 1) != 0 && (int)s[len >> 1] != 4) {
		s[len >> 1] = (char)((int)s[len >> 1] ^ 3);
	}
}

/**
 * Return the reverse-complement of s.
 */
template<typename TStr>
static inline TStr reverseCopy(const TStr& s) {
	TStr s_rc;
	size_t slen = length(s);
	resize(s_rc, slen);
	for(size_t i = 0; i < slen; i++) {
		s_rc[i] = (char)((int)s[slen-i-1]);
	}
	return s_rc;
}

/**
 * Return true iff the first string is dollar-less-than the second.
 * This means that we pretend that a 'dollar sign' character,
 * lexicographically larger than all other characters, exists at the
 * end of both strings.
 */
template <typename TStr>
static inline bool
dollarLt(const TStr& l, const TStr& r) {
	return isPrefix(r, l) || (l < r && !isPrefix(l, r));
}

/**
 * Return true iff the first string is dollar-greater-than the second.
 * This means that we pretend that a 'dollar sign' character,
 * lexicographically larger than all other characters, exists at the
 * end of both strings.
 */
template <typename TStr>
static inline bool
dollarGt(const TStr& l, const TStr& r) {
	return !dollarLt(l, r);
}

/**
 * Return a copy of the suffix of l starting at 'off'.
 */
template <typename TStr>
static inline std::string
suffixStr(const TStr& l, size_t off) {
	std::string ret;
	size_t len = l.length();
	for(size_t i = off; i < len; i++) {
		ret.push_back(l[i]);
	}
	return ret;
}

/**
 * Return the DNA complement of the given ASCII char.
 */
static inline char comp(char c) {
	switch(c) {
	case 'a': return 't';
	case 'A': return 'T';
	case 'c': return 'g';
	case 'C': return 'G';
	case 'g': return 'c';
	case 'G': return 'C';
	case 't': return 'a';
	case 'T': return 'A';
	default: return c;
	}
}

extern uint8_t dna4Cat[];
extern uint8_t charToDna5[];

/// Convert an ascii char to a 2-bit base: 0=A, 1=C, 2=G, 3=T, 4=N
extern uint8_t asc2dna[];

/// Convert an ascii char to a 2-bit base: 0=A, 1=C, 2=G, 3=T, 4=N
extern uint8_t asc2col[];
extern uint8_t rcCharToDna5[];

/// Convert an ascii char to a DNA category.  Categories are:
/// 0 -> invalid
/// 1 -> unambiguous a, c, g or t
/// 2 -> ambiguous
/// 3 -> unmatchable
extern uint8_t asc2dnacat[];

/// Convert an ascii char to a color category.  Categories are:
/// 0 -> invalid
/// 1 -> unambiguous 0, 1, 2 or 3
/// 2 -> ambiguous (not applicable for colors)
/// 3 -> unmatchable
extern uint8_t asc2colcat[];

/// Convert a 2-bit nucleotide (and 4=N) and a color to the
/// corresponding 2-bit nucleotide
extern uint8_t nuccol2nuc[5][5];

/// Convert ambiguous ASCII nuceleotide to mask
extern uint8_t asc2dnamask[];

/// Convert a 4-bit mask into an IUPAC code
extern signed char mask2iupac[16];

/**
 * Return true iff c is an unambiguous Dna character.
 */
static inline bool isUnambigDna(char c) {
	return asc2dnacat[(int)c] == 1;
}

/**
 * Return true iff c is a Dna character.
 */
static inline bool isDna(char c) {
	return asc2dnacat[(int)c] > 0;
}

/**
 * Return true iff c is an unambiguous color character (0,1,2,3).
 */
static inline bool isUnambigColor(char c) {
	return asc2colcat[(int)c] == 1;
}

/**
 * Return true iff c is a color character.
 */
static inline bool isColor(char c) {
	return asc2colcat[(int)c] > 0;
}

/// Convert bit encoded DNA char to its complement
extern int dnacomp[5];

/// String of all DNA and IUPAC characters
extern const char *iupacs;

/**
 * Return the reverse complement of a bit-encoded nucleotide.
 */
static inline int compDna(int c) {
	assert_leq(c, 4);
	return dnacomp[c];
}

/// Convert a pair of 2-bit (and 4=N) encoded DNA bases to a color
extern uint8_t dinuc2color[5][5];

/// Map from masks to their reverse-complement masks
extern int maskcomp[16];

#endif /*ALPHABETS_H_*/
