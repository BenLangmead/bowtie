#ifndef ALPHABETS_H_
#define ALPHABETS_H_

#include <stdexcept>
#include <string>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <sstream>
#include "assert_helpers.h"

using namespace std;
using namespace seqan;

/// Reverse a string in-place
template <typename TStr>
static inline void reverseInPlace(TStr& s) {
	typedef typename Value<TStr>::Type TVal;
	size_t len = length(s);
	for(size_t i = 0; i < (len>>1); i++) {
		TVal tmp = s[i];
		s[i] = s[len-i-1];
		s[len-i-1] = tmp;
	}
}

/**
 * Return a new TStr containing the reverse-complement of s.  Ns go to
 * Ns.
 */
template<typename TStr>
static inline TStr reverseComplement(const TStr& s, bool color) {
	typedef typename Value<TStr>::Type TVal;
	TStr s_rc;
	size_t slen = length(s);
	resize(s_rc, slen);
	if(color) {
		for(size_t i = 0; i < slen; i++) {
			s_rc[i] = s[slen-i-1];
		}
	} else {
		for(size_t i = 0; i < slen; i++) {
			int sv = (int)s[slen-i-1];
			if(sv == 4) {
				s_rc[i] = (TVal)4;
			} else {
				s_rc[i] = (TVal)(sv ^ 3);
			}
		}
	}
	return s_rc;
}

/**
 * Reverse-complement s in-place.  Ns go to Ns.
 */
template<typename TStr>
static inline void reverseComplementInPlace(TStr& s, bool color) {
	typedef typename Value<TStr>::Type TVal;
	if(color) {
		reverseInPlace(s);
		return;
	}
	size_t len = length(s);
	size_t i;
	for(i = 0; i < (len>>1); i++) {
		int sv = (int)s[len-i-1];
		int sf = (int)s[i];
		if(sv == 4) {
			s[i] = (TVal)4;
		} else {
			s[i] = (TVal)(sv ^ 3);
		}
		if(sf == 4)  {
			s[len-i-1] = (TVal)4;
		} else {
			s[len-i-1] = (TVal)(sf ^ 3);
		}
	}
	if((len & 1) != 0 && (int)s[len >> 1] != 4) {
		s[len >> 1] = (TVal)((int)s[len >> 1] ^ 3);
	}
}

/**
 * Return the reverse-complement of s.
 */
template<typename TStr>
static inline TStr reverseCopy(const TStr& s) {
	typedef typename Value<TStr>::Type TVal;
	TStr s_rc;
	size_t slen = length(s);
	resize(s_rc, slen);
	for(size_t i = 0; i < slen; i++) {
		s_rc[i] = (TVal)((int)s[slen-i-1]);
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
	typedef typename Value<TStr>::Type TVal;
	std::string ret;
	size_t len = seqan::length(l);
	for(size_t i = off; i < len; i++) {
		ret.push_back((char)(TVal)l[i]);
	}
	return ret;
}

/**
 * Calculate the entropy of the given read.  Handle Ns by charging them
 * to the most frequent non-N character.
 */
static inline float entropyDna5(const String<Dna5>& read) {
	size_t cs[5] = {0, 0, 0, 0, 0};
	size_t readLen = seqan::length(read);
	for(size_t i = 0; i < readLen; i++) {
		int c = (int)read[i];
		assert_lt(c, 5);
		assert_geq(c, 0);
		cs[c]++;
	}
	if(cs[4] > 0) {
		// Charge the Ns to the non-N character with maximal count and
		// then exclude them from the entropy calculation (i.e.,
		// penalize Ns as much as possible)
		if(cs[0] >= cs[1] && cs[0] >= cs[2] && cs[0] >= cs[3]) {
			// Charge Ns to As
			cs[0] += cs[4];
		} else if(cs[1] >= cs[2] && cs[1] >= cs[3]) {
			// Charge Ns to Cs
			cs[1] += cs[4];
		} else if(cs[2] >= cs[3]) {
			// Charge Ns to Gs
			cs[2] += cs[4];
		} else {
			// Charge Ns to Ts
			cs[3] += cs[4];
		}
	}
	float ent = 0.0;
	for(int i = 0; i < 4; i++) {
		if(cs[i] > 0) {
			float frac = (float)cs[i] / (float)readLen;
			ent += (frac * log(frac));
		}
	}
	ent = -ent;
	assert_geq(ent, 0.0);
	return ent;
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

/// Convert a pair of 2-bit (and 4=N) encoded DNA bases to a color
extern uint8_t dinuc2color[5][5];

#endif /*ALPHABETS_H_*/
