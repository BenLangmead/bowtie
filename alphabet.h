#ifndef ALPHABETS_H_
#define ALPHABETS_H_

#include <stdexcept>
#include <string>
#include <seqan/sequence.h>
#include <sstream>
#include "assert_helpers.h"

using namespace std;
using namespace seqan;

/**
 * Helper function to print a uint32_t as a DNA string where each 2-bit
 * stretch is a character and more significiant bits appear to the left
 * of less singificant bits.
 */
static inline std::string u32ToDna(uint32_t a, int len) {
	char buf[17]; // TODO: return a new string; by value I guess
	assert_leq(len, 16);
	for(int i = 0; i < len; i++) {
		buf[len-i-1] = "ACGT"[a & 3];
		a >>= 2;
	}
	buf[len] = '\0';
	return std::string(buf);
}

/**
 * Return a new TStr containing the reverse-complement of s.  Ns go to
 * Ns.
 */
template<typename TStr>
static inline TStr reverseComplement(const TStr& s) {
	typedef typename Value<TStr>::Type TVal;
	TStr s_rc;
	size_t slen = length(s);
	resize(s_rc, slen);
	for(size_t i = 0; i < slen; i++) {
		int sv = (int)s[slen-i-1];
		if(sv == 4) {
			s_rc[i] = (TVal)4;
		} else {
			s_rc[i] = (TVal)(sv ^ 3);
		}
	}
	return s_rc;
}

/**
 * Reverse-complement s in-place.  Ns go to Ns.
 */
template<typename TStr>
static inline void reverseComplementInPlace(TStr& s) {
	typedef typename Value<TStr>::Type TVal;
	size_t len = length(s);
	for(size_t i = 0; i < (len>>1); i++) {
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
}

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
 * Return the reverse-complement of s.
 */
static inline bool isReverseComplement(const String<Dna5>& s1,
                                       const String<Dna5>& s2)
{
	if(length(s1) != length(s2)) return false;
	size_t slen = length(s1);
	for(size_t i = 0; i < slen; i++) {
		int i1 = (int)s1[i];
		int i2 = (int)s2[slen - i - 1];
		if(i1 == 4) {
			if(i2 != 4) return false;
		}
		else if(i1 != (i2 ^ 3)) return false;
	}
	return true;
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


extern uint8_t dna4Cat[];
extern uint8_t charToDna5[];
extern uint8_t rcCharToDna5[];

#endif /*ALPHABETS_H_*/
