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

extern uint8_t dna4Cat[];
extern uint8_t charToDna5[];
extern uint8_t rcCharToDna5[];

#endif /*ALPHABETS_H_*/
