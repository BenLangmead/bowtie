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
 * Exception to throw when a Fasta file is malformed.
 */
class AlphabetException : public runtime_error {
public:
	AlphabetException(const string& msg = "") : runtime_error(msg) {}
};

/**
 * Convert an int to an 'A', 'C', 'G' or 'T'
 */
static inline char toDna(int i) {
	switch(i) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: {
			stringstream ss; ss << "No such DNA character: " << i;
			throw AlphabetException(ss.str());
		}
	}
	throw;
}

/**
 * Convert an int to an 'A', 'C', 'G', 'T' or '$'
 */
static inline char toDna5(int i) {
	switch(i) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		case 4: return '$';
		default: {
			stringstream ss; ss << "No such DNA+$ character: " << i;
			throw AlphabetException(ss.str());
		}
	}
	throw;
}

/**
 * Helper function to print a uint32_t as a DNA string where each 2-bit
 * stretch is a character and more significiant bits appear to the left
 * of less singificant bits. 
 */
static inline char * u32ToDna(uint32_t a, int len) {
	static char buf[17]; // TODO: return a new string; by value I guess
	assert_leq(len, 16);
	for(int i = 0; i < len; i++) {
		buf[len-i-1] = toDna(a & 3);
		a >>= 2;
	}
	buf[len] = '\0';
	return buf;
}

/**
 * Return the complement of Dna (or Rna) char c.
 */
static inline int complement(int c) {
	return (c < 4)? c ^ 3 : c;
}

/**
 * Return the reverse-complement of s.
 */
template<typename TStr>
static inline TStr reverseComplement(const TStr& s) {
	typedef typename Value<TStr>::Type TVal;
	TStr s_rc;
	size_t slen = length(s);
	resize(s_rc, slen);
	for(size_t i = 0; i < slen; i++) {
		s_rc[i] = (TVal)complement((int)s[slen-i-1]);
	}
	return s_rc;
}

/**
 * Return the reverse-complement of s.
 */
template<typename TStr>
static inline TStr reverse(const TStr& s) {
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
static inline bool isReverseComplement(const String<Dna>& s1,
                                       const String<Dna>& s2)
{
	if(length(s1) != length(s2)) return false;
	size_t slen = length(s1);
	for(size_t i = 0; i < slen; i++) {
		int i1 = (int)s1[i];
		int i2 = (int)s2[slen - i - 1];
		if(i1 != (i2 ^ 3)) return false;
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

#endif /*ALPHABETS_H_*/
