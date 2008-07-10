#ifndef MULTIKEY_QSORT_H_
#define MULTIKEY_QSORT_H_

#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include "sequence_io.h"
#include "alphabet.h"
#include "assert_helpers.h"
#include "diff_sample.h"

using namespace std;
using namespace seqan;

/**
 * Swap elements a and b in seqan::String s
 */
template <typename TStr, typename TPos>
static inline void swap(TStr& s, size_t slen, TPos a, TPos b) {
	typedef typename Value<TStr>::Type TAlphabet;
	assert_lt(a, slen);
	assert_lt(b, slen);
	TAlphabet tmp = s[a];
	s[a] = s[b];
	s[b] = tmp;
}

/**
 * Swap elements a and b in array s
 */
template <typename TVal, typename TPos>
static inline void swap(TVal* s, size_t slen, TPos a, TPos b) {
	assert_lt(a, slen);
	assert_lt(b, slen);
	TVal tmp = s[a];
	s[a] = s[b];
	s[b] = tmp;
}

/**
 * Helper macro for swapping elements a and b in seqan::String s.  Does
 * some additional sainty checking w/r/t begin and end (which are
 * parameters to the sorting routines below).
 */
#define SWAP(s, a, b) { \
	assert_geq(a, begin); \
	assert_geq(b, begin); \
	assert_lt(a, end); \
	assert_lt(b, end); \
	swap(s, slen, a, b); \
}

/**
 * Helper macro for swapping the same pair of elements a and b in
 * two different seqan::Strings s and s2.  This is a helpful variant
 * if, for example, the caller would like to see how their input was
 * permuted by the sort routine (in that case, the caller would let s2
 * be an array s2[] where s2 is the same length as s and s2[i] = i). 
 */
#define SWAP2(s, s2, a, b) { \
	SWAP(s, a, b); \
	swap(s2, slen, a, b); \
}

#define SWAP1(s, s2, a, b) { \
	SWAP(s, a, b); \
}

/**
 * Helper macro that swaps a range of elements [i, i+n) with another
 * range [j, j+n) in seqan::String s.
 */
#define VECSWAP(s, i, j, n) { \
	if(n > 0) { vecswap(s, slen, i, j, n, begin, end); } \
}

/**
 * Helper macro that swaps a range of elements [i, i+n) with another
 * range [j, j+n) both in seqan::String s and seqan::String s2.
 */
#define VECSWAP2(s, s2, i, j, n) { \
	if(n > 0) { vecswap2(s, slen, s2, i, j, n, begin, end); } \
}

/**
 * Helper function that swaps a range of elements [i, i+n) with another
 * range [j, j+n) in seqan::String s.  begin and end represent the
 * current range under consideration by the caller (one of the
 * recursive multikey_quicksort routines below).
 */
template <typename TStr, typename TPos>
static inline void vecswap(TStr& s, size_t slen, TPos i, TPos j, TPos n, TPos begin, TPos end) {
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
	}
}

template <typename TVal, typename TPos>
static inline void vecswap(TVal *s, size_t slen, TPos i, TPos j, TPos n, TPos begin, TPos end) {
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
	}
}

/**
 * Helper function that swaps a range of elements [i, i+n) with another
 * range [j, j+n) both in seqan::String s and seqan::String s2.  begin
 * and end represent the current range under consideration by the
 * caller (one of the recursive multikey_quicksort routines below).
 */
template <typename TStr, typename TPos>
static inline void vecswap2(TStr& s, size_t slen, TStr& s2, TPos i, TPos j, TPos n, TPos begin, TPos end) {
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
		swap(s2, slen, a, b);
	}
}

template <typename TVal, typename TPos>
static inline void vecswap2(TVal* s, size_t slen, TVal* s2, TPos i, TPos j, TPos n, TPos begin, TPos end) {
	assert_geq(i, begin);
	assert_geq(j, begin);
	assert_lt(i, end);
	assert_lt(j, end);
	while(n-- > 0) {
		assert_geq(n, 0);
		TPos a = i+n;
		TPos b = j+n;
		assert_geq(a, begin);
		assert_geq(b, begin);
		assert_lt(a, end);
		assert_lt(b, end);
		swap(s, slen, a, b);
		swap(s2, slen, a, b);
	}
}

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.
#define CHAR_AT(ss, aa) ((length(s[ss]) > aa) ? (int)(Dna)(s[ss][aa]) : hi)

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.
#define CHAR_AT_SUF(si, off) (((off+s[si]) < hlen) ? ((int)(Dna)(host[off+s[si]])) : (hi))

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.
#ifndef PACKED_STRINGS
#define CHAR_AT_SUF_U8(si, off) (((off+s[si]) < hlen) ? (int)host[off+s[si]] : (hi))
#else
#define CHAR_AT_SUF_U8(si, off) (((off+s[si]) < hlen) ? (int)(Dna)host[off+s[si]] : (hi))
#endif

#define CHOOSE_AND_SWAP_RANDOM_PIVOT(sw, ch) {                            \
	/* Note: rand() didn't really cut it here; it seemed to run out of */ \
	/* randomness and, after a time, returned the same thing over and */  \
	/* over again */                                                      \
	a = (random() % n) + begin; /* choose pivot between begin and end */  \
	assert_lt(a, end); assert_geq(a, begin);                              \
	sw(s, s2, begin, a); /* move pivot to beginning */                    \
}

/**
 * Ad-hoc DNA-centric way of choose a pretty good pivot without using
 * the pseudo-random number generator.  We try to get a 1 or 2 if
 * possible, since they'll split things more evenly than a 0 or 4.  We
 * also avoid swapping in the event that we choose the first element.
 */
#define CHOOSE_AND_SWAP_SMART_PIVOT(sw, ch) {                                  \
	a = begin; /* choose first elt */                                          \
	/* now try to find a better elt */                                         \
	if(n >= 5) {                                                               \
		if     (ch(begin+1, depth) == 1 || ch(begin+1, depth) == 2) a = begin+1; \
		else if(ch(begin+2, depth) == 1 || ch(begin+2, depth) == 2) a = begin+2; \
		else if(ch(begin+3, depth) == 1 || ch(begin+3, depth) == 2) a = begin+3; \
		else if(ch(begin+4, depth) == 1 || ch(begin+4, depth) == 2) a = begin+4; \
		if(a != begin) sw(s, s2, begin, a); /* move pivot to beginning */      \
	}                                                                          \
}

#define CHOOSE_AND_SWAP_PIVOT CHOOSE_AND_SWAP_SMART_PIVOT

#ifndef NDEBUG
/**
 * Assert that the range of chars at depth 'depth' in strings 'begin'
 * to 'end' in string-of-strings s is parititioned properly according
 * to the ternary paritioning strategy of Bentley and McIlroy (*prior
 * to* swapping the = regions to the center)
 */
template<typename TStr>
bool assertPartitioned(TStr& s, int hi, int pivot, size_t begin, size_t end, size_t depth) {
	int state = 0; // 0 -> 1st = section, 1 -> < section, 2 -> > section, 3 -> 2nd = section
	for(size_t i = begin; i < end; i++) {
		switch(state) {
		case 0:
			if       (CHAR_AT(i, depth) < pivot)  { state = 1; break; }
			else if  (CHAR_AT(i, depth) > pivot)  { state = 2; break; }
			assert_eq(CHAR_AT(i, depth), pivot);  break;
		case 1:
			if       (CHAR_AT(i, depth) > pivot)  { state = 2; break; }
			else if  (CHAR_AT(i, depth) == pivot) { state = 3; break; }
			assert_lt(CHAR_AT(i, depth), pivot);  break;
		case 2:
			if       (CHAR_AT(i, depth) == pivot) { state = 3; break; }
			assert_gt(CHAR_AT(i, depth), pivot);	 break;
		case 3:
			assert_eq(CHAR_AT(i, depth), pivot);	 break;
		}
	}
	return true;
}

/**
 * Assert that the range of chars at depth 'depth'�in strings 'begin'
 * to 'end' in string-of-strings s is parititioned properly according
 * to the ternary paritioning strategy of Bentley and McIlroy (*after*
 * swapping the = regions to the center)
 */
template<typename TStr>
bool assertPartitioned2(TStr& s, int hi, int pivot, size_t begin, size_t end, size_t depth) {
	int state = 0; // 0 -> < section, 1 -> = section, 2 -> > section
	for(size_t i = begin; i < end; i++) {
		switch(state) {
		case 0:
			if       (CHAR_AT(i, depth) == pivot) { state = 1; break; }
			else if  (CHAR_AT(i, depth) > pivot)  { state = 2; break; }
			assert_lt(CHAR_AT(i, depth), pivot);  break;
		case 1:
			if       (CHAR_AT(i, depth) > pivot)  { state = 2; break; }
			assert_eq(CHAR_AT(i, depth), pivot);  break;
		case 2:
			assert_gt(CHAR_AT(i, depth), pivot);  break;
		}
	}
	return true;
}
#endif

/**
 * Simple helper to print a list of Strings.
 */
template <typename T>
void printStringList(const T& ss, ostream& out) {
	for(size_t i = 0; i < length(ss); i++) {
		out << "  " << ss[i] << endl;
	}
}

/**
 * Driver for the 4-arg version of multikey_qsort.  Optionally sanity-
 * checks the sorted result before returning.
 */
template<typename TStr>
void mkeyQSort(TStr& s, int hi, bool verbose = false, bool sanityCheck = false) {
	mkeyQSort(s, hi, 0, length(s), 0);
	if(sanityCheck) {
		for(size_t i = 0; i < length(s)-1; i++) {
			if(s[i] > s[i+1]) {
				// operator > treats shorter strings as
				// lexicographically smaller, but we want to opposite
				assert(isPrefix(s[i+1], s[i]));
			}
		}
	}
}

/**
 * Main multikey quicksort function.  Based on Bentley & Sedgewick's
 * algorithm on p.5 of their paper "Fast Algorithms for Sorting and
 * Searching Strings".
 */
template<typename TStr>
void mkeyQSort(TStr& s, int hi, size_t begin, size_t end, size_t depth) {
	size_t slen = length(s);
	typedef typename Value<TStr>::Type TAlphabet;
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		mkeyQSort(s, hi, nbegin, nend, ndepth); \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t a, b, c, d, e, r;
	size_t n = end - begin;
	if(n <= 1) return;                // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT);
	int v = CHAR_AT(begin, depth); // v <- randomly-selected pivot value
	a = b = begin;
	c = d = e = end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc;
		while(b <= c && v >= (bc = CHAR_AT(b, depth))) {
			if(v == bc) {
				SWAP(s, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc;
		while(b <= c && v <= (cc = CHAR_AT(c, depth))) {
			if(v == cc) {
				SWAP(s, c, d); d--; e--;
			}
			else if(c == e && v == hi) e--;
			c--;
		}
		if(b > c) break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(e-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	assert(assertPartitioned(s, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	assert(assertPartitioned2(s, hi, v, begin, end, depth)); // check post-=-swap invariant
	r = b-a; // r <- # of <'s
	MQS_RECURSE(begin, begin + r, depth); // recurse on <'s
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted 
	if(v != hi) {
		MQS_RECURSE(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = e-c;                        // r <- # of >'s excluding those exhausted
	MQS_RECURSE(end-r, end, depth); // recurse on >'s
}

/**
 * Simple helper to print a list of suffixes with indices.
 */
template<typename THost, typename TElt1, typename TElt2>
void printSuffixList(const THost& host,
                     const String<TElt1> s,
                     const String<TElt2> idxs,
                     ostream& out)
{
	for(size_t i = 0; i < length(s); i++) {
		out << "  " << suffix(host, s[i]) << " (" << idxs[i] << ")" << endl;
	}
}

/**
 * Simple helper to print a list of suffixes with indices.
 */
template<typename THost>
void printSuffixList(const THost& host,
                     uint32_t *s,
                     size_t slen,
                     uint32_t *idxs,
                     ostream& out)
{
	for(size_t i = 0; i < slen; i++) {
		out << "  " << suffix(host, s[i]) << " (" << idxs[i] << ")" << endl;
	}
}

/**
 * Simple helper to print a list of suffixes.
 */
template<typename THost, typename TElt1>
void printSuffixList(const THost& host,
                     const String<TElt1> s,
                     ostream& out)
{
	for(size_t i = 0; i < length(s); i++) {
		out << "  " << suffix(host, s[i]) << endl;
	}
}

/**
 * Simple helper to print a list of suffixes.
 */
template<typename THost>
void printSuffixList(const THost& host,
                     uint32_t *s,
                     size_t slen,
                     ostream& out)
{
	for(size_t i = 0; i < slen; i++) {
		out << "  " << suffix(host, s[i]) << endl;
	}
}

#ifndef NDEBUG
/**
 * Assert that the range of chars at depth 'depth'�in strings 'begin'
 * to 'end' in string-of-suffix-offsets s is parititioned properly
 * according to the ternary paritioning strategy of Bentley and McIlroy
 * (*prior to* swapping the = regions to the center)
 */
template<typename THost>
bool assertPartitionedSuf(const THost& host,
                          uint32_t *s,
                          size_t slen,
                          int hi,
                          int pivot,
                          size_t begin,
                          size_t end,
                          size_t depth)
{
	size_t hlen = length(host);
	int state = 0; // 0 -> 1st = section, 1 -> < section, 2 -> > section, 3 -> 2nd = section
	for(size_t i = begin; i < end; i++) {
		switch(state) {
		case 0:
			if       (CHAR_AT_SUF(i, depth) < pivot)  { state = 1; break; }
			else if  (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			assert_eq(CHAR_AT_SUF(i, depth), pivot);  break;
		case 1:
			if       (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			else if  (CHAR_AT_SUF(i, depth) == pivot) { state = 3; break; }
			assert_lt(CHAR_AT_SUF(i, depth), pivot);  break;
		case 2:
			if       (CHAR_AT_SUF(i, depth) == pivot) { state = 3; break; }
			assert_gt(CHAR_AT_SUF(i, depth), pivot);	 break;
		case 3:
			assert_eq(CHAR_AT_SUF(i, depth), pivot);	 break;
		}
	}
	return true;
}

/**
 * Assert that the range of chars at depth 'depth'�in strings 'begin'
 * to 'end' in string-of-suffix-offsets s is parititioned properly
 * according to the ternary paritioning strategy of Bentley and McIlroy
 * (*after* swapping the = regions to the center)
 */
template<typename THost>
bool assertPartitionedSuf2(const THost& host,
                           uint32_t *s,
                           size_t slen,
                           int hi,
                           int pivot,
                           size_t begin,
                           size_t end,
                           size_t depth)
{
	size_t hlen = length(host);
	int state = 0; // 0 -> < section, 1 -> = section, 2 -> > section
	for(size_t i = begin; i < end; i++) {
		switch(state) {
		case 0:
			if       (CHAR_AT_SUF(i, depth) == pivot) { state = 1; break; }
			else if  (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			assert_lt(CHAR_AT_SUF(i, depth), pivot);  break;
		case 1:
			if       (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			assert_eq(CHAR_AT_SUF(i, depth), pivot);  break;
		case 2:
			assert_gt(CHAR_AT_SUF(i, depth), pivot);  break;
		}
	}
	return true;
}
#endif

/**
 * Assert that the seqan::String s of suffix offsets into seqan::String
 * 'host' is a seemingly legitimate suffix-offset list (at this time,
 * we just check that it doesn't list any suffix twice).
 */
static void sanityCheckInputSufs(uint32_t *s, size_t slen) {
	assert_gt(slen, 0);
	for(size_t i = 0; i < slen; i++) {
		// Actually, it's convenient to allow the caller to provide
		// suffix offsets thare are off the end of the host string.
		// See, e.g., build() in diff_sample.cpp.
		//assert_lt(s[i], length(host));
		for(size_t j = i+1; j < slen; j++) {
			assert_neq(s[i], s[j]);
		}
	}
}

/**
 * Assert that the seqan::String s of suffix offsets into seqan::String
 * 'host' really are in lexicographical order up to depth 'upto'.
 */
template <typename T>
void sanityCheckOrderedSufs(const T& host,
                            size_t hlen,
                            uint32_t *s,
                            size_t slen,
                            size_t upto,
                            uint32_t lower = 0,
                            uint32_t upper = 0xffffffff)
{
	assert_lt(s[0], hlen);
	upper = min<uint32_t>(upper, slen-1);
	for(size_t i = lower; i < upper; i++) {
		// Allow s[i+t] to point off the end of the string; this is
		// convenient for some callers
		if(s[i+1] >= hlen) continue;
		if(upto == 0xffffffff) {
			if(!(dollarLt(suffix(host, s[i]), suffix(host, s[i+1])))) {
				assert(false);
			}
		} else {
			if(prefix(suffix(host, s[i]), upto) > prefix(suffix(host, s[i+1]), upto)) {
				// operator > treats shorter strings as
				// lexicographically smaller, but we want to opposite
				if(!isPrefix(suffix(host, s[i+1]), suffix(host, s[i]))) {
					//cout << "  " << suffix(host, s[i]) << endl;
					//cout << "  " << suffix(host, s[i+1]) << endl;
				}
				assert(isPrefix(suffix(host, s[i+1]), suffix(host, s[i])));
			}
		}
	}
}

/**
 * Toplevel function for multikey quicksort over suffixes.
 */
template<typename T>
void mkeyQSortSuf(const T& host,
                  uint32_t *s,
                  size_t slen,
                  int hi,
                  bool verbose = false,
                  bool sanityCheck = false,
                  size_t upto = 0xffffffff)
{
	size_t hlen = length(host);
	assert(!empty(s));
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	mkeyQSortSuf(host, hlen, s, slen, hi, 0, slen, 0, upto);
	if(sanityCheck) sanityCheckOrderedSufs(host, hlen, s, slen, upto);
}

/**
 * Toplevel function for multikey quicksort over suffixes with double
 * swapping.
 */
template<typename T>
void mkeyQSortSuf2(const T& host,
                   uint32_t *s,
                   size_t slen,
                   uint32_t *s2,
                   int hi,
                   bool verbose = false,
                   bool sanityCheck = false,
                   size_t upto = 0xffffffff)
{
	size_t hlen = length(host);
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	uint32_t *sOrig = NULL;
	if(sanityCheck) {
		sOrig = new uint32_t[slen];
		memcpy(sOrig, s, 4 * slen);
	}
	mkeyQSortSuf2(host, hlen, s, slen, s2, hi, 0, slen, 0, upto);
	if(sanityCheck) {
		sanityCheckOrderedSufs(host, hlen, s, slen, upto);
		for(size_t i = 0; i < slen; i++) {
			assert_eq(s[i], sOrig[s2[i]]);
		}
	}
}

/**
 * Main multikey quicksort function for suffixes.  Based on Bentley &
 * Sedgewick's algorithm on p.5 of their paper "Fast Algorithms for
 * Sorting and Searching Strings".  That algorithm has been extended in
 * three ways: 
 * 
 *  1. Deal with keys of different lengths by checking bounds and
 *     considering off-the-end values to be 'hi' (b/c our goal is the
 *     BWT transform, we're biased toward considring prefixes as
 *     lexicographically *greater* than their extensions).
 *  2. The multikey_qsort_suffixes version takes a single host string
 *     and a list of suffix offsets as input.  This reduces memory
 *     footprint compared to an approach that treats its input
 *     generically as a set of strings (not necessarily suffixes), thus
 *     requiring that we store at least two integers worth of
 *     information for each string.
 *  3. Sorting functions take an extra "upto" parameter that upper-
 *     bounds the depth to which the function sorts.
 * 
 * TODO: Consult a tie-breaker (like a difference cover sample) if two
 * keys share a long prefix.
 */
template<typename T>
void mkeyQSortSuf(const T& host,
                  size_t hlen,
                  uint32_t *s,
                  size_t slen,
                  int hi,
                  size_t begin,
                  size_t end,
                  size_t depth,
                  size_t upto = 0xffffffff)
{
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		if(ndepth < upto) { /* don't exceed depth of 'upto' */ \
			mkeyQSortSuf(host, hlen, s, slen, hi, nbegin, nend, ndepth, upto); \
		} \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t a, b, c, d, /*e,*/ r;
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF);
	int v = CHAR_AT_SUF(begin, depth); // v <- randomly-selected pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (hlen-s[i])) {
				stillInBounds = true;
				break;
			} else { /* already fell off this suffix */ }
		}
		assert(stillInBounds); // >=1 suffix must still be in bounds
	}
	#endif
	a = b = begin;
	c = d = end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF(b, depth))) {
			if(v == bc) {
				SWAP(s, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v <= (cc = CHAR_AT_SUF(c, depth))) {
			if(v == cc) {
				SWAP(s, c, d); d--;
			}
			c--;
		}
		if(b > c) break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(d-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth)); // check post-=-swap invariant
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted 
	if(v != hi) {
		MQS_RECURSE_SUF(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c; // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF(end-r, end, depth); // recurse on >'s
	}
}

/**
 * Just like mkeyQSortSuf but all swaps are applied to s2 as well as s.
 * This is a helpful variant if, for example, the caller would like to
 * see how their input was permuted by the sort routine (in that case,
 * the caller would let s2 be an array s2[] where s2 is the same length
 * as s and s2[i] = i). 
 */
template<typename T>
void mkeyQSortSuf2(const T& host,
                   size_t hlen,
                   uint32_t *s,
                   size_t slen,
                   uint32_t *s2,
                   int hi,
                   size_t begin,
                   size_t end,
                   size_t depth,
                   size_t upto = 0xffffffff)
{
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF_DS(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		if(ndepth < upto) { /* don't exceed depth of 'upto' */ \
			mkeyQSortSuf2(host, hlen, s, slen, s2, hi, nbegin, nend, ndepth, upto); \
		} \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t a, b, c, d, /*e,*/ r;
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	CHOOSE_AND_SWAP_PIVOT(SWAP2, CHAR_AT_SUF);
	int v = CHAR_AT_SUF(begin, depth); // v <- randomly-selected pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (hlen-s[i])) {
				stillInBounds = true;
				break;
			} else { /* already fell off this suffix */ }
		}
		assert(stillInBounds); // >=1 suffix must still be in bounds
	}
	#endif
	a = b = begin;
	c = d = /*e =*/ end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF(b, depth))) {
			if(v == bc) {
				SWAP2(s, s2, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v <= (cc = CHAR_AT_SUF(c, depth))) {
			if(v == cc) {
				SWAP2(s, s2, c, d); d--; /*e--;*/
			}
			//else if(c == e && v == hi) e--;
			c--;
		}
		if(b > c) break;
		SWAP2(s, s2, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(/*e*/d-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP2(s, s2, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP2(s, s2, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth)); // check post-=-swap invariant
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF_DS(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted 
	if(v != hi) {
		MQS_RECURSE_SUF_DS(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c;   // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF_DS(end-r, end, depth); // recurse on >'s
	}
}

// Ugly but necessary; otherwise the compiler chokes dramatically on
// the DifferenceCoverSample<> template args to the next few functions
template <typename T>
class DifferenceCoverSample;

/**
 * Toplevel function for multikey quicksort over suffixes.
 */
template<typename T>
void mkeyQSortSufDc(const T& host,
                    uint32_t* s,
                    size_t slen,
                    const DifferenceCoverSample<T>& dc,
                    int hi,
                    bool verbose = false,
                    bool sanityCheck = false)
{
	size_t hlen = length(host);
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	mkeyQSortSufDc(host, hlen, s, slen, dc, hi, 0, slen, 0, sanityCheck);
	if(sanityCheck) sanityCheckOrderedSufs(host, hlen, s, slen, 0xffffffff);
}

/**
 * Constant time
 */
template<typename T1, typename T2> inline
bool sufDcLt(const T1& host,
             const T2& s1,
             const T2& s2,
             const DifferenceCoverSample<T1>& dc,
             bool sanityCheck = false)
{
	uint32_t diff = dc.tieBreakOff(s1, s2);
	assert_lt(diff, dc.v());
	assert_lt(diff, length(host)-s1);
	assert_lt(diff, length(host)-s2);
	if(sanityCheck) {
		for(uint32_t i = 0; i < diff; i++) {
			assert_eq(host[s1+i], host[s2+i]);
		}
	}
	bool ret = dc.breakTie(s1+diff, s2+diff) < 0;
	if(sanityCheck && ret != dollarLt(suffix(host, s1), suffix(host, s2))) {
		assert(false);
	}
	return ret;
}

/**
 * k log(k)
 */
template<typename T> inline
void qsortSufDc(const T& host,
                size_t hlen,
                uint32_t* s,
                size_t slen,
                const DifferenceCoverSample<T>& dc,
                size_t begin,
                size_t end,
                bool sanityCheck = false)
{
	assert_leq(end, slen);
	assert_lt(begin, slen);
	assert_gt(end, begin);
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	size_t a = (random() % n) + begin; // choose pivot between begin and end
	assert_lt(a, end);
	assert_geq(a, begin);
	SWAP(s, end-1, a); // move pivot to end
	size_t cur = 0;
	for(size_t i = begin; i < end-1; i++) {
		if(sufDcLt(host, s[i], s[end-1], dc, sanityCheck)) {
			if(sanityCheck)
				assert(dollarLt(suffix(host, s[i]), suffix(host, s[end-1])));
			assert_lt(begin + cur, end-1);
			SWAP(s, i, begin + cur);
			cur++;
		}
	}
	// Put pivot into place
	assert_lt(cur, end-begin);
	SWAP(s, end-1, begin+cur);
	if(begin+cur > begin) qsortSufDc(host, hlen, s, slen, dc, begin, begin+cur);
	if(end > begin+cur+1) qsortSufDc(host, hlen, s, slen, dc, begin+cur+1, end);
}

/**
 * Main multikey quicksort function for suffixes.  Based on Bentley &
 * Sedgewick's algorithm on p.5 of their paper "Fast Algorithms for
 * Sorting and Searching Strings".  That algorithm has been extended in
 * three ways: 
 * 
 *  1. Deal with keys of different lengths by checking bounds and
 *     considering off-the-end values to be 'hi' (b/c our goal is the
 *     BWT transform, we're biased toward considring prefixes as
 *     lexicographically *greater* than their extensions).
 *  2. The multikey_qsort_suffixes version takes a single host string
 *     and a list of suffix offsets as input.  This reduces memory
 *     footprint compared to an approach that treats its input
 *     generically as a set of strings (not necessarily suffixes), thus
 *     requiring that we store at least two integers worth of
 *     information for each string.
 *  3. Sorting functions take an extra "upto" parameter that upper-
 *     bounds the depth to which the function sorts.
 * 
 * TODO: Consult a tie-breaker (like a difference cover sample) if two
 * keys share a long prefix.
 */
template<typename T>
void mkeyQSortSufDc(const T& host,
                    size_t hlen,
                    uint32_t* s,
                    size_t slen,
                    const DifferenceCoverSample<T>& dc,
                    int hi,
                    size_t begin,
                    size_t end,
                    size_t depth,
                    bool sanityCheck = false)
{
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF_DC(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		mkeyQSortSufDc(host, hlen, s, slen, dc, hi, nbegin, nend, ndepth, sanityCheck); \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t n = end - begin;
	if(n <= 1) return; // 1-element list already sorted
	if(depth > dc.v()) {
		// Quicksort the remaining suffixes using difference cover
		// for constant-time comparisons; this is O(k*log(k)) where
		// k=(end-begin)
		qsortSufDc<T>(host, hlen, s, slen, dc, begin, end, sanityCheck);
		return;
	}
	size_t a, b, c, d, r;
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF); // choose pivot, swap to begin
	int v = CHAR_AT_SUF(begin, depth); // v <- pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (hlen-s[i])) {
				stillInBounds = true;
				break;
			} else { /* already fell off this suffix */ }
		}
		assert(stillInBounds); // >=1 suffix must still be in bounds
	}
	#endif
	a = b = begin;
	c = d = end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF(b, depth))) {
			if(v == bc) {
				SWAP(s, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v <= (cc = CHAR_AT_SUF(c, depth))) {
			if(v == cc) {
				SWAP(s, c, d); d--;
			}
			c--;
		}
		if(b > c) break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(d-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth)); // check post-=-swap invariant
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF_DC(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted 
	if(v != hi) {
		MQS_RECURSE_SUF_DC(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c; // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF_DC(end-r, end, depth); // recurse on >'s
	}
}

/**
 * Toplevel function for multikey quicksort over suffixes.
 */
template<typename T1, typename T2>
void mkeyQSortSufDcU8(const T1& seqanHost,
                      const T2& host,
                      size_t hlen,
                      uint32_t* s,
                      size_t slen,
                      const DifferenceCoverSample<T1>& dc,
                      int hi,
                      bool verbose = false,
                      bool sanityCheck = false)
{
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	mkeyQSortSufDcU8(seqanHost, host, hlen, s, slen, dc, hi, 0, slen, 0, sanityCheck);
	if(sanityCheck) sanityCheckOrderedSufs(seqanHost, hlen, s, slen, 0xffffffff);
}

/**
 * Constant time
 */
template<typename T1, typename T2> inline
bool sufDcLtU8(const T1& seqanHost,
               const T2& host,
               size_t hlen,
               uint32_t s1,
               uint32_t s2,
               const DifferenceCoverSample<T1>& dc,
               bool sanityCheck = false)
{
	hlen += 0;
	uint32_t diff = dc.tieBreakOff(s1, s2);
	assert_lt(diff, dc.v());
	assert_lt(diff, hlen-s1);
	assert_lt(diff, hlen-s2);
	if(sanityCheck) {
		for(uint32_t i = 0; i < diff; i++) {
			assert_eq(host[s1+i], host[s2+i]);
		}
	}
	bool ret = dc.breakTie(s1+diff, s2+diff) < 0;
	if(sanityCheck && ret != dollarLt(suffix(seqanHost, s1), suffix(seqanHost, s2))) {
		assert(false);
	}
	return ret;
}

/**
 * k log(k)
 */
template<typename T1, typename T2> inline
void qsortSufDcU8(const T1& seqanHost,
                  const T2& host,
                  size_t hlen,
                  uint32_t* s,
                  size_t slen,
                  const DifferenceCoverSample<T1>& dc,
                  size_t begin,
                  size_t end,
                  bool sanityCheck = false)
{
	assert_leq(end, slen);
	assert_lt(begin, slen);
	assert_gt(end, begin);
	size_t n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	size_t a = (random() % n) + begin; // choose pivot between begin and end
	assert_lt(a, end);
	assert_geq(a, begin);
	SWAP(s, end-1, a); // move pivot to end
	size_t cur = 0;
	for(size_t i = begin; i < end-1; i++) {
		if(sufDcLtU8(seqanHost, host, hlen, s[i], s[end-1], dc, sanityCheck)) {
			if(sanityCheck)
				assert(dollarLt(suffix(seqanHost, s[i]), suffix(seqanHost, s[end-1])));
			assert_lt(begin + cur, end-1);
			SWAP(s, i, begin + cur);
			cur++;
		}
	}
	// Put pivot into place
	assert_lt(cur, end-begin);
	SWAP(s, end-1, begin+cur);
	if(begin+cur > begin) qsortSufDcU8(seqanHost, host, hlen, s, slen, dc, begin, begin+cur);
	if(end > begin+cur+1) qsortSufDcU8(seqanHost, host, hlen, s, slen, dc, begin+cur+1, end);
}

#define BUCKET_SORT_CUTOFF (4 * 1024 * 1024)
#define SELECTION_SORT_CUTOFF 4

// 5 64-element buckets for bucket-sorting A, C, G, T, $
static uint32_t bkts[4][4 * 1024 * 1024];

template<typename T1, typename T2>
static void selectionSortSufDcU8(
		const T1& seqanHost,
		const T2& host,
        size_t hlen,
        uint32_t* s,
        size_t slen,
        const DifferenceCoverSample<T1>& dc,
        uint8_t hi,
        size_t begin,
        size_t end,
        size_t depth,
        bool sanityCheck = false)
{
	assert_gt(end, begin+1);
	assert_leq(end-begin, SELECTION_SORT_CUTOFF);
	assert_eq(hi, 4);
	uint32_t v = dc.v();
	if(end == begin+2) {
		uint32_t off = dc.tieBreakOff(s[begin], s[begin+1]);
		if(off != 0xffffffff) {
			if(off < depth) {
				qsortSufDcU8<T1,T2>(seqanHost, host, hlen, s, slen, dc, begin, end, sanityCheck);
				// It's helpful for debugging if we call this here
				if(sanityCheck) {
					sanityCheckOrderedSufs(seqanHost, hlen, s, slen,
					                       0xffffffff, begin, end);
				}
				return;
			}
			v = off - depth;
		}
	}
	assert_leq(v, dc.v());
	uint32_t lim = v;
	assert_geq(lim, 0);
	bool atLeastOneTie = false;
	for(size_t i = begin; i < end-1; i++) {
		uint32_t targ = i;
		uint32_t targoff = depth + s[i];
		for(size_t j = i+1; j < end; j++) {
			uint32_t joff = depth + s[j];
			uint32_t k;
			for(k = 0; k <= lim; k++) {
				#ifndef PACKED_STRINGS
				uint8_t jc = (k + joff < hlen) ? host[k + joff] : hi;
				uint8_t tc = (k + targoff < hlen) ? host[k + targoff] : hi;
				#else
				uint8_t jc = (k + joff < hlen) ? (uint8_t)(Dna)host[k + joff] : hi;
				uint8_t tc = (k + targoff < hlen) ? (uint8_t)(Dna)host[k + targoff] : hi;
				#endif
				if(jc > tc) {
					// the jth suffix is greater than the current
					// smallest suffix
					break;
				} else if(jc < tc) {
					// the jth suffix is less than the current smallest
					// suffix, so update smallest to be j
					targ = j;
					targoff = joff;
					break;
				} else {
					// They're equal so far, keep going
				}
			}
			// The jth suffix was equal to the current smallest suffix
			// up to the difference-cover period, so disambiguation
			// with difference cover will be necessary.  Note that we
			// never change the relative order of two equal keys.
			if(k == lim+1) {
				atLeastOneTie = true;
			}
		}
		if(i != targ) {
			// swap i and targ
			uint32_t tmp = s[i];
			s[i] = s[targ];
			s[targ] = tmp;
		}
	}
	if(atLeastOneTie) {
		qsortSufDcU8<T1,T2>(seqanHost, host, hlen, s, slen, dc, begin, end, sanityCheck);
		// It's helpful for debugging if we call this here
		if(sanityCheck) {
			sanityCheckOrderedSufs(seqanHost, hlen, s, slen,
			                       0xffffffff, begin, end);
		}
	} else {
		// It's helpful for debugging if we call this here
		if(sanityCheck) {
			sanityCheckOrderedSufs(seqanHost, hlen, s, slen,
			                       0xffffffff, begin, end);
		}
	}
}

template<typename T1, typename T2>
static void bucketSortSufDcU8(
		const T1& seqanHost,
		const T2& host,
        size_t hlen,
        uint32_t* s,
        size_t slen,
        const DifferenceCoverSample<T1>& dc,
        uint8_t hi,
        size_t begin,
        size_t end,
        size_t depth,
        bool sanityCheck = false)
{
	uint32_t cnts[] = { 0, 0, 0, 0, 0 };
	#define BKT_RECURSE_SUF_DC_U8(nbegin, nend) { \
		bucketSortSufDcU8<T1,T2>(seqanHost, host, hlen, s, slen, dc, hi, \
		                         (nbegin), (nend), depth+1, sanityCheck); \
	}
	assert_gt(end, begin);
	assert_leq(end-begin, BUCKET_SORT_CUTOFF);
	assert_eq(hi, 4);
	if(end == begin+1) return; // 1-element list already sorted
	if(depth > dc.v()) {
		// Quicksort the remaining suffixes using difference cover
		// for constant-time comparisons; this is O(k*log(k)) where
		// k=(end-begin)
		qsortSufDcU8<T1,T2>(seqanHost, host, hlen, s, slen, dc, begin, end, sanityCheck);
		return;
	}
	if(end-begin <= SELECTION_SORT_CUTOFF) {
		// Bucket sort remaining items
		selectionSortSufDcU8(seqanHost, host, hlen, s, slen, dc, hi,
		                     begin, end, depth, sanityCheck);
		if(sanityCheck) {
			sanityCheckOrderedSufs(seqanHost, hlen, s, slen,
			                       0xffffffff, begin, end);
		}
		return;
	}
	for(size_t i = begin; i < end; i++) {
		uint32_t off = depth+s[i];
		// Following line is a huge chunk of the gprof profile
		#ifndef PACKED_STRINGS
		uint8_t c = (off < hlen) ? host[off] : hi;
		#else
		uint8_t c = (off < hlen) ? (uint8_t)(Dna)host[off] : hi;
		#endif
		assert_leq(c, 4);
		if(c == 0) {
			s[begin + cnts[0]++] = s[i];
		} else {
			bkts[c-1][cnts[c]++] = s[i];
		}
	}
	assert_eq(cnts[0] + cnts[1] + cnts[2] + cnts[3] + cnts[4], end - begin);
	uint32_t cur = begin + cnts[0];
	// TODO: are straight copies faster than memcpys?
#if 1
	if(cnts[1] > 0) { memcpy(&s[cur], bkts[0], cnts[1] << 2); cur += cnts[1]; }
	if(cnts[2] > 0) { memcpy(&s[cur], bkts[1], cnts[2] << 2); cur += cnts[2]; }
	if(cnts[3] > 0) { memcpy(&s[cur], bkts[2], cnts[3] << 2); cur += cnts[3]; }
	if(cnts[4] > 0) { memcpy(&s[cur], bkts[3], cnts[4] << 2); }
#else
	if(cnts[1] > 0) {
		size_t j = 0, i = cur;
		for(; j < cnts[1]; i++, j++) s[i] = bkts[0][j];
		cur += cnts[1];
	}
	if(cnts[2] > 0) {
		size_t j = 0, i = cur;
		for(; j < cnts[2]; i++, j++) s[i] = bkts[1][j];
		cur += cnts[2];
	}
	if(cnts[3] > 0) {
		size_t j = 0, i = cur;
		for(; j < cnts[3]; i++, j++) s[i] = bkts[2][j];
		cur += cnts[3];
	}
	if(cnts[4] > 0) {
		size_t j = 0, i = cur;
		for(; j < cnts[4]; i++, j++) s[i] = bkts[3][j];
	}
#endif
	// This frame is now totally finished with bkts[][], so recursive
	// callees can safely clobber it; we're not done with cnts[], but
	// that's local to the stack frame.
	cur = begin;
	if(cnts[0] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[0]); cur += cnts[0];
	}
	if(cnts[1] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[1]); cur += cnts[1];
	}
	if(cnts[2] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[2]); cur += cnts[2];
	}
	if(cnts[3] > 0) {
		BKT_RECURSE_SUF_DC_U8(cur, cur + cnts[3]);
	}
	// Done
}

/**
 * Main multikey quicksort function for suffixes.  Based on Bentley &
 * Sedgewick's algorithm on p.5 of their paper "Fast Algorithms for
 * Sorting and Searching Strings".  That algorithm has been extended in
 * three ways: 
 * 
 *  1. Deal with keys of different lengths by checking bounds and
 *     considering off-the-end values to be 'hi' (b/c our goal is the
 *     BWT transform, we're biased toward considring prefixes as
 *     lexicographically *greater* than their extensions).
 *  2. The multikey_qsort_suffixes version takes a single host string
 *     and a list of suffix offsets as input.  This reduces memory
 *     footprint compared to an approach that treats its input
 *     generically as a set of strings (not necessarily suffixes), thus
 *     requiring that we store at least two integers worth of
 *     information for each string.
 *  3. Sorting functions take an extra "upto" parameter that upper-
 *     bounds the depth to which the function sorts.
 * 
 * TODO: Consult a tie-breaker (like a difference cover sample) if two
 * keys share a long prefix.
 */
template<typename T1, typename T2>
void mkeyQSortSufDcU8(const T1& seqanHost,
                      const T2& host,
                      size_t hlen,
                      uint32_t* s,
                      size_t slen,
                      const DifferenceCoverSample<T1>& dc,
                      int hi,
                      size_t begin,
                      size_t end,
                      size_t depth,
                      bool sanityCheck = false)
{
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF_DC_U8(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		mkeyQSortSufDcU8(seqanHost, host, hlen, s, slen, dc, hi, nbegin, nend, ndepth, sanityCheck); \
	}
	assert_leq(begin, slen);
	assert_leq(end, slen);
	size_t n = end - begin;
	if(n <= 1) return; // 1-element list already sorted
	if(depth > dc.v()) {
		// Quicksort the remaining suffixes using difference cover
		// for constant-time comparisons; this is O(k*log(k)) where
		// k=(end-begin)
		qsortSufDcU8<T1,T2>(seqanHost, host, hlen, s, slen, dc, begin, end, sanityCheck);
		if(sanityCheck) {
			sanityCheckOrderedSufs(seqanHost, hlen, s, slen, 0xffffffff, begin, end);
		}
		return;
	}
	if(n <= BUCKET_SORT_CUTOFF) {
		// Bucket sort remaining items
		bucketSortSufDcU8(seqanHost, host, hlen, s, slen, dc,
		                  (uint8_t)hi, begin, end, depth, sanityCheck);
		if(sanityCheck) {
			sanityCheckOrderedSufs(seqanHost, hlen, s, slen, 0xffffffff, begin, end);
		}
		return;
	}
	size_t a, b, c, d, r;
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF_U8); // choose pivot, swap to begin
	int v = CHAR_AT_SUF_U8(begin, depth); // v <- pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (hlen-s[i])) {
				stillInBounds = true;
				break;
			} else { /* already fell off this suffix */ }
		}
		assert(stillInBounds); // >=1 suffix must still be in bounds
	}
	#endif
	a = b = begin;
	c = d = end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF_U8(b, depth))) {
			if(v == bc) {
				SWAP(s, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		bool hiLatch = true;
		while(b <= c && v <= (cc = CHAR_AT_SUF_U8(c, depth))) {
			if(v == cc) {
				SWAP(s, c, d); d--;
			}
			else if(hiLatch && cc == hi) {
				
			}
			c--;
		}
		if(b > c) break;
		SWAP(s, b, c);
		b++;
		c--;
	}
	assert(a > begin || c < end-1);                      // there was at least one =s
	assert_lt(d-c, n); // they can't all have been > pivot
	assert_lt(b-a, n); // they can't all have been < pivot
	//assert(assertPartitionedSuf(host, s, slen, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	//assert(assertPartitionedSuf2(host, s, slen, hi, v, begin, end, depth)); // check post-=-swap invariant
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF_DC_U8(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted 
	if(v != hi) {
		MQS_RECURSE_SUF_DC_U8(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c; // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF_DC_U8(end-r, end, depth); // recurse on >'s
	}
}


#endif /*MULTIKEY_QSORT_H_*/
