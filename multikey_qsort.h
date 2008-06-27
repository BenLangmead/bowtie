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
static inline void swap(TStr& s, TPos a, TPos b) {
	typedef typename Value<TStr>::Type TAlphabet;
	assert_lt(a, length(s));
	assert_lt(b, length(s));
	TAlphabet tmp = s[a];
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
	swap(s, a, b); \
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
	swap(s2, a, b); \
}

#define SWAP1(s, s2, a, b) { \
	SWAP(s, a, b); \
}

/**
 * Helper macro that swaps a range of elements [i, i+n) with another
 * range [j, j+n) in seqan::String s.
 */
#define VECSWAP(s, i, j, n) { \
	if(n > 0) { vecswap(s, i, j, n, begin, end); } \
}

/**
 * Helper macro that swaps a range of elements [i, i+n) with another
 * range [j, j+n) both in seqan::String s and seqan::String s2.
 */
#define VECSWAP2(s, s2, i, j, n) { \
	if(n > 0) { vecswap2(s, s2, i, j, n, begin, end); } \
}

/**
 * Helper function that swaps a range of elements [i, i+n) with another
 * range [j, j+n) in seqan::String s.  begin and end represent the
 * current range under consideration by the caller (one of the
 * recursive multikey_quicksort routines below).
 */
template <typename TStr, typename TPos>
static inline void vecswap(TStr& s, TPos i, TPos j, TPos n, TPos begin, TPos end) {
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
		swap(s, a, b);
	}
}

/**
 * Helper function that swaps a range of elements [i, i+n) with another
 * range [j, j+n) both in seqan::String s and seqan::String s2.  begin
 * and end represent the current range under consideration by the
 * caller (one of the recursive multikey_quicksort routines below).
 */
template <typename TStr, typename TPos>
static inline void vecswap2(TStr& s, TStr& s2, TPos i, TPos j, TPos n, TPos begin, TPos end) {
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
		swap(s, a, b);
		swap(s2, a, b);
	}
}

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.
#define CHAR_AT(ss, aa) ((length(s[ss]) > aa) ? (int)(Dna)(s[ss][aa]) : hi)

/// Retrieve an int-ized version of the ath character of string s, or,
/// if a goes off the end of s, return a (user-specified) int greater
/// than any TAlphabet character - 'hi'.
#define CHAR_AT_SUF(si, off) (((off+s[si]) < length(host)) ? ((int)(Dna)(host[off+s[si]])) : (hi))

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
 * Assert that the range of chars at depth 'depth'�in strings 'begin'
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
	typedef typename Value<TStr>::Type TAlphabet;
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		mkeyQSort(s, hi, nbegin, nend, ndepth); \
	}
	assert_leq(begin, length(s));
	assert_leq(end, length(s));
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
 * Simple helper to print a list of suffixes.
 */
template<typename THost, typename TElt1, typename TElt2>
void printSuffixList(const THost& host,
                     const String<TElt1> s,
                     ostream& out)
{
	for(size_t i = 0; i < length(s); i++) {
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
template<typename THost, typename TStr>
bool assertPartitionedSuf(const THost& host,
                          TStr& s,
                          int hi,
                          int pivot,
                          size_t begin,
                          size_t end,
                          size_t depth)
{
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
template<typename THost, typename TStr>
bool assertPartitionedSuf2(const THost& host,
                           TStr& s,
                           int hi,
                           int pivot,
                           size_t begin,
                           size_t end,
                           size_t depth)
{
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
 * 'host' is a legitimate suffix-offset list (doesn't list any suffix
 * twice and doesn't fall off of host).
 */
template <typename THost, typename TStr>
void sanityCheckInputSufs(const THost& host, TStr& s) {
	assert_gt(length(host), 0);
	assert_gt(length(s), 0);
	for(size_t i = 0; i < length(s); i++) {
		// Actually, it's convenient to allow the caller to provide
		// suffix offsets thare are off the end of the host string.
		// See, e.g., build() in diff_sample.cpp.
		//assert_lt(s[i], length(host));
		for(size_t j = i+1; j < length(s); j++) {
			assert_neq(s[i], s[j]);
		}
	}
}

/**
 * Assert that the seqan::String s of suffix offsets into seqan::String
 * 'host' really are in lexicographical order up to depth 'upto'.
 */
template <typename THost, typename TStr>
void sanityCheckOrderedSufs(const THost& host, TStr& s, size_t upto) {
	assert_lt(s[0], length(host));
	for(size_t i = 0; i < length(s)-1; i++) {
		// Allow s[i+t] to point off the end of the string; this is
		// convenient for some callers
		if(s[i+1] >= length(host)) continue;
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
template<typename THost, typename TStr>
void mkeyQSortSuf(const THost& host,
                  TStr& s,
                  int hi,
                  bool verbose = false,
                  bool sanityCheck = false,
                  size_t upto = 0xffffffff)
{
	assert(!empty(s));
	if(sanityCheck) sanityCheckInputSufs(host, s);
	mkeyQSortSuf(host, s, hi, 0, length(s), 0, upto);
	if(sanityCheck) sanityCheckOrderedSufs(host, s, upto);
}

/**
 * Toplevel function for multikey quicksort over suffixes with double
 * swapping.
 */
template<typename THost, typename TStr>
void mkeyQSortSuf2(const THost& host,
                   TStr& s,
                   TStr& s2,
                   int hi,
                   bool verbose = false,
                   bool sanityCheck = false,
                   size_t upto = 0xffffffff)
{
	assert_eq(length(s), length(s2));
	if(sanityCheck) sanityCheckInputSufs(host, s);
	TStr sOrig = s;
	mkeyQSortSuf2(host, s, s2, hi, 0, length(s), 0, upto);
	if(sanityCheck) {
		sanityCheckOrderedSufs(host, s, upto);
		for(size_t i = 0; i < length(s); i++) {
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
template<typename THost, typename TStr>
void mkeyQSortSuf(const THost& host,
                  TStr& s,
                  int hi,
                  size_t begin,
                  size_t end,
                  size_t depth,
                  size_t upto = 0xffffffff)
{
	typedef typename Value<TStr>::Type TAlphabet;
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		if(ndepth < upto) { /* don't exceed depth of 'upto' */ \
			mkeyQSortSuf(host, s, hi, nbegin, nend, ndepth, upto); \
		} \
	}
	assert_leq(begin, length(s));
	assert_leq(end, length(s));
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
			if(depth < (length(host)-s[i])) {
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
	assert(assertPartitionedSuf(host, s, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, hi, v, begin, end, depth)); // check post-=-swap invariant
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
template<typename THost, typename TStr>
void mkeyQSortSuf2(const THost& host,
                   TStr& s,
                   TStr& s2,
                   int hi,
                   size_t begin,
                   size_t end,
                   size_t depth,
                   size_t upto = 0xffffffff)
{
	typedef typename Value<TStr>::Type TAlphabet;
	// Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF_DS(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		if(ndepth < upto) { /* don't exceed depth of 'upto' */ \
			mkeyQSortSuf2(host, s, s2, hi, nbegin, nend, ndepth, upto); \
		} \
	}
	assert_leq(begin, length(s));
	assert_leq(end, length(s));
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
			if(depth < (length(host)-s[i])) {
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
	assert(assertPartitionedSuf(host, s, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP2(s, s2, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP2(s, s2, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, hi, v, begin, end, depth)); // check post-=-swap invariant
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
template<typename T1, typename T2>
void mkeyQSortSufDc(const T1& host,
                    String<T2>& s,
                    const DifferenceCoverSample<T1>& dc,
                    int hi,
                    bool verbose = false,
                    bool sanityCheck = false)
{
	if(sanityCheck) sanityCheckInputSufs(host, s);
	mkeyQSortSufDc(host, s, dc, hi, 0, length(s), 0, sanityCheck);
	if(sanityCheck) sanityCheckOrderedSufs(host, s, 0xffffffff);
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
template<typename T1, typename TAlphabet> inline
void qsortSufDc(const T1& host,
                String<TAlphabet>& s,
                const DifferenceCoverSample<T1>& dc,
                size_t begin,
                size_t end,
                bool sanityCheck = false)
{
	assert_leq(end, length(s));
	assert_lt(begin, length(s));
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
	if(begin+cur > begin) qsortSufDc(host, s, dc, begin, begin+cur);
	if(end > begin+cur+1) qsortSufDc(host, s, dc, begin+cur+1, end);
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
template<typename T1, typename TAlphabet>
void mkeyQSortSufDc(const T1& host,
                    String<TAlphabet>& s,
                    const DifferenceCoverSample<T1>& dc,
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
		mkeyQSortSufDc(host, s, dc, hi, nbegin, nend, ndepth, sanityCheck); \
	}
	assert_leq(begin, length(s));
	assert_leq(end, length(s));
	size_t n = end - begin;
	if(n <= 1) return; // 1-element list already sorted
	if(depth > dc.v()) {
		// Quicksort the remaining suffixes using difference cover
		// for constant-time comparisons; this is O(k*log(k)) where
		// k=(end-begin)
		qsortSufDc<T1, TAlphabet>(host, s, dc, begin, end, sanityCheck);
		return;
	}
	size_t a, b, c, d, r;
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF); // choose pivot, swap to begin
	int v = CHAR_AT_SUF(begin, depth); // v <- pivot value
	#ifndef NDEBUG
	{
		bool stillInBounds = false;
		for(size_t i = begin; i < end; i++) {
			if(depth < (length(host)-s[i])) {
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
	assert(assertPartitionedSuf(host, s, hi, v, begin, end, depth));  // check pre-=-swap invariant
	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
	assert(assertPartitionedSuf2(host, s, hi, v, begin, end, depth)); // check post-=-swap invariant
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

#endif /*MULTIKEY_QSORT_H_*/
