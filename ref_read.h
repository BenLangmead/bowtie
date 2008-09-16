#ifndef REF_READ_H_
#define REF_READ_H_

#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <ctype.h>
#include <fstream>
#include <seqan/sequence.h>
#include "assert_helpers.h"

using namespace std;
using namespace seqan;

static uint8_t dna4Cat[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0,
	       /*                                        - */
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
	       /*    A  B  C  D        G  H        K     M  N */
	/*  80 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
	       /*       R  S  T     V  W  X  Y */
	/*  96 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
           /*    a  b  c  d        g  h        k     m  n */
	/* 112 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
           /*       r  s  t     v  w  x  y */
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

/// Skip to the end of the current line; return the first character
/// of the next line
static inline int skipWhitespace(istream& in) {
	int c;
	while(isspace(c = in.get()));
	return c;
}

struct RefRecord {
	RefRecord(uint32_t _off, uint32_t _len, bool _first) :
		off(_off), len(_len), first(_first)
	{
	}
	size_t off;
	size_t len;
	bool   first;
};

/**
 * Parameters governing treatment of references as they're read in.
 */
struct RefReadInParams {
	RefReadInParams(int64_t bc, int64_t sc, bool r) :
		baseCutoff(bc), numSeqCutoff(sc), reverse(r) { }
	// stop reading references once we've finished reading this many
	// total reference bases
	int64_t baseCutoff;
	// stop reading references once we've finished reading this many
	// distinct sequences
	int64_t numSeqCutoff;
	// reverse each reference sequence before passing it along
	bool reverse;
};

extern RefRecord fastaRefReadSize(istream& in,
                                  const RefReadInParams& refparams,
                                  bool first);
extern size_t fastaRefReadSizes(vector<istream*>& in,
                                vector<RefRecord>& recs,
                                const RefReadInParams& refparams);

/**
 * For given filehandle, read to the end of the current line and return
 * the first character on the next line, or -1 if eof is encountered at
 * any time.
 */
static inline int skipLine(istream& in) {
	while(true) {
		int c = in.get(); if(in.eof()) return -1;
		if(c == '\n' || c == '\r') {
			while(c == '\n' || c == '\r') {
				c = in.get(); if(in.eof()) return -1;
			}
			// c now holds first character of next line
			return c;
		}
	}
}

/**
 * Reads the next sequence from the given FASTA file and appends it to
 * the end of dst, optionally reversing it.
 */
template <typename TStr>
static RefRecord fastaRefReadAppend(istream& in,
                                    bool first,
                                    TStr& dst,
                                    RefReadInParams& refparams,
                                    string* name = NULL)
{
	typedef typename Value<TStr>::Type TVal;
	int c;
	static int lastc = '>';
	if(first) {
		c = skipWhitespace(in);
		if(c != '>') {
			cerr << "Reference file does not seem to be a FASTA file" << endl;
			exit(1);
		}
		lastc = c;
	}
	assert_neq(-1, lastc);

	assert_neq(refparams.baseCutoff, 0);
	assert_neq(refparams.numSeqCutoff, 0);

	// RefRecord params
	size_t seqCharsRead = 0;
	size_t seqOff = 0;
	bool seqFirst = true;

	size_t ilen = length(dst);
	bool found_space_in_name = false;

	// Chew up the id line; if the next line is either
	// another id line or a comment line, keep chewing
	c = lastc;
	if(c == '>' || c == '#') {
		do {
			while (c == '#') {
				if((c = skipLine(in)) == -1) {
					lastc = -1;
					goto bail;
				}
			}
			assert_eq('>', c);
			while(true) {
				c = in.get();
				if(c == -1) {
					lastc = -1;
					goto bail;
				}
				if (name)
				{
					if (isspace(c))
					{
						if (name->length())
							found_space_in_name = true;
					}
					else if (!found_space_in_name)
					{
						name->push_back(c);
					}
				}

				if(c == '\n' || c == '\r') {
					while(c == '\n' || c == '\r') {
						c = in.get();
						if(c == -1) {
							lastc = -1;
							goto bail;
						}
					}
					// c now holds first character of next line
					break;
				}
			}
		} while (c == '>' || c == '#');
	} else {
		ASSERT_ONLY(int cc = toupper(c));
		assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
		seqFirst = false;
	}

	// Skip over gaps
	while(true) {
		int cat = dna4Cat[c];
		if(cat == 1) {
			// This is a DNA character
			break; // to read-in loop
		} else if(cat == 2) {
			ASSERT_ONLY(int cc = toupper(c));
			assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
			seqOff++; // skip it
		} else if(c == '>') {
			lastc = '>';
			goto bail;
		}
		c = in.get();
		if(c == -1) {
			lastc = -1;
			goto bail;
		}
	}
	assert_eq(1, dna4Cat[c]);

	// in now points just past the first character of a sequence
	// line, and c holds the first character
	while(true) {
		// Note: can't have a comment in the middle of a sequence,
		// though a comment can end a sequence
		int cat = dna4Cat[c];
		assert_neq(2, cat);
		if(cat == 1) {
			appendValue(dst, (Dna)(char)c);
			assert_lt((uint8_t)(Dna)dst[length(dst)-1], 4);
			seqCharsRead++;
			if((int64_t)seqCharsRead >= refparams.baseCutoff) {
				lastc = -1;
				goto bail;
			}
		}
		c = in.get();
		if (c == -1 || c == '>' || c == '#' || dna4Cat[c] == 2) {
			lastc = c;
			break;
		}
	}

  bail:
	// Optionally reverse the portion that we just appended
	if(refparams.reverse) {
		// Find limits of the portion we just appended
		size_t nlen = length(dst);
		assert_eq(nlen - ilen, seqCharsRead);
		if(seqCharsRead > 0) {
			size_t halfway =  ilen + (seqCharsRead>>1);
			// Reverse it in-place
			for(size_t i = ilen; i < halfway; i++) {
				size_t diff = i-ilen;
				size_t j = nlen-diff-1;
				TVal tmp = dst[i];
				dst[i] = dst[j];
				dst[j] = tmp;
			}
		}
	}
	return RefRecord(seqOff, seqCharsRead, seqFirst);
}

#endif /*ndef REF_READ_H_*/
