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

extern size_t fastaRefReadSize(istream& in,
                               const RefReadInParams& refparams, 
                               bool first);
extern size_t fastaRefReadSizes(vector<istream*>& in,
                                vector<uint32_t>& szs,
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
static size_t fastaRefReadAppend(istream& in,
                                 TStr& dst,
                                 RefReadInParams& refparams, 
								 string* name = NULL)
{
	typedef typename Value<TStr>::Type TVal;
	int c;

	assert_neq(refparams.baseCutoff, 0);
	assert_neq(refparams.numSeqCutoff, 0);
	size_t seqCharsRead = 0;
	size_t ilen = length(dst);
	bool found_space_in_name = false;
	// Pick off the first carat
	
	c = in.get(); if(in.eof()) goto bail;
	assert(c == '>' || c == '#');
	
	// Chew up the id line; if the next line is either
	// another id line or a comment line, keep chewing
	do {
		if (c == '#')
			if((c = skipLine(in)) == -1) goto bail;
		
		while(true) {
			c = in.get(); if(in.eof()) goto bail;
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
					c = in.get(); if(in.eof()) goto bail;
				}
				// c now holds first character of next line
				break;
			}
		}
		//if((c = skipLine(in)) == -1) goto bail;
	} while (c == '>' || c == '#');
	
	// in now points just past the first character of a sequence
	// line, and c holds the first character
	while(c != '>' && c != '#') {
		// Note: can't have a comment in the middle of a sequence,
		// though a comment can end a sequence
		
		if(isalpha(c)) {
			appendValue(dst, (Dna)(char)c);
			assert_lt((uint8_t)(Dna)dst[length(dst)-1], 4);
			seqCharsRead++;
			if((int64_t)seqCharsRead >= refparams.baseCutoff) {
				return seqCharsRead;
			}
		}
		if (in.peek() == '>' || in.peek() == '#')
			break;
		c = in.get();
		if(in.eof()) break;
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
	return seqCharsRead;
}

#endif /*ndef REF_READ_H_*/
