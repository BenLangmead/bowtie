#ifndef REF_READ_H_
#define REF_READ_H_

#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <ctype.h>
#include <fstream>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "filebuf.h"
#include "word_io.h"

using namespace std;
using namespace seqan;

/// Skip to the end of the current line; return the first character
/// of the next line
static inline int skipWhitespace(FileBuf& in) {
	int c;
	while(isspace(c = in.get())) {
		if(in.eof()) return -1;
	}
	return c;
}

/**
 * Encapsulates a stretch of the reference containing only unambiguous
 * characters.  From an ordered list of RefRecords, one can (almost)
 * deduce the "shape" of the reference sequences (almost because we
 * lose information about stretches of ambiguous characters at the end
 * of reference sequences).
 */
struct RefRecord {
	RefRecord(uint32_t _off, uint32_t _len, bool _first) :
		off(_off), len(_len), first(_first)
	{ }

	RefRecord(FILE *in, bool swap) {
		assert(in != NULL);
		if(!fread(&off, 4, 1, in)) {
			cerr << "Error reading RefRecord offset from FILE" << endl;
			exit(1);
		}
		if(swap) off = endianSwapU32(off);
		if(!fread(&len, 4, 1, in)) {
			cerr << "Error reading RefRecord offset from FILE" << endl;
			exit(1);
		}
		if(swap) len = endianSwapU32(len);
		first = fgetc(in) ? true : false;
	}

	RefRecord(int in, bool swap) {
		off = readU32(in, swap);
		len = readU32(in, swap);
		char c;
		if(!read(in, &c, 1)) {
			cerr << "Error reading RefRecord 'first' flag" << endl;
			exit(1);
		}
		first = (c ? true : false);
	}

	void write(std::ostream& out, bool be) {
		writeU32(out, off, be);
		writeU32(out, len, be);
		out.put(first ? 1 : 0);
	}

	uint32_t off; /// Offset of the first character in the record
	uint32_t len; /// Length of the record
	bool   first; /// Whether this record is the first for a reference sequence
};

/**
 * Parameters governing treatment of references as they're read in.
 */
struct RefReadInParams {
	RefReadInParams(int64_t bc, int64_t sc, bool r, bool nsToA) :
		baseCutoff(bc), numSeqCutoff(sc), reverse(r), nsToAs(nsToA) { }
	// stop reading references once we've finished reading this many
	// total reference bases
	int64_t baseCutoff;
	// stop reading references once we've finished reading this many
	// distinct sequences
	int64_t numSeqCutoff;
	// reverse each reference sequence before passing it along
	bool reverse;
	// convert ambiguous characters to As
	bool nsToAs;
};

extern RefRecord fastaRefReadSize(FileBuf& in,
                                  const RefReadInParams& refparams,
                                  bool first,
                                  BitpairOutFileBuf* bpout = NULL);
extern size_t fastaRefReadSizes(vector<FileBuf*>& in,
                                vector<RefRecord>& recs,
                                const RefReadInParams& refparams,
                                BitpairOutFileBuf* bpout = NULL);

/**
 * For given filehandle, read to the end of the current line and return
 * the first character on the next line, or -1 if eof is encountered at
 * any time.
 */
static inline int skipLine(FileBuf& in) {
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
static RefRecord fastaRefReadAppend(FileBuf& in,
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
		if(refparams.nsToAs && dna4Cat[c] == 2) {
			c = 'A';
		}
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
			if(refparams.baseCutoff > 0 &&
			   (int64_t)seqCharsRead >= refparams.baseCutoff)
			{
				lastc = -1;
				goto bail;
			}
		}
		c = in.get();
		if(refparams.nsToAs && dna4Cat[c] == 2) {
			c = 'A';
		}
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
