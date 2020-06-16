#ifndef REF_READ_H_
#define REF_READ_H_

#include <cassert>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "alphabet.h"
#include "assert_helpers.h"
#include "ds.h"
#include "endian_swap.h"
#include "filebuf.h"
#include "word_io.h"
#include "ds.h"

using namespace std;

class RefTooLongException : public exception {

public:
	RefTooLongException() {
#ifdef BOWTIE_64BIT_INDEX
		// This should never happen!
		msg = "Error: Reference sequence has more than 2^64-1 characters!  "
		      "Please divide the reference into smaller chunks and index each "
			  "independently.";
#else
		msg = "Error: Reference sequence has more than 2^32-1 characters!  "
		      "Please build a large index by passing the --large-index option "
			  "to bowtie2-build";
#endif
	}

	~RefTooLongException() throw() {}

	const char* what() const throw() {
		return msg.c_str();
	}

protected:

	string msg;

};


/**
 * Encapsulates a stretch of the reference containing only unambiguous
 * characters.  From an ordered list of RefRecords, one can (almost)
 * deduce the "shape" of the reference sequences (almost because we
 * lose information about stretches of ambiguous characters at the end
 * of reference sequences).
 */
struct RefRecord {
	RefRecord() : off(), len(), first() { }
	RefRecord(TIndexOffU _off, TIndexOffU _len, bool _first) :
		off(_off), len(_len), first(_first)
	{ }

	RefRecord(FILE *in, bool swap) {
		assert(in != NULL);
		if(!fread(&off, OFF_SIZE, 1, in)) {
			cerr << "Error reading RefRecord offset from FILE" << endl;
			throw 1;
		}
		if(swap) off = endianSwapU(off);
		if(!fread(&len, OFF_SIZE, 1, in)) {
			cerr << "Error reading RefRecord offset from FILE" << endl;
			throw 1;
		}
		if(swap) len = endianSwapU(len);
		first = fgetc(in) ? true : false;
	}

	void write(std::ostream& out, bool be) {
		writeU<TIndexOffU>(out, off, be);
		writeU<TIndexOffU>(out, len, be);
		out.put(first ? 1 : 0);
	}

	TIndexOffU off; /// Offset of the first character in the record
	TIndexOffU len; /// Length of the record
	bool   first; /// Whether this record is the first for a reference sequence
};

enum {
	REF_READ_FORWARD = 0, // don't reverse reference sequence
	REF_READ_REVERSE,     // reverse entire reference sequence
	REF_READ_REVERSE_EACH // reverse each unambiguous stretch of reference
};

/**
 * Parameters governing treatment of references as they're read in.
 */
struct RefReadInParams {
	RefReadInParams(int r, bool nsToA, bool bisulf) :
		reverse(r), nsToAs(nsToA), bisulfite(bisulf) { }
	// reverse each reference sequence before passing it along
	int reverse;
	// convert ambiguous characters to As
	bool nsToAs;
	// bisulfite-convert the reference
	bool bisulfite;
};

extern RefRecord
fastaRefReadSize(
	FileBuf& in,
	const RefReadInParams& rparms,
	bool first,
	BitpairOutFileBuf* bpout = NULL);

extern std::pair<size_t, size_t>
fastaRefReadSizes(
	EList<FileBuf*>& in,
	EList<RefRecord>& recs,
	EList<uint32_t>& plens,
	const RefReadInParams& rparms,
	BitpairOutFileBuf* bpout,
	TIndexOff& numSeqs);

extern void
reverseRefRecords(
	const EList<RefRecord>& src,
	EList<RefRecord>& dst,
	bool recursive = false,
	bool verbose = false);

/**
 * Reads the next sequence from the given FASTA file and appends it to
 * the end of dst, optionally reversing it.
 */
template <typename TStr>
static RefRecord fastaRefReadAppend(FileBuf& in,
                                    bool first,
                                    TStr& dst,
				    TIndexOffU& dstoff,
                                    RefReadInParams& rparms,
                                    string* name = NULL)
{
	int c;
	static int lastc = '>';
	if(first) {
		c = in.getPastWhitespace();
		if(c != '>') {
			cerr << "Reference file does not seem to be a FASTA file" << endl;
			throw 1;
		}
		lastc = c;
	}
	assert_neq(-1, lastc);

	// RefRecord params
	size_t len = 0;
	size_t off = 0;
	first = true;

	size_t ilen = dstoff;

	// Chew up the id line; if the next line is either
	// another id line or a comment line, keep chewing
	c = lastc;
	if(c == '>' || c == '#') {
		do {
			while (c == '#') {
				if((c = in.getPastNewline()) == -1) {
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
				if(c == '\n' || c == '\r') {
					while(c == '\r' || c == '\n') c = in.get();
					if(c == -1) {
						lastc = -1;
						goto bail;
					}
					break;
				}
				if (name) name->push_back(c);
			}
			// c holds the first character on the line after the name
			// line
		} while (c == '>' || c == '#');
	} else {
		ASSERT_ONLY(int cc = toupper(c));
		assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
		first = false;
	}

	// Skip over an initial stretch of gaps or ambiguous characters.
	while(true) {
		int cat = dna4Cat[c];
		if(rparms.nsToAs && cat == 2) {
			c = 'A';
		}
		int cc = toupper(c);
		if(rparms.bisulfite && cc == 'C') c = cc = 'T';
		if(cat == 1) {
			break; // to read-in loop
		} else if(cat == 2) {
			off++; // skip it
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
		int cat = asc2dnacat[c];
		assert_neq(2, cat);
		if(cat == 1) {
			// Consume it
			len++;
			// Add it to referenece buffer
			dst.set(asc2dna[c], dstoff++);
			assert_lt((int)dst[dstoff-1], 4);
		}
		c = in.get();
		if(rparms.nsToAs && dna4Cat[c] == 2) c = 'A';
		if (c == -1 || c == '>' || c == '#' || dna4Cat[c] == 2) {
			lastc = c;
			break;
		}
		if(rparms.bisulfite && toupper(c) == 'C') c = 'T';
	}

  bail:
	// Optionally reverse the portion that we just appended.
	// ilen = length of buffer before this last sequence was appended.
	if(rparms.reverse == REF_READ_REVERSE_EACH) {
		// Find limits of the portion we just appended
		dst.reverseWindow(ilen, len);
	}
	return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
}

#endif /*ndef REF_READ_H_*/
