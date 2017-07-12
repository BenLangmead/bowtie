//
//  read.h
//  bowtie
//
//  Created by Ben Langmead on 7/6/17.
//  Copyright © 2017 jhu. All rights reserved.
//

#ifndef READ_H_
#define READ_H_

#include <stdint.h>
#include <sys/time.h>
#include "sstring.h"
#include "filebuf.h"

typedef uint64_t TReadId;
typedef size_t TReadOff;
typedef int64_t TAlScore;

class HitSet;

/**
 * A buffer for keeping all relevant information about a single read.
 * Each search thread has one.
 */
struct Read {

	typedef SStringExpandable<char, 1024, 2, 1024> TBuf;

	Read() { reset(); }

	~Read() {
		clearAll(); reset();
		// Prevent seqan from trying to free buffers
		_setBegin(patFw, NULL);
		_setBegin(patRc, NULL);
		_setBegin(qual, NULL);
		_setBegin(patFwRev, NULL);
		_setBegin(patRcRev, NULL);
		_setBegin(qualRev, NULL);
		_setBegin(name, NULL);
	}

#define RESET_BUF(str, buf, typ) _setBegin(str, (typ*)buf); _setLength(str, 0); _setCapacity(str, BUF_SIZE);
#define RESET_BUF_LEN(str, buf, len, typ) _setBegin(str, (typ*)buf); _setLength(str, len); _setCapacity(str, BUF_SIZE);

	/// Point all Strings to the beginning of their respective buffers
	/// and set all lengths to 0
	void reset() {
		patid = 0;
		trimmed5 = trimmed3 = 0;
		color = false;
		primer = '?';
		trimc = '?';
		seed = 0;
		parsed = false;
		RESET_BUF(patFw, patBufFw, Dna5);
		RESET_BUF(patRc, patBufRc, Dna5);
		RESET_BUF(qual, qualBuf, char);
		RESET_BUF(patFwRev, patBufFwRev, Dna5);
		RESET_BUF(patRcRev, patBufRcRev, Dna5);
		RESET_BUF(qualRev, qualBufRev, char);
		RESET_BUF(name, nameBuf, char);
		readOrigBuf.clear();
		qualOrigBuf.clear();
	}

	void clearAll() {
		seqan::clear(patFw);
		seqan::clear(patRc);
		seqan::clear(qual);
		seqan::clear(patFwRev);
		seqan::clear(patRcRev);
		seqan::clear(qualRev);
		seqan::clear(name);
		parsed = false;
		trimmed5 = trimmed3 = 0;
		color = false;
		primer = '?';
		trimc = '?';
		seed = 0;
		readOrigBuf.clear();
		qualOrigBuf.clear();
	}

	/// Return true iff the read (pair) is empty
	bool empty() const {
		return seqan::empty(patFw);
	}

	/// Return length of the read in the buffer
	uint32_t length() const {
		return (uint32_t)seqan::length(patFw);
	}

	/**
	 * Construct reverse complement of the pattern.  If read is in
	 * colorspace, reverse color string.
	 */
	void constructRevComps() {
		uint32_t len = length();
		RESET_BUF_LEN(patRc, patBufRc, len, Dna5);
		if(color) {
			for(uint32_t i = 0; i < len; i++) {
				// Reverse the sequence
				patBufRc[i]  = patBufFw[len-i-1];
			}
		} else {
			for(uint32_t i = 0; i < len; i++) {
				// Reverse-complement the sequence
				patBufRc[i]  = (patBufFw[len-i-1] == 4) ? 4 : (patBufFw[len-i-1] ^ 3);
			}
		}
	}

	/**
	 * Given patFw, patRc, and qual, construct the *Rev versions in
	 * place.  Assumes constructRevComps() was called previously.
	 */
	void constructReverses() {
		uint32_t len = length();
		RESET_BUF_LEN(patFwRev, patBufFwRev, len, Dna5);
		RESET_BUF_LEN(patRcRev, patBufRcRev, len, Dna5);
		RESET_BUF_LEN(qualRev, qualBufRev, len, char);
		for(uint32_t i = 0; i < len; i++) {
			patFwRev[i]  = patFw[len-i-1];
			patRcRev[i]  = patRc[len-i-1];
			qualRev[i]   = qual[len-i-1];
		}
	}

	/**
	 * Append a "/1" or "/2" string onto the end of the name buf if
	 * it's not already there.
	 */
	void fixMateName(int i) {
		assert(i == 1 || i == 2);
		size_t namelen = seqan::length(name);
		bool append = false;
		if(namelen < 2) {
			// Name is too short to possibly have /1 or /2 on the end
			append = true;
		} else {
			if(i == 1) {
				// append = true iff mate name does not already end in /1
				append =
					nameBuf[namelen-2] != '/' ||
					nameBuf[namelen-1] != '1';
			} else {
				// append = true iff mate name does not already end in /2
				append =
					nameBuf[namelen-2] != '/' ||
					nameBuf[namelen-1] != '2';
			}
		}
		if(append) {
			assert_leq(namelen, BUF_SIZE-2);
			_setLength(name, namelen + 2);
			nameBuf[namelen] = '/';
			nameBuf[namelen+1] = "012"[i];
		}
	}

	/**
	 * Dump basic information about this read to the given ostream.
	 */
	void dump(std::ostream& os) const {
		os << name << ' ';
		if(color) {
			for(size_t i = 0; i < seqan::length(patFw); i++) {
				os << "0123."[(int)patFw[i]];
			}
		} else {
			os << patFw;
		}
		os << qual << " ";
	}

	static const int BUF_SIZE = 1024;

	String<Dna5>  patFw;               // forward-strand sequence
	uint8_t       patBufFw[BUF_SIZE];  // forward-strand sequence buffer
	String<Dna5>  patRc;               // reverse-complement sequence
	uint8_t       patBufRc[BUF_SIZE];  // reverse-complement sequence buffer
	String<char>  qual;                // quality values
	char          qualBuf[BUF_SIZE];   // quality value buffer

	String<Dna5>  patFwRev;               // forward-strand sequence reversed
	uint8_t       patBufFwRev[BUF_SIZE];  // forward-strand sequence buffer reversed
	String<Dna5>  patRcRev;               // reverse-complement sequence reversed
	uint8_t       patBufRcRev[BUF_SIZE];  // reverse-complement sequence buffer reversed
	String<char>  qualRev;                // quality values reversed
	char          qualBufRev[BUF_SIZE];   // quality value buffer reversed

	// For remembering the exact input text used to define a read
	TBuf readOrigBuf;
	TBuf qualOrigBuf;

	String<char>  name;                // read name
	char          nameBuf[BUF_SIZE];   // read name buffer
	bool          parsed;              // whether read has been fully parsed
	uint32_t      patid;               // unique 0-based id based on order in read file(s)
	int           mate;                // 0 = single-end, 1 = mate1, 2 = mate2
	uint32_t      seed;                // random seed
	bool          color;               // whether read is in color space
	char          primer;              // primer base, for csfasta files
	char          trimc;               // trimmed color, for csfasta files
	int           trimmed5;            // amount actually trimmed off 5' end
	int           trimmed3;            // amount actually trimmed off 3' end
	HitSet        hitset;              // holds previously-found hits; for chaining
};

/**
 * Pass this around when you want to be able to parse a buffer with data for
 * many reads in it.  Alows it to be consumed bit by bit with a cursor.
 */
struct ParsingCursor {
	
	ParsingCursor() : buf(NULL), off(0) { }
	
	ParsingCursor(const Read::TBuf* buf_, size_t off_) :
		buf(buf_),
		off(off_) { }
	
	const Read::TBuf* buf;
	size_t off;
};

#endif /*READ_H_*/
