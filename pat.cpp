#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "pat.h"
#include "filebuf.h"

using namespace std;
using namespace seqan;

/**
 * Parse a single quality string from fb and store qualities in r.
 * Assume the next character obtained via fb.get() is the first
 * character of the quality string.  When returning, the next
 * character returned by fb.peek() or fb.get() should be the first
 * character of the following line.
 */
int parseQuals(ReadBuf& r,
               FileBuf& fb,
               int readLen,
               int trim3,
               int trim5,
               bool intQuals,
               bool phred64,
               bool solexa64)
{
	int qualsRead = 0;
	int c = 0;
	assert(fb.peek() != '\n' && fb.peek() != '\r');
	_setBegin (r.qual, (char*)r.qualBuf);
	_setLength(r.qual, 0);
	if (intQuals) {
		while (c != '\r' && c != '\n' && c != -1) {
			bool neg = false;
			int num = 0;
			while(!isspace(c = fb.peek()) && !fb.eof()) {
				if(c == '-') {
					neg = true;
					assert_eq(num, 0);
				} else {
					if(!isdigit(c)) {
						char buf[2048];
						cerr << "Warning: could not parse quality line:" << endl;
						fb.getPastNewline();
						cerr << fb.copyLastN(buf);
						buf[2047] = '\0';
						cerr << buf;
						throw 1;
					}
					assert(isdigit(c));
					num *= 10;
					num += (c - '0');
				}
				fb.get();
			}
			if(neg) num = 0;
			// Phred-33 ASCII encode it and add it to the back of the
			// quality string
			r.qualBuf[qualsRead++] = ('!' + num);
			// Skip over next stretch of whitespace
			c = fb.peek();
			while(c != '\r' && c != '\n' && isspace(c) && !fb.eof()) {
				fb.get();
				c = fb.peek();
			}
		}
	} else {
		while (c != '\r' && c != '\n' && c != -1) {
			c = fb.get();
			r.qualBuf[qualsRead++] = charToPhred33(c, solexa64, phred64);
			c = fb.peek();
			while(c != '\r' && c != '\n' && isspace(c) && !fb.eof()) {
				fb.get();
				c = fb.peek();
			}
		}
	}
	if (qualsRead < readLen-1 ||
	    (qualsRead < readLen && !r.color))
	{
		tooFewQualities(r.name);
	}
	qualsRead -= trim3;
	if(qualsRead <= 0) return 0;
	int trimmedReadLen = readLen-trim3-trim5;
	if(trimmedReadLen < 0) trimmedReadLen = 0;
	if(qualsRead > trimmedReadLen) {
		// Shift everybody left
		for(int i = 0; i < readLen; i++) {
			r.qualBuf[i] = r.qualBuf[i+qualsRead-trimmedReadLen];
		}
	}
	_setLength(r.qual, trimmedReadLen);
	while(fb.peek() == '\n' || fb.peek() == '\r') fb.get();
	return qualsRead;
}

/**
 * "Light" parser.  This is inside the critical section, so the key is to do
 * just enough parsing so that another function downstream (finalize()) can do
 * the rest of the parsing.  Really this function's only job is to stick every
 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
 * then parses the contents of r.readOrigBuf later.
 */
pair<bool, bool> FastqPatternSource::readLight(ReadBuf& r) {
	int c = 0;
	size_t& len = r.readOrigBufLen;
	len = 0;
	if(first_) {
		c = getc_unlocked(fp_);
		while(c == '\r' || c == '\n') {
			c = getc_unlocked(fp_);
		}
		if(c != '@') {
			cerr << "Error: reads file does not look like a FASTQ file" << endl;
			throw 1;
		}
		r.readOrigBuf[len++] = c;
		assert_eq('@', c);
		first_ = false;
	}
	// Note: to reduce the number of times we have to enter the critical
	// section (each entrance has some assocaited overhead), we could populate
	// the buffer with several reads worth of data here, instead of just one.
	int newlines = 4;
	while(newlines) {
		c = getc_unlocked(fp_);
		if(c == '\n' || (c < 0 && newlines == 1)) {
			newlines--;
			c = '\n';
		} else if(c < 0) {
			return make_pair(false, true);
		}
		r.readOrigBuf[len++] = c;
	}
	readCnt_++;
	r.patid = (uint32_t)(readCnt_-1);
	return make_pair(true, c < 0);
}

/// Read another pattern from a FASTQ input file
void FastqPatternSource::finalize(ReadBuf &r) const {
	int c;
	size_t name_len = 0;
	size_t cur = 1;
	r.color = color_;
	r.primer = -1;

	// Parse read name
	while(true) {
		assert(cur < r.readOrigBufLen);
		c = r.readOrigBuf[cur++];
		if(c == '\n' || c == '\r') {
			do {
				c = r.readOrigBuf[cur++];
			} while(c == '\n' || c == '\r');
			break;
		}
		r.nameBuf[name_len++] = c;
	}
	_setBegin(r.name, r.nameBuf);
	_setLength(r.name, name_len);

	if(color_) {
		// May be a primer character.  If so, keep it 'primer' field
		// of read buf and parse the rest of the read without it.
//		c = toupper(c);
//		if(asc2dnacat[c] > 0) {
//			// First char is a DNA char
//			int c2 = toupper(fb_.peek());
//			// Second char is a color char
//			if(asc2colcat[c2] > 0) {
//				r.primer = c;
//				r.trimc = c2;
//				mytrim5 += 2; // trim primer and first color
//			}
//		}
//		if(c < 0) { bail(r); return; }
		throw 1;
	}
	
	// Parse sequence
	int nchar = 0;
	uint8_t *seqbuf = r.patBufFw;
	while(c != '+') {
		if(c == '.') {
			c = 'N';
		}
		if(color_ && c >= '0' && c <= '4') {
			c = "ACGTN"[(int)c - '0'];
		}
		if(isalpha(c)) {
			// If it's past the 5'-end trim point
			if(nchar++ >= trim5_) {
				*seqbuf++ = charToDna5[c];
			}
		}
		assert(cur < r.readOrigBufLen);
		c = r.readOrigBuf[cur++];
	}
	int seq_len = (int)(seqbuf - r.patBufFw);
	r.trimmed5 = (int)(nchar - seq_len);
	r.trimmed3 = min(seq_len, trim3_);
	seq_len = max(seq_len - trim3_, 0);
	_setBegin(r.patFw, (Dna5*)r.patBufFw);
	_setLength(r.patFw, seq_len);
	
	assert_eq('+', c);
	do {
		assert(cur < r.readOrigBufLen);
		c = r.readOrigBuf[cur++];
	} while(c != '\n' && c != '\r');
	do {
		assert(cur < r.readOrigBufLen);
		c = r.readOrigBuf[cur++];
	} while(c == '\n' || c == '\r');
	
	// Now we're on the next non-blank line after the + line
	if(seq_len == 0) {
		return; // done parsing empty read
	}

	int nqual = 0;
	char *qualbuf = r.qualBuf;
	if (intQuals_) {
		// TODO: must implement this for compatibility with other Bowtie
		throw 1; // not yet implemented
	} else {
		while(c != '\r' && c != '\n') {
			c = charToPhred33(c, solQuals_, phred64Quals_);
			if(c == ' ') {
				wrongQualityFormat(r.name);
			}
			if(nqual++ >= r.trimmed5) {
				*qualbuf++ = c;
			}
			c = r.readOrigBuf[cur++];
		}
		int qual_len = (int)(qualbuf - r.qualBuf);
		qual_len = max(0, qual_len - r.trimmed3);
		if(qual_len < seq_len) {
			tooFewQualities(r.name);
		} else if(qual_len > seq_len+1) {
			tooManyQualities(r.name);
		}
		_setBegin(r.qual, (char*)r.qualBuf);
		_setLength(r.qual, seq_len);
	}
	// Set up a default name if one hasn't been set
	if(name_len == 0) {
		itoa10((int)readCnt_, r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		name_len = (int)strlen(r.nameBuf);
		_setLength(r.name, name_len);
	}
}

void wrongQualityFormat(const String<char>& read_name) {
	cerr << "Encountered a space parsing the quality string for read " << read_name << endl
	     << "If this is a FASTQ file with integer (non-ASCII-encoded) qualities, please" << endl
	     << "re-run Bowtie with the --integer-quals option." << endl;
	throw 1;
}

void tooFewQualities(const String<char>& read_name) {
	cerr << "Too few quality values for read: " << read_name << endl
		 << "\tare you sure this is a FASTQ-int file?" << endl;
	throw 1;
}

void tooManyQualities(const String<char>& read_name) {
	cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
		 << "Please truncate reads and quality values and and re-run Bowtie" << endl;
	throw 1;
}

void tooManySeqChars(const String<char>& read_name) {
	cerr << "Reads file contained a pattern with more than 1024 sequence characters." << endl
		 << "Please truncate reads and quality values and and re-run Bowtie." << endl
		 << "Offending read: " << read_name << endl;
	throw 1;
}

/**
 * C++ version char* style "itoa":
 */
char* itoa10(int value, char* result) {
	// Check that base is valid
	char* out = result;
	int quotient = value;
	do {
		*out = "0123456789"[ std::abs( quotient % 10 ) ];
		++out;
		quotient /= 10;
	} while ( quotient );

	// Only apply negative sign for base 10
	if (value < 0) *out++ = '-';
	std::reverse( result, out );

	*out = 0; // terminator
	return out;
}
