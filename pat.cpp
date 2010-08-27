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

void wrongQualityFormat(const String<char>& read_name) {
	cerr << "Encountered a space parsing the quality string for read " << read_name << endl
	     << "If this is a FASTQ file with integer (non-ASCII-encoded) qualities, please" << endl
	     << "re-run Bowtie with the --integer-quals option.  If this is a FASTQ file with" << endl
	     << "alternate basecall information, please re-run Bowtie with the --fuzzy option." << endl;
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
