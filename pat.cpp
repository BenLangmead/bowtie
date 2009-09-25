#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace std;
using namespace seqan;

void wrongQualityFormat(const String<char>& read_name) {
	cerr << "Encounterd a space parsing the quality string for read " << read_name << endl
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
		 << "Please truncate reads and quality values and and re-run Bowtie" << endl;
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
