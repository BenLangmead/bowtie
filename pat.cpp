#include <cmath>
#include <iostream>
#include <string>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace std;
using namespace seqan;

void wrongQualityFormat() {
	cerr << "Encounterd space-separated qualities"<<endl
	     <<  "This appears to be an FASTQ-int file" << endl
	     << "Please re-run Bowtie with the --integer-quals option.";
}

void tooFewQualities(const String<char>& read_name) {
	cerr << "Too few quality values for read: " << read_name << endl
		 << "\tare you sure this is a FASTQ-int file?" << endl;
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
	return result;
}
