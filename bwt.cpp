#include "blockwise_sa.h"
#include <seqan/file.h>

/**
 * Example application showing how to use the KarkkainenBlockwiseSA
 * class to obtain a stream ordered suffixes of an input text.  The
 * KarkkainenBlockwiseSA class is working behind-the-scenes to divide
 * the suffix array into blocks and calculate each block as needed.
 */
int main(int argc, char** argv) {
	typedef String<Dna> TStr;
	TStr input;
	if(argc < 2) {
		cerr << "Must specify Fasta input file as first argument" << endl;
		return 1;
	}
	
	// Read a DNA sequence in from a Fasta file
	std::ifstream in(argv[1], std::ios_base::in);
	read(in, input, Fasta());
	if(empty(input)) {
		cerr << "Bad or empty input file" << endl;
		return 1;
	}

	// Construct a blockwise SA builder.  First argument is the input
	// text, second argument is the difference-cover periodicity, and
	// third argument is the maximum block size.
	KarkkainenBlockwiseSA<TStr> sa(input, length(input)/5, 1024);
	while(sa.hasMoreSuffixes()) {
		// Get the offset of the next greatest suffix
		uint32_t suf = sa.nextSuffix();
		// Output BWT char
		cout << ((suf == 0)? input[length(input)-1] : input[suf-1]);
	}
	cout << endl;
	return 0;
}
