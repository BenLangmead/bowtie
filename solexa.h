#ifndef SOLEXA_H_
#define SOLEXA_H_

#include <string>
#include <vector>
#include <seqan/sequence.h>

/**
 * Read all sequences from a Solexa _seq.txt-style sequence file and
 * add the corresponding TStrs to vector ss.  Skips sequences with .'s.
 * TStr may be any Seqan String type, but in practice we expect it to
 * be a Dna string.
 */
template <typename TStr>
static void readSolexaSeq(const std::string& infile,
                          std::vector<TStr>& ss,
                          int upto = -1)
{
	typedef typename seqan::Value<TStr>::Type TVal;
	static char buf[256 * 1024]; // fairly large input buffer
	FILE *in = fopen(infile.c_str(), "r");
	if(in == NULL) {
		throw runtime_error("Could not open sequence file");
	}
	// Associate large input buffer with FILE *in
	if(setvbuf(in, buf, _IOFBF, 256 * 1024) != 0) {
		throw runtime_error("Could not create input buffer for sequence file");
	}
	TStr s; reserve(s, 64);
	bool skip = false;
	while(!feof(in)) {
		// We don't try to parse the line; we simply append any and all
		// DNA characters encountered.   Everything else on a line is 
		// either whitespace or a number.
		int c = fgetc(in);
		if(c == '.') {
			// Skip sequences containing one or more .'s
			skip = true;
		}
		else if(c == 'A' || c == 'a' ||
				c == 'T' || c == 't' ||
				c == 'C' || c == 'c' ||
				c == 'G' || c == 'g')
		{
			// Append DNA character to s
			seqan::append(s, (TVal)c);
		}
		else if(c == '\n') {
			if(!seqan::empty(s) && !skip) {
				ss.push_back(s);
				if(upto != -1 && ss.size() >= (size_t)upto) {
					fclose(in);
					return;
				}
			}
			seqan::clear(s);
			skip = false;
		}
	}
	fclose(in);
}

/**
 * Read all sequences from a vector of Solexa _seq.txt-style sequence
 * files and add the corresponding TStrs to vector ss.  Skips sequences
 * with .'s.
 */
template <typename TStr>
static void readSolexaSeqs(const std::vector<std::string>& infiles,
                           std::vector<TStr>& ss,
                           int upto = -1)
{
	for(size_t i = 0; i < infiles.size(); i++) {
		readSolexaSeq<TStr>(infiles[i], ss, upto);
	}
}

#endif /*SOLEXA_H_*/
