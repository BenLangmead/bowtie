#ifndef READ_BFQ_H_
#define READ_BFQ_H_

#include <stdint.h>
// The following is necessary for zlib.h on the Mac; unsure why
#define Byte uint8_t
#include <zlib.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <seqan/sequence.h>

/**
 * Load one named sequence from a gzipped file.  At this point,
 * the character and quality value for a single character are
 * packed into a single byte, with the character occupying the
 * top 2 bits and the quality value occupying the bottom 6.
 * 
 * Note that .bfq files are only meant to handle DNA data.
 * 
 * Adapted from ma_load_1read in maq's read.h
 */
template <typename TStr>
static bool ma_load_1read(gzFile fp, TStr& ret) {
	typedef typename seqan::Value<TStr>::Type TVal;
	static char name[2048]; // buffer for .bfq read names
	static unsigned char seq[2048];  // buffer 
	int len;
	assert(seqan::empty(ret));
	// Read name length
	if (gzread(fp, &len, sizeof(int)) != sizeof(int)) return false;
	if(len > 2047) {
		throw std::runtime_error(
			"One or more .bfq read names are longer than 2047 characters");
	}
	// Read name
	gzread(fp, name, sizeof(char) * len);
	name[len] = '\0'; // TODO: do something with name (we just ignore it)
	// Read sequence length
	if(gzread(fp, &len, sizeof(int)) != sizeof(int)) return false;
	if(len > 2048) {
		throw std::runtime_error(
			"One or more .bfq read sequences are longer than 2048 bases");
	}
	seqan::reserve(ret, len, seqan::Exact());
	gzread(fp, seq, sizeof(char) * len);
	for(int i = 0; i < len; i++) {
		appendValue(ret, (TVal)((int)seq[i] >> 6));
	}
	return true;
}

/**
 * Read a sequence file of the given format and alphabet type.  Store
 * all of the extracted sequences in vector ss.  Note that SeqAn's
 * policy for when it encounters characters not from the specified
 * alphabet is to convert them to the lexicographically smallest
 * character in the alphabet.
 */
template <typename TStr>
static void readBfq(const std::string& infile,
                    std::vector<TStr>& ss,
                    int upto = -1)
{
	gzFile bfq;
	bfq = gzopen(infile.c_str(), "r");
	while(true) {
		TStr s;
		if(ma_load_1read(bfq, s)) {
			ss.push_back(s);
			if(upto != -1 && ss.size() >= (size_t)upto) {
				gzclose(bfq);
				return;
			}
		} else {
			break;
		}
	}
	gzclose(bfq);
}

/**
 * Read a set of sequence files of the given format and alphabet type.
 * Store all of the extracted sequences in vector ss.
 */
template <typename TStr>
static void readBfqs(const std::vector<std::string>& infiles,
                     std::vector<TStr>& ss,
                     int upto = -1)
{
	for(size_t i = 0; i < infiles.size(); i++) {
		readBfq<TStr>(infiles[i], ss, upto);
	}
}

#endif /*READ_BFQ_H_*/
