#ifndef SEQUENCE_IO_H_
#define SEQUENCE_IO_H_

#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include "assert_helpers.h"
#include "packed_io.h"
#include "maq/read_bfq.h"
#include "solexa.h"
#include "pat.h"

/**
 * Encapsualtes a source of patterns where reading is handled by a
 * SeqAn read() function. 
 */
template <typename TStr, typename TFile>
class SeqanPatternSource : public PatternSource<TStr> {
public:
	SeqanPatternSource(const std::string& infile) :
		PatternSource<TStr>(), _sz(64), _cur(), _next(), _in(NULL)
	{
		// Open input file
		_in = fopen(infile.c_str(), "r");
		if(_in == NULL) {
			throw runtime_error("Could not open sequence file");
		}
		// Associate large input buffer with FILE *in
		if(setvbuf(_in, _buf, _IOFBF, 256 * 1024) != 0) {
			throw runtime_error("Could not create input buffer for sequence file");
		}
		read();
		assert(!seqan::empty(_next));
	}
	virtual ~SeqanPatternSource() {
		if(_in != NULL) fclose(_in);
	}
	virtual const TStr& nextPattern() {
		assert(!seqan::empty(_next));
		_cur = _next;
		read();
		return _cur;
	}
	virtual bool hasMorePatterns() {
		return !seqan::empty(_next);
	}
private:
	void read() {
		seqan::clear(_next);
		seqan::resize(_next, _sz, Exact());
		while(true) {
			seqan::read2(_in, _next, TFile());
			if(seqan::length(_next) > _sz) {
				_sz = seqan::length(_next) + 4;
			} else {
				break;
			}
		}
	}
	size_t _sz;  // default pattern size (expandable)
	TStr _cur;   // current pattern; outsider might hold a reference to it
	TStr _next;  // next pattern; might be empty
	FILE *_in;   // file to read patterns from
	char _buf[256 * 1024]; // (large) input buffer
};

/**
 * Read a sequence file of the given format and alphabet type.  Store
 * all of the extracted sequences in vector ss.  Note that SeqAn's
 * policy for when it encounters characters not from the specified
 * alphabet is to convert them to the lexicographically smallest
 * character in the alphabet.
 */
template <typename TStr, typename TFile>
static void readSequenceFile(const std::string& infile,
                             std::vector<TStr>& ss,
                             int64_t& baseCutoff, // limit for total bases
                             int seqCutoff = -1,  // limit for sequences
                             bool reverse = false)
{
	typedef typename Value<TStr>::Type TVal;
	static char buf[256 * 1024]; // fairly large input buffer
	if(baseCutoff <= 0) return;
	FILE *in = fopen(infile.c_str(), "r");
	if(in == NULL) {
		throw runtime_error("Could not open sequence file");
	}
	// Associate large input buffer with FILE *in
	if(setvbuf(in, buf, _IOFBF, 256 * 1024) != 0) {
		throw runtime_error("Could not create input buffer for sequence file");
	}
	// Read entries using SeqAn
	int cnt = 0;
	while(!feof(in)) {
		TStr tmp;
		while(true) {
			seqan::read(in, tmp, TFile());
			if(seqan::empty(tmp)) {
				break;
			}
			if((int64_t)length(tmp) > baseCutoff) {
				resize(tmp, baseCutoff);
				baseCutoff = 0;
			} else {
				baseCutoff -= length(tmp);
			}
			if(reverse) {
				size_t len = length(tmp);
				for(size_t i = 0; i < len/2; i++) {
					TVal t = tmp[i];
					tmp[i] = tmp[len-i-1];
					tmp[len-i-1] = t;
				}
			}
			#ifndef NDEBUG
			for(size_t i = 0; i < length(tmp); i++) {
				assert_lt(tmp[i], 4);
				assert_geq(tmp[i], 0);
			}
			#endif
			ss.push_back(tmp);
			cnt++;
			if(seqCutoff != -1 && cnt >= seqCutoff) {
				fclose(in);
				return;
			}
		}
	}
	fclose(in);
}

/**
 * Read a set of sequence files of the given format and alphabet type.
 * Store all of the extracted sequences in vector ss.
 */
template <typename TStr, typename TFile>
static void readSequenceFiles(const std::vector<std::string>& infiles,
                              std::vector<TStr>& ss,
                              int64_t& baseCutoff,
                              int seqCutoff = -1,
                              bool reverse = false)
{
	for(size_t i = 0; i < infiles.size() && baseCutoff > 0; i++) {
		readSequenceFile<TStr,TFile>(infiles[i], ss, baseCutoff, seqCutoff, reverse);
		if(baseCutoff <= 0) break;
	}
}

/**
 * Read a set of sequence files of the given format and alphabet type.
 * Store all of the extracted sequences in vector ss.
 */
template <typename TStr, typename TFile>
static void readSequenceFiles(const std::vector<std::string>& infiles,
                              std::vector<TStr>& ss,
                              int seqCutoff = -1,
                              bool reverse = false)
{
	int64_t i = 0xffffffffll;
	readSequenceFiles<TStr,TFile>(infiles, ss, i, seqCutoff, reverse);
}

/**
 * Parse a comma-delimited list of strings of type T into a vector. 
 */
template <typename T>
void readSequenceString(const std::string& s,
                        std::vector<T>& ss,
                        int64_t& baseCutoff,
                        int seqCutoff = -1,
                        bool reverse = false)
{
	// Split string s using comma as a delimiter.  Borrowed from C++
	// Programming HOWTO 7.3
	std::string::size_type lastPos = s.find_first_not_of(",", 0);
	std::string::size_type pos = s.find_first_of(",", lastPos);
	while (baseCutoff > 0 && (std::string::npos != pos || std::string::npos != lastPos)) {
		string stmp = s.substr(lastPos, pos - lastPos);
		if((int64_t)stmp.length() < baseCutoff) {
			baseCutoff -= stmp.length();
		} else {
			stmp = stmp.substr(0, baseCutoff);
			baseCutoff = 0;
		}
		if(reverse) {
			size_t len = stmp.length();
			for(size_t i = 0; i < len/2; i++) {
				char tmp = stmp[i];
				stmp[i] = stmp[len-i-1];
				stmp[len-i-1] = tmp;
			}
			ss.push_back(T(stmp.c_str()));
		} else {
			ss.push_back(T(stmp.c_str()));
		}
		if(seqCutoff != -1 && ss.size() >= (size_t)seqCutoff) {
			return;
		}
	    lastPos = s.find_first_not_of(",", pos);
	    pos = s.find_first_of(",", lastPos);
	}
}

/**
 * Parse a comma-delimited list of strings of type T into a vector. 
 * Doesn't require callee to supply a baseCutoff.
 */
template <typename T>
void readSequenceString(const std::string& s,
                        std::vector<T>& ss,
                        int seqCutoff = -1,
                        bool reverse = false)
{
	int64_t i = 0xffffffffll;
	readSequenceString(s, ss, i, seqCutoff, reverse);
}

#endif /*SEQUENCE_IO_H_*/
