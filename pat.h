#ifndef PAT_H_
#define PAT_H_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <ctype.h>
#include <zlib.h>
#include <fstream>
#include <seqan/sequence.h>
#include "assert_helpers.h"
#include "tokenize.h"

using namespace std;
using namespace seqan;

/// Reverse a string
template <typename TStr>
static inline void reverse(TStr& s) {
	typedef typename Value<TStr>::Type TVal;
	size_t len = length(s);
	for(size_t i = 0; i < (len>>1); i++) {
		TVal tmp = s[i];
		s[i] = s[len-i-1];
		s[len-i-1] = tmp;
	}
}

static bool nRcA = false;

/**
 * Install a character in forward-oriented and (optionally) a reverse-
 * complemented buffers where the size of the buffer is known
 * beforehand.
 */
template <typename TStr>
static inline void appendCharResized(TStr *dst,
                                     TStr *rcDst,
                                     uint32_t i,
                                     uint32_t len,
                                     char c)
{
	assert_geq(length(*dst), len);
	if(rcDst != NULL) {
		assert_geq(length(*rcDst), len);
	}
	if(c == 'A' || c == 'a') {
		(*dst)[i] = 'A';
		if(rcDst != NULL) (*rcDst)[len-i-1] = 'T';
	} else if(c == 'C' || c == 'c') {
		(*dst)[i] = 'C';
		if(rcDst != NULL) (*rcDst)[len-i-1] = 'G';
	} else if(c == 'G' || c == 'g') {
		(*dst)[i] = 'G';
		if(rcDst != NULL) (*rcDst)[len-i-1] = 'C';
	} else if(c == 'T' || c == 't') {
		(*dst)[i] = 'T';
		if(rcDst != NULL) (*rcDst)[len-i-1] = 'A';
	} else {
		// Wildcard; give forward and reverse complements 'A's
		(*dst)[i] = 'A';
		if(rcDst != NULL) (*rcDst)[len-i-1] = (nRcA ? 'A' : 'T');
	}
}

/**
 * Install a character in forward-oriented and (optionally) a reverse-
 * complemented buffers where the size of the buffer is not known
 * beforehand.  If rcDst is non-NULL, the caller must reverse it before
 * using it. 
 */
template <typename TStr>
static inline void appendChar(TStr *dst, TStr *rcDst, char c) {
	if(c == 'A' || c == 'a') {
		append(*dst, 'A');
		if(rcDst != NULL) append(*rcDst, 'T');
	} else if(c == 'C' || c == 'c') {
		append(*dst, 'C');
		if(rcDst != NULL) append(*rcDst, 'G');
	} else if(c == 'G' || c == 'g') {
		append(*dst, 'G');
		if(rcDst != NULL) append(*rcDst, 'C');
	} else if(c == 'T' || c == 't') {
		append(*dst, 'T');
		if(rcDst != NULL) append(*rcDst, 'A');
	} else {
		// Wildcard; give forward and reverse complements 'A's
		append(*dst, 'A');
		if(rcDst != NULL) append(*rcDst, (nRcA ? 'A' : 'T'));
	}
}

/**
 * Encapsualtes a source of patterns; usually a file.  Handles dumping
 * patterns to a logfile (useful for debugging) and optionally
 * reversing them before returning them.
 */
template <typename TStr>
class PatternSource {
public:
	PatternSource(bool __reverse = false, const char *__dumpfile = NULL) :
		_reverse(__reverse), 
		_def_qual("EDCCCBAAAA@@@@?>===<;;9:9998777666655444433333333333333333333"),
		_def_rqual("333333333333333333334444336666778999:9;;<===>?@@@@AAAABCCCDE"),
		_dumpfile(__dumpfile), _out()
	{
		if(_dumpfile != NULL) {
			_out.open(_dumpfile, ios_base::out);
			if(!_out.good()) {
				throw runtime_error("Could not open pattern dump file for writing");
			}
		}
	}
	virtual ~PatternSource() { }
	virtual void nextPattern(TStr& s, string& qual, string& name) {
		nextPatternImpl(s, qual, name);
		// Reverse it, if desired
		if(_reverse) {
			size_t len = length(s);
			resize(_revtmp, len, Exact());
			_revqual.resize(len);
			// tmp = reverse of s
			for(size_t i = 0; i < len; i++)
			{
				_revtmp[i] = s[len-i-1];
				_revqual[i] = qual[len-i-1];
			}
			// Output it, if desired
			if(_dumpfile != NULL) {
				_out << _revtmp << "\t" << _revqual << endl;
			}
			s = _revtmp;
			qual = _revqual;
			return;
		}
		// Output it, if desired
		if(_dumpfile != NULL) {
			_out << s << endl;
		}
	}
	virtual void nextPatternImpl(TStr& , string&, string&) = 0;
	virtual bool hasMorePatterns() = 0;
	virtual bool nextIsReverseComplement() = 0;
	virtual void reset() { }
	const char *dumpfile() const { return _dumpfile; }
	bool reverse() const { return _reverse; }
	void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
private:
	bool _reverse;         // reverse patterns before returning them
	TStr _revtmp;          // temporary buffer for reversed patterns
	string _revqual;
	
	const string _def_qual;
	const string _def_rqual;
	const char *_dumpfile; // dump patterns to this file before returning them
	ofstream _out;         // output stream for dumpfile
};

/**
 * Encapsualtes a source of patterns where each raw pattern is trimmed
 * by some user-defined amount on the 3' and 5' ends.  Doesn't
 * implement the actual trimming - that's up to the concrete
 * descendants.
 */
template <typename TStr>
class TrimmingPatternSource : public PatternSource<TStr> {
public:
	TrimmingPatternSource(bool __reverse = false,
	                      const char *__dumpfile = NULL,
	                      int __trim3 = 0,
	                      int __trim5 = 0) :
		PatternSource<TStr>(__reverse, __dumpfile),
		_trim3(__trim3),
		_trim5(__trim5) { }
protected:
	int _trim3;
	int _trim5;
};

/**
 * Encapsualtes a source of patterns which is an in-memory vector. 
 */
template <typename TStr>
class VectorPatternSource : public TrimmingPatternSource<TStr> {
public:
	VectorPatternSource(vector<TStr> v,
	                    bool __revcomp = true,
	                    bool __reverse = false,
	                    const char *__dumpfile = NULL,
	                    size_t cur = 0,
	                    int __trim3 = 0,
	                    int __trim5 = 0) :
		TrimmingPatternSource<TStr>(__reverse, __dumpfile, __trim3, __trim5),
		_revcomp(__revcomp), _cur(cur), _v()
	{
		for(size_t i = 0; i < v.size(); i++) {
			TStr s(v[i]);
			_v.push_back(s);
			if(_revcomp) {
				_v.push_back(reverseComplement(s));
			}
		}
		assert(!_v.empty());
	}
	VectorPatternSource(vector<string> v,
	                    bool __revcomp = true,
	                    bool __reverse = false,
	                    const char *__dumpfile = NULL,
	                    size_t cur = 0,
	                    int __trim3 = 0,
	                    int __trim5 = 0) :
		TrimmingPatternSource<TStr>(__reverse, __dumpfile, __trim3, __trim5),
		_revcomp(__revcomp), _cur(cur), _v()
	{
		for(size_t i = 0; i < v.size(); i++) {
			TStr s(v[i]);
			_v.push_back(s);
			if(_revcomp) {
				_v.push_back(reverseComplement(s));
			}
		}
		assert(!_v.empty());
	}
	virtual bool nextIsReverseComplement() {
		return _revcomp && (_cur & 1) == 1;
	}
	virtual ~VectorPatternSource() { }
	virtual void nextPatternImpl(TStr& s, string& qual, string& name) {
		s = _v[_cur++];
	}
	virtual bool hasMorePatterns() {
		return _cur < _v.size();
	}
	virtual void reset() {
		TrimmingPatternSource<TStr>::reset();
		_cur = 0;
	}
private:
	bool _revcomp;
	size_t _cur;
	vector<TStr> _v;
};

/**
 * Supports reversing all strings as they're read in.  Also supports
 * returning reverse complements interspersed with forward versions of
 * patterns.
 */
template <typename TStr>
class BufferedFilePatternSource : public TrimmingPatternSource<TStr> {
public:
	BufferedFilePatternSource(const vector<string>& infiles,
	                          bool __gz = false,
	                          bool __revcomp = true,
	                          bool __reverse = false,
	                          const char *__dumpfile = NULL,
	                          int __trim3 = 0,
	                          int __trim5 = 0) :
		TrimmingPatternSource<TStr>(__reverse, __dumpfile, __trim3, __trim5),
		_infiles(infiles),
		_gz(__gz),
		_revcomp(__revcomp),
		_filecur(0),
		_aCur(true),
		_a(),
		_aRc(),
		_b(),
		_bRc(),
		_tmp(),
		_tmpRc(),
		_fw(true),
		_in(NULL)
	{
		assert_gt(infiles.size(), 0);
		open(); _filecur++;
	}
	virtual ~BufferedFilePatternSource() {
		if(_gz) gzclose(_gzin);
		else fclose(_in);
	}
	/// Return the next pattern from the file
	virtual void nextPatternImpl(TStr& s, string& qual, string& name) {
		assert(hasMorePatterns());
		readNext(s, qual, name);
	}
	/// Return true iff the next call to nextPattern() will succeed
	virtual bool hasMorePatterns() {
		if(_revcomp && !_fw) return true; // still a reverse comp to dish out
		if(!empty(_tmp)) return true;    // still another pattern to dish out
		this->read(&_tmp, 
				   (_revcomp? &_tmpRc : NULL), 
				   &_tmpQual,
				   (_revcomp? &_tmpRQual : NULL),
				   &_tmpName);
		return !empty(_tmp) || _filecur < _infiles.size(); // still another pattern to dish out
	}
	/// Return true iff the next will be a reverse complement
	virtual bool nextIsReverseComplement() {
		return _revcomp && !_fw;
	}
	/// Reset state so that we read start reading again from the
	/// beginning of the first file
	virtual void reset() {
		TrimmingPatternSource<TStr>::reset();
		// Close current file
		if(_gz) gzclose(_gzin);
		else fclose(_in);
		_filecur = 0,
		_fw = true; // dish forward next
		_aCur = true; // currently dishing a
		clear(_a);   clear(_aRc); _aQual.clear(); _aRQual.clear(); _aName.clear();
		clear(_b);   clear(_bRc); _bQual.clear(); _bRQual.clear(); _bName.clear();
		clear(_tmp); clear(_tmpRc); _tmpQual.clear(); _tmpRQual.clear(); _tmpName.clear();
		open(); _filecur++;
	}
protected:
	/// Read another pattern from the input file; this is overridden
	/// to deal with specific file formats
	virtual void read(TStr* dst, TStr* rcDst, string* qual, string* rqual, string* name) = 0;
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() = 0;
	/// Swap "current" status between _a/_aRc and _b/_bRc
	void swap() {
		_aCur = !_aCur;
	}
	/// Swap "current" with tmp
	void swapCur() {
		TStr tmp;
		tmp = (_aCur? _a : _b);
		if(_aCur) _a = _tmp; else _b = _tmp;
		_tmp = tmp;
		tmp = (_aCur? _aRc : _bRc);
		if(_aCur) _aRc = _tmpRc; else _bRc = _tmpRc;
		_tmpRc = tmp;
	}
	
	void swapCurQual() {
		string tmp;
		tmp = (_aCur? _aQual : _bQual);
		if(_aCur) _aQual = _tmpQual; else _bQual = _tmpQual;
		_tmpQual = tmp;
		tmp = (_aCur? _aRQual : _bRQual);
		if(_aCur) _aRQual = _tmpRQual; else _bRQual = _tmpRQual;
		_tmpRQual = tmp;
	}
	
	void swapCurName() {
		string tmp;
		tmp = (_aCur ? _aName : _bName);
		if (_aCur) _aName = _tmpName; else _bName = _tmpName;
		_tmpName = tmp;
	}
	
	/// Return "current" forward-oriented pattern
	TStr& cur() {
		return _aCur? _a : _b;
	}
	/// Return "current" reverse-complemented pattern
	TStr& curRc() {
		return _aCur? _aRc : _bRc;
	}
	
	string& curQual() {
		return _aCur? _aQual: _bQual;
	}
	
	string& curRQual() {
		return _aCur? _aRQual: _bRQual;
	}
	
	string& curName() {
		return _aCur? _aName: _bName;
	}
	/// Read in next pattern
	virtual void readNext(TStr& s, string& qual, string& name) {
		s = cur();
		qual = curQual();
		name = curName();
		if(_fw) {
			clear(cur());
			clear(curRc());
			curQual().clear();
			curRQual().clear();
			curName().clear();
			if(!empty(_tmp)) {
				// Just swap
				swapCur();
				swapCurQual();
				swapCurName();
				s = cur();
				qual = curQual();
				name = curName();
			} else {
				// Read in a fresh pair
				this->read(&cur(), 
						   (_revcomp? &curRc() : NULL), 
						   &curQual(), 
						   (_revcomp? &curRQual() : NULL),
						   &curName());
			}
			assert(empty(_tmp));
			assert(empty(_tmpRc));
			assert(empty(_tmpQual));
			assert(empty(_tmpName));
		} else {
			assert(empty(_tmp));
			assert(empty(_tmpRc));
			assert(empty(_tmpQual));
			assert(empty(_tmpName));
			s = curRc();
			qual = curRQual();
			name = curName();
			assert(!empty(s));
			swap();
		}
		assert(empty(_tmp));
		assert(empty(_tmpRc));
		assert(empty(_tmpQual));
		assert(empty(_tmpName));
		// If we exhausted all of the patterns in this file, go on to
		// the next one
		assert(!seqan::empty(s) || _fw);
		if(seqan::empty(s) && _filecur < _infiles.size()) {
			assert(_fw);
			// Close current file
			if(_gz) gzclose(_gzin);
			else fclose(_in);
			// Open next file
			open(); _filecur++;
			this->resetForNextFile(); // reset state to handle a fresh file
			this->read(&cur(), 
					   (_revcomp? &curRc() : NULL), 
					   &curQual(), 
					   (_revcomp? &curRQual() : NULL),
					   &curName());
			
			s = cur();
			qual = curQual();
			name = curName();
			assert(!empty(s));
		}
		// if !_revcomp, _fw always remains true
		if(_revcomp) _fw = !_fw;
		return;
	}
	void open() {
		if(_gz) {
			// Open gzipped input file
			if((_gzin = gzopen(_infiles[_filecur].c_str(), "r")) == NULL) {
				throw runtime_error("Could not open gzipped sequence file");
			}
		} else {
			// Open input file
			if((_in = fopen(_infiles[_filecur].c_str(), "r")) == NULL) {
				throw runtime_error("Could not open sequence file");
			}
			// Associate large input buffer with FILE *in
			if(setvbuf(_in, _buf, _IOFBF, 256 * 1024) != 0) {
				throw runtime_error("Could not create input buffer for sequence file");
			}
		}
	}
	const vector<string>& _infiles; // input filenames
	bool _gz;     // whether input file/files are gzipped
	bool _revcomp;   // whether to calculate reverse complements
	size_t _filecur; // index into _infiles of next file to read
	bool _aCur;   // true iff _a is the "current" pattern (one that an
	              // outsider might hold a reference to); otherwise _b
	              // is the current pattern and _a will be the target
	              // of the next read()
	TStr _a;      // pattern (might be either current or next)
	TStr _aRc;    // pattern (reverse-comp)
	string _aQual;// pattern (qualities)
	string _aRQual; //pattern (reversed qualities)
	string _aName; //pattern (sequence name)
	TStr _b;      // next-up buffer
	TStr _bRc;    // next-up buffer (reverse-comp)
	string _bQual;// next-up buffer (qualities)
	string _bRQual; //next-up buffer (reversed qualities)
	string _bName; //pattern (sequence name)
	TStr _tmp;      // next-up buffer
	TStr _tmpRc;    // next-up buffer (reverse-comp)
	string _tmpQual; // next-up buffer (qualities)
	string _tmpRQual; //next-up buffer (reversed qualities)
	string _tmpName; //pattern (sequence name)
	bool _fw;     // whether reverse complement should be returned next
	FILE *_in;    // file to read patterns from
	gzFile _gzin; // file to read patterns from
private:
	char _buf[256 * 1024]; // (large) input buffer
};

/**
 * 
 */
template <typename TStr>
class FastaPatternSource : public BufferedFilePatternSource<TStr> {
public:
	typedef typename Value<TStr>::Type TVal;
	FastaPatternSource(const vector<string>& infiles,
	                   bool __revcomp = true,
	                   bool __reverse = false,
	                   const char *__dumpfile = NULL,
	                   int __trim3 = 0,
	                   int __trim5 = 0) :
		BufferedFilePatternSource<TStr>(infiles, false, __revcomp, __reverse, __dumpfile, __trim3, __trim5),
		_first(true)
	{
		assert(this->hasMorePatterns());
	}
	virtual void reset() {
		_first = true;
		BufferedFilePatternSource<TStr>::reset();
	}
protected:
	int skipToEnd() {
		while(true) {
			int c = fgetc(this->_in); if(feof(this->_in)) return -1;
			if(c == '\n' || c == '\r') {
				while(c == '\n' || c == '\r') {
					c = fgetc(this->_in); if(feof(this->_in)) return -1;
				}
				// c now holds first character of next line
				return c;
			}
		}
	}
	/// Read another pattern from a FASTA input file
	// TODO: store default qualities in qual
	virtual void read(TStr *dst, TStr *rcDst, string* qual, string* rqual, string* name) {
		int c;
		assert(dst != NULL);
		assert(empty(*dst)); // caller should have cleared this
		// Pick off the first carat
		if(_first) {
			c = fgetc(this->_in); if(feof(this->_in)) return;
			assert(c == '>' || c == '#');
			_first = false;
		}
		// Skip to the end of the id line; if the next line is either
		// another id line or a comment line, keep skipping
		// TODO: grab name here, stick in *name
		do {
			if((c = skipToEnd()) == -1) return;
		} while (c == '>' || c == '#');
		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		int begin = 0;
		while(c != '>' && c != '#') {
			// Note: can't have a comment in the middle of a sequence,
			// though a comment can end a sequence
			if(isalpha(c)) {
				// 5' trim
				if(begin++ >= this->_trim5) {
					appendChar(dst, rcDst, (char)c);
					assert_lt((int)(*dst)[length(*dst)-1], 4);
				}
			}
			c = fgetc(this->_in);
			if(feof(this->_in)) break;
		}
		resize(*dst, length(*dst) - this->_trim3);
		if(rcDst != NULL) {
			resize(*rcDst, length(*rcDst) - this->_trim3);
			::reverse(*rcDst);
		}
	}
	virtual void resetForNextFile() {
		_first = true;
	}
private:
	bool _first;
};

/**
 * Read a FASTQ-format file.
 * See: http://maq.sourceforge.net/fastq.shtml
 */
template <typename TStr>
class FastqPatternSource : public BufferedFilePatternSource<TStr> {
public:
	typedef typename Value<TStr>::Type TVal;
	FastqPatternSource(const vector<string>& infiles,
	                   bool __revcomp = true,
	                   bool __reverse = false,
	                   const char *__dumpfile = NULL,
	                   int __trim3 = 0,
	                   int __trim5 = 0,
					   bool solexa_quals = false) :
		BufferedFilePatternSource<TStr>(infiles, false, __revcomp, __reverse, __dumpfile, __trim3, __trim5),
		_first(true), _solexa_quals(solexa_quals)
	{
		assert(this->hasMorePatterns());
		
		for (int l = 0; l != 128; ++l) {
			_table[l] = (int)(10.0 * log(1.0 + pow(10.0, (l - 64) / 10.0)) / log(10.0) + .499);
			if (_table[l] >= 63) _table[l] = 63;
			if (_table[l] == 0) _table[l] = 1;
		}
	}
	virtual void reset() {
		_first = true;
		BufferedFilePatternSource<TStr>::reset();
	}
protected:
	/// Read another pattern from a FASTQ input file
	virtual void read(TStr *dst, TStr *rcDst, string* qual, string* rqual, string* name) {
		int c;
		assert(dst != NULL);
		assert(empty(*dst)); // caller should have cleared this
		assert(qual->length() == 0);
		// Pick off the first at
		if(_first) {
			c = fgetc(this->_in); if(feof(this->_in)) return;
			assert_eq('@', c);
			_first = false;
		}
		// Read to the end of the id line, sticking everything after the '@'
		// into *name
		while(true) {
			c = fgetc(this->_in); if(feof(this->_in)) return;
			if(c == '\n' || c == '\r') {
				while(c == '\n' || c == '\r') {
					c = fgetc(this->_in); if(feof(this->_in)) return;
				}
				break;
			}
			name->push_back(c);
		}
		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		int begin = 0;
		while(c != '+') {
			if(isalpha(c)) {
				if(begin++ >= this->_trim5) {
					appendChar(dst, rcDst, (char)c);
					assert_lt((int)(*dst)[length(*dst)-1], 4);
				}
			}
			c = fgetc(this->_in);
			if(feof(this->_in)) break;
		}
		// Trim dst and rcDst on 3' end and reverse rcDst
		resize(*dst, length(*dst) - this->_trim3);
		if(rcDst != NULL) {
			resize(*rcDst, length(*rcDst) - this->_trim3);
			::reverse(*rcDst);
		}
		
		// Chew up the optional name on the '+' line
		while(c == '+' || (c != '\n' && c != '\r'))
		{
			c = fgetc(this->_in); if(feof(this->_in)) return;
		}
		
		// Move to the line after the + line, which has the qualities
		do {
			c = fgetc(this->_in); if(feof(this->_in)) return;
		} while(c == '\n' && c == '\r');
		
		// Now read the qualities
		begin = 0;
		
		while(true) {
			if(begin >= length(*dst)|| feof(this->_in)) {
				
				if (_solexa_quals)
				{
					// TODO: refactor!
					vector<string> s_quals;
					tokenize(*qual, " ", s_quals);
					qual->clear();
					for (unsigned int j = 0; j < s_quals.size(); ++j)
					{
						int sQ = atoi(s_quals[j].c_str());
						int pQ = (int)(10.0 * log(1.0 + pow(10.0, sQ / 10.0)) / log(10.0) + .499);
						qual->push_back((char)(pQ + 33));
					}
					s_quals.clear();
					if (rqual != NULL)
					{
						tokenize(*rqual, " ", s_quals);
						rqual->clear();
						for (unsigned int j = 0; j < s_quals.size(); ++j)
						{
							int sQ = atoi(s_quals[j].c_str());
							int pQ = (int)(10.0 * log(1.0 + pow(10.0, sQ / 10.0)) / log(10.0) + .499);
							rqual->push_back((char)(pQ + 33));
						}
					}
				}
				
				qual->resize(qual->length() - this->_trim3);
				if (rqual != NULL) {
					rqual->resize(rqual->length() - this->_trim3);
					std::reverse(rqual->begin(), rqual->end());
				}
				
				// chew up any qualities beyond length(*dst)
				while((c != '\n' && c != '\r'))
				{
					c = fgetc(this->_in); if(feof(this->_in)) return;
				}
				
				// Skip additional linebreak chars
				while(c == '\n' || c == '\r') {
					c = fgetc(this->_in); if(feof(this->_in)) return;
				}
				break;
			}
			else if (c != '\r' && c != '\n')
			{
				if (begin++ >= this->_trim5){
					qual->push_back(c);
					if (rqual != NULL)
						rqual->push_back(c);
				}	
			}
			c = fgetc(this->_in); if(feof(this->_in)) return;
		}
	}
	virtual void resetForNextFile() {
		_first = true;
	}
private:
	bool _first;
	bool _solexa_quals;
	int _table[128];
};

/**
 * 
 */
template <typename TStr>
class SolexaPatternSource : public BufferedFilePatternSource<TStr> {
public:
	typedef typename Value<TStr>::Type TVal;
	SolexaPatternSource(const vector<string>& infiles,
	                    bool __revcomp = true,
	                    bool __reverse = false,
	                    const char * __dumpfile = NULL,
	                    int __trim3 = 0,
	                    int __trim5 = 0) :
		BufferedFilePatternSource<TStr>(infiles, false, __revcomp, __reverse, __dumpfile, __trim3, __trim5),
		_savedC(-1)
	{
		assert(this->hasMorePatterns());
	}
	virtual void reset() {
		_savedC = -1;
		BufferedFilePatternSource<TStr>::reset();
	}
protected:
	//TODO: stick default quals in qual, rqual
	virtual void read(TStr *dst, TStr *rcDst, string* qual, string* rqual, string* name) {
		assert(dst != NULL);
		assert(empty(*dst)); // caller should have cleared this
		int c = (_savedC == -1) ? fgetc(this->_in) : _savedC;  
		if(feof(this->_in)) return;
		int begin = 0;
		while(c != '\n' && c != '\r') {
			if(isalpha(c)) {
				if(begin++ >= this->_trim5) {
					appendChar(dst, rcDst, (char)c);
					assert_lt((int)(*dst)[length(*dst)-1], 4);
				}
			}
			c = fgetc(this->_in);
			if(feof(this->_in)) break;
		}
		// Trim dst, rcDst; reverse rcDst
		resize(*dst, length(*dst) - this->_trim3);
		if(rcDst != NULL) {
			resize(*rcDst, length(*rcDst) - this->_trim3);
			::reverse(*rcDst);
		}
		// Skip to next line
		while(c == '\n' || c == '\r') {
			c = fgetc(this->_in); if(feof(this->_in)) return;
		}
		_savedC = c;
	}
	virtual void resetForNextFile() {
		_savedC = -1;
	}
private:
	int _savedC;
};

/**
 * 
 */
template <typename TStr>
class BfqPatternSource : public BufferedFilePatternSource<TStr> {
public:
	typedef typename Value<TStr>::Type TVal;
	BfqPatternSource(const vector<string>& infiles,
	                 bool __revcomp = true,
	                 bool __reverse = false,
	                 const char *__dumpfile = NULL,
	                 int __trim3 = 0,
	                 int __trim5 = 0) :
	BufferedFilePatternSource<TStr>(infiles, true, __revcomp, __reverse, __dumpfile, __trim3, __trim5)
	{
		assert(this->hasMorePatterns());
	}
	virtual void reset() {
		BufferedFilePatternSource<TStr>::reset();
	}
protected:
	/// Read another pattern from a FASTA input file
	virtual void read(TStr *dst, TStr *rcDst, string* qual, string* rqual, string* name) {
		typedef typename seqan::Value<TStr>::Type TVal;
		assert(dst != NULL);
		static char _name[2048]; // buffer for .bfq read names
		static unsigned char seq[2048];  // buffer 
		int len;
		assert(empty(*dst));
		// Read name length
		if(gzread(this->_gzin, &len, sizeof(int)) != sizeof(int)) return;
		if(len > 2047) {
			throw std::runtime_error(
				"One or more .bfq read names are longer than 2047 characters");
		}
		// Read name
		if(gzread(this->_gzin, _name, len) != len) return;
		_name[len] = '\0'; // TODO: do something with name (we just ignore it)
		*name = _name;
		// Read sequence length
		if(gzread(this->_gzin, &len, 4) != 4) return;
		if(len > 2048) {
			throw std::runtime_error(
				"One or more .bfq read sequences are longer than 2048 bases");
		}
		uint32_t tlen = len - this->_trim5 - this->_trim3;
		resize(*dst, tlen, Exact());
		if(rcDst != NULL) {
			resize(*rcDst, tlen, Exact());
		}
		if(gzread(this->_gzin, seq, len) != len) return;
		for(int i = this->_trim5; i < len - this->_trim3; i++) {
			char c = (char)(seq[i] >> 6);
			appendCharResized(dst, rcDst, i - this->_trim5, tlen, c);
		}
	}
	virtual void resetForNextFile() { }
};

#endif /*PAT_H_*/
