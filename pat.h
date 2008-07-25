#ifndef PAT_H_
#define PAT_H_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstring>
#include <ctype.h>
#include <zlib.h>
#include <fstream>
#include <seqan/sequence.h>
#include "alphabet.h"
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
	virtual void nextPattern(TStr** s, String<char>** qual, String<char>** name) {
		nextPatternImpl(s, qual, name);
		// Reverse it, if desired
		if(_reverse) {
			size_t len = length(**s);
			resize(_revtmp, len, Exact());
			if(*qual != NULL) {
				resize(_revqual, len);
			}
			// tmp = reverse of s
			for(size_t i = 0; i < len; i++)
			{
				_revtmp[i] = (**s)[len-i-1];
				if(*qual != NULL) {
					_revqual[i] = (**qual)[len-i-1];
				}
			}
			// Output it, if desired
			if(_dumpfile != NULL) {
				dump(_out, _revtmp,
				     ((*qual) == NULL) ? String<char>("NULL") : _revqual,
				     ((*name) == NULL) ? String<char>("NULL") : (**name));
			}
			(*s) = &_revtmp;
			if(*qual != NULL) {
				(*qual) = &_revqual;
			}
			return;
		}
		// Output it, if desired
		if(_dumpfile != NULL) {
			dump(_out, (**s),
			     ((*qual) == NULL) ? String<char>("NULL") : (**qual),
			     ((*name) == NULL) ? String<char>("NULL") : (**name));
		}
	}
	virtual void skipPattern() {
		TStr         *tmp1;
		String<char> *tmp2;
		String<char> *tmp3;
		nextPattern(&tmp1, &tmp2, &tmp3);
	}
	virtual void nextPatternImpl(TStr**, String<char>**, String<char>**) = 0;
	virtual bool hasMorePatterns() = 0;
	virtual bool nextIsReverseComplement() = 0;
	virtual void reset() { }
	const char *dumpfile() const { return _dumpfile; }
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
protected:
	virtual void dump(ostream& out,
	                  const TStr& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << name << ": " << seq << " " << qual << endl;
	}
private:
	bool _reverse;         // reverse patterns before returning them
	TStr _revtmp;          // temporary buffer for reversed patterns
	String<char> _revqual;
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
	VectorPatternSource(const vector<string>& v,
	                    bool __revcomp = true,
	                    bool __reverse = false,
	                    const char *__dumpfile = NULL,
	                    size_t cur = 0,
	                    int __trim3 = 0,
	                    int __trim5 = 0) :
		TrimmingPatternSource<TStr>(false, __dumpfile, __trim3, __trim5),
		_revcomp(__revcomp), _reverse(__reverse), _cur(cur), _v(),
		_quals(), _vrev(), _qualsrev()
	{
		for(size_t i = 0; i < v.size(); i++) {
			vector<string> ss;
			tokenize(v[i], ":", ss);
			assert_gt(ss.size(), 0);
			assert_leq(ss.size(), 2);
			// Initialize s
			string s = ss[0];
			//  Initialize vq
			string vq;
			if(ss.size() == 2) {
				vq = ss[1];
			}
			while(vq.length() < length(s)) {
				vq.push_back('I');
			}
			if(vq.length() > length(s)) {
				vq.erase(length(s));
			}
			assert_eq(vq.length(), length(s));
			_v.push_back(s);
			_quals.push_back(vq);
			{
				_vrev.push_back(s);
				::reverse(_vrev.back());
				_qualsrev.push_back(vq);
				::reverse(_qualsrev.back());
			}
			if(_revcomp) {
				_v.push_back(reverseComplement(TStr(s)));
				_quals.push_back(vq);
				::reverse(_quals.back());
				{
					_vrev.push_back(reverseComplement(TStr(s)));
					::reverse(_vrev.back());
					_qualsrev.push_back(vq);
				}
			}
		}
		assert(!_v.empty());
		assert_eq(_v.size(), _vrev.size());
		assert_eq(_v.size(), _quals.size());
		assert_eq(_v.size(), _qualsrev.size());
	}
	virtual bool nextIsReverseComplement() {
		return _revcomp && (_cur & 1) == 1;
	}
	virtual ~VectorPatternSource() { }
	virtual void nextPatternImpl(TStr** s, String<char>** qual, String<char>** name) {
		assert(hasMorePatterns());
		if(!_reverse) {
			(*s) = &(_v[_cur]);
			(*qual) = &(_quals[_cur]);
		} else {
			(*s) = &(_vrev[_cur]);
			(*qual) = &(_qualsrev[_cur]);
		}
		_cur++;
	}
	virtual bool hasMorePatterns() {
		return _cur < _v.size();
	}
	virtual void reset() {
		TrimmingPatternSource<TStr>::reset();
		_cur = 0;
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
private:
	bool _revcomp;
	bool _reverse;
	size_t _cur;
	vector<TStr>          _v;        // forward/rev-comp sequences
	vector<String<char> > _quals;    // quality values parallel to _v
	vector<TStr>          _vrev;     // reversed forward and rev-comp sequences 
	vector<String<char> > _qualsrev; // quality values parallel to _vrev
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
		_fw(true),
		_in(NULL),
		_gzin(),
		_aLen(0),
		_aTStr(),
		_aRcTStr(),
		_aQualStr(),
		_aRcQualStr(),
		_aNameLen(0),
		_aNameStr(),
		_bLen(0),
		_bTStr(),
		_bRcTStr(),
		_bQualStr(),
		_bRcQualStr(),
		_bNameLen(0),
		_bNameStr(),
		_tmpLen(0),
		_tmpTStr(),
		_tmpRcTStr(),
		_tmpQualStr(),
		_tmpRcQualStr(),
		_tmpNameLen(0),
		_tmpNameStr(),
		_a(_buf1),
		_aRc(_rcBuf1),
		_aQual(_qualBuf1),
		_aRcQual(_rcQualBuf1),
		_aName(_nameBuf1),
		_b(_buf2),
		_bRc(_rcBuf2),
		_bQual(_qualBuf2),
		_bRcQual(_rcQualBuf2),
		_bName(_nameBuf2),
		_tmp(_buf3),
		_tmpRc(_rcBuf3),
		_tmpQual(_qualBuf3),
		_tmpRcQual(_rcQualBuf3),
		_tmpName(_nameBuf3)
	{
		_setBegin(_aTStr, (Dna*)_a); _setLength(_aTStr, 0);
		_setBegin(_aRcTStr, (Dna*)_aRc); _setLength(_aRcTStr, 0);
		_setBegin(_aQualStr, (char*)_aQual); _setLength(_aQualStr, 0);
		_setBegin(_aRcQualStr, (char*)_aRcQual); _setLength(_aRcQualStr, 0);
		_setBegin(_aNameStr, (char*)_aName); _setLength(_aNameStr, 0);
		_setBegin(_bTStr, (Dna*)_b); _setLength(_bTStr, 0);
		_setBegin(_bRcTStr, (Dna*)_bRc); _setLength(_bRcTStr, 0);
		_setBegin(_bQualStr, (char*)_bQual); _setLength(_bQualStr, 0);
		_setBegin(_bRcQualStr, (char*)_bRcQual); _setLength(_bRcQualStr, 0);
		_setBegin(_bNameStr, (char*)_bName); _setLength(_bNameStr, 0);
		_setBegin(_tmpTStr, (Dna*)_tmp); _setLength(_tmpTStr, 0);
		_setBegin(_tmpRcTStr, (Dna*)_tmpRc); _setLength(_tmpRcTStr, 0);
		_setBegin(_tmpQualStr, (char*)_tmpQual); _setLength(_tmpQualStr, 0);
		_setBegin(_tmpRcQualStr, (char*)_tmpRcQual); _setLength(_tmpRcQualStr, 0);
		_setBegin(_tmpNameStr, (char*)_tmpName); _setLength(_tmpNameStr, 0);
		assert_gt(infiles.size(), 0);
		open(); _filecur++;
	}
	virtual ~BufferedFilePatternSource() {
		if(_gz) gzclose(_gzin);
		else fclose(_in);
	}
	/// Return the next pattern from the file
	virtual void nextPatternImpl(TStr** s,
	                             String<char>** qual,
	                             String<char>** name)
	{
		assert(hasMorePatterns());
		readNext(s, qual, name);
	}
	/// Return true iff the next call to nextPattern() will succeed
	virtual bool hasMorePatterns() {
		if(_revcomp && !_fw) return true; // still a reverse comp to dish out
		if(_tmpLen > 0) return true;    // still another pattern to dish out
		this->read(_tmp,
		           _tmpTStr,
				   _tmpRc,
				   _tmpRcTStr,
				   _tmpQual,
				   _tmpQualStr,
				   _tmpRcQual,
				   _tmpRcQualStr,
				   &_tmpLen,
				   _tmpName,
				   _tmpNameStr,
				   &_tmpNameLen);
		return _tmpLen > 0 || _filecur < _infiles.size(); // still another pattern to dish out
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
		_aLen = 0;
		_aNameLen = 0;
		_bLen = 0;
		_bNameLen = 0;
		_tmpLen = 0;
		_tmpNameLen = 0;
		open(); _filecur++;
	}
protected:
	/// Read another pattern from the input file; this is overridden
	/// to deal with specific file formats
	virtual void read(char* dst,      // destination buf for sequence
	                  TStr& dstTStr,  // destination TStr for sequence
	                  char* rcDst,    // destination buf for reverse-comp of dst
	                  TStr& rcDstTStr,// destination TStr for reverse-comp of dst
	                  char* qual,     // destination buf for qualities
	                  String<char>& qualStr, // destination String for qualities
	                  char* rcQual,   // destination buf for reverse-comp qualities
	                  String<char>& rcQualStr, // destination String for reverse-comp quals
	                  size_t* dstLen, // length of sequence installed in dst
	                  char* name,     // destination buf for read name
	                  String<char>& nameStr, // destination String for name
	                  size_t* nameLen) = 0; // length of name
	/// Skip to the end of the current line; return the first character
	/// of the next line
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
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() = 0;
	/// Swap "current" status between _a/_aRc and _b/_bRc
	void swap() {
		_aCur = !_aCur;
	}
	/// Swap "current" with tmp
	void swapCur() {
		char   *tmp;
		char   *tmpRc;
		char   *tmpQual;
		char   *tmpRcQual;
		char   *tmpName;
		if(_aCur) {
			tmp       = _a;
			tmpRc     = _aRc;
			tmpQual   = _aQual;
			tmpRcQual = _aRcQual;
			tmpName   = _aName; 
			_a       = _tmp;       _setBegin(_aTStr,      (Dna*)begin(_tmpTStr));
			_aRc     = _tmpRc;     _setBegin(_aRcTStr,    (Dna*)begin(_tmpRcTStr));
			_aQual   = _tmpQual;   _setBegin(_aQualStr,   (char*)begin(_tmpQualStr));
			_aRcQual = _tmpRcQual; _setBegin(_aRcQualStr, (char*)begin(_tmpRcQualStr));
			_aName   = _tmpName;   _setBegin(_aNameStr,   (char*)begin(_tmpNameStr));
			_aLen = _tmpLen;
			_setLength(_aTStr, _aLen);
			_setLength(_aRcTStr, _aLen);
			_setLength(_aQualStr, _aLen);
			_setLength(_aRcQualStr, _aLen);
			_aNameLen = _tmpNameLen;
			_setLength(_aNameStr, _aNameLen);
		} else {
			tmp       = _b;
			tmpRc     = _bRc;
			tmpQual   = _bQual;
			tmpRcQual = _bRcQual;
			tmpName   = _bName;
			_b       = _tmp;       _setBegin(_bTStr,      (Dna*)begin(_tmpTStr));
			_bRc     = _tmpRc;     _setBegin(_bRcTStr,    (Dna*)begin(_tmpRcTStr));
			_bQual   = _tmpQual;   _setBegin(_bQualStr,   (char*)begin(_tmpQualStr));
			_bRcQual = _tmpRcQual; _setBegin(_bRcQualStr, (char*)begin(_tmpRcQualStr));
			_bName   = _tmpName;   _setBegin(_bNameStr,   (char*)begin(_tmpNameStr));
			_bLen = _tmpLen;
			_setLength(_bTStr, _bLen);
			_setLength(_bRcTStr, _bLen);
			_setLength(_bQualStr, _bLen);
			_setLength(_bRcQualStr, _bLen);
			_bNameLen = _tmpNameLen;
			_setLength(_bNameStr, _bNameLen);
		}
		//_tmp = tmp; _setBegin(_tmpTStr, (Dna*)_tmp);
		//_tmpRc = tmpRc; _setBegin(_tmpRcTStr, (Dna*)_tmpRc);
		//_tmpQual = tmpQual; _setBegin(_tmpQualStr, (char*)_tmpQual);
		//_tmpRcQual = tmpRcQual; _setBegin(_tmpRcQualStr, (char*)_tmpRcQual);
		_tmpLen = 0; // clear
		//_tmpName = tmpName; _setBegin(_tmpNameStr, (char*)_tmpName);
		_tmpNameLen = 0; // clear
	}
	
	/// Return "current" forward-oriented pattern
	char *cur() {
		return _aCur? _a : _b;
	}
	TStr *curStr() {
		return _aCur? &_aTStr : &_bTStr;
	}
	/// Return "current" reverse-complemented pattern
	char *curRc() {
		if(_revcomp) {
			return _aCur? _aRc : _bRc;
		} else {
			return NULL;
		}
	}
	TStr *curRcStr() {
		if(_revcomp) {
			return _aCur? &_aRcTStr : &_bRcTStr;
		} else {
			return NULL;
		}
	}
	char *curQual() {
		return _aCur? _aQual: _bQual;
	}
	String<char> *curQualStr() {
		return _aCur? &_aQualStr: &_bQualStr;
	}
	char *curRcQual() {
		if(_revcomp) {
			return _aCur? _aRcQual: _bRcQual;
		} else {
			return NULL;
		}
	}
	String<char> *curRcQualStr() {
		return _aCur? &_aRcQualStr: &_bRcQualStr;
	}
	size_t *curLen() {
		return _aCur? &_aLen : &_bLen;
	}
	char *curName() {
		return _aCur? _aName : _bName;
	}
	String<char> *curNameStr() {
		return _aCur? &_aNameStr : &_bNameStr;
	}
	size_t *curNameLen() {
		return _aCur? &_aNameLen : &_bNameLen;
	}
	void readCur() {
		if(_aCur) {
			this->read(_a,
			           _aTStr,
			           _aRc,
			           _aRcTStr,
			           _aQual,
			           _aQualStr,
			           _aRcQual,
			           _aRcQualStr,
			           &_aLen,
			           _aName,
			           _aNameStr,
			           &_aNameLen);
		} else {
			this->read(_b,
			           _bTStr,
			           _bRc,
			           _bRcTStr,
			           _bQual,
			           _bQualStr,
			           _bRcQual,
			           _bRcQualStr,
			           &_bLen,
			           _bName,
			           _bNameStr,
			           &_bNameLen);
		}
	}
	/// Read in next pattern
	virtual void readNext(TStr** s, String<char>** qual, String<char>** name) {
		(*s)    = curStr();
		(*qual) = curQualStr();
		(*name) = curNameStr();
		if(_fw) {
			if(_tmpLen > 0) {
				// Just swap
				swapCur();
				(*s)    = curStr();
				(*qual) = curQualStr();
				(*name) = curNameStr();
			} else {
				// Read in a fresh pair
				readCur();
			}
		} else {
			assert(curRc() != NULL);
			(*s)    = curRcStr();
			(*qual) = curRcQualStr();
			(*name) = curNameStr();
			assert(!empty(**s));
			swap();
		}
		assert_eq(0, _tmpLen);
		// If we exhausted all of the patterns in this file, go on to
		// the next one
		assert(!seqan::empty(**s) || _fw);
		if(seqan::empty(**s) && _filecur < _infiles.size()) {
			assert(_fw);
			// Close current file
			if(_gz) gzclose(_gzin);
			else fclose(_in);
			// Open next file
			open(); _filecur++;
			this->resetForNextFile(); // reset state to handle a fresh file
			// Read in a fresh pair
			readCur();
			(*s)    = curStr();
			(*qual) = curQualStr();
			(*name) = curNameStr();
			assert(!empty(**s));
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

	bool _fw;     // whether reverse complement should be returned next
	FILE *_in;    // file to read patterns from
	gzFile _gzin; // file to read patterns from

	// Pattern a
	size_t _aLen;
	TStr   _aTStr;         // host is set to _a
	TStr   _aRcTStr;       // host is set to _aRc
	String<char> _aQualStr;
	String<char> _aRcQualStr;
	size_t       _aNameLen;      // length of read name
	String<char> _aNameStr;      // if _aName changes, must assign to this
	
	// Pattern b
	size_t _bLen;
	TStr   _bTStr;         // host is set to _b
	TStr   _bRcTStr;       // host is set to _bRc
	String<char> _bQualStr;
	String<char> _bRcQualStr;
	size_t       _bNameLen;      // length of read name
	String<char> _bNameStr;
	
	// Pattern tmp (for hasMorePatterns())
	size_t _tmpLen;
	TStr   _tmpTStr;
	TStr   _tmpRcTStr;
	String<char> _tmpQualStr;
	String<char> _tmpRcQualStr;
	size_t _tmpNameLen;
	String<char> _tmpNameStr;
	
private:
	char* _a;
	char* _aRc;
	char* _aQual;   // quality values for forward-strand version
	char* _aRcQual; // quality values for reverse-comp version
	char* _aName;   // read name
	
	char* _b;
	char* _bRc;
	char* _bQual;   // quality values for forward-strand version
	char* _bRcQual; // quality values for reverse-comp version
	char* _bName;   // read name

	char* _tmp;
	char* _tmpRc;
	char* _tmpQual;   // next-up buffer (qualities)
	char* _tmpRcQual; // next-up buffer (reversed qualities)
	char* _tmpName;

	char _buf1[1024];
	char _buf2[1024];
	char _buf3[1024];
	char _rcBuf1[1024];
	char _rcBuf2[1024];
	char _rcBuf3[1024];
	char _qualBuf1[1024];
	char _qualBuf2[1024];
	char _qualBuf3[1024];
	char _rcQualBuf1[1024];
	char _rcQualBuf2[1024];
	char _rcQualBuf3[1024];
	char _nameBuf1[1024];
	char _nameBuf2[1024];
	char _nameBuf3[1024];

	char _buf[256 * 1024]; // (large) input buffer
};


static uint8_t charToDna[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	       /*    A     C           G */              
	/*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	       /*             T */
	/*  96 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 112 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

static uint8_t rcCharToDna[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
	       /*    A     C           G */              
	/*  80 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	       /*             T */
	/*  96 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 112 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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
		BufferedFilePatternSource<TStr>(infiles, false, __revcomp, false, __dumpfile, __trim3, __trim5),
		_first(true),
		_reverse(__reverse)
	{
		assert(this->hasMorePatterns());
	}
	virtual void reset() {
		_first = true;
		BufferedFilePatternSource<TStr>::reset();
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
protected:
	
	/// Read another pattern from a FASTA input file
	// TODO: store default qualities in qual
	virtual void read(char* dst,
	                  TStr& dstTStr,
	                  char *rcDst,
	                  TStr& rcDstTStr,
	                  char* qual,
	                  String<char>& qualTStr,
	                  char* rcQual,
	                  String<char>& rcQualTStr,
	                  size_t* dstLen,
	                  char* name,
	                  String<char>& nameStr,
	                  size_t* nameLen)
	{
		assert(dst != NULL);
		assert(qual != NULL);
		assert(dstLen != NULL);
		assert(name != NULL);
		assert(nameLen != NULL);
		
		int c;
		*dstLen = 0;
		*nameLen = 0;
		
		// Pick off the first carat
		if(_first) {
			c = fgetc(this->_in); if(c < 0) return;
			assert(c == '>' || c == '#');
			_first = false;
		}
		
		// Read to the end of the id line, sticking everything after the '>'
		// into *name
		while(true) {
			c = fgetc(this->_in); if(c < 0) return;
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					c = fgetc(this->_in); if(c < 0) return;
				}
				break;
			}
			name[(*nameLen)++] = c;
		}
		_setBegin(nameStr, (char*)name);
		_setLength(nameStr, *nameLen);
		
		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		int begin = 0;
		if(!_reverse) {
			while(c != '>' && c != '#') {
				// Note: can't have a comment in the middle of a sequence,
				// though a comment can end a sequence
				if(isalpha(c) && begin++ >= this->_trim5) {
					if(c == 'N' || c == 'n') c = 'A';
					dst[(*dstLen)] = charToDna[c];
					rcDst[1024-(*dstLen)-1] = rcCharToDna[c];
					(*dstLen)++;
				}
				if((c = fgetc(this->_in)) < 0) break;
			}
			(*dstLen) -= this->_trim3;
			_setBegin(dstTStr, (Dna*)dst);
			_setLength(dstTStr, (*dstLen));
			if(rcDst != NULL) {
				_setBegin(rcDstTStr, (Dna*)&rcDst[1024-(*dstLen)]);
				_setLength(rcDstTStr, (*dstLen));
			}
		} else {
			while(c != '>' && c != '#') {
				// Note: can't have a comment in the middle of a sequence,
				// though a comment can end a sequence
				if(isalpha(c) && begin++ >= this->_trim5) {
					if(c == 'N' || c == 'n') c = 'A';
					dst[1024-(*dstLen)-1] = charToDna[c];
					rcDst[(*dstLen)] = rcCharToDna[c];
					(*dstLen)++;
				}
				if((c = fgetc(this->_in)) < 0) break;
			}
			(*dstLen) -= this->_trim3;
			_setBegin(dstTStr, (Dna*)&dst[1024-(*dstLen)]);
			_setLength(dstTStr, (*dstLen));
			if(rcDst != NULL) {
				_setBegin(rcDstTStr, (Dna*)rcDst);
				_setLength(rcDstTStr, (*dstLen));
			}
		}
	}
	virtual void resetForNextFile() {
		_first = true;
	}
	virtual void dump(ostream& out,
	                  const TStr& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << ">" << name << endl << seq << endl;
	}
private:
	bool _first;
	bool _reverse;
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
		BufferedFilePatternSource<TStr>(infiles, false, __revcomp, false, __dumpfile, __trim3, __trim5),
		_first(true), _reverse(__reverse), _solexa_quals(solexa_quals)
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
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
protected:
	/// Read another pattern from a FASTQ input file
	virtual void read(char* dst,
	                  TStr& dstTStr,
	                  char *rcDst,
	                  TStr& rcDstTStr,
	                  char* qual,
	                  String<char>& qualTStr,
	                  char* rcQual,
	                  String<char>& rcQualTStr,
	                  size_t* dstLen,
	                  char* name,
	                  String<char>& nameStr,
	                  size_t* nameLen)
	{
		assert(dst != NULL);
		assert(qual != NULL);
		assert(dstLen != NULL);
		assert(name != NULL);
		assert(nameLen != NULL);

		int c;
		*dstLen = 0;
		*nameLen = 0;
		
		// Pick off the first at
		if(_first) {
			c = fgetc(this->_in); if(c < 0) return;
			assert_eq('@', c);
			_first = false;
		}
		
		// Read to the end of the id line, sticking everything after the '@'
		// into *name
		while(true) {
			c = fgetc(this->_in); if(c < 0) return;
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					c = fgetc(this->_in); if(c < 0) return;
				}
				break;
			}
			name[(*nameLen)++] = c;
		}
		_setBegin(nameStr, (char*)name);
		_setLength(nameStr, *nameLen);
		
		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		int charsRead = 0;
		if(!_reverse) {
			while(c != '+') {
				if(isalpha(c) && charsRead >= this->_trim5) {
					if(c == 'N' || c == 'n') c = 'A';
					dst[(*dstLen)] = charToDna[c];
					rcDst[1024-(*dstLen)-1] = rcCharToDna[c];
					charsRead++; (*dstLen)++;
				}
				c = fgetc(this->_in);
				if(c < 0) break;
			}
			(*dstLen) -= this->_trim3;
			_setBegin(dstTStr, (Dna*)dst);
			_setLength(dstTStr, (*dstLen));
			if(rcDst != NULL) {
				_setBegin(rcDstTStr, (Dna*)&rcDst[1024-(*dstLen)]);
				_setLength(rcDstTStr, (*dstLen));
			}
		} else {
			while(c != '+') {
				if(isalpha(c) && charsRead >= this->_trim5) {
					if(c == 'N' || c == 'n') c = 'A';
					dst[1024-(*dstLen)-1] = charToDna[c];
					rcDst[(*dstLen)] = rcCharToDna[c];
					charsRead++; (*dstLen)++;
				}
				c = fgetc(this->_in);
				if(c < 0) break;
			}
			(*dstLen) -= this->_trim3;
			_setBegin(dstTStr, (Dna*)&dst[1024-(*dstLen)]);
			_setLength(dstTStr, (*dstLen));
			if(rcDst != NULL) {
				_setBegin(rcDstTStr, (Dna*)rcDst);
				_setLength(rcDstTStr, (*dstLen));
			}
		}
		
		// Chew up the optional name on the '+' line
		while(c == '+' || (c != '\n' && c != '\r'))
		{
			c = fgetc(this->_in); if(c < 0) return;
		}

		// Now read the qualities
		size_t qualsRead = 0;
		
		//Move to the line after the + line, which has the qualities
//		while(c == '\n' && c == '\r'){
//			c = fgetc(this->_in); if(feof(this->_in)) return;
//		} 
		
		//qual->clear();
		if (_solexa_quals)
		{
			char buf[1024];
				
			while (fgets(buf, sizeof(buf), this->_in) && qualsRead < (*dstLen) + this->_trim5)
			{
				char* nl = strrchr(buf, '\n');
				if (nl) *nl = 0;
				
				vector<string> s_quals;
				tokenize(string(buf), " ", s_quals);
				
				if(!_reverse) {
					for (unsigned int j = 0; j < s_quals.size(); ++j)
					{
						int sQ = atoi(s_quals[j].c_str());
						int pQ = (int)(10.0 * log(1.0 + pow(10.0, sQ / 10.0)) / log(10.0) + .499);
						if (qualsRead >= (size_t)this->_trim5)
						{
							size_t off = qualsRead - this->_trim5;
							c = (char)(pQ + 33);
							qual[off] = c;
							if (rcQual != NULL) {
								rcQual[1024 - off - 1] = c;
							}
						}
						++qualsRead;
					}
				} else {
					for (unsigned int j = 0; j < s_quals.size(); ++j)
					{
						int sQ = atoi(s_quals[j].c_str());
						int pQ = (int)(10.0 * log(1.0 + pow(10.0, sQ / 10.0)) / log(10.0) + .499);
						if (qualsRead >= (size_t)this->_trim5)
						{
							size_t off = qualsRead - this->_trim5;
							c = (char)(pQ + 33);
							qual[1024 - off - 1] = c;
							if (rcQual != NULL) {
								rcQual[off] = c;
							}
						}
						++qualsRead;
					}
				}
			} // done reading Solexa quality lines
			if(!_reverse) {
				_setBegin(qualTStr, (char*)qual);
				_setLength(qualTStr, (*dstLen));
				if(rcQual != NULL) {
					_setBegin(rcQualTStr, (char*)&rcQual[1024-(*dstLen)]);
					_setLength(rcQualTStr, (*dstLen));
				}
			} else {
				_setBegin(qualTStr, (char*)&qual[1024-(*dstLen)]);
				_setLength(qualTStr, (*dstLen));
				if(rcQual != NULL) {
					_setBegin(rcQualTStr, (char*)&rcQual[(*dstLen)]);
					_setLength(rcQualTStr, (*dstLen));
				}
			}
		}
		else
		{
			if(!_reverse) {
				while(qualsRead < (*dstLen) + this->_trim5 && c >= 0) {
					c = fgetc(this->_in);
					if (c != '\r' && c != '\n')
					{
						if (qualsRead >= (size_t)this->_trim5) {
							size_t off = qualsRead - this->_trim5;
							assert_geq(c, 33);
							assert_leq(c, 73);
							qual[off] = c;
							if (rcQual != NULL) {
								rcQual[1024 - off - 1] = c;
							}
							qualsRead++;
						}	
					}
				}
				assert_eq(qualsRead, (*dstLen) + this->_trim5);
				_setBegin(qualTStr, (char*)qual);
				_setLength(qualTStr, (*dstLen));
				if(rcQual != NULL) {
					_setBegin(rcQualTStr, (char*)&rcQual[1024-(*dstLen)]);
					_setLength(rcQualTStr, (*dstLen));
				}
			} else {
				while(qualsRead < (*dstLen) + this->_trim5 && c >= 0) {
					c = fgetc(this->_in);
					if (c != '\r' && c != '\n')
					{
						if (qualsRead >= (size_t)this->_trim5) {
							size_t off = qualsRead - this->_trim5;
							assert_geq(c, 33);
							assert_leq(c, 73);
							qual[1024 - off - 1] = c;
							if (rcQual != NULL) {
								rcQual[off] = c;
							}
							qualsRead++;
						}	
					}
				}
				assert_eq(qualsRead, (*dstLen) + this->_trim5);
				_setBegin(qualTStr, (char*)&qual[1024-(*dstLen)]);
				_setLength(qualTStr, (*dstLen));
				if(rcQual != NULL) {
					_setBegin(rcQualTStr, (char*)rcQual);
					_setLength(rcQualTStr, (*dstLen));
				}
			}
		}
		
		// Trim the 3' end of the quality string
		//qual->resize(qual->length() - this->_trim3);
		
		// Now make the reversed qual string if necessary
		//if (rqual != NULL) {
		//	rqual->resize(rqual->length() - this->_trim3);
		//	std::reverse(rqual->begin(), rqual->end());
		//}
		
		if (feof(this->_in))
			return;
		else
		{
			do {
				c = fgetc(this->_in);
			} while(c != '@' && c >= 0);
		}
//		
//		
//		while(true) {
//			if(begin >= length(*dst)|| feof(this->_in)) {
//				
//				if (_solexa_quals)
//				{
//					// TODO: refactor!
//					vector<string> s_quals;
//					tokenize(*qual, " ", s_quals);
//					qual->clear();
//					for (unsigned int j = 0; j < s_quals.size(); ++j)
//					{
//						int sQ = atoi(s_quals[j].c_str());
//						int pQ = (int)(10.0 * log(1.0 + pow(10.0, sQ / 10.0)) / log(10.0) + .499);
//						qual->push_back((char)(pQ + 33));
//					}
//					s_quals.clear();
//					if (rqual != NULL)
//					{
//						tokenize(*rqual, " ", s_quals);
//						rqual->clear();
//						for (unsigned int j = 0; j < s_quals.size(); ++j)
//						{
//							int sQ = atoi(s_quals[j].c_str());
//							int pQ = (int)(10.0 * log(1.0 + pow(10.0, sQ / 10.0)) / log(10.0) + .499);
//							rqual->push_back((char)(pQ + 33));
//						}
//					}
//				}
//				
//				qual->resize(qual->length() - this->_trim3);
//				if (rqual != NULL) {
//					rqual->resize(rqual->length() - this->_trim3);
//					std::reverse(rqual->begin(), rqual->end());
//				}
//				
//				// chew up any qualities beyond length(*dst)
//				//while((c != '\n' && c != '\r'))
////				{
////					c = fgetc(this->_in); if(feof(this->_in)) return;
////				}
////				
////				// Skip additional linebreak chars
////				while(c == '\n' || c == '\r') {
////					c = fgetc(this->_in); if(feof(this->_in)) return;
////				}
//				
//				// Chew up anything left over until we get to the next read
//				while(c != '@' && !feof(this->_in))
//				{
//					c = fgetc(this->_in); if(feof(this->_in)) return;
//				}
//				break;
//			}
//			else if (c != '\r' && c != '\n')
//			{
//				if (begin++ >= this->_trim5){
//					qual->push_back(c);
//					if (rqual != NULL)
//						rqual->push_back(c);
//				}	
//			}
//			c = fgetc(this->_in); if(feof(this->_in)) return;
//		}
	}
	virtual void resetForNextFile() {
		_first = true;
	}
	virtual void dump(ostream& out,
	                  const TStr& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
	}
private:
	bool _first;
	bool _reverse;
	bool _solexa_quals;
	int _table[128];
};

/**
 * 
 */
//template <typename TStr>
//class SolexaPatternSource : public BufferedFilePatternSource<TStr> {
//public:
//	typedef typename Value<TStr>::Type TVal;
//	SolexaPatternSource(const vector<string>& infiles,
//	                    bool __revcomp = true,
//	                    bool __reverse = false,
//	                    const char * __dumpfile = NULL,
//	                    int __trim3 = 0,
//	                    int __trim5 = 0) :
//		BufferedFilePatternSource<TStr>(infiles, false, __revcomp, __reverse, __dumpfile, __trim3, __trim5),
//		_savedC(-1)
//	{
//		assert(this->hasMorePatterns());
//	}
//	virtual void reset() {
//		_savedC = -1;
//		BufferedFilePatternSource<TStr>::reset();
//	}
//protected:
//	//TODO: stick default quals in qual, rqual
//	virtual void read(TStr *dst, TStr *rcDst, string* qual, string* rqual, string* name) {
//		assert(dst != NULL);
//		assert(empty(*dst)); // caller should have cleared this
//		int c = (_savedC == -1) ? fgetc(this->_in) : _savedC;  
//		if(feof(this->_in)) return;
//		int begin = 0;
//		while(c != '\n' && c != '\r') {
//			if(isalpha(c)) {
//				if(begin++ >= this->_trim5) {
//					appendChar(dst, rcDst, (char)c);
//					assert_lt((int)(*dst)[length(*dst)-1], 4);
//				}
//			}
//			c = fgetc(this->_in);
//			if(feof(this->_in)) break;
//		}
//		// Trim dst, rcDst; reverse rcDst
//		resize(*dst, length(*dst) - this->_trim3);
//		if(rcDst != NULL) {
//			resize(*rcDst, length(*rcDst) - this->_trim3);
//			::reverse(*rcDst);
//		}
//		// Skip to next line
//		while(c == '\n' || c == '\r') {
//			c = fgetc(this->_in); if(feof(this->_in)) return;
//		}
//		_savedC = c;
//	}
//	virtual void resetForNextFile() {
//		_savedC = -1;
//	}
//private:
//	int _savedC;
//};

/**
 * 
 */
//template <typename TStr>
//class BfqPatternSource : public BufferedFilePatternSource<TStr> {
//public:
//	typedef typename Value<TStr>::Type TVal;
//	BfqPatternSource(const vector<string>& infiles,
//	                 bool __revcomp = true,
//	                 bool __reverse = false,
//	                 const char *__dumpfile = NULL,
//	                 int __trim3 = 0,
//	                 int __trim5 = 0) :
//	BufferedFilePatternSource<TStr>(infiles, true, __revcomp, __reverse, __dumpfile, __trim3, __trim5)
//	{
//		assert(this->hasMorePatterns());
//	}
//	virtual void reset() {
//		BufferedFilePatternSource<TStr>::reset();
//	}
//protected:
//	/// Read another pattern from a FASTA input file
//	virtual void read(TStr *dst, TStr *rcDst, string* qual, string* rqual, string* name) {
//		typedef typename seqan::Value<TStr>::Type TVal;
//		assert(dst != NULL);
//		static char _name[2048]; // buffer for .bfq read names
//		static unsigned char seq[2048];  // buffer 
//		int len;
//		assert(empty(*dst));
//		// Read name length
//		if(gzread(this->_gzin, &len, sizeof(int)) != sizeof(int)) return;
//		if(len > 2047) {
//			throw std::runtime_error(
//				"One or more .bfq read names are longer than 2047 characters");
//		}
//		// Read name
//		if(gzread(this->_gzin, _name, len) != len) return;
//		_name[len] = '\0'; // TODO: do something with name (we just ignore it)
//		*name = _name;
//		// Read sequence length
//		if(gzread(this->_gzin, &len, 4) != 4) return;
//		if(len > 2048) {
//			throw std::runtime_error(
//				"One or more .bfq read sequences are longer than 2048 bases");
//		}
//		uint32_t tlen = len - this->_trim5 - this->_trim3;
//		resize(*dst, tlen, Exact());
//		if(rcDst != NULL) {
//			resize(*rcDst, tlen, Exact());
//		}
//		if(gzread(this->_gzin, seq, len) != len) return;
//		for(int i = this->_trim5; i < len - this->_trim3; i++) {
//			char c = (char)(seq[i] >> 6);
//			appendCharResized(dst, rcDst, i - this->_trim5, tlen, c);
//		}
//	}
//	virtual void resetForNextFile() { }
//};

#endif /*PAT_H_*/
