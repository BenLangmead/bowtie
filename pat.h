#ifndef PAT_H_
#define PAT_H_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstring>
#include <ctype.h>
#include <fstream>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "tokenize.h"
#include "random_source.h"
#include "spinlock.h"
#include "threading.h"

using namespace std;
using namespace seqan;

/// Wildcard policies
enum {
	NS_TO_NS    = 1, // Ns stay Ns and don't match anything
	NS_TO_AS    = 2, // Ns become As
	NS_TO_RANDS = 3  // Ns become random characters
};

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

/// Reverse a string in-place
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
 * A buffer for keeping all relevant information about a single read.
 * Each search thread has one.
 */
struct ReadBuf {
	~ReadBuf() {
		clearAll(); reset();
		_setBegin(patFw, NULL);
		_setBegin(patRc, NULL);
		_setBegin(qualFw, NULL);
		_setBegin(qualRc, NULL);
		_setBegin(name, NULL);
	}
	/// Point all Strings to the beginning of their respective buffers
	/// and set all lengths to 0
	void reset() {
		_setBegin(patFw,  (Dna5*)patBufFw);  _setLength(patFw, 0);  _setCapacity(patFw, 1024);
		_setBegin(patRc,  (Dna5*)patBufRc);  _setLength(patRc, 0);  _setCapacity(patRc, 1024);
		_setBegin(qualFw, (char*)qualBufFw); _setLength(qualFw, 0); _setCapacity(qualFw, 1024);
		_setBegin(qualRc, (char*)qualBufRc); _setLength(qualRc, 0); _setCapacity(qualRc, 1024);
		_setBegin(name,   (char*)nameBuf);   _setLength(name, 0);   _setCapacity(name, 1024);
	}
	void clearAll() {
		seqan::clear(patFw);
		seqan::clear(patRc);
		seqan::clear(qualFw);
		seqan::clear(qualRc);
		seqan::clear(name);
	}
	static const int BUF_SIZE = 1024;
	String<Dna5>  patFw;               // forward-strand sequence
	uint8_t       patBufFw[BUF_SIZE];  // forward-strand sequence buffer
	String<Dna5>  patRc;               // reverse-complement sequence
	uint8_t       patBufRc[BUF_SIZE];  // reverse-complement sequence buffer
	String<char>  qualFw;              // quality values
	char          qualBufFw[BUF_SIZE]; // quality value buffer
	String<char>  qualRc;              // reverse quality values
	char          qualBufRc[BUF_SIZE]; // reverse quality value buffer
	String<char>  name;                // read name
	char          nameBuf[BUF_SIZE];   // read name buffer
};

/**
 * Encapsualtes a source of patterns; usually a file.  Handles dumping
 * patterns to a logfile (useful for debugging) and optionally
 * reversing them before returning them.
 */
class PatternSource {
public:
	PatternSource(bool __reverse = false,
	              const char *__dumpfile = NULL) :
	    _readCnt(0),
	    _reverse(__reverse),
		_dumpfile(__dumpfile),
		_lock()
	{
		// Open dumpfile, if specified
		if(_dumpfile != NULL) {
			_out.open(_dumpfile, ios_base::out);
			if(!_out.good()) {
				cerr << "Could not open pattern dump file \"" << _dumpfile << "\" for writing" << endl;
				exit(1);
			}
		}
#ifdef USE_SPINLOCK
		// No initialization
#else
		MUTEX_INIT(_lock);
#endif
	}
	virtual ~PatternSource() { }
	/**
	 * The main member function for dispensing patterns.
	 */
	virtual void nextRead(ReadBuf& r, uint32_t& patid) {
		// nextPatternImpl does the reading from the ultimate source;
		// it is implemented in concrete subclasses
		nextReadImpl(r, patid);
		// If it's this class's responsibility to reverse the pattern,
		// do so here.  Usually it's the responsibility of one of the
		// concrete subclasses, since they can usually do it more
		// efficiently.
		if(_reverse) {
			::reverse(r.patFw);
			::reverse(r.patRc);
			::reverse(r.qualFw);
			::reverse(r.qualRc);
		}
		// Output it, if desired
		if(_dumpfile != NULL) {
			dump(_out, r.patFw,
			     empty(r.qualFw) ? String<char>("(empty)") : r.qualFw,
			     empty(r.name)   ? String<char>("(empty)") : r.name);
			dump(_out, r.patRc,
			     empty(r.qualRc) ? String<char>("(empty)") : r.qualRc,
			     empty(r.name)   ? String<char>("(empty)") : r.name);
		}
	}
	/// Implementation to be provided by concrete subclasses
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) = 0;
	/// Reset state to start over again with the first read
	virtual void reset() { _readCnt = 0; }
	/**
	 * Whether to reverse reads as they're read in (useful when using
	 * the mirror index)
	 */
	virtual bool reverse() const { return _reverse; }
	/**
	 * Set whether to reverse reads as they're read in (useful when
	 * using the mirror index.
	 */
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
	uint32_t patid() { return _readCnt; }
protected:
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << name << ": " << seq << " " << qual << endl;
	}
	uint32_t _readCnt;
	/**
	 * Concrete subclasses call lock() to enter a critical region.
	 * What constitutes a critical region depends on the subclass.
	 */
	void mylock() {
#ifdef USE_SPINLOCK
		_lock.Enter();
#else
		MUTEX_LOCK(_lock);
#endif
	}
	/**
	 * Concrete subclasses call unlock() to exit a critical region
	 * What constitutes a critical region depends on the subclass.
	 */
	void myunlock() {
#ifdef USE_SPINLOCK
		_lock.Leave();
#else
		MUTEX_UNLOCK(_lock);
#endif
	}
private:
	bool _reverse;         /// reverse patterns before returning them
	const char *_dumpfile; /// dump patterns to this file before returning them
	ofstream _out;         /// output stream for dumpfile
#ifdef USE_SPINLOCK
	SpinLock _lock;
#else
	MUTEX_T _lock; /// mutex for locking critical regions
#endif
};

/**
 * Encapsulates a single thread's interaction with the PatternSource.
 * Most notably, this class holds the buffers into which the
 * PatterSource will write sequences.  This class is *not* threadsafe
 * - it doesn't need to be since there's one per thread.  PatternSource
 * is thread-safe.
 */
class PatternSourcePerThread {
public:
	PatternSourcePerThread(PatternSource& __patsrc) :
		_patsrc(__patsrc), _buf(), _patid(0xffffffff) { }

	void nextRead() {
		ASSERT_ONLY(uint32_t lastPatid = _patid);
		_buf.clearAll();
		_patsrc.nextRead(_buf, _patid);
		assert(empty() || _patid != lastPatid);
	}

	String<Dna5>& patFw()  { return _buf.patFw;               }
	String<Dna5>& patRc()  { return _buf.patRc;               }
	String<char>& qualFw() { return _buf.qualFw;              }
	String<char>& qualRc() { return _buf.qualRc;              }
	String<char>& name()   { return _buf.name;                }
	uint32_t      patid()  { return _patid;                   }
	bool          empty()  { return seqan::empty(_buf.patFw); }
	void          reset()  { _patid = 0xffffffff;             }
private:
	PatternSource& _patsrc;
	ReadBuf  _buf;   // read buffer
	uint32_t _patid; // index of read just read
};

/**
 * Encapsualtes a source of patterns where each raw pattern is trimmed
 * by some user-defined amount on the 3' and 5' ends.  Doesn't
 * implement the actual trimming - that's up to the concrete
 * descendants.
 */
class TrimmingPatternSource : public PatternSource {
public:
	TrimmingPatternSource(bool __reverse = false,
	                      const char *__dumpfile = NULL,
	                      int __trim3 = 0,
	                      int __trim5 = 0) :
		PatternSource(__reverse, __dumpfile),
		_trim3(__trim3),
		_trim5(__trim5) { }
protected:
	int _trim3;
	int _trim5;
};

class RandomPatternSource : public PatternSource {
public:
	RandomPatternSource(uint32_t numReads = 2000000,
	                    int length = 35,
	                    const char *__dumpfile = NULL,
	                    uint32_t seed = 0) :
		PatternSource(false, __dumpfile),
		_numReads(numReads),
		_length(length),
		_seed(seed),
		_rand(seed),
		_reverse(false) { }

	/// Implementation to be provided by concrete subclasses
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		mylock();
		if(_readCnt >= _numReads) {
			r.clearAll();
			myunlock();
			return;
		}
		uint32_t ra = _rand.nextU32() & 3;
		patid = _readCnt;
		_readCnt++;
		myunlock();
		if(!_reverse) {
			for(int i = 0; i < _length; i++) {
				ra = RandomSource::nextU32(ra) & 3;
				r.patBufFw[i]            = ra;
				r.patBufRc[_length-i-1]  = ra ^ 3;
				char c                   = 'I' - ((ra >> 2) & 31);
				r.qualBufFw[i]           = c;
				r.qualBufRc[_length-i-1] = c;
			}
		} else {
			for(int i = 0; i < _length; i++) {
				ra = RandomSource::nextU32(ra) & 3;
				r.patBufFw[_length-i-1]  = ra;
				r.patBufRc[i]            = ra ^ 3;
				char c                   = 'I' - ((ra >> 2) & 31);
				r.qualBufFw[_length-i-1] = c;
				r.qualBufRc[i]           = c;
			}
		}
		_setBegin (r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, _length);
		_setBegin (r.patRc, (Dna5*)r.patBufRc);
		_setLength(r.patRc, _length);
		_setBegin (r.qualFw, r.qualBufFw);
		_setLength(r.qualFw, _length);
		_setBegin (r.qualRc, r.qualBufRc);
		_setLength(r.qualRc, _length);

		itoa10(patid, r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, strlen(r.nameBuf));
	}
	virtual void reset() {
		PatternSource::reset();
		_rand.init(_seed);
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
private:
	uint32_t _numReads;
	int _length;
	uint32_t _seed;
	RandomSource _rand;
	bool _reverse;
};

class FileBuf {
public:
	FileBuf(FILE *__in) : _in(__in), _cur(BUF_SZ), _buf_sz(BUF_SZ), _done(false) {
		assert(_in != NULL);
	}

	int get() {
		assert(_in != NULL);
		int c = peek();
		if(c != -1) _cur++;
		return c;
	}

	int peek() {
		assert(_in != NULL);
		assert_leq(_cur, _buf_sz);
		if(_cur == _buf_sz) {
			if(_done) { return -1; }
			else {
				// Get the next chunk
				_buf_sz = fread(_buf, 1, BUF_SZ, _in);
				_cur = 0;
				if(_buf_sz == 0) {
					_done = true;
					return -1;
				} else if(_buf_sz < BUF_SZ) {
					_done = true;
				}
			}
		}
		return (int)_buf[_cur];
	}

	size_t gets(char *buf, size_t len) {
		size_t stored = 0;
		while(true) {
			int c = get();
			if(c == -1) {
				buf[stored] = '\0';
				return stored;
			}
			if(stored == len-1 || c == '\n' || c == '\r') {
				buf[stored] = '\0';
				return stored;
			}
			buf[stored++] = (char)c;
		}
	}
private:
	static const size_t BUF_SZ = 256 * 1024;
	FILE   *_in;
	size_t  _cur;
	size_t  _buf_sz;
	bool    _done;
	char    _buf[BUF_SZ]; // (large) input buffer
};

/// Skip to the end of the current string of newline chars and return
/// the first character after the newline chars, or -1 for EOF
static inline int getOverNewline(FileBuf *in) {
	int c;
	while(isspace(c = in->get()));
	return c;
}

/// Skip to the end of the current string of newline chars such that
/// the next call to get() returns the first character after the
/// whitespace
static inline void peekOverNewline(FileBuf *in) {
	while(true) {
		int c = in->peek();
		if(c != '\r' && c != '\n') {
			return;
		}
		in->get();
	}
}

/// Skip to the end of the current line; return the first character
/// of the next line or -1 for EOF
static inline int getToEndOfLine(FileBuf *in) {
	while(true) {
		int c = in->get(); if(c < 0) return -1;
		if(c == '\n' || c == '\r') {
			while(c == '\n' || c == '\r') {
				c = in->get(); if(c < 0) return -1;
			}
			// c now holds first character of next line
			return c;
		}
	}
}

/// Skip to the end of the current line such that the next call to
/// get() returns the first character on the next line
static inline void peekToEndOfLine(FileBuf *in) {
	while(true) {
		int c = in->get(); if(c < 0) return;
		if(c == '\n' || c == '\r') {
			c = in->peek();
			while(c == '\n' || c == '\r') {
				in->get(); if(c < 0) return; // consume \r or \n
				c = in->peek();
			}
			// next get() gets first character of next line
			return;
		}
	}
}

/**
 * Encapsualtes a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public TrimmingPatternSource {
public:
	VectorPatternSource(const vector<string>& v,
	                    bool __reverse = false,
	                    const char *__dumpfile = NULL,
	                    int __trim3 = 0,
	                    int __trim5 = 0,
		                int __policy = NS_TO_NS,
	                    int __maxNs = 9999,
	                    uint32_t seed = 0) :
		TrimmingPatternSource(false, __dumpfile, __trim3, __trim5),
		_reverse(__reverse), _cur(0), _maxNs(__maxNs),
		_v(), _vrev(), _vrc(), _vrcrev(), _quals(), _qualsrev(), _rand(seed)
	{
		int nrejects = 0;
		for(size_t i = 0; i < v.size(); i++) {
			vector<string> ss;
			tokenize(v[i], ":", ss);
			assert_gt(ss.size(), 0);
			assert_leq(ss.size(), 2);
			// Initialize s
			string s = ss[0];
			if(s.length() <= (size_t)(_trim3 + _trim5)) {
				// Entire read is trimmed away
				continue;
			} else {
				// Trim on 5' (high-quality) end
				if(_trim5 > 0) {
					s.erase(0, _trim5);
				}
				// Trim on 3' (low-quality) end
				if(_trim3 > 0) {
					s.erase(s.length()-_trim3);
				}
			}
			// Count Ns and possibly reject
			int ns = 0;
			for(size_t j = 0; j < s.length(); j++) {
				if(s[j] == 'N' || s[j] == 'n') {
					ns++;
					if(__policy == NS_TO_NS) {
						// Leave s[j] == 'N'
					} else if(__policy == NS_TO_RANDS) {
						s[j] = "ACGT"[_rand.nextU32() & 3];
					} else {
						assert_eq(NS_TO_AS, __policy);
						s[j] = 'A';
					}
				}
			}
			// Enforce upper limit on Ns; note that this limit is in
			// effect even if the N policy is not NS_TO_NS
			if(ns > _maxNs) {
				nrejects++;
				continue;
			}
			//  Initialize vq
			string vq;
			if(ss.size() == 2) {
				vq = ss[1];
			}
			// Trim qualities
			if(vq.length() > (size_t)(_trim3 + _trim5)) {
				// Trim on 5' (high-quality) end
				if(_trim5 > 0) {
					vq.erase(0, _trim5);
				}
				// Trim on 3' (low-quality) end
				if(_trim3 > 0) {
					vq.erase(vq.length()-_trim3);
				}
			}
			// Pad quals with Is if necessary; this shouldn't happen
			while(vq.length() < length(s)) {
				vq.push_back('I');
			}
			// Truncate quals to match length of read if necessary;
			// this shouldn't happen
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
			_vrc.push_back(reverseComplement(String<Dna5>(s)));
			{
				_vrcrev.push_back(reverseComplement(String<Dna5>(s)));
				::reverse(_vrcrev.back());
			}
			ostringstream os;
			os << (_names.size());
			_names.push_back(os.str());
		}
		assert_gt(_v.size() + nrejects, 0);
		assert_eq(_v.size(), _vrev.size());
		assert_eq(_v.size(), _quals.size());
		assert_eq(_v.size(), _qualsrev.size());
	}
	virtual ~VectorPatternSource() { }
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		// Let Strings begin at the beginning of the respective bufs
		r.reset();
		mylock();
		if(_cur >= _v.size()) {
			myunlock();
			// Clear all the Strings, as a signal to the caller that
			// we're out of reads
			clear(r.patFw);
			clear(r.patRc);
			clear(r.qualFw);
			clear(r.qualRc);
			clear(r.name);
			return;
		}
		// Copy _v*, _quals* strings into the respective Strings
		if(!_reverse) {
			// not reversed
			r.patFw  = _v[_cur];
			r.patRc  = _vrc[_cur];
			r.qualFw = _quals[_cur];
			r.qualRc = _qualsrev[_cur];
		} else {
			// reversed
			r.patFw  = _vrev[_cur];
			r.patRc  = _vrcrev[_cur];
			r.qualFw = _qualsrev[_cur];
			r.qualRc = _quals[_cur];
		}
		ostringstream os;
		os << _cur;
		r.name = os.str();
		_cur++;
		_readCnt++;
		patid = _readCnt;
		myunlock();
	}
	virtual void reset() {
		TrimmingPatternSource::reset();
		_cur = 0;
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
private:
	bool   _reverse;
	size_t _cur;
	int    _maxNs;
	vector<String<Dna5> > _v;        /// forward sequences
	vector<String<Dna5> > _vrev;     /// reversed forward sequences
	vector<String<Dna5> > _vrc;      /// rev-comp sequences
	vector<String<Dna5> > _vrcrev;   /// reversed rev-comp sequences
	vector<String<char> > _quals;    /// quality values parallel to _v
	vector<String<char> > _qualsrev; /// quality values parallel to _vrev
	vector<String<char> > _names;    /// names
	RandomSource _rand;
};

/**
 * Supports reversing all strings as they're read in.  Also supports
 * returning reverse complements interspersed with forward versions of
 * patterns.
 */
class BufferedFilePatternSource : public TrimmingPatternSource {
public:
	BufferedFilePatternSource(const vector<string>& infiles,
	                          bool __reverse = false,
	                          const char *__dumpfile = NULL,
	                          int __trim3 = 0,
	                          int __trim5 = 0) :
		TrimmingPatternSource(__reverse, __dumpfile, __trim3, __trim5),
		_infiles(infiles),
		_filecur(0),
		_filebuf(NULL),
		_first(true)
	{
		assert_gt(infiles.size(), 0);
		open();
		_filecur++;
	}
	virtual ~BufferedFilePatternSource() {
		// close currently-open file
		if(_filebuf != NULL) {
			delete _filebuf;
			_filebuf = NULL;
		}
	}
	/**
	 * Fill ReadBuf with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 */
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		// We are entering a critical region, because we're
		// manipulating our file handle and _filecur state
		mylock();
		read(r, patid);
		if(_first && seqan::empty(r.patFw)) {
			// No reads could be extracted from the first _infile
			cerr << "Warning: Could not find any reads in \"" << _infiles[0] << "\"" << endl;
		}
		_first = false;
		while(seqan::empty(r.patFw) && _filecur < _infiles.size()) {
			// Close current file
			assert(_filebuf != NULL);
			delete _filebuf;
			_filebuf = NULL;
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			read(r, patid);
			if(seqan::empty(r.patFw)) {
				// No reads could be extracted from this _infile
				cerr << "Warning: Could not find any reads in \"" << _infiles[_filecur] << "\"" << endl;
			}
			_filecur++;
		}
		// Leaving critical region
		myunlock();
		// If r.patFw is empty, then the caller knows that we are
		// finished with the reads
	}
	/**
	 * Reset state so that we read start reading again from the
	 * beginning of the first file.  Should only be called by the
	 * master thread.
	 */
	virtual void reset() {
		TrimmingPatternSource::reset();
		delete _filebuf;
		_filebuf = NULL;
		_filecur = 0,
		open();
		_filecur++;
	}
protected:
	/// Read another pattern from the input file; this is overridden
	/// to deal with specific file formats
	virtual void read(ReadBuf& r, uint32_t& patid) = 0; // length of name
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() = 0;
	void open() {
		assert(_filebuf == NULL);
		while(true) {
			// Open read
			FILE *in;
			if((in = fopen(_infiles[_filecur].c_str(), "r")) == NULL) {
				cerr << "Warning: Could not open file \"" << _infiles[_filecur] << "\" for reading" << endl;
				_filecur++;
				continue;
			}
			_filebuf = new FileBuf(in);
			break;
		}
	}
	const vector<string>& _infiles; // filenames for read files
	size_t _filecur;   // index into _infiles of next file to read
	FileBuf *_filebuf; // read file currently being read from
	bool _first;
};

/// For converting from ASCII to the Dna5 code where A=0, C=1, G=2,
/// T=3, N=4
static uint8_t charToDna5[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
	       /*    A     C           G                    N */
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

/// For converting from ASCII to the reverse-complement Dna5 code where
/// A=3, C=2, G=1, T=0, N=4
static uint8_t rcCharToDna5[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 4, 0,
	       /*    A     C           G                    N */
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

/// Default quality values for use with the FASTA pattern source
const char* qualDefault = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

/**
 *
 */
class FastaPatternSource : public BufferedFilePatternSource {
public:
	FastaPatternSource(const vector<string>& infiles,
	                   bool __reverse = false,
	                   const char *__dumpfile = NULL,
	                   int __trim3 = 0,
	                   int __trim5 = 0,
	                   int __policy = NS_TO_NS,
	                   int __maxNs = 9999,
	                   uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, __dumpfile, __trim3, __trim5),
		_first(true), _reverse(__reverse), _policy(__policy), _maxNs(__maxNs), _rand(seed)
	{ }
	virtual void reset() {
		_first = true;
		BufferedFilePatternSource::reset();
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
protected:

	/// Read another pattern from a FASTA input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		const int bufSz = ReadBuf::BUF_SIZE;
		int ns;
		do {
			int c;
			int dstLen = 0;
			int nameLen = 0;
			ns = 0; // reset 'N' count
			// Pick off the first carat
			if(_first) {
				c = getOverNewline(_filebuf); if(c < 0) {
					r.clearAll(); return;
				}
				if(c != '>') {
					int cc = toupper(c);
					cerr << "Error: reads file does not look like a FASTA file" << endl;
					if(c == '@') {
						cerr << "Reads file looks like a FASTQ file; please use -f" << endl;
					} else if(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T') {
						cerr << "Reads file may be a raw file; please use -r" << endl;
					}
					exit(1);
				}
				assert(c == '>' || c == '#');
				_first = false;
			}

			// Read to the end of the id line, sticking everything after the '>'
			// into *name
			while(true) {
				c = _filebuf->get(); if(c < 0) {
					r.clearAll(); return;
				}
				if(c == '\n' || c == '\r') {
					// Break at end of line, after consuming all \r's, \n's
					while(c == '\n' || c == '\r') {
						c = _filebuf->get(); if(c < 0) {
							r.clearAll(); return;
						}
					}
					break;
				}
				r.nameBuf[nameLen++] = c;
			}
			_setBegin(r.name, r.nameBuf);
			_setLength(r.name, nameLen);

			// _in now points just past the first character of a sequence
			// line, and c holds the first character
			int begin = 0;
			if(!_reverse) {
				while(c != '>' && c != '#') {
					// Note: can't have a comment in the middle of a sequence,
					// though a comment can end a sequence
					if(isalpha(c) && begin++ >= this->_trim5) {
						if(c == 'N' || c == 'n') {
							if(_policy == NS_TO_NS) {
								// Leave c = 'N'
							} else if(_policy == NS_TO_RANDS) {
								c = "ACGT"[_rand.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, _policy);
								c = 'A';
							}
						}
						r.patBufFw[dstLen] = charToDna5[c];
						r.patBufRc[bufSz-dstLen-1] = rcCharToDna5[c];
						dstLen++;
					}
					if((c = _filebuf->get()) < 0) break;
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[i] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)r.patBufFw);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, const_cast<char*>(qualDefault));
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)&r.patBufRc[bufSz-dstLen]);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, const_cast<char*>(qualDefault));
				_setLength(r.qualRc, dstLen);
			} else {
				while(c != '>' && c != '#') {
					// Note: can't have a comment in the middle of a sequence,
					// though a comment can end a sequence
					if(isalpha(c) && begin++ >= this->_trim5) {
						if(c == 'N' || c == 'n') {
							if(_policy == NS_TO_NS) {
								// Leave c = 'N'
							} else if(_policy == NS_TO_RANDS) {
								c = "ACGT"[_rand.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, _policy);
								c = 'A';
							}
						}
						r.patBufFw[bufSz-dstLen-1] = charToDna5[c];
						r.patBufRc[dstLen] = rcCharToDna5[c];
						dstLen++;
					}
					if((c = _filebuf->get()) < 0) break;
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[bufSz-i-1] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)&r.patBufFw[bufSz-dstLen]);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, const_cast<char*>(qualDefault));
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)r.patBufRc);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, const_cast<char*>(qualDefault));
				_setLength(r.qualRc, dstLen);
			}

			// Set up a default name if one hasn't been set
			if(nameLen == 0) {
				itoa10(_readCnt, r.nameBuf);
				_setBegin(r.name, r.nameBuf);
				nameLen = strlen(r.nameBuf);
				_setLength(r.name, nameLen);
			}
			assert_gt(nameLen, 0);
			_readCnt++;

		} while(ns > _maxNs);
		patid = _readCnt-1;
	}
	virtual void resetForNextFile() {
		_first = true;
	}
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << ">" << name << endl << seq << endl;
	}
private:
	bool _first;
	bool _reverse;
	int _policy;
	int _maxNs;
	RandomSource _rand;
};

/**
 * Read a FASTQ-format file.
 * See: http://maq.sourceforge.net/fastq.shtml
 */
class FastqPatternSource : public BufferedFilePatternSource {
public:
	FastqPatternSource(const vector<string>& infiles,
	                   bool __reverse = false,
	                   const char *__dumpfile = NULL,
	                   int __trim3 = 0,
	                   int __trim5 = 0,
	                   int __policy = NS_TO_NS,
					   bool solexa_quals = false,
					   int __maxNs = 9999,
	                   uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, __dumpfile, __trim3, __trim5),
		_first(true), _reverse(__reverse), _solexa_quals(solexa_quals),
		_policy(__policy), _maxNs(__maxNs), _rand(seed)
	{
		for (int l = 0; l != 128; ++l) {
			_table[l] = (int)(10.0 * log(1.0 + pow(10.0, (l - 64) / 10.0)) / log(10.0) + .499);
			if (_table[l] >= 63) _table[l] = 63;
			if (_table[l] == 0) _table[l] = 1;
		}
	}
	virtual void reset() {
		_first = true;
		BufferedFilePatternSource::reset();
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
protected:
	/// Read another pattern from a FASTQ input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		int ns;
		const int bufSz = ReadBuf::BUF_SIZE;
		do {
			int c;
			ns = 0;
			int dstLen = 0;
			int nameLen = 0;

			// Pick off the first at
			if(_first) {
				c = getOverNewline(_filebuf); if(c < 0) {
					seqan::clear(r.patFw);
					return;
				}
				if(c != '@') {
					cerr << "Error: reads file does not look like a FASTQ file" << endl;
					if(c == '>') {
						cerr << "Reads file looks like a FASTA file; please use -f" << endl;
					}
					exit(1);
				}
				assert_eq('@', c);
				_first = false;
			}

			// Read to the end of the id line, sticking everything after the '@'
			// into *name
			while(true) {
				c = _filebuf->get(); if(c < 0) {
					seqan::clear(r.patFw);
					return;
				}
				if(c == '\n' || c == '\r') {
					// Break at end of line, after consuming all \r's, \n's
					while(c == '\n' || c == '\r') {
						c = _filebuf->get();
						if(c < 0) {
							seqan::clear(r.patFw);
							return;
						}
					}
					break;
				}
				r.nameBuf[nameLen++] = c;
			}
			_setBegin(r.name, r.nameBuf);
			_setLength(r.name, nameLen);
			// c now holds the first character on the line after the
			// @name line

			// _filebuf now points just past the first character of a
			// sequence line, and c holds the first character
			if(!_reverse) {
				while(c != '+') {
					if(isalpha(c) && dstLen >= this->_trim5) {
						if(c == 'N' || c == 'n') {
							if(_policy == NS_TO_NS) {
								// Leave c = 'N'
							} else if(_policy == NS_TO_RANDS) {
								c = "ACGT"[_rand.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, _policy);
								c = 'A';
							}
						}
						r.patBufFw[dstLen] = charToDna5[c];
						r.patBufRc[bufSz-dstLen-1] = rcCharToDna5[c];
						dstLen++;
					}
					c = _filebuf->get();
					if(c < 0) break;
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[i] == 4) ns++;
				}
				_setBegin(r.patFw, (Dna5*)r.patBufFw);
				_setLength(r.patFw, dstLen);
				_setBegin(r.patRc, (Dna5*)&r.patBufRc[bufSz-dstLen]);
				_setLength(r.patRc, dstLen);
			} else {
				while(c != '+') {
					if(isalpha(c) && dstLen >= this->_trim5) {
						if(c == 'N' || c == 'n') {
							if(_policy == NS_TO_NS) {
								// Leave c = 'N'
							} else if(_policy == NS_TO_RANDS) {
								c = "ACGT"[_rand.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, _policy);
								c = 'A';
							}
						}
						r.patBufFw[bufSz-dstLen-1] = charToDna5[c];
						r.patBufRc[dstLen] = rcCharToDna5[c];
						dstLen++;
					}
					c = _filebuf->get();
					if(c < 0) break;
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[bufSz-i-1] == 4) ns++;
				}
				_setBegin(r.patFw, (Dna5*)&r.patBufFw[bufSz-dstLen]);
				_setLength(r.patFw, dstLen);
				_setBegin(r.patRc, (Dna5*)r.patBufRc);
				_setLength(r.patRc, dstLen);
			}
			assert_eq('+', c);

			// Chew up the optional name on the '+' line
			peekToEndOfLine(_filebuf);

			// Now read the qualities
			int qualsRead = 0;
			if (_solexa_quals) {
				char buf[1024];
				while (_filebuf->gets(buf, sizeof(buf)) && qualsRead < dstLen + this->_trim5)
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
							if (qualsRead >= _trim5)
							{
								size_t off = qualsRead - _trim5;
								c = (char)(pQ + 33);
								assert_geq(c, 33);
								assert_leq(c, 73);
								r.qualBufFw[off] = c;
								r.qualBufRc[bufSz - off - 1] = c;
							}
							++qualsRead;
						}
					} else {
						for (unsigned int j = 0; j < s_quals.size(); ++j)
						{
							int sQ = atoi(s_quals[j].c_str());
							int pQ = (int)(10.0 * log(1.0 + pow(10.0, sQ / 10.0)) / log(10.0) + .499);
							if (qualsRead >= _trim5)
							{
								size_t off = qualsRead - _trim5;
								c = (char)(pQ + 33);
								assert_geq(c, 33);
								assert_leq(c, 73);
								r.qualBufFw[bufSz - off - 1] = c;
								r.qualBufRc[off] = c;
							}
							++qualsRead;
						}
					}
				} // done reading Solexa quality lines
				if(!_reverse) {
					_setBegin(r.qualFw, (char*)r.qualBufFw);
					_setLength(r.qualFw, dstLen);
					_setBegin(r.qualRc, (char*)&r.qualBufRc[bufSz-dstLen]);
					_setLength(r.qualRc, dstLen);
				} else {
					_setBegin(r.qualFw, (char*)&r.qualBufFw[bufSz-dstLen]);
					_setLength(r.qualFw, dstLen);
					_setBegin(r.qualRc, (char*)r.qualBufRc);
					_setLength(r.qualRc, dstLen);
				}
				c = getOverNewline(_filebuf);
			}
			else
			{
				// Non-solexa qualities
				if(!_reverse) {
					while((qualsRead < dstLen + this->_trim5) && c >= 0) {
						c = _filebuf->get();
						if (c != '\r' && c != '\n') {
							if (qualsRead >= _trim5) {
								size_t off = qualsRead - _trim5;
								assert_geq(c, 33);
								assert_leq(c, 73);
								r.qualBufFw[off] = c;
								r.qualBufRc[bufSz - off - 1] = c;
							}
							qualsRead++;
						} else {
							break;
						}
					}
					assert_eq(qualsRead, dstLen + this->_trim5);
					_setBegin (r.qualFw, (char*)r.qualBufFw);
					_setLength(r.qualFw, dstLen);
					_setBegin (r.qualRc, (char*)&r.qualBufRc[bufSz-dstLen]);
					_setLength(r.qualRc, dstLen);
				} else {
					while((qualsRead < dstLen + this->_trim5) && c >= 0) {
						c = _filebuf->get();
						if (c != '\r' && c != '\n') {
							if (qualsRead >= _trim5) {
								size_t off = qualsRead - _trim5;
								assert_geq(c, 33);
								assert_leq(c, 73);
								r.qualBufFw[bufSz - off - 1] = c;
								r.qualBufRc[off] = c;
							}
							qualsRead++;
						} else {
							break;
						}
					}
					assert_eq(qualsRead, dstLen + this->_trim5);
					_setBegin (r.qualFw, (char*)&r.qualBufFw[bufSz-dstLen]);
					_setLength(r.qualFw, dstLen);
					_setBegin (r.qualRc, (char*)r.qualBufRc);
					_setLength(r.qualRc, dstLen);
				}
				if(c == '\r' || c == '\n') {
					c = getOverNewline(_filebuf);
				} else {
					c = getToEndOfLine(_filebuf);
				}
			}
			assert(c == -1 || c == '@');

			// Set up a default name if one hasn't been set
			if(nameLen == 0) {
				itoa10(_readCnt, r.nameBuf);
				_setBegin(r.name, r.nameBuf);
				nameLen = strlen(r.nameBuf);
				_setLength(r.name, nameLen);
			}
			assert_gt(nameLen, 0);
			_readCnt++;
		} while(ns > _maxNs);
		patid = _readCnt-1;
	}
	virtual void resetForNextFile() {
		_first = true;
	}
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
	}
private:
	bool _first;
	bool _reverse;
	bool _solexa_quals;
	int _policy;
	int _maxNs;
	int _table[128];
	RandomSource _rand;
};

/**
 * Read a Raw-format file (one sequence per line).
 */
class RawPatternSource : public BufferedFilePatternSource {
public:
	RawPatternSource(const vector<string>& infiles,
	                 bool __reverse = false,
	                 const char *__dumpfile = NULL,
	                 int __trim3 = 0,
	                 int __trim5 = 0,
	                 int __policy = NS_TO_NS,
	                 int __maxNs = 9999,
	                 uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, __dumpfile, __trim3, __trim5),
		_first(true), _reverse(__reverse), _policy(__policy), _maxNs(__maxNs), _rand(seed)
	{ }
	virtual void reset() {
		_first = true;
		BufferedFilePatternSource::reset();
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
protected:
	/// Read another pattern from a Raw input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		const int bufSz = ReadBuf::BUF_SIZE;
		int ns;
		do {
			int c;
			int dstLen = 0;
			int nameLen = 0;
			ns = 0; // reset 'N' count
			c = getOverNewline(this->_filebuf);
			assert(!isspace(c));
			if(_first) {
				int cc = toupper(c);
				if(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T') {
					cerr << "Error: reads file does not look like a Raw file" << endl;
					if(c == '>') {
						cerr << "Reads file looks like a FASTA file; please use -f" << endl;
					}
					if(c == '@') {
						cerr << "Reads file looks like a FASTQ file; please use -q" << endl;
					}
					exit(1);
				}
				_first = false;
			}

			// _in now points just past the first character of a sequence
			// line, and c holds the first character
			if(!_reverse) {
				while(!isspace(c) && c >= 0) {
					if(isalpha(c) && dstLen >= this->_trim5) {
						if(c == 'N' || c == 'n') {
							if(_policy == NS_TO_NS) {
								// Leave c = 'N'
							} else if(_policy == NS_TO_RANDS) {
								c = "ACGT"[_rand.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, _policy);
								c = 'A';
							}
						}
						r.patBufFw[dstLen] = charToDna5[c];
						r.patBufRc[bufSz-dstLen-1] = rcCharToDna5[c];
						dstLen++;
					}
					c = _filebuf->get();
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[i] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)r.patBufFw);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, const_cast<char*>(qualDefault));
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)&r.patBufRc[bufSz-dstLen]);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, const_cast<char*>(qualDefault));
				_setLength(r.qualRc, dstLen);
			} else {
				while(!isspace(c) && c >= 0) {
					if(isalpha(c) && dstLen >= this->_trim5) {
						if(c == 'N' || c == 'n') {
							if(_policy == NS_TO_NS) {
								// Leave c = 'N'
							} else if(_policy == NS_TO_RANDS) {
								c = "ACGT"[_rand.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, _policy);
								c = 'A';
							}
						}
						r.patBufFw[bufSz-dstLen-1] = charToDna5[c];
						r.patBufRc[dstLen] = rcCharToDna5[c];
						dstLen++;
					}
					c = _filebuf->get();
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[bufSz-i-1] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)&r.patBufFw[bufSz-dstLen]);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, const_cast<char*>(qualDefault));
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)r.patBufRc);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, const_cast<char*>(qualDefault));
				_setLength(r.qualRc, dstLen);
			}

			// Set up name
			itoa10(_readCnt, r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			nameLen = strlen(r.nameBuf);
			_setLength(r.name, nameLen);
			_readCnt++;

			if(c == -1) {
				if(ns > _maxNs) {
					// This read violates the Ns constraint and we're
					// about to return it; make sure that caller
					// doesn't see this read
					r.clearAll();
				}
				return;
			} else {
				assert(isspace(c));
			}
		} while(ns > _maxNs);
		patid = _readCnt-1;
	}
	virtual void resetForNextFile() {
		_first = true;
	}
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << seq << endl;
	}
private:
	bool _first;
	bool _reverse;
	int _policy;
	int _maxNs;
	RandomSource _rand;
};

#endif /*PAT_H_*/
