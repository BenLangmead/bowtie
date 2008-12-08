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
#include "filebuf.h"

/**
 * Classes and routines for reading reads from various input sources.
 */

using namespace std;
using namespace seqan;

/// Constructs string base-10 representation of integer 'value'
extern char* itoa10(int value, char* result);

/// Wildcard policies
enum {
	NS_TO_NS    = 1, // Ns stay Ns and don't match anything
	NS_TO_AS    = 2, // Ns become As
	NS_TO_RANDS = 3  // Ns become random characters
};

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
	void reverseAll() {
		::reverseInPlace(patFw);
		::reverseInPlace(patRc);
		::reverseInPlace(qualFw);
		::reverseInPlace(qualRc);
	}
	bool empty() {
		return seqan::empty(patFw);
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
 * Struct encapsulating a pair of mates.
 */
struct ReadBufPair {
	ReadBuf a; // first mate (from -1 arg)
	ReadBuf b; // second mate (from -2 arg)
};

/**
 * Encapsualtes a synchronized source of patterns; usually a file.
 * Handles dumping patterns to a logfile (useful for debugging).  Also
 * optionally reverses reads and quality strings before returning them,
 * though that is usually more efficiently done by the concrete
 * subclass.  Concrete subclasses should delimit critical sections with
 * calls to lock() and unlock().
 */
class PatternSource {
public:
	PatternSource(bool __reverse = false,
	              bool __useSpinlock = true,
	              const char *__dumpfile = NULL) :
	    _readCnt(0),
	    _reverse(__reverse),
		_dumpfile(__dumpfile),
		_numWrappers(0),
		_doLocking(true),
		_useSpinlock(__useSpinlock),
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
		MUTEX_INIT(_lock);
	}

	virtual ~PatternSource() { }

	/**
	 * Call this whenever this PatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks will be contended.
	 */
	void addWrapper() {
		_numWrappers++;
	}

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
			::reverseInPlace(r.patFw);
			::reverseInPlace(r.patRc);
			::reverseInPlace(r.qualFw);
			::reverseInPlace(r.qualRc);
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

	/**
	 * Default format for dumping a read to an output stream.  Concrete
	 * subclasses might want to do something fancier.
	 */
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << name << ": " << seq << " " << qual << endl;
	}

	/// The number of reads read by this PatternSource
	uint32_t _readCnt;

	/**
	 * Concrete subclasses call lock() to enter a critical region.
	 * What constitutes a critical region depends on the subclass.
	 */
	void lock() {
		if(!_doLocking || _numWrappers < 2) return; // no contention
#ifdef USE_SPINLOCK
		if(_useSpinlock) {
			// User can ask to use the normal pthreads lock even if
			// spinlocks are compiled in.
			_spinlock.Enter();
		} else {
#endif
			MUTEX_LOCK(_lock);
#ifdef USE_SPINLOCK
		}
#endif
	}

	/**
	 * Concrete subclasses call unlock() to exit a critical region
	 * What constitutes a critical region depends on the subclass.
	 */
	void unlock() {
		if(!_doLocking || _numWrappers < 2) return; // no contention
#ifdef USE_SPINLOCK
		if(_useSpinlock) {
			// User can ask to use the normal pthreads lock even if
			// spinlocks are compiled in.
			_spinlock.Leave();
		} else {
#endif
			MUTEX_UNLOCK(_lock);
#ifdef USE_SPINLOCK
		}
#endif
	}

	bool _reverse;         /// reverse patterns before returning them
	const char *_dumpfile; /// dump patterns to this file before returning them
	ofstream _out;         /// output stream for dumpfile
	int _numWrappers;      /// # threads that own a wrapper for this PatternSource
	bool _doLocking;       ///
	/// User can ask to use the normal pthreads-style lock even if
	/// spinlocks is enabled and compiled in.  This is sometimes better
	/// if we expect bad I/O latency on some reads.
	bool _useSpinlock;
#ifdef USE_SPINLOCK
	SpinLock _spinlock;
#endif
	MUTEX_T _lock; /// mutex for locking critical regions
};

/**
 * Encapsulates a synchronized source of both paired-end patterns and
 * unpaired patterns.
 */
class PairedPatternSource {

public:

	PairedPatternSource(const vector<PatternSource*>& srca,
	                    const vector<PatternSource*>& srcb) :
		cur_(0), basePatid_(0), srca_(srca), srcb_(srcb)
	{
	    // srca_ and srcb_ must be parallel
	    assert_eq(srca_.size(), srcb_.size());
	    for(size_t i = 0; i < srca_.size(); i++) {
	    	// Can't have NULL first-mate sources.  Second-mate sources
	    	// can be NULL, in the case when the corresponding first-
	    	// mate source is unpaired.
	    	assert(srca_[i] != NULL);
		    for(size_t j = 0; j < srcb_.size(); j++) {
		    	assert_neq(srca_[i], srcb_[j]);
		    }
	    }
	}

	/**
	 * Call this whenever this PairedPatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks within PatternSources will be contended.
	 */
	void addWrapper() {
		for(size_t i = 0; i < srca_.size(); i++) {
			srca_[i]->addWrapper();
			if(srcb_[i] != NULL) {
				srcb_[i]->addWrapper();
			}
		}
	}

	/**
	 * Reset this object and all the PatternSources under it so that
	 * the next call to nextReadPair gets the very first read pair.
	 */
	void reset() {
		for(size_t i = 0; i < srca_.size(); i++) {
			srca_[i]->reset();
			if(srcb_[i] != NULL) {
				srcb_[i]->reset();
			}
		}
		cur_ = 0;
		basePatid_ = 0;
	}

	/**
	 * Set whether to reverse reads as they're read in (useful when
	 * using the mirror index).
	 */
	void setReverse(bool r) {
		for(size_t i = 0; i < srca_.size(); i++) {
			srca_[i]->setReverse(r);
			if(srcb_[i] != NULL) {
				srcb_[i]->setReverse(r);
			}
		}
	}

	/**
	 * Return true iff the contained PatternSources are currently
	 * reversing their output.
	 */
	bool reverse() {
		return srca_[0]->reverse();
	}

	/**
	 * Return the patid of the most recently read pair.
	 */
	uint32_t patid() {
		uint32_t ret = basePatid_;
		if(cur_ < srca_.size()) {
			ret += srca_[cur_]->patid();
		}
		return ret;
	}

	/**
	 * The main member function for dispensing pairs of reads or
	 * singleton reads.  Returns true iff ra and rb contain a new
	 * pair; returns false if ra contains a new unpaired read.
	 */
	bool nextReadPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		while(cur_ < srca_.size()) {
			if(srcb_[cur_] == NULL) {
				// Patterns from srca_[cur_] are unpaired
				srca_[cur_]->nextRead(ra, patid);
				if(seqan::empty(ra.patFw)) {
					// If patFw is empty, that's our signal that the
					// input dried up
					basePatid_ += srca_[cur_]->patid();
					cur_++;
					continue; // on to next pair of PatternSources
				}
				return false; // unpaired
			} else {
				// Patterns from srca_[cur_] and srcb_[cur_] are paired
				uint32_t patid_a = 0;
				uint32_t patid_b = 0;
				srca_[cur_]->nextRead(ra, patid_a);
				srcb_[cur_]->nextRead(rb, patid_b);
				bool cont = false;
				while(patid_a != patid_b) {
					// Is either input exhausted?  If so, bail.
					if(seqan::empty(ra.patFw) || seqan::empty(rb.patFw)) {
						seqan::clear(ra.patFw);
						basePatid_ += srca_[cur_]->patid();
						cur_++; // on to next pair of PatternSources
						cont = true;
						break;
					}
					if(patid_a < patid_b) {
						srca_[cur_]->nextRead(ra, patid_a);
					} else {
						srcb_[cur_]->nextRead(rb, patid_b);
					}
				}
				if(cont) continue; // on to next pair of PatternSources
				assert_eq(patid_a, patid_b);
				patid = patid_a;
				return true; // paired
			}
		}
		return false;
	}

protected:
	uint32_t cur_; // current element in parallel srca_, srcb_ vectors
	uint32_t basePatid_;
	vector<PatternSource*> srca_; /// PatternSources for 1st mates and/or unpaired reads
	vector<PatternSource*> srcb_; /// PatternSources for 2nd mates
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
	PatternSourcePerThread() :
		_bufa(), _bufb(), _patid(0xffffffff) { }

	virtual ~PatternSourcePerThread() { }

	virtual void nextReadPair()  = 0;

	ReadBuf& bufa()        { return _bufa;         }
	ReadBuf& bufb()        { return _bufb;         }

	uint32_t      patid()  { return _patid;        }
	virtual void  reset()  { _patid = 0xffffffff;  }
	void    reverseRead()  { _bufa.reverseAll();
	                         _bufb.reverseAll();   }
	bool          empty()  { return _bufa.empty(); }

protected:
	ReadBuf  _bufa;  // read buffer
	ReadBuf  _bufb;  // read buffer
	uint32_t _patid; // index of read just read
};

/**
 * A per-thread wrapper for a PairedPatternSource.
 */
class WrappedPatternSourcePerThread : public PatternSourcePerThread {
public:
	WrappedPatternSourcePerThread(PairedPatternSource& __patsrc) :
		_patsrc(__patsrc)
	{
		_patsrc.addWrapper();
	}

	/**
	 * Get the next paired or unpaired read from the wrapped
	 * PairedPatternSource.
	 */
	virtual void nextReadPair() {
		ASSERT_ONLY(uint32_t lastPatid = _patid);
		_bufa.clearAll();
		_bufb.clearAll();
		_patsrc.nextReadPair(_bufa, _bufb, _patid);
		assert(_bufa.empty() || _patid != lastPatid);
	}
private:
	/// Container for obtaining paired reads from PatternSources
	PairedPatternSource& _patsrc;
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
	                      bool __useSpinlock = true,
	                      const char *__dumpfile = NULL,
	                      int __trim3 = 0,
	                      int __trim5 = 0) :
		PatternSource(__reverse, __useSpinlock, __dumpfile),
		_trim3(__trim3),
		_trim5(__trim5) { }
protected:
	int _trim3;
	int _trim5;
};

/**
 * A synchronized pattern source that simply returns random reads
 * without reading from the disk or storing lists of reads in memory.
 * Reads are generated with a RandomSource.
 */
class RandomPatternSource : public PatternSource {
public:
	RandomPatternSource(uint32_t numReads = 2000000,
	                    int length = 35,
	                    bool __useSpinlock = true,
	                    const char *__dumpfile = NULL,
	                    uint32_t seed = 0) :
		PatternSource(false, __useSpinlock, __dumpfile),
		_numReads(numReads),
		_length(length),
		_seed(seed),
		_rand(seed),
		_reverse(false)
	{
		if(_length > 1024) {
			cerr << "Read length for RandomPatternSource may not exceed 1024; got " << _length << endl;
			exit(1);
		}
	}

	/** Get the next random read and set patid */
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		// Begin critical section
		lock();
		if(_readCnt >= _numReads) {
			r.clearAll();
			unlock();
			return;
		}
		uint32_t ra = _rand.nextU32();
		patid = _readCnt;
		_readCnt++;
		unlock();
		fillRandomRead(r, ra, _length, patid, _reverse);
	}
	/** */
	static void fillRandomRead(ReadBuf& r,
	                           uint32_t ra,
	                           int length,
	                           uint32_t patid,
	                           bool reverse)
	{
		// End critical section
		if(!reverse) {
			for(int i = 0; i < length; i++) {
				ra = RandomSource::nextU32(ra) >> 8;
				r.patBufFw[i]           = (ra & 3);
				r.patBufRc[length-i-1]  = (ra & 3) ^ 3;
				char c                  = 'I' - ((ra >> 2) & 31);
				r.qualBufFw[i]          = c;
				r.qualBufRc[length-i-1] = c;
			}
		} else {
			for(int i = 0; i < length; i++) {
				ra = RandomSource::nextU32(ra) >> 8;
				r.patBufFw[length-i-1]  = (ra & 3);
				r.patBufRc[i]           = (ra & 3) ^ 3;
				char c                  = 'I' - ((ra >> 2) & 31);
				r.qualBufFw[length-i-1] = c;
				r.qualBufRc[i]          = c;
			}
		}
		_setBegin (r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, length);
		_setBegin (r.patRc, (Dna5*)r.patBufRc);
		_setLength(r.patRc, length);
		_setBegin (r.qualFw, r.qualBufFw);
		_setLength(r.qualFw, length);
		_setBegin (r.qualRc, r.qualBufRc);
		_setLength(r.qualRc, length);
		itoa10(patid, r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, strlen(r.nameBuf));
	}

	/** Reset the pattern source to the beginning */
	virtual void reset() {
		PatternSource::reset();
		// reset pseudo-random generator; next string of calls to
		// nextU32() will return same pseudo-randoms as the last
		_rand.init(_seed);
	}
	virtual bool reverse() const { return _reverse; }
	virtual void setReverse(bool __reverse) {
		_reverse = __reverse;
	}
private:
	uint32_t     _numReads; /// number of reads to dish out
	int          _length;   /// length of reads
	uint32_t     _seed;     /// seed for pseudo-randoms
	RandomSource _rand;     /// pseudo-random generator
	bool         _reverse;  /// whether to reverse reads
};

/**
 * A version of PatternSourcePerThread that dishes out random patterns
 * without any synchronization.
 */
class RandomPatternSourcePerThread : public PatternSourcePerThread {
public:
	RandomPatternSourcePerThread(uint32_t __numreads,
	                             int __length,
	                             int __numthreads,
	                             int thread,
	                             bool __reverse) :
		PatternSourcePerThread(),
		_numreads(__numreads),
		_length(__length),
		_numthreads(__numthreads),
		_thread(thread),
		_reverse(__reverse),
		_rand(_thread)
	{
		_patid = _thread;
		if(_length > 1024) {
			cerr << "Read length for RandomPatternSourcePerThread may not exceed 1024; got " << _length << endl;
			exit(1);
		}
	}

	virtual void nextReadPair() {
		if(_patid >= _numreads) {
			_bufa.clearAll();
			_bufb.clearAll();
			return;
		}
		RandomPatternSource::fillRandomRead(
			_bufa, _rand.nextU32(), _length, _patid, _reverse);
		RandomPatternSource::fillRandomRead(
			_bufb, _rand.nextU32(), _length, _patid, _reverse);
		_patid += _numthreads;
	}

	virtual void reset() {
		PatternSourcePerThread::reset();
		_patid = _thread;
		_rand.init(_thread);
	}

private:
	uint32_t     _numreads;
	int          _length;
	int          _numthreads;
	int          _thread;
	bool         _reverse;
	RandomSource _rand;
};

/**
 * Scan to the next FASTQ record (starting with @) and return the first
 * character of the record (which will always be @).  Since the quality
 * line may start with @, we keep scanning until we've seen a line
 * beginning with @ where the line two lines back began with +.
 */
static inline int getToFirstRecord(FileBuf& in, bool sawPlus) {
	int line = 0;
	int plusLine = -1;
	int c = in.get();
	int firstc = c;
	while(true) {
		if(line > 20) {
			// If we couldn't find our desired '@' in the first 20
			// lines, it's time to give up
			if(firstc == '>') {
				// That firstc is '>' may be a hint that this is
				// actually a FASTA file, so return it intact
				return '>';
			}
			// Return an error
			return -1;
		}
		if(c == -1) return -1;
		if(c == '\n') {
			c = in.get();
			if(c == '@' && sawPlus && plusLine == (line-2)) {
				return '@';
			}
			else if(c == '+') {
				// Saw a '+' at the beginning of a line; remember where
				// we saw it
				sawPlus = true;
				plusLine = line;
			}
			else if(c == -1) {
				return -1;
			}
			line++;
		}
		c = in.get();
	}
}

/// Skip to the end of the current string of newline chars and return
/// the first character after the newline chars, or -1 for EOF
static inline int getOverNewline(FileBuf& in) {
	int c;
	while(isspace(c = in.get()));
	return c;
}

/// Skip to the end of the current string of newline chars such that
/// the next call to get() returns the first character after the
/// whitespace
static inline void peekOverNewline(FileBuf& in) {
	while(true) {
		int c = in.peek();
		if(c != '\r' && c != '\n') {
			return;
		}
		in.get();
	}
}

/// Skip to the end of the current line; return the first character
/// of the next line or -1 for EOF
static inline int getToEndOfLine(FileBuf& in) {
	while(true) {
		int c = in.get(); if(c < 0) return -1;
		if(c == '\n' || c == '\r') {
			while(c == '\n' || c == '\r') {
				c = in.get(); if(c < 0) return -1;
			}
			// c now holds first character of next line
			return c;
		}
	}
}

/// Skip to the end of the current line such that the next call to
/// get() returns the first character on the next line
static inline void peekToEndOfLine(FileBuf& in) {
	while(true) {
		int c = in.get(); if(c < 0) return;
		if(c == '\n' || c == '\r') {
			c = in.peek();
			while(c == '\n' || c == '\r') {
				in.get(); if(c < 0) return; // consume \r or \n
				c = in.peek();
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
	                    bool __useSpinlock = true,
	                    const char *__dumpfile = NULL,
	                    int __trim3 = 0,
	                    int __trim5 = 0,
		                int __policy = NS_TO_NS,
	                    int __maxNs = 9999,
	                    uint32_t seed = 0) :
		TrimmingPatternSource(false, __useSpinlock, __dumpfile, __trim3, __trim5),
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
				::reverseInPlace(_vrev.back());
				_qualsrev.push_back(vq);
				::reverseInPlace(_qualsrev.back());
			}
			_vrc.push_back(reverseComplement(String<Dna5>(s)));
			{
				_vrcrev.push_back(reverseComplement(String<Dna5>(s)));
				::reverseInPlace(_vrcrev.back());
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
		lock();
		if(_cur >= _v.size()) {
			unlock();
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
		unlock();
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
	                          bool __useSpinlock = true,
	                          bool __forgiveInput = false,
	                          const char *__dumpfile = NULL,
	                          int __trim3 = 0,
	                          int __trim5 = 0) :
		TrimmingPatternSource(__reverse, __useSpinlock, __dumpfile, __trim3, __trim5),
		_infiles(infiles),
		_errs(),
		_filecur(0),
		_filebuf(),
		_forgiveInput(__forgiveInput),
		_first(true)
	{
		assert_gt(infiles.size(), 0);
		_errs.resize(_infiles.size(), false);
		open(); // open first file in the list
		_filecur++;
	}

	virtual ~BufferedFilePatternSource() {
		if(_filebuf.isOpen()) _filebuf.close();
	}

	/**
	 * Fill ReadBuf with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 */
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		// We are entering a critical region, because we're
		// manipulating our file handle and _filecur state
		lock();
		read(r, patid);
		if(_first && seqan::empty(r.patFw) && !_forgiveInput) {
			// No reads could be extracted from the first _infile
			cerr << "Warning: Could not find any reads in \"" << _infiles[0] << "\"" << endl;
		}
		_first = false;
		while(seqan::empty(r.patFw) && _filecur < _infiles.size()) {
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			read(r, patid);
			if(seqan::empty(r.patFw) && !_forgiveInput) {
				// No reads could be extracted from this _infile
				cerr << "Warning: Could not find any reads in \"" << _infiles[_filecur] << "\"" << endl;
			}
			_filecur++;
		}
		// Leaving critical region
		unlock();
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
		if(_filebuf.isOpen()) _filebuf.close();
		while(_filecur < _infiles.size()) {
			// Open read
			FILE *in;
			if(_infiles[_filecur] == "-") {
				in = stdin;
			} else if((in = fopen(_infiles[_filecur].c_str(), "r")) == NULL) {
				if(!_errs[_filecur]) {
					cerr << "Warning: Could not open file \"" << _infiles[_filecur] << "\" for reading" << endl;
					_errs[_filecur] = true;
				}
				_filecur++;
				continue;
			}
			_filebuf.newFile(in);
			return;
		}
		exit(1);
	}
	vector<string> _infiles; /// filenames for read files
	vector<bool> _errs; /// whether we've already printed an error for each file
	size_t _filecur;   /// index into _infiles of next file to read
	FileBuf _filebuf;  /// read file currently being read from
	bool _forgiveInput; /// try hard to parse input even if it's malformed
	bool _first;
};

/**
 *
 */
class FastaPatternSource : public BufferedFilePatternSource {
public:
	FastaPatternSource(const vector<string>& infiles,
	                   bool __reverse = false,
	                   bool __useSpinlock = true,
	                   const char *__dumpfile = NULL,
	                   int __trim3 = 0,
	                   int __trim5 = 0,
	                   int __policy = NS_TO_NS,
	                   int __maxNs = 9999,
	                   uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, false, __useSpinlock, __dumpfile, __trim3, __trim5),
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
				c = _filebuf.get(); if(c < 0) {
					r.clearAll(); return;
				}
				if(c == '\n' || c == '\r') {
					// Break at end of line, after consuming all \r's, \n's
					while(c == '\n' || c == '\r') {
						c = _filebuf.get(); if(c < 0) {
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
						if(dstLen + 1 > 1024) {
							cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
							     << "reads and re-run Bowtie";
							exit(1);
						}
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
						r.patBufFw [dstLen] = charToDna5[c];
						r.qualBufFw[dstLen] = 'I';
						r.patBufRc [bufSz-dstLen-1] = rcCharToDna5[c];
						r.qualBufRc[bufSz-dstLen-1] = 'I';
						dstLen++;
					}
					if((c = _filebuf.get()) < 0) break;
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[i] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)r.patBufFw);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, r.qualBufFw);
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)&r.patBufRc[bufSz-dstLen]);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, &r.qualBufRc[bufSz-dstLen]);
				_setLength(r.qualRc, dstLen);
			} else {
				while(c != '>' && c != '#') {
					// Note: can't have a comment in the middle of a sequence,
					// though a comment can end a sequence
					if(isalpha(c) && begin++ >= this->_trim5) {
						if(dstLen + 1 > 1024) {
							cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
							     << "reads and re-run Bowtie";
							exit(1);
						}
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
						r.patBufFw [bufSz-dstLen-1] = charToDna5[c];
						r.qualBufFw[bufSz-dstLen-1] = 'I';
						r.patBufRc [dstLen] = rcCharToDna5[c];
						r.qualBufRc[dstLen] = 'I';
						dstLen++;
					}
					if((c = _filebuf.get()) < 0) break;
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[bufSz-i-1] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)&r.patBufFw[bufSz-dstLen]);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, &r.qualBufFw[bufSz-dstLen]);
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)r.patBufRc);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, r.qualBufRc);
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

extern void wrongQualityScale();
extern void wrongQualityFormat();
extern void tooFewQualities(const String<char>& read_name);

/**
 * Read a FASTQ-format file.
 * See: http://maq.sourceforge.net/fastq.shtml
 */
class FastqPatternSource : public BufferedFilePatternSource {
public:
	FastqPatternSource(const vector<string>& infiles,
	                   bool __reverse = false,
	                   bool __useSpinlock = true,
	                   const char *__dumpfile = NULL,
	                   int __trim3 = 0,
	                   int __trim5 = 0,
	                   int __policy = NS_TO_NS,
	                   bool __forgiveInput = false,
					   bool solexa_quals = false,
					   bool integer_quals = true,
					   int __maxNs = 9999,
	                   uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, __useSpinlock, __forgiveInput, __dumpfile, __trim3, __trim5),
		_first(true), _reverse(__reverse),
		_solexa_quals(solexa_quals),
		_integer_quals(integer_quals),
		_policy(__policy), _maxNs(__maxNs),
		_rand(seed)
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
				c = _filebuf.get();
				if(c != '@') {
					if(_forgiveInput) {
						c = getToFirstRecord(_filebuf, c == '+'); if(c < 0) {
							seqan::clear(r.patFw);
							return;
						}
					} else {
						c = getOverNewline(_filebuf); if(c < 0) {
							seqan::clear(r.patFw);
							return;
						}
					}
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
				c = _filebuf.get(); if(c < 0) {
					seqan::clear(r.patFw);
					return;
				}
				if(c == '\n' || c == '\r') {
					// Break at end of line, after consuming all \r's, \n's
					while(c == '\n' || c == '\r') {
						c = _filebuf.get();
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
			int charsRead = 0;
			if(!_reverse) {
				while(c != '+') {
					if(isalpha(c)) {
						// If it's past the 5'-end trim point
						if(charsRead >= this->_trim5) {
							if(dstLen + 1 > 1024) {
								cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
								     << "reads and re-run Bowtie";
								exit(1);
							}
							// Add it to the read buffer
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
						charsRead++;
					}
					c = _filebuf.get();
					if(c < 0) {
						// EOF occurred in the middle of a read - abort
						seqan::clear(r.patFw);
						return;
					}
				}
				// Trim from 3' end
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[i] == 4) ns++;
				}
				// Set trimmed bounds of buffers
				_setBegin(r.patFw, (Dna5*)r.patBufFw);
				_setLength(r.patFw, dstLen);
				_setBegin(r.patRc, (Dna5*)&r.patBufRc[bufSz-dstLen]);
				_setLength(r.patRc, dstLen);
			} else {
				while(c != '+') {
					if(isalpha(c)) {
						// If it's past the 5'-end trim point
						if(charsRead >= this->_trim5) {
							if(dstLen + 1 > 1024) {
								cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
								     << "reads and re-run Bowtie";
								exit(1);
							}
							// Add it to the read buffer
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
						charsRead++;
					}
					c = _filebuf.get();
					if(c < 0) {
						// EOF occurred in the middle of a read - abort
						seqan::clear(r.patFw);
						return;
					}
				}
				// Trim from 3' end
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[bufSz-i-1] == 4) ns++;
				}
				// Set trimmed bounds of buffers
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
			if (_integer_quals) {
				char buf[4096];
				while (qualsRead < charsRead) {
					size_t rd = _filebuf.gets(buf, sizeof(buf));
					if(rd == 0) break;
					assert(NULL == strrchr(buf, '\n'));
					vector<string> s_quals;
					tokenize(string(buf), " ", s_quals);
					if(!_reverse) {
						for (unsigned int j = 0; j < s_quals.size(); ++j)
						{
							int iQ = atoi(s_quals[j].c_str());
							int pQ;
							if (_solexa_quals)
							{
								// Convert from solexa quality to phred
								// quality and translate to ASCII
								// http://maq.sourceforge.net/qual.shtml
								pQ = (int)(10.0 * log(1.0 + pow(10.0, (iQ) / 10.0)) / log(10.0) + .499) + 33;
							}
							else
							{
								// Keep the phred quality and translate
								// to ASCII
								pQ = (iQ <= 93 ? iQ : 93) + 33;
								if (pQ < 33)
								{
									cerr << "Saw ASCII character " << ((int)pQ) << "." << endl;
									wrongQualityScale();
									exit(1);
								}
							}

							if (qualsRead >= _trim5)
							{
								size_t off = qualsRead - _trim5;
								if(off + 1 > 1024) {
									cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
									     << "Please truncate reads and quality values and and re-run Bowtie";
									exit(1);
								}
								c = (char)(pQ);
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
							int iQ = atoi(s_quals[j].c_str());
							int pQ;
							if (_solexa_quals)
							{
								// Convert from solexa quality to phred
								// quality and translate to ASCII
								// http://maq.sourceforge.net/qual.shtml
								pQ = (int)(10.0 * log(1.0 + pow(10.0, (iQ) / 10.0)) / log(10.0) + .499) + 33;
							}
							else
							{
								// Keep the phred quality and translate
								// to ASCII
								pQ = (iQ <= 93 ? iQ : 93) + 33;
								if (pQ < 33)
								{
									cerr << "Saw ASCII character " << ((int)pQ) << "." << endl;
									wrongQualityScale();
									exit(1);
								}
							}
							if (qualsRead >= _trim5)
							{
								size_t off = qualsRead - _trim5;
								if(off + 1 > 1024) {
									cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
									     << "Please truncate reads and quality values and and re-run Bowtie";
									exit(1);
								}
								c = (char)(pQ);
								assert_geq(c, 33);
								assert_leq(c, 73);
								r.qualBufFw[bufSz - off - 1] = c;
								r.qualBufRc[off] = c;
							}
							++qualsRead;
						}
					}
				} // done reading integer quality lines
				//assert_eq(charsRead, qualsRead);
				if (charsRead > qualsRead)
				{
					tooFewQualities(r.name);
					exit(1);
				}
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
				c = _filebuf.get();
			}
			else
			{
				// Non-integer qualities
				if(!_reverse) {
					while((qualsRead < dstLen + this->_trim5) && c >= 0) {
						c = _filebuf.get();
						if (c == ' ') {
							wrongQualityFormat();
							exit(1);
						}
						if(c < 0) {
							// EOF occurred in the middle of a read - abort
							seqan::clear(r.patFw);
							return;
						}
						if (c != '\r' && c != '\n') {
							if (qualsRead >= _trim5) {
								size_t off = qualsRead - _trim5;
								if(off + 1 > 1024) {
									cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
									     << "Please truncate reads and quality values and and re-run Bowtie";
									exit(1);
								}

								if (_solexa_quals)
								{
									// Convert solexa-scaled chars to phred
									// http://maq.sourceforge.net/fastq.shtml
									int pQ = (int)(10.0 * log(1.0 + pow(10.0, ((int)c - 64) / 10.0)) / log(10.0) + .499) + 33;
									c = (char)(pQ);
								}
								else
								{
									// Keep the phred quality
									if (c < 33)
									{
										cerr << "Saw ASCII character " << ((int)c) << "." << endl;
										wrongQualityScale();
										exit(1);
									}
								}

								assert_geq(c, 33);
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
						c = _filebuf.get();
						if (c == ' ') {
							wrongQualityFormat();
							exit(1);
						}
						if(c < 0) {
							// EOF occurred in the middle of a read - abort
							seqan::clear(r.patFw);
							return;
						}
						if (c != '\r' && c != '\n') {
							if (qualsRead >= _trim5) {
								size_t off = qualsRead - _trim5;
								if(off + 1 > 1024) {
									cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
									     << "Please truncate reads and quality values and and re-run Bowtie";
									exit(1);
								}

								if (_solexa_quals)
								{
									// Convert solexa-scaled chars to phred
									// http://maq.sourceforge.net/fastq.shtml
									int pQ = (int)(10.0 * log(1.0 + pow(10.0, ((int)c - 64) / 10.0)) / log(10.0) + .499) + 33;
									c = (char)(pQ);
								}
								else
								{
									// Keep the phred quality
									if (c < 33)
									{
										cerr << "Saw ASCII character " << ((int)c) << "." << endl;
										wrongQualityScale();
										exit(1);
									}
								}

								assert_geq(c, 33);
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
	bool _integer_quals;
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
	                 bool __useSpinlock = true,
	                 const char *__dumpfile = NULL,
	                 int __trim3 = 0,
	                 int __trim5 = 0,
	                 int __policy = NS_TO_NS,
	                 int __maxNs = 9999,
	                 uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, false, __useSpinlock, __dumpfile, __trim3, __trim5),
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
				if(dna4Cat[c] == 0) {
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
						if(dstLen + 1 > 1024) {
							cerr << "Reads file contained a pattern with more than 1024 characters." << endl
							     << "Please truncate reads and and re-run Bowtie";
							exit(1);
						}
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
						r.patBufFw [dstLen] = charToDna5[c];
						r.qualBufFw[dstLen] = 'I';
						r.patBufRc [bufSz-dstLen-1] = rcCharToDna5[c];
						r.qualBufRc[bufSz-dstLen-1] = 'I';
						dstLen++;
					}
					c = _filebuf.get();
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[i] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)r.patBufFw);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, r.qualBufFw);
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)&r.patBufRc[bufSz-dstLen]);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, &r.qualBufRc[bufSz-dstLen]);
				_setLength(r.qualRc, dstLen);
			} else {
				while(!isspace(c) && c >= 0) {
					if(isalpha(c) && dstLen >= this->_trim5) {
						if(dstLen + 1 > 1024) {
							cerr << "Reads file contained a pattern with more than 1024 characters.." << endl
							     << "Please truncate reads and and re-run Bowtie";
							exit(1);
						}
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
						r.patBufFw [bufSz-dstLen-1] = charToDna5[c];
						r.qualBufFw[bufSz-dstLen-1] = 'I';
						r.patBufRc [dstLen] = rcCharToDna5[c];
						r.qualBufRc[dstLen] = 'I';
						dstLen++;
					}
					c = _filebuf.get();
				}
				dstLen -= this->_trim3;
				// Now that we've trimmed on both ends, count the Ns
				for(int i = 0; i < dstLen; i++) {
					if(r.patBufFw[bufSz-i-1] == 4) ns++;
				}
				_setBegin (r.patFw,  (Dna5*)&r.patBufFw[bufSz-dstLen]);
				_setLength(r.patFw,  dstLen);
				_setBegin (r.qualFw, &r.qualBufFw[bufSz-dstLen]);
				_setLength(r.qualFw, dstLen);
				_setBegin (r.patRc,  (Dna5*)r.patBufRc);
				_setLength(r.patRc,  dstLen);
				_setBegin (r.qualRc, r.qualBufRc);
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
