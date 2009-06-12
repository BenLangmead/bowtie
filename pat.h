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
#include "qual.h"

/**
 * Classes and routines for reading reads from various input sources.
 */

using namespace std;
using namespace seqan;

/// Constructs string base-10 representation of integer 'value'
extern char* itoa10(int value, char* result);

/**
 * Calculate a per-read random seed based on a combination of
 * the read data (incl. sequence, name, quals) and the global
 * seed in '_randSeed'.
 */
static inline uint32_t genRandSeed(const String<Dna5>& qry,
                                   const String<char>& qual,
                                   const String<char>& name)
{
	// Calculate a per-read random seed based on a combination of
	// the read data (incl. sequence, name, quals) and the global
	// seed
	uint32_t rseed = 0;
	size_t qlen = seqan::length(qry);
	// Throw all the characters of the read into the random seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qry[i];
		assert_leq(p, 4);
		size_t off = ((i & 15) << 1);
		rseed |= (p << off);
	}
	// Throw all the quality values for the read into the random
	// seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qual[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed |= (p << off);
	}
	// Throw all the characters in the read name into the random
	// seed
	size_t namelen = seqan::length(name);
	for(size_t i = 0; i < namelen; i++) {
		int p = (int)name[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed |= (p << off);
	}
	return rseed;
}

/**
 * A buffer for keeping all relevant information about a single read.
 * Each search thread has one.
 */
struct ReadBuf {
	~ReadBuf() {
		clearAll(); reset();
		// Prevent seqan from trying to free buffers
		_setBegin(patFw, NULL);
		_setBegin(patRc, NULL);
		_setBegin(qualFw, NULL);
		_setBegin(qualRc, NULL);
		_setBegin(patFwRev, NULL);
		_setBegin(patRcRev, NULL);
		_setBegin(qualFwRev, NULL);
		_setBegin(qualRcRev, NULL);
		_setBegin(name, NULL);
	}

	/// Point all Strings to the beginning of their respective buffers
	/// and set all lengths to 0
	void reset() {
		patid = 0;
		readOrigBufLen = 0;
		_setBegin(patFw,     (Dna5*)patBufFw);     _setLength(patFw, 0);     _setCapacity(patFw, BUF_SIZE);
		_setBegin(patRc,     (Dna5*)patBufRc);     _setLength(patRc, 0);     _setCapacity(patRc, BUF_SIZE);
		_setBegin(qualFw,    (char*)qualBufFw);    _setLength(qualFw, 0);    _setCapacity(qualFw, BUF_SIZE);
		_setBegin(qualRc,    (char*)qualBufRc);    _setLength(qualRc, 0);    _setCapacity(qualRc, BUF_SIZE);
		_setBegin(patFwRev,  (Dna5*)patBufFwRev);  _setLength(patFwRev, 0);  _setCapacity(patFwRev, BUF_SIZE);
		_setBegin(patRcRev,  (Dna5*)patBufRcRev);  _setLength(patRcRev, 0);  _setCapacity(patRcRev, BUF_SIZE);
		_setBegin(qualFwRev, (char*)qualBufFwRev); _setLength(qualFwRev, 0); _setCapacity(qualFwRev, BUF_SIZE);
		_setBegin(qualRcRev, (char*)qualBufRcRev); _setLength(qualRcRev, 0); _setCapacity(qualRcRev, BUF_SIZE);
		_setBegin(name,      (char*)nameBuf);      _setLength(name, 0);      _setCapacity(name, BUF_SIZE);
	}

	void clearAll() {
		seqan::clear(patFw);
		seqan::clear(patRc);
		seqan::clear(qualFw);
		seqan::clear(qualRc);
		seqan::clear(patFwRev);
		seqan::clear(patRcRev);
		seqan::clear(qualFwRev);
		seqan::clear(qualRcRev);
		seqan::clear(name);
		readOrigBufLen = 0;
	}

	/// Return true iff the read (pair) is empty
	bool empty() {
		return seqan::empty(patFw);
	}

	/// Return length of the read in the buffer
	uint32_t length() {
		return seqan::length(patFw);
	}

	/**
	 * Given patFw and qualFw, construct patRc and qualRc in place.
	 */
	void constructRevComps() {
		uint32_t len = length();
		assert_gt(len, 0);
		_setBegin(patRc,  (Dna5*)patBufRc);
		_setBegin(qualRc, (char*)qualBufRc);
		_setLength(patRc,  len);
		_setLength(qualRc, len);
		_setCapacity(patRc,  BUF_SIZE);
		_setCapacity(qualRc, BUF_SIZE);
		for(uint32_t i = 0; i < len; i++) {
			// Reverse-complement the sequence
			patBufRc[i]  = (patBufFw[len-i-1] == 4) ? 4 : (patBufFw[len-i-1] ^ 3);
			// Reverse the quality
			qualBufRc[i] = qualBufFw[len-i-1];
		}
	}

	/**
	 * Given patFw, patRc, qualFw and qualRc, construct the *Rev
	 * versions in place.  Assumes constructRevComps() was called
	 * previously.
	 */
	void constructReverses() {
		uint32_t len = length();
		assert_gt(len, 0);
		_setBegin(patFwRev,  (Dna5*)patBufFwRev);
		_setBegin(patRcRev,  (Dna5*)patBufRcRev);
		_setBegin(qualFwRev, (char*)qualBufFwRev);
		_setBegin(qualRcRev, (char*)qualBufRcRev);
		_setLength(patFwRev,  len);
		_setLength(patRcRev,  len);
		_setLength(qualFwRev, len);
		_setLength(qualRcRev, len);
		_setCapacity(patFwRev,  BUF_SIZE);
		_setCapacity(patRcRev,  BUF_SIZE);
		_setCapacity(qualFwRev, BUF_SIZE);
		_setCapacity(qualRcRev, BUF_SIZE);
		for(uint32_t i = 0; i < len; i++) {
			patFwRev[i]  = patFw[len-i-1];
			patRcRev[i]  = patRc[len-i-1];
			qualFwRev[i] = qualFw[len-i-1];
			qualRcRev[i] = qualRc[len-i-1];
		}
	}

	/**
	 * Append a "/1" or "/2" string onto the end of the name buf if
	 * it's not already there.
	 */
	void fixMateName(int i) {
		assert(i == 1 || i == 2);
		size_t namelen = seqan::length(name);
		bool append = false;
		if(namelen < 2) {
			// Name is too short to possibly have /1 or /2 on the end
			append = true;
		} else {
			if(i == 1) {
				// append = true iff mate name does not already end in /1
				append =
					nameBuf[namelen-2] != '/' ||
					nameBuf[namelen-1] != '1';
			} else {
				// append = true iff mate name does not already end in /2
				append =
					nameBuf[namelen-2] != '/' ||
					nameBuf[namelen-1] != '2';
			}
		}
		if(append) {
			assert_leq(namelen, BUF_SIZE-2);
			_setLength(name, namelen + 2);
			nameBuf[namelen] = '/';
			nameBuf[namelen+1] = "012"[i];
		}
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

	String<Dna5>  patFwRev;               // forward-strand sequence reversed
	uint8_t       patBufFwRev[BUF_SIZE];  // forward-strand sequence buffer reversed
	String<Dna5>  patRcRev;               // reverse-complement sequence reversed
	uint8_t       patBufRcRev[BUF_SIZE];  // reverse-complement sequence buffer reversed
	String<char>  qualFwRev;              // quality values reversed
	char          qualBufFwRev[BUF_SIZE]; // quality value buffer reversed
	String<char>  qualRcRev;              // reverse quality values reversed
	char          qualBufRcRev[BUF_SIZE]; // reverse quality value buffer reversed

	char          readOrigBuf[FileBuf::LASTN_BUF_SZ];
	size_t        readOrigBufLen;

	String<char>  name;                // read name
	char          nameBuf[BUF_SIZE];   // read name buffer
	uint32_t      patid;               // unique 0-based id based on order in read file(s)
	int           mate;                // 0 = single-end, 1 = mate1, 2 = mate2
	uint32_t      seed;                // random seed
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
	PatternSource(bool randomizeQuals = false,
	              bool useSpinlock = true,
	              const char *dumpfile = NULL) :
	    readCnt_(0),
		dumpfile_(dumpfile),
		numWrappers_(0),
		doLocking_(true),
		useSpinlock_(useSpinlock),
		randomizeQuals_(randomizeQuals),
		lock_()
	{
		// Open dumpfile, if specified
		if(dumpfile_ != NULL) {
			out_.open(dumpfile_, ios_base::out);
			if(!out_.good()) {
				cerr << "Could not open pattern dump file \"" << dumpfile_ << "\" for writing" << endl;
				exit(1);
			}
		}
		MUTEX_INIT(lock_);
	}

	virtual ~PatternSource() { }

	/**
	 * Call this whenever this PatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks will be contended.
	 */
	void addWrapper() {
		numWrappers_++;
	}

	/**
	 * The main member function for dispensing patterns.
	 */
	virtual void nextReadPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// nextPatternImpl does the reading from the ultimate source;
		// it is implemented in concrete subclasses
		nextReadPairImpl(ra, rb, patid);
		if(!ra.empty()) {
			// Possibly randomize the qualities so that they're more
			// scattered throughout the range of possible values
			if(randomizeQuals_) {
				randomizeQuals(ra);
				if(!rb.empty()) {
					randomizeQuals(rb);
				}
			}
			// TODO: Perhaps bundle all of the following up into a
			// finalize() member in the ReadBuf class?

			// Construct the reversed versions of the fw and rc seqs
			// and quals
			ra.constructRevComps();
			ra.constructReverses();
			if(!rb.empty()) {
				rb.constructRevComps();
				rb.constructReverses();
			}
			// Fill in the random-seed field using a combination of
			// information from the user-specified seed and the read
			// sequence, qualities, and name
			ra.seed = genRandSeed(ra.patFw, ra.qualFw, ra.name);
			if(!rb.empty()) {
				rb.seed = genRandSeed(rb.patFw, rb.qualFw, rb.name);
			}
			// Output it, if desired
			if(dumpfile_ != NULL) {
				dumpBuf(ra);
				if(!rb.empty()) {
					dumpBuf(rb);
				}
			}
		}
	}

	/**
	 * The main member function for dispensing patterns.
	 */
	virtual void nextRead(ReadBuf& r, uint32_t& patid) {
		// nextPatternImpl does the reading from the ultimate source;
		// it is implemented in concrete subclasses
		nextReadImpl(r, patid);
		if(!r.empty()) {
			// Possibly randomize the qualities so that they're more
			// scattered throughout the range of possible values
			if(randomizeQuals_) {
				randomizeQuals(r);
			}
			// Construct the reversed versions of the fw and rc seqs
			// and quals
			r.constructRevComps();
			r.constructReverses();
			// Fill in the random-seed field using a combination of
			// information from the user-specified seed and the read
			// sequence, qualities, and name
			r.seed = genRandSeed(r.patFw, r.qualFw, r.name);
			// Output it, if desired
			if(dumpfile_ != NULL) {
				dumpBuf(r);
			}
		}
	}

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats that
	 * can read in a pair of reads in a single transaction with a
	 * single input source.  If paired-end input is given as a pair of
	 * parallel files, this member should throw an error and exit.
	 */
	virtual void nextReadPairImpl(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) = 0;

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats
	 * where individual input sources look like single-end-read
	 * sources, e.g., formats where paired-end reads are specified in
	 * parallel read files.
	 */
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) = 0;

	/// Reset state to start over again with the first read
	virtual void reset() { readCnt_ = 0; }

	/**
	 * Concrete subclasses call lock() to enter a critical region.
	 * What constitutes a critical region depends on the subclass.
	 */
	void lock() {
		if(!doLocking_) return; // no contention
#ifdef USE_SPINLOCK
		if(useSpinlock_) {
			// User can ask to use the normal pthreads lock even if
			// spinlocks are compiled in.
			spinlock_.Enter();
		} else {
#endif
			MUTEX_LOCK(lock_);
#ifdef USE_SPINLOCK
		}
#endif
	}

	/**
	 * Concrete subclasses call unlock() to exit a critical region
	 * What constitutes a critical region depends on the subclass.
	 */
	void unlock() {
		if(!doLocking_) return; // no contention
#ifdef USE_SPINLOCK
		if(useSpinlock_) {
			// User can ask to use the normal pthreads lock even if
			// spinlocks are compiled in.
			spinlock_.Leave();
		} else {
#endif
			MUTEX_UNLOCK(lock_);
#ifdef USE_SPINLOCK
		}
#endif
	}

protected:

	/**
	 * Mix up the quality values for ReadBuf r.  There's probably a
	 * more (pseudo-)randomly rigorous way to do this; the output looks
	 * pretty cyclic.
	 */
	void randomizeQuals(ReadBuf& r) {
		const size_t rlen = r.length();
		for(size_t i = 0; i < rlen; i++) {
			if(i < rlen-1) {
				r.qualFw[i] *= (r.qualFw[i+1] + 7);
			}
			if(i > 0) {
				r.qualFw[i] *= (r.qualFw[i-1] + 11);
			}
			// A user says that g++ complained here about "comparison
			// is always false due to limited range of data type", but
			// I can't see why.  I added the (int) cast to try to avoid
			// the warning.
			if((int)r.qualFw[i] < 0) r.qualFw[i] = -(r.qualFw[i]+1);
			r.qualFw[i] %= 41;
			assert_leq(r.qualFw[i], 40);
			r.qualFw[i] += 33;
		}
	}

	/**
	 * Dump the contents of the ReadBuf to the dump file.
	 */
	void dumpBuf(const ReadBuf& r) {
		assert(dumpfile_ != NULL);
		dump(out_, r.patFw,
		     empty(r.qualFw) ? String<char>("(empty)") : r.qualFw,
		     empty(r.name)   ? String<char>("(empty)") : r.name);
		dump(out_, r.patRc,
		     empty(r.qualRc) ? String<char>("(empty)") : r.qualRc,
		     empty(r.name)   ? String<char>("(empty)") : r.name);
	}

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
	uint32_t readCnt_;

	const char *dumpfile_; /// dump patterns to this file before returning them
	ofstream out_;         /// output stream for dumpfile
	int numWrappers_;      /// # threads that own a wrapper for this PatternSource
	bool doLocking_;       /// override whether to lock (true = don't override)
	/// User can ask to use the normal pthreads-style lock even if
	/// spinlocks is enabled and compiled in.  This is sometimes better
	/// if we expect bad I/O latency on some reads.
	bool useSpinlock_;
	bool randomizeQuals_;  /// true -> mess up qualities in a random way
#ifdef USE_SPINLOCK
	SpinLock spinlock_;
#endif
	MUTEX_T lock_; /// mutex for locking critical regions
};

/**
 * Abstract parent class for synhconized sources of paired-end reads
 * (and possibly also single-end reads).
 */
class PairedPatternSource {
public:
	PairedPatternSource() {
		MUTEX_INIT(lock_);
	}
	virtual ~PairedPatternSource() { }

	virtual void addWrapper() = 0;
	virtual void reset() = 0;
	virtual bool nextReadPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) = 0;

	/**
	 * Lock this PairedPatternSource, usually because one of its shared
	 * fields is being updated.
	 */
	void lock() {
#ifdef USE_SPINLOCK
		spinlock_.Enter();
#else
		MUTEX_LOCK(lock_);
#endif
	}

	/**
	 * Unlock this PairedPatternSource.
	 */
	void unlock() {
#ifdef USE_SPINLOCK
		spinlock_.Leave();
#else
		MUTEX_UNLOCK(lock_);
#endif
	}

protected:

#ifdef USE_SPINLOCK
	SpinLock spinlock_;
#endif
	MUTEX_T lock_; /// mutex for locking critical regions
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class PairedSoloPatternSource : public PairedPatternSource {

public:

	PairedSoloPatternSource(const vector<PatternSource*>& src) :
		cur_(0), src_(src)
	{
	    for(size_t i = 0; i < src_.size(); i++) {
	    	assert(src_[i] != NULL);
	    }
	}

	virtual ~PairedSoloPatternSource() { }

	/**
	 * Call this whenever this PairedPatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks within PatternSources will be contended.
	 */
	virtual void addWrapper() {
		for(size_t i = 0; i < src_.size(); i++) {
			src_[i]->addWrapper();
		}
	}

	/**
	 * Reset this object and all the PatternSources under it so that
	 * the next call to nextReadPair gets the very first read pair.
	 */
	virtual void reset() {
		for(size_t i = 0; i < src_.size(); i++) {
			src_[i]->reset();
		}
		cur_ = 0;
	}

	/**
	 * The main member function for dispensing pairs of reads or
	 * singleton reads.  Returns true iff ra and rb contain a new
	 * pair; returns false if ra contains a new unpaired read.
	 */
	virtual bool nextReadPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		uint32_t cur = cur_;
		while(cur < src_.size()) {
			// Patterns from srca_[cur_] are unpaired
			src_[cur]->nextReadPair(ra, rb, patid);
			if(seqan::empty(ra.patFw)) {
				// If patFw is empty, that's our signal that the
				// input dried up
				lock();
				if(cur + 1 > cur_) cur_++;
				cur = cur_;
				unlock();
				continue; // on to next pair of PatternSources
			}
			if(!rb.empty()) {
				ra.fixMateName(1);
				rb.fixMateName(2);
			}
			ra.patid = patid;
			ra.mate  = 1;
			rb.mate  = 2;
			return true; // paired
		}
		return false;
	}

protected:

	volatile uint32_t cur_; // current element in parallel srca_, srcb_ vectors
	vector<PatternSource*> src_; /// PatternSources for paired-end reads
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class PairedDualPatternSource : public PairedPatternSource {

public:

	PairedDualPatternSource(const vector<PatternSource*>& srca,
	                        const vector<PatternSource*>& srcb) :
		cur_(0), srca_(srca), srcb_(srcb)
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

	virtual ~PairedDualPatternSource() { }

	/**
	 * Call this whenever this PairedPatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks within PatternSources will be contended.
	 */
	virtual void addWrapper() {
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
	virtual void reset() {
		for(size_t i = 0; i < srca_.size(); i++) {
			srca_[i]->reset();
			if(srcb_[i] != NULL) {
				srcb_[i]->reset();
			}
		}
		cur_ = 0;
	}

	/**
	 * The main member function for dispensing pairs of reads or
	 * singleton reads.  Returns true iff ra and rb contain a new
	 * pair; returns false if ra contains a new unpaired read.
	 */
	virtual bool nextReadPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		uint32_t cur = cur_;
		while(cur < srca_.size()) {
			if(srcb_[cur] == NULL) {
				// Patterns from srca_[cur_] are unpaired
				srca_[cur]->nextRead(ra, patid);
				if(seqan::empty(ra.patFw)) {
					// If patFw is empty, that's our signal that the
					// input dried up
					lock();
					if(cur + 1 > cur_) cur_++;
					cur = cur_;
					unlock();
					continue; // on to next pair of PatternSources
				}
				ra.patid = patid;
				ra.mate  = 0;
				return false; // unpaired
			} else {
				// Patterns from srca_[cur_] and srcb_[cur_] are paired
				uint32_t patid_a = 0;
				uint32_t patid_b = 0;
				// Lock to ensure that this thread gets parallel reads
				// in the two mate files
				lock();
				srca_[cur]->nextRead(ra, patid_a);
				srcb_[cur]->nextRead(rb, patid_b);
				bool cont = false;
				// Did the pair obtained fail to match up?
				while(patid_a != patid_b) {
					// Is either input exhausted?  If so, bail.
					if(seqan::empty(ra.patFw) || seqan::empty(rb.patFw)) {
						seqan::clear(ra.patFw);
						lock();
						if(cur + 1 > cur_) cur_++;
						cur = cur_;
						unlock();
						cont = true;
						break;
					}
					if(patid_a < patid_b) {
						srca_[cur]->nextRead(ra, patid_a);
						ra.fixMateName(1);
					} else {
						srcb_[cur]->nextRead(rb, patid_b);
						rb.fixMateName(2);
					}
				}
				unlock();
				if(cont) continue; // on to next pair of PatternSources
				ra.fixMateName(1);
				rb.fixMateName(2);
				if(seqan::empty(ra.patFw)) {
					// If patFw is empty, that's our signal that the
					// input dried up
					lock();
					if(cur + 1 > cur_) cur_++;
					cur = cur_;
					unlock();
					continue; // on to next pair of PatternSources
				}
				assert_eq(patid_a, patid_b);
				patid = patid_a;
				ra.patid = patid;
				rb.patid = patid;
				ra.mate  = 1;
				rb.mate  = 2;
				return true; // paired
			}
		}
		return false;
	}

protected:

	volatile uint32_t cur_; // current element in parallel srca_, srcb_ vectors
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
		buf1_(), buf2_(), patid_(0xffffffff) { }

	virtual ~PatternSourcePerThread() { }

	/**
	 * Read the next read pair.
	 */
	virtual void nextReadPair() { }

	ReadBuf& bufa()        { return buf1_;         }
	ReadBuf& bufb()        { return buf2_;         }

	uint32_t      patid()  { return patid_;        }
	virtual void  reset()  { patid_ = 0xffffffff;  }
	bool          empty()  { return buf1_.empty(); }

	/**
	 * Return true iff the buffers jointly contain a paired-end read.
	 */
	bool paired() {
		bool ret = !buf2_.empty();
	    assert(!ret || !empty());
	    return ret;
	}

protected:
	ReadBuf  buf1_;    // read buffer for mate a
	ReadBuf  buf2_;    // read buffer for mate b
	uint32_t patid_;   // index of read just read
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class PatternSourcePerThreadFactory {
public:
	virtual ~PatternSourcePerThreadFactory() { }
	virtual PatternSourcePerThread* create() const = 0;
	virtual std::vector<PatternSourcePerThread*>* create(uint32_t n) const = 0;

	/// Free memory associated with a pattern source
	virtual void destroy(PatternSourcePerThread* patsrc) const {
		assert(patsrc != NULL);
		// Free the PatternSourcePerThread
		delete patsrc;
	}

	/// Free memory associated with a pattern source list
	virtual void destroy(std::vector<PatternSourcePerThread*>* patsrcs) const {
		assert(patsrcs != NULL);
		// Free all of the PatternSourcePerThreads
		for(size_t i = 0; i < patsrcs->size(); i++) {
			if((*patsrcs)[i] != NULL) {
				delete (*patsrcs)[i];
				(*patsrcs)[i] = NULL;
			}
		}
		// Free the vector
		delete patsrcs;
	}
};

/**
 * A per-thread wrapper for a PairedPatternSource.
 */
class WrappedPatternSourcePerThread : public PatternSourcePerThread {
public:
	WrappedPatternSourcePerThread(PairedPatternSource& __patsrc) :
		patsrc_(__patsrc)
	{
		patsrc_.addWrapper();
	}

	/**
	 * Get the next paired or unpaired read from the wrapped
	 * PairedPatternSource.
	 */
	virtual void nextReadPair() {
		PatternSourcePerThread::nextReadPair();
		ASSERT_ONLY(uint32_t lastPatid = patid_);
		buf1_.clearAll();
		buf2_.clearAll();
		patsrc_.nextReadPair(buf1_, buf2_, patid_);
		assert(buf1_.empty() || patid_ != lastPatid);
	}

private:

	/// Container for obtaining paired reads from PatternSources
	PairedPatternSource& patsrc_;
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class WrappedPatternSourcePerThreadFactory : public PatternSourcePerThreadFactory {
public:
	WrappedPatternSourcePerThreadFactory(PairedPatternSource& patsrc) :
		patsrc_(patsrc) { }

	/**
	 * Create a new heap-allocated WrappedPatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new WrappedPatternSourcePerThread(patsrc_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * WrappedPatternSourcePerThreads.
	 */
	virtual std::vector<PatternSourcePerThread*>* create(uint32_t n) const {
		std::vector<PatternSourcePerThread*>* v = new std::vector<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new WrappedPatternSourcePerThread(patsrc_));
			assert(v->back() != NULL);
		}
		return v;
	}

private:
	/// Container for obtaining paired reads from PatternSources
	PairedPatternSource& patsrc_;
};

/**
 * Encapsualtes a source of patterns where each raw pattern is trimmed
 * by some user-defined amount on the 3' and 5' ends.  Doesn't
 * implement the actual trimming - that's up to the concrete
 * descendants.
 */
class TrimmingPatternSource : public PatternSource {
public:
	TrimmingPatternSource(bool randomizeQuals = false,
	                      bool useSpinlock = true,
	                      const char *dumpfile = NULL,
	                      int trim3 = 0,
	                      int trim5 = 0) :
		PatternSource(randomizeQuals, useSpinlock, dumpfile),
		trim3_(trim3),
		trim5_(trim5) { }
protected:
	int trim3_;
	int trim5_;
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
	                    bool useSpinlock = true,
	                    const char *dumpfile = NULL,
	                    uint32_t seed = 0) :
		PatternSource(false, useSpinlock, dumpfile),
		numReads_(numReads),
		length_(length),
		seed_(seed)
	{
		if(length_ > 1024) {
			cerr << "Read length for RandomPatternSource may not exceed 1024; got " << length_ << endl;
			exit(1);
		}
		rand_.init(seed_);
	}

	/** Get the next random read and set patid */
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		// Begin critical section
		lock();
		if(readCnt_ >= numReads_) {
			r.clearAll();
			unlock();
			return;
		}
		uint32_t ra = rand_.nextU32();
		patid = readCnt_;
		readCnt_++;
		unlock();
		fillRandomRead(r, ra, length_, patid);
	}

	/** Get the next random read and set patid */
	virtual void nextReadPairImpl(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// Begin critical section
		lock();
		if(readCnt_ >= numReads_) {
			ra.clearAll();
			rb.clearAll();
			unlock();
			return;
		}
		uint32_t rna = rand_.nextU32();
		uint32_t rnb = rand_.nextU32();
		patid = readCnt_;
		readCnt_++;
		unlock();
		fillRandomRead(ra, rna, length_, patid);
		fillRandomRead(rb, rnb, length_, patid);
	}

	/** */
	static void fillRandomRead(ReadBuf& r,
	                           uint32_t ra,
	                           int length,
	                           uint32_t patid)
	{
		// End critical section
		for(int i = 0; i < length; i++) {
			ra = RandomSource::nextU32(ra) >> 8;
			r.patBufFw[i]           = (ra & 3);
			char c                  = 'I' - ((ra >> 2) & 31);
			r.qualBufFw[i]          = c;
		}
		_setBegin (r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, length);
		_setBegin (r.qualFw, r.qualBufFw);
		_setLength(r.qualFw, length);
		itoa10(patid, r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, strlen(r.nameBuf));
	}

	/** Reset the pattern source to the beginning */
	virtual void reset() {
		PatternSource::reset();
		// reset pseudo-random generator; next string of calls to
		// nextU32() will return same pseudo-randoms as the last
		rand_.init(seed_);
	}
private:
	uint32_t     numReads_; /// number of reads to dish out
	int          length_;   /// length of reads
	uint32_t     seed_;     /// seed for pseudo-randoms
	RandomSource rand_;     /// pseudo-random generator
};

/**
 * A version of PatternSourcePerThread that dishes out random patterns
 * without any synchronization.
 */
class RandomPatternSourcePerThread : public PatternSourcePerThread {
public:
	RandomPatternSourcePerThread(uint32_t numreads,
	                             int length,
	                             int numthreads,
	                             int thread) :
		PatternSourcePerThread(),
		numreads_(numreads),
		length_(length),
		numthreads_(numthreads),
		thread_(thread)
	{
		patid_ = thread_;
		if(length_ > 1024) {
			cerr << "Read length for RandomPatternSourcePerThread may not exceed 1024; got " << length_ << endl;
			exit(1);
		}
		rand_.init(thread_);
	}

	virtual void nextReadPair() {
		PatternSourcePerThread::nextReadPair();
		if(patid_ >= numreads_) {
			buf1_.clearAll();
			buf2_.clearAll();
			return;
		}
		RandomPatternSource::fillRandomRead(
			buf1_, rand_.nextU32(), length_, patid_);
		RandomPatternSource::fillRandomRead(
			buf2_, rand_.nextU32(), length_, patid_);
		patid_ += numthreads_;
	}

	virtual void reset() {
		PatternSourcePerThread::reset();
		patid_ = thread_;
		rand_.init(thread_);
	}

private:
	uint32_t     numreads_;
	int          length_;
	int          numthreads_;
	int          thread_;
	RandomSource rand_;
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class RandomPatternSourcePerThreadFactory : public PatternSourcePerThreadFactory {
public:
	RandomPatternSourcePerThreadFactory(
			uint32_t numreads,
            int length,
            int numthreads,
            int thread) :
            numreads_(numreads),
            length_(length),
            numthreads_(numthreads),
            thread_(thread) { }

	/**
	 * Create a new heap-allocated WrappedPatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new RandomPatternSourcePerThread(
			numreads_, length_, numthreads_, thread_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * WrappedPatternSourcePerThreads.
	 */
	virtual std::vector<PatternSourcePerThread*>* create(uint32_t n) const {
		std::vector<PatternSourcePerThread*>* v = new std::vector<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new RandomPatternSourcePerThread(
				numreads_, length_, numthreads_, thread_));
			assert(v->back() != NULL);
		}
		return v;
	}

private:
	uint32_t numreads_;
    int length_;
    int numthreads_;
    int thread_;
};

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
static inline int peekOverNewline(FileBuf& in) {
	while(true) {
		int c = in.peek();
		if(c != '\r' && c != '\n') {
			return c;
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
static inline int peekToEndOfLine(FileBuf& in) {
	while(true) {
		int c = in.get(); if(c < 0) return c;
		if(c == '\n' || c == '\r') {
			c = in.peek();
			while(c == '\n' || c == '\r') {
				in.get(); if(c < 0) return c; // consume \r or \n
				c = in.peek();
			}
			// next get() gets first character of next line
			return c;
		}
	}
}

/**
 * Encapsualtes a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public TrimmingPatternSource {
public:
	VectorPatternSource(const vector<string>& v,
	                    bool randomizeQuals = false,
	                    bool useSpinlock = true,
	                    const char *dumpfile = NULL,
	                    int trim3 = 0,
	                    int trim5 = 0,
		                uint32_t skip = 0) :
		TrimmingPatternSource(randomizeQuals, useSpinlock, dumpfile, trim3, trim5),
		cur_(skip), skip_(skip), paired_(false), v_(), quals_()
	{
		for(size_t i = 0; i < v.size(); i++) {
			vector<string> ss;
			tokenize(v[i], ":", ss, 2);
			assert_gt(ss.size(), 0);
			assert_leq(ss.size(), 2);
			// Initialize s
			string s = ss[0];
			if(s.length() <= (size_t)(trim3_ + trim5_)) {
				// Entire read is trimmed away
				continue;
			} else {
				// Trim on 5' (high-quality) end
				if(trim5_ > 0) {
					s.erase(0, trim5_);
				}
				// Trim on 3' (low-quality) end
				if(trim3_ > 0) {
					s.erase(s.length()-trim3_);
				}
			}
			//  Initialize vq
			string vq;
			if(ss.size() == 2) {
				vq = ss[1];
			}
			// Trim qualities
			if(vq.length() > (size_t)(trim3_ + trim5_)) {
				// Trim on 5' (high-quality) end
				if(trim5_ > 0) {
					vq.erase(0, trim5_);
				}
				// Trim on 3' (low-quality) end
				if(trim3_ > 0) {
					vq.erase(vq.length()-trim3_);
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
			v_.push_back(s);
			quals_.push_back(vq);
			ostringstream os;
			os << (names_.size());
			names_.push_back(os.str());
		}
		assert_eq(v_.size(), quals_.size());
	}
	virtual ~VectorPatternSource() { }
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		// Let Strings begin at the beginning of the respective bufs
		r.reset();
		lock();
		if(cur_ >= v_.size()) {
			unlock();
			// Clear all the Strings, as a signal to the caller that
			// we're out of reads
			r.clearAll();
			assert(r.empty());
			return;
		}
		// Copy v_*, quals_* strings into the respective Strings
		r.patFw  = v_[cur_];
		r.qualFw = quals_[cur_];
		ostringstream os;
		os << cur_;
		r.name = os.str();
		cur_++;
		readCnt_++;
		patid = readCnt_;
		unlock();
	}
	/**
	 * This is unused, but implementation is given for completeness.
	 */
	virtual void nextReadPairImpl(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// Let Strings begin at the beginning of the respective bufs
		ra.reset();
		rb.reset();
		if(!paired_) {
			paired_ = true;
			cur_ <<= 1;
		}
		lock();
		if(cur_ >= v_.size()-1) {
			unlock();
			// Clear all the Strings, as a signal to the caller that
			// we're out of reads
			ra.clearAll();
			rb.clearAll();
			assert(ra.empty());
			assert(rb.empty());
			return;
		}
		// Copy v_*, quals_* strings into the respective Strings
		ra.patFw  = v_[cur_];
		ra.qualFw = quals_[cur_];
		cur_++;
		rb.patFw  = v_[cur_];
		rb.qualFw = quals_[cur_];
		ostringstream os;
		os << readCnt_;
		ra.name = os.str();
		rb.name = os.str();
		cur_++;
		readCnt_++;
		patid = readCnt_;
		unlock();
	}
	virtual void reset() {
		TrimmingPatternSource::reset();
		cur_ = skip_;
		paired_ = false;
	}
private:
	size_t cur_;
	uint32_t skip_;
	bool paired_;
	vector<String<Dna5> > v_;        /// forward sequences
	vector<String<char> > quals_;    /// quality values parallel to v_
	vector<String<char> > names_;    /// names
};

/**
 *
 */
class BufferedFilePatternSource : public TrimmingPatternSource {
public:
	BufferedFilePatternSource(const vector<string>& infiles,
	                          bool randomizeQuals = false,
	                          bool useSpinlock = true,
	                          bool __forgiveInput = false,
	                          const char *dumpfile = NULL,
	                          int trim3 = 0,
	                          int trim5 = 0,
	                          uint32_t skip = 0) :
		TrimmingPatternSource(randomizeQuals, useSpinlock, dumpfile, trim3, trim5),
		infiles_(infiles),
		filecur_(0),
		filebuf_(),
		forgiveInput_(__forgiveInput),
		skip_(skip),
		first_(true)
	{
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size(), false);
		open(); // open first file in the list
		filecur_++;
	}

	virtual ~BufferedFilePatternSource() {
		if(filebuf_.isOpen()) filebuf_.close();
	}

	/**
	 * Fill ReadBuf with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 */
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		lock();
		bool notDone = true;
		do {
			read(r, patid);
			// Try again if r is empty (indicating an error) and input
			// is not yet exhausted, OR if we have more reads to skip
			// over
			notDone = seqan::empty(r.patFw) && !filebuf_.eof();
		} while(notDone || (!filebuf_.eof() && patid < skip_));
		if(patid < skip_) {
			unlock();
			r.clearAll();
			assert(seqan::empty(r.patFw));
			return;
		}
		if(first_ && seqan::empty(r.patFw) && !forgiveInput_) {
			// No reads could be extracted from the first _infile
			cerr << "Warning: Could not find any reads in \"" << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(seqan::empty(r.patFw) && filecur_ < infiles_.size()) {
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			do {
				read(r, patid);
			} while((seqan::empty(r.patFw) && !filebuf_.eof()));
			assert_geq(patid, skip_);
			if(seqan::empty(r.patFw) && !forgiveInput_) {
				// No reads could be extracted from this _infile
				cerr << "Warning: Could not find any reads in \"" << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		// Leaving critical region
		unlock();
		// If r.patFw is empty, then the caller knows that we are
		// finished with the reads
	}
	/**
	 *
	 */
	virtual void nextReadPairImpl(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		lock();
		bool notDone = true;
		do {
			readPair(ra, rb, patid);
			// Try again if ra is empty (indicating an error) and input
			// is not yet exhausted, OR if we have more reads to skip
			// over
			notDone = seqan::empty(ra.patFw) && !filebuf_.eof();
		} while(notDone || (!filebuf_.eof() && patid < skip_));
		if(patid < skip_) {
			unlock();
			ra.clearAll();
			rb.clearAll();
			assert(seqan::empty(ra.patFw));
			return;
		}
		if(first_ && seqan::empty(ra.patFw) && !forgiveInput_) {
			// No reads could be extracted from the first _infile
			cerr << "Warning: Could not find any read pairs in \"" << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(seqan::empty(ra.patFw) && filecur_ < infiles_.size()) {
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			do {
				readPair(ra, rb, patid);
			} while((seqan::empty(ra.patFw) && !filebuf_.eof()));
			assert_geq(patid, skip_);
			if(seqan::empty(ra.patFw) && !forgiveInput_) {
				// No reads could be extracted from this _infile
				cerr << "Warning: Could not find any reads in \"" << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		// Leaving critical region
		unlock();
		// If ra.patFw is empty, then the caller knows that we are
		// finished with the reads
	}
	/**
	 * Reset state so that we read start reading again from the
	 * beginning of the first file.  Should only be called by the
	 * master thread.
	 */
	virtual void reset() {
		TrimmingPatternSource::reset();
		filecur_ = 0,
		open();
		filecur_++;
	}
protected:
	/// Read another pattern from the input file; this is overridden
	/// to deal with specific file formats
	virtual void read(ReadBuf& r, uint32_t& patid) = 0;
	/// Read another pattern pair from the input file; this is
	/// overridden to deal with specific file formats
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) = 0;
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() { }
	void open() {
		if(filebuf_.isOpen()) filebuf_.close();
		while(filecur_ < infiles_.size()) {
			// Open read
			FILE *in;
			if(infiles_[filecur_] == "-") {
				in = stdin;
			} else if((in = fopen(infiles_[filecur_].c_str(), "r")) == NULL) {
				if(!errs_[filecur_]) {
					cerr << "Warning: Could not open file \"" << infiles_[filecur_] << "\" for reading" << endl;
					errs_[filecur_] = true;
				}
				filecur_++;
				continue;
			}
			filebuf_.newFile(in);
			return;
		}
		exit(1);
	}
	vector<string> infiles_; /// filenames for read files
	vector<bool> errs_; /// whether we've already printed an error for each file
	size_t filecur_;   /// index into infiles_ of next file to read
	FileBuf filebuf_;  /// read file currently being read from
	bool forgiveInput_; /// try hard to parse input even if it's malformed
	uint32_t skip_;     /// number of reads to skip
	bool first_;
};

/**
 * Synchronized concrete pattern source for a list of FASTA files.
 */
class FastaPatternSource : public BufferedFilePatternSource {
public:
	FastaPatternSource(const vector<string>& infiles,
	                   bool randomizeQuals = false,
	                   bool useSpinlock = true,
	                   const char *dumpfile = NULL,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool __forgiveInput = false,
	                   uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, useSpinlock,
		                          __forgiveInput, dumpfile,
		                          trim3, trim5, skip),
		first_(true)
	{ }
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}
protected:
	/**
	 * Scan to the next FASTA record (starting with >) and return the first
	 * character of the record (which will always be >).
	 */
	static int skipToNextFastaRecord(FileBuf& in) {
		int c;
		while((c = in.get()) != '>') {
			if(in.eof()) return -1;
		}
		return c;
	}
	/// Read another pattern from a FASTA input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		int c;
		int dstLen = 0;
		int nameLen = 0;
		// Pick off the first carat
		c = filebuf_.get();
		if(c == -1) {
			r.clearAll();
			filebuf_.resetLastN();
			return;
		}
		assert_eq('>', c);
		if(first_) {
			if(c != '>') {
				if(forgiveInput_) {
					c = FastaPatternSource::skipToNextFastaRecord(filebuf_);
					if(c < 0) {
						r.clearAll();
						filebuf_.resetLastN();
						return;
					}
				} else {
					c = getOverNewline(filebuf_); if(c < 0) {
						r.clearAll();
						filebuf_.resetLastN();
						return;
					}
				}
			}
			if(c != '>') {
				cerr << "Error: reads file does not look like a FASTA file" << endl;
				exit(1);
			}
			assert(c == '>' || c == '#');
			first_ = false;
		}

		// Read to the end of the id line, sticking everything after the '>'
		// into *name
		while(true) {
			c = filebuf_.get(); if(c < 0) {
				r.clearAll();
				filebuf_.resetLastN();
				return;
			}
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					c = filebuf_.get(); if(c < 0) {
						r.clearAll();
						filebuf_.resetLastN();
						return;
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
		while(c != '>') {
			if(isalpha(c) && begin++ >= this->trim5_) {
				if(dstLen + 1 > 1024) {
					cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
						 << "reads and re-run Bowtie";
					exit(1);
				}
				r.patBufFw [dstLen] = charToDna5[c];
				r.qualBufFw[dstLen] = 'I';
				dstLen++;
			}
			c = filebuf_.get();
			if(c == '\r' || c == '\n' || c == -1) {
				// Either we hit EOL or EOF; time to copy the buffered
				// read data into the 'orig' buffer.
				c = peekOverNewline(filebuf_);
				r.readOrigBufLen = filebuf_.copyLastN(r.readOrigBuf);
				filebuf_.resetLastN();
			}
			if(c < 0) break;
		}
		dstLen -= this->trim3_;
		_setBegin (r.patFw,  (Dna5*)r.patBufFw);
		_setLength(r.patFw,  dstLen);
		_setBegin (r.qualFw, r.qualBufFw);
		_setLength(r.qualFw, dstLen);

		// Set up a default name if one hasn't been set
		if(nameLen == 0) {
			itoa10(readCnt_, r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			nameLen = strlen(r.nameBuf);
			_setLength(r.name, nameLen);
		}
		assert_gt(nameLen, 0);
		readCnt_++;
		patid = readCnt_-1;
	}
	/// Read another pair of patterns from a FASTA input file
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In FastaPatternSource.readPair()" << endl;
		exit(1);
		read(ra, patid);
		readCnt_--;
		read(rb, patid);
	}
	virtual void resetForNextFile() {
		first_ = true;
	}
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << ">" << name << endl << seq << endl;
	}
private:
	bool first_;
	int policy_;
};

extern void wrongQualityFormat();
extern void tooFewQualities(const String<char>& read_name);

/**
 * Synchronized concrete pattern source for a list of files with tab-
 * delimited name, seq, qual fields (or, for paired-end reads,
 * basename, seq1, qual1, seq2, qual2).
 */
class TabbedPatternSource : public BufferedFilePatternSource {
public:
	TabbedPatternSource(const vector<string>& infiles,
	                   bool randomizeQuals = false,
	                   bool useSpinlock = true,
	                   const char *dumpfile = NULL,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool forgiveInput = false,
					   bool solQuals = false,
					   bool phred64Quals = false,
					   bool intQuals = false,
	                   uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, useSpinlock,
		                          forgiveInput, dumpfile,
		                          trim3, trim5, skip),
		solQuals_(solQuals),
		phred64Quals_(phred64Quals),
		intQuals_(intQuals)
	{ }

protected:

	/// Read another pattern from a FASTA input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		// filebuf_ is about to dish out the first character of the
		// name field
		if(parseName(r, NULL, '\t') == -1) {
			peekOverNewline(filebuf_); // skip rest of line
			r.clearAll();
			return;
		}
		assert_neq('\t', filebuf_.peek());

		// filebuf_ is about to dish out the first character of the
		// sequence field
		int charsRead = 0;
		int dstLen = parseSeq(r, charsRead, '\t');
		assert_neq('\t', filebuf_.peek());
		if(dstLen <= 0) {
			peekOverNewline(filebuf_); // skip rest of line
			r.clearAll();
			return;
		}

		// filebuf_ is about to dish out the first character of the
		// quality-string field
		char ct = 0;
		if(parseQuals(r, charsRead, dstLen, ct, '\n') <= 0) {
			peekOverNewline(filebuf_); // skip rest of line
			r.clearAll();
			return;
		}
		assert_eq(ct, '\n');
		assert_neq('\n', filebuf_.peek());
		r.readOrigBufLen = filebuf_.copyLastN(r.readOrigBuf);
		filebuf_.resetLastN();
		// The last character read in parseQuals should have been a
		// '\n'

		readCnt_++;
		patid = readCnt_-1;
	}

	/// Read another pair of patterns from a FASTA input file
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// filebuf_ is about to dish out the first character of the
		// name field
		if(parseName(ra, &rb, '\t') == -1) {
			peekOverNewline(filebuf_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			filebuf_.resetLastN();
			return;
		}
		assert_neq('\t', filebuf_.peek());

		// filebuf_ is about to dish out the first character of the
		// sequence field for the first mate
		int charsRead1 = 0;
		int dstLen1 = parseSeq(ra, charsRead1, '\t');
		if(dstLen1 <= -1) {
			peekOverNewline(filebuf_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			filebuf_.resetLastN();
			return;
		}
		assert_neq('\t', filebuf_.peek());

		// filebuf_ is about to dish out the first character of the
		// quality-string field
		char ct = 0;
		if(parseQuals(ra, charsRead1, dstLen1, ct, '\t', '\n') <= 0) {
			peekOverNewline(filebuf_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			filebuf_.resetLastN();
			return;
		}
		assert(ct == '\t' || ct == '\n');
		if(ct == '\n') {
			rb.clearAll();
			peekOverNewline(filebuf_);
			ra.readOrigBufLen = filebuf_.copyLastN(ra.readOrigBuf);
			filebuf_.resetLastN();
			readCnt_++;
			patid = readCnt_-1;
			return;
		}
		assert_neq('\t', filebuf_.peek());

		// filebuf_ is about to dish out the first character of the
		// sequence field for the second mate
		int charsRead2 = 0;
		int dstLen2 = parseSeq(rb, charsRead2, '\t');
		if(dstLen2 <= 0) {
			peekOverNewline(filebuf_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			filebuf_.resetLastN();
			return;
		}
		assert_neq('\t', filebuf_.peek());

		// filebuf_ is about to dish out the first character of the
		// quality-string field
		if(parseQuals(rb, charsRead2, dstLen2, ct, '\n') <= 0) {
			peekOverNewline(filebuf_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			filebuf_.resetLastN();
			return;
		}
		assert_eq('\n', ct);
		if(filebuf_.peek() == '\n') {
			if(!forgiveInput_) {
				assert(false);
			}
		}
		peekOverNewline(filebuf_);
		ra.readOrigBufLen = filebuf_.copyLastN(ra.readOrigBuf);
		filebuf_.resetLastN();

		// The last character read in parseQuals should have been a
		// '\n'

		readCnt_++;
		patid = readCnt_-1;
	}

	/**
	 * Dump a FASTQ-style record for the read.
	 */
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << "@" << name << endl << seq << endl
		    << "+" << endl << qual << endl;
	}
private:

	/**
	 * Parse a name from filebuf_ and store in r.  Assume that the next
	 * character obtained via filebuf_.get() is the first character of
	 * the sequence and the string stops at the next char upto (could
	 * be tab, newline, etc.).
	 */
	int parseName(ReadBuf& r, ReadBuf* r2, char upto = '\t') {
		// Read the name out of the first field
		int c = 0;
		int nameLen = 0;
		while(true) {
			if((c = filebuf_.get()) < 0) {
				return -1;
			}
			if(c == upto) {
				// Finished with first field
				break;
			}
			if(c == '\n' || c == '\r') {
				return -1;
			}
			if(r2 != NULL) (*r2).nameBuf[nameLen] = c;
			r.nameBuf[nameLen++] = c;
		}
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, nameLen);
		if(r2 != NULL) {
			_setBegin((*r2).name, (*r2).nameBuf);
			_setLength((*r2).name, nameLen);
		}
		// Set up a default name if one hasn't been set
		if(nameLen == 0) {
			itoa10(readCnt_, r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			nameLen = strlen(r.nameBuf);
			_setLength(r.name, nameLen);
			if(r2 != NULL) {
				itoa10(readCnt_, (*r2).nameBuf);
				_setBegin((*r2).name, (*r2).nameBuf);
				_setLength((*r2).name, nameLen);
			}
		}
		assert_gt(nameLen, 0);
		return nameLen;
	}

	/**
	 * Parse a single sequence from filebuf_ and store in r.  Assume
	 * that the next character obtained via filebuf_.get() is the first
	 * character of the sequence and the sequence stops at the next
	 * char upto (could be tab, newline, etc.).
	 */
	int parseSeq(ReadBuf& r, int& charsRead, char upto = '\t') {
		int begin = 0;
		int dstLen = 0;
		int c = filebuf_.get();
		assert(c != upto);
		while(c != upto) {
			// Note: can't have a comment in the middle of a sequence,
			// though a comment can end a sequence
			if(isalpha(c)) {
				if(begin++ >= this->trim5_) {
					if(dna4Cat[c] == 0) {
						if(forgiveInput_) {
							return -1;
						} else {
							assert(false);
						}
					}
					if(dstLen + 1 > 1024) {
						cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
							 << "reads and re-run Bowtie";
						exit(1);
					}
					r.patBufFw[dstLen] = charToDna5[c];
					dstLen++;
				}
				charsRead++;
			}
			if((c = filebuf_.get()) < 0) {
				return -1;
			}
		}
		dstLen -= this->trim3_;
		_setBegin (r.patFw,  (Dna5*)r.patBufFw);
		_setLength(r.patFw,  dstLen);
		return dstLen;
	}

	/**
	 * Parse a single quality string from filebuf_ and store in r.
	 * Assume that the next character obtained via filebuf_.get() is
	 * the first character of the quality string and the string stops
	 * at the next char upto (could be tab, newline, etc.).
	 */
	int parseQuals(ReadBuf& r, int charsRead, int dstLen, char& c2,
	               char upto = '\t', char upto2 = -1)
	{
		int qualsRead = 0;
		int c = 0;
		if (intQuals_) {
			char buf[4096];
			while (qualsRead < charsRead) {
				// TODO: fixme
				size_t rd = filebuf_.gets(buf, sizeof(buf));
				if(rd == 0) break;
				assert(NULL == strrchr(buf, '\n'));
				vector<string> s_quals;
				tokenize(string(buf), " ", s_quals);
				for (unsigned int j = 0; j < s_quals.size(); ++j) {
					int iQ = atoi(s_quals[j].c_str());
					char c = intToPhred33(iQ, solQuals_);
					if (qualsRead >= trim5_) {
						size_t off = qualsRead - trim5_;
						if(off + 1 > 1024) {
							cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
								 << "Please truncate reads and quality values and and re-run Bowtie";
							exit(1);
						}
						assert_geq(c, 33);
						assert_leq(c, 73);
						r.qualBufFw[off] = c;
					}
					++qualsRead;
				}
			} // done reading integer quality lines
			if (charsRead > qualsRead) {
				tooFewQualities(r.name);
				exit(1);
			}
		} else {
			// Non-integer qualities
			while((qualsRead < dstLen + this->trim5_) && c >= 0) {
				c = filebuf_.get();
				c2 = c;
				if (c == ' ') {
					wrongQualityFormat();
					exit(1);
				}
				if(c < 0) {
					// EOF occurred in the middle of a read - abort
					return -1;
				}
				if(!isspace(c) && c != upto && (upto2 == -1 || c != upto2)) {
					if (qualsRead >= trim5_) {
						size_t off = qualsRead - trim5_;
						if(off + 1 > 1024) {
							cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
								 << "Please truncate reads and quality values and and re-run Bowtie";
							exit(1);
						}
						c = charToPhred33(c, solQuals_, phred64Quals_);
						assert_geq(c, 33);
						r.qualBufFw[off] = c;
					}
					qualsRead++;
				} else {
					break;
				}
			}
			if(qualsRead != dstLen + this->trim5_) {
				if(forgiveInput_) {
					return -1;
				} else {
					assert(false);
				}
			}
		}
		_setBegin (r.qualFw, (char*)r.qualBufFw);
		_setLength(r.qualFw, dstLen);
		while(c != upto && (upto2 == -1 || c != upto2)) {
			c = filebuf_.get();
			c2 = c;
		}
		return qualsRead;
	}

	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	int policy_;
};

/**
 * Synchronized concrete pattern source for a list of FASTA files where
 * reads need to be extracted from long continuous sequences.
 */
class FastaContinuousPatternSource : public BufferedFilePatternSource {
public:
	FastaContinuousPatternSource(
			const vector<string>& infiles,
			size_t length,
			size_t freq,
	        bool useSpinlock = true,
	        const char *dumpfile = NULL,
	        uint32_t skip = 0,
	        uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, useSpinlock,
		                          false, dumpfile, 0, 0, skip),
		length_(length), freq_(freq),
		eat_(length_), bufCur_(0)
	{
		assert_lt(length_, (size_t)ReadBuf::BUF_SIZE);
	}

	virtual void reset() {
		BufferedFilePatternSource::reset();
	}
protected:
	/// Read another pattern from a FASTA input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		while(true) {
			int c = filebuf_.get();
			if(c < 0) {
				seqan::clear(r.patFw);
				return;
			} if(c == '>') {
				resetForNextFile();
			} else {
				int cat = dna4Cat[c];
				if(cat == 2) {
					eat_ = length_;
				} else if(cat == 0) {
					continue;
				} else {
					buf_[bufCur_++] = c;
					if(bufCur_ == 1024) bufCur_ = 0;
					if(eat_ > 0) eat_--;
					if(eat_ == 0) {
						for(size_t i = 0; i < length_; i++) {
							if(length_ + i <= bufCur_) {
								c = buf_[bufCur_ - length_ + i];
							} else {
								c = buf_[bufCur_ - length_ + i + 1024];
							}
							r.patBufFw [i] = charToDna5[c];
							r.qualBufFw[i] = 'I';
						}
						_setBegin (r.patFw,  (Dna5*)r.patBufFw);
						_setLength(r.patFw,  length_);
						_setBegin (r.qualFw, r.qualBufFw);
						_setLength(r.qualFw, length_);
						// Set up a default name if one hasn't been set
						itoa10(readCnt_, r.nameBuf);
						_setBegin(r.name, r.nameBuf);
						_setLength(r.name, strlen(r.nameBuf));
						readCnt_++;
						patid = readCnt_-1;
						break;
					}
				}
			}
		}
	}
	/// Shouldn't ever be here; it's not sensible to obtain read pairs
	// from a continuous input.
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		cerr << "In FastaContinuousPatternSource.readPair()" << endl;
		exit(1);
	}
	/**
	 * Reset state to be read for the next file.
	 */
	virtual void resetForNextFile() {
		eat_ = length_;
		bufCur_ = 0;
	}
private:
	size_t length_;     /// length of reads to generate
	size_t freq_;       /// frequency to sample reads
	int policy_;        /// policy for handling Ns

	size_t eat_;        /// number of characters we need to skip before
	                    /// we have flushed all of the ambiguous or
	                    /// non-existent characters out of our read
	                    /// window
	char buf_[1024];    /// read buffer
	size_t bufCur_;     /// buffer cursor; points to where we should
	                    /// insert the next character
};

/**
 * Read a FASTQ-format file.
 * See: http://maq.sourceforge.net/fastq.shtml
 */
class FastqPatternSource : public BufferedFilePatternSource {
public:
	FastqPatternSource(const vector<string>& infiles,
	                   bool randomizeQuals = false,
	                   bool useSpinlock = true,
	                   const char *dumpfile = NULL,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool __forgiveInput = false,
					   bool solexa_quals = false,
					   bool phred64Quals = false,
					   bool integer_quals = false,
					   uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, useSpinlock,
		                          __forgiveInput, dumpfile,
		                          trim3, trim5, skip),
		first_(true),
		solQuals_(solexa_quals),
		phred64Quals_(phred64Quals),
		intQuals_(integer_quals)
	{ }
	virtual void reset() {
		first_ = true;
		filebuf_.resetLastN();
		BufferedFilePatternSource::reset();
	}
protected:
	/**
	 * Scan to the next FASTQ record (starting with @) and return the first
	 * character of the record (which will always be @).  Since the quality
	 * line may start with @, we keep scanning until we've seen a line
	 * beginning with @ where the line two lines back began with +.
	 */
	static int skipToNextFastqRecord(FileBuf& in, bool sawPlus) {
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

	/// Read another pattern from a FASTQ input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		const int bufSz = ReadBuf::BUF_SIZE;
		int c;
		int dstLen = 0;
		int nameLen = 0;
		// Pick off the first at
		if(first_) {
			c = filebuf_.get();
			if(c != '@') {
				if(forgiveInput_) {
					c = FastqPatternSource::skipToNextFastqRecord(filebuf_, c == '+');
					if(c < 0) {
						filebuf_.resetLastN();
						seqan::clear(r.patFw);
						return;
					}
				} else {
					c = getOverNewline(filebuf_);
					if(c < 0) {
						filebuf_.resetLastN();
						seqan::clear(r.patFw);
						return;
					}
				}
			}
			if(c != '@') {
				cerr << "Error: reads file does not look like a FASTQ file" << endl;
				exit(1);
			}
			assert_eq('@', c);
			first_ = false;
		}

		// Read to the end of the id line, sticking everything after the '@'
		// into *name
		while(true) {
			c = filebuf_.get();
			if(c < 0) {
				seqan::clear(r.patFw);
				filebuf_.resetLastN();
				return;
			}
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					c = filebuf_.get();
					if(c < 0) {
						seqan::clear(r.patFw);
						filebuf_.resetLastN();
						return;
					}
				}
				break;
			}
			r.nameBuf[nameLen++] = c;
			if(nameLen > bufSz-2) {
				// Too many chars in read name; print friendly error message
				_setBegin(r.name, r.nameBuf);
				_setLength(r.name, nameLen);
				cerr << "FASTQ read name is too long; read names must be " << (bufSz-2) << " characters or fewer." << endl;
				cerr << "Beginning of bad read name: " << r.name << endl;
				exit(1);
			}
		}
		_setBegin(r.name, r.nameBuf);
		assert_leq(nameLen, bufSz-2);
		_setLength(r.name, nameLen);
		// c now holds the first character on the line after the
		// @name line

		// filebuf_ now points just past the first character of a
		// sequence line, and c holds the first character
		int charsRead = 0;
		while(c != '+') {
			if(isalpha(c)) {
				// If it's past the 5'-end trim point
				if(charsRead >= this->trim5_) {
					if(dstLen + 1 > 1024) {
						cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
							 << "reads and re-run Bowtie";
						exit(1);
					}
					r.patBufFw[dstLen] = charToDna5[c];
					dstLen++;
				}
				charsRead++;
			}
			c = filebuf_.get();
			if(c < 0) {
				// EOF occurred in the middle of a read - abort
				seqan::clear(r.patFw);
				filebuf_.resetLastN();
				return;
			}
		}
		// Trim from 3' end
		dstLen -= this->trim3_;
		// Set trimmed bounds of buffers
		_setBegin(r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, dstLen);
		assert_eq('+', c);

		// Chew up the optional name on the '+' line
		peekToEndOfLine(filebuf_);

		// Now read the qualities
		int qualsRead = 0;
		if (intQuals_) {
			char buf[4096];
			while (qualsRead < charsRead) {
				size_t rd = filebuf_.gets(buf, sizeof(buf));
				if(rd == 0) break;
				assert(NULL == strrchr(buf, '\n'));
				vector<string> s_quals;
				tokenize(string(buf), " ", s_quals);
				for (unsigned int j = 0; j < s_quals.size(); ++j)
				{
					int iQ = atoi(s_quals[j].c_str());
					char c = intToPhred33(iQ, solQuals_);
					if (qualsRead >= trim5_)
					{
						size_t off = qualsRead - trim5_;
						if(off + 1 > 1024) {
							cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
								 << "Please truncate reads and quality values and and re-run Bowtie";
							exit(1);
						}
						assert_geq(c, 33);
						assert_leq(c, 73);
						r.qualBufFw[off] = c;
					}
					++qualsRead;
				}
			} // done reading integer quality lines
			if (charsRead > qualsRead) {
				tooFewQualities(r.name);
				exit(1);
			}
			_setBegin(r.qualFw, (char*)r.qualBufFw);
			_setLength(r.qualFw, dstLen);
			peekOverNewline(filebuf_);
		} else {
			// Non-integer qualities
			while((qualsRead < dstLen + this->trim5_) && c >= 0) {
				c = filebuf_.get();
				if (c == ' ') {
					wrongQualityFormat();
					exit(1);
				}
				if(c < 0) {
					// EOF occurred in the middle of a read - abort
					seqan::clear(r.patFw);
					filebuf_.resetLastN();
					return;
				}
				if (c != '\r' && c != '\n') {
					if (qualsRead >= trim5_) {
						size_t off = qualsRead - trim5_;
						if(off + 1 > 1024) {
							cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
								 << "Please truncate reads and quality values and and re-run Bowtie";
							exit(1);
						}
						c = charToPhred33(c, solQuals_, phred64Quals_);
						assert_geq(c, 33);
						r.qualBufFw[off] = c;
					}
					qualsRead++;
				} else {
					break;
				}
			}
			assert_eq(qualsRead, dstLen + this->trim5_);
			_setBegin (r.qualFw, (char*)r.qualBufFw);
			_setLength(r.qualFw, dstLen);
			if(c == '\r' || c == '\n') {
				c = peekOverNewline(filebuf_);
			} else {
				c = peekToEndOfLine(filebuf_);
			}
		}
		r.readOrigBufLen = filebuf_.copyLastN(r.readOrigBuf);
		filebuf_.resetLastN();

		c = filebuf_.get();
		assert(c == -1 || c == '@');

		// Set up a default name if one hasn't been set
		if(nameLen == 0) {
			itoa10(readCnt_, r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			nameLen = strlen(r.nameBuf);
			_setLength(r.name, nameLen);
		}
		assert_gt(nameLen, 0);
		readCnt_++;
		patid = readCnt_-1;
	}
	/// Read another read pair from a FASTQ input file
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In FastqPatternSource.readPair()" << endl;
		exit(1);
		read(ra, patid);
		readCnt_--;
		read(rb, patid);
	}
	virtual void resetForNextFile() {
		first_ = true;
	}
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
	}
private:
	bool first_;
	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	int policy_;
};

/**
 * Read a Raw-format file (one sequence per line).  No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class RawPatternSource : public BufferedFilePatternSource {
public:
	RawPatternSource(const vector<string>& infiles,
	                 bool randomizeQuals = false,
	                 bool useSpinlock = true,
	                 const char *dumpfile = NULL,
	                 int trim3 = 0,
	                 int trim5 = 0,
	                 uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, false, useSpinlock,
		                          dumpfile, trim3, trim5, skip),
		first_(true)
	{ }
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}
protected:
	/// Read another pattern from a Raw input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		int c;
		int dstLen = 0;
		int nameLen = 0;
		c = peekOverNewline(this->filebuf_);
		if(c < 0) {
			seqan::clear(r.patFw);
			return;
		}
		assert(!isspace(c));
		if(first_) {
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
			first_ = false;
		}

		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		while(!isspace(c) && c >= 0) {
			c = filebuf_.get();
			if(isalpha(c) && dstLen >= this->trim5_) {
				size_t len = dstLen - this->trim5_;
				if(len + 1 > 1024) {
					cerr << "Reads file contained a pattern with more than 1024 characters." << endl
						 << "Please truncate reads and and re-run Bowtie";
					exit(1);
				}
				r.patBufFw [len] = charToDna5[c];
				r.qualBufFw[len] = 'I';
				dstLen++;
			} else if(isalpha(c)) dstLen++;
			c = filebuf_.peek();
		}
		if(dstLen >= (this->trim3_ + this->trim5_)) {
			dstLen -= (this->trim3_ + this->trim5_);
		} else {
			dstLen = 0;
		}
		_setBegin (r.patFw,  (Dna5*)r.patBufFw);
		_setLength(r.patFw,  dstLen);
		_setBegin (r.qualFw, r.qualBufFw);
		_setLength(r.qualFw, dstLen);

		c = peekToEndOfLine(filebuf_);
		r.readOrigBufLen = filebuf_.copyLastN(r.readOrigBuf);
		filebuf_.resetLastN();

		// Set up name
		itoa10(readCnt_, r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		nameLen = strlen(r.nameBuf);
		_setLength(r.name, nameLen);
		readCnt_++;

		patid = readCnt_-1;
	}
	/// Read another read pair from a FASTQ input file
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In RawPatternSource.readPair()" << endl;
		exit(1);
		read(ra, patid);
		readCnt_--;
		read(rb, patid);
	}
	virtual void resetForNextFile() {
		first_ = true;
	}
	virtual void dump(ostream& out,
	                  const String<Dna5>& seq,
	                  const String<char>& qual,
	                  const String<char>& name)
	{
		out << seq << endl;
	}
private:
	bool first_;
	int policy_;
};

#endif /*PAT_H_*/
