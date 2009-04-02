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
	}

	void reverseAll() {
		::reverseInPlace(patFw);
		::reverseInPlace(patRc);
		::reverseInPlace(qualFw);
		::reverseInPlace(qualRc);
		::reverseInPlace(patFwRev);
		::reverseInPlace(patRcRev);
		::reverseInPlace(qualFwRev);
		::reverseInPlace(qualRcRev);
	}

	/// Return true iff the read (pair) is empty
	bool empty() {
		return seqan::empty(patFw);
	}

	/// Return length of the read in the buffer
	uint32_t length() {
		return seqan::length(patFw);
	}

	void constructReverse() {
		uint32_t len = length();

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
		if(namelen < 2) append = true;
		else {
			if(i == 1) {
				append =
					nameBuf[namelen-2] != '/' ||
					nameBuf[namelen-1] != '1';
			} else {
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

	bool          reversed;

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
	PatternSource(bool reverse = false,
	              bool randomizeQuals = false,
	              bool useSpinlock = true,
	              const char *dumpfile = NULL) :
	    readCnt_(0),
	    reverse_(reverse),
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
	virtual void nextRead(ReadBuf& r, uint32_t& patid) {
		// nextPatternImpl does the reading from the ultimate source;
		// it is implemented in concrete subclasses
		nextReadImpl(r, patid);
		if(!r.empty()) {
			// Possibly randomize the qualities so that they're more
			// scattered throughout the range of possible values
			if(randomizeQuals_) {
				const size_t rlen = r.length();
				for(size_t i = 0; i < rlen; i++) {
					if(i < rlen-1) {
						r.qualFw[i] *= (r.qualFw[i+1] + 7);
					}
					if(i > 0) {
						r.qualFw[i] *= (r.qualFw[i-1] + 11);
					}
					if(r.qualFw[i] < 0) r.qualFw[i] = -(r.qualFw[i]+1);
					r.qualFw[i] %= 41;
					assert_leq(r.qualFw[i], 40);
					r.qualFw[i] += 33;
					r.qualRc[rlen-i-1] = r.qualFw[i];
				}
			}
			// Construct the reversed versions of the fw and rc seqs
			// and quals
			r.constructReverse();
			r.seed = genRandSeed(r.patFw, r.qualFw, r.name);
			// If it's this class's responsibility to reverse the pattern,
			// do so here.  Usually it's the responsibility of one of the
			// concrete subclasses, since they can usually do it more
			// efficiently.
			if(reverse_) {
				::reverseInPlace(r.patFw);
				::reverseInPlace(r.patRc);
				::reverseInPlace(r.qualFw);
				::reverseInPlace(r.qualRc);
			}
		}
		// Output it, if desired
		if(dumpfile_ != NULL) {
			dump(out_, r.patFw,
			     empty(r.qualFw) ? String<char>("(empty)") : r.qualFw,
			     empty(r.name)   ? String<char>("(empty)") : r.name);
			dump(out_, r.patRc,
			     empty(r.qualRc) ? String<char>("(empty)") : r.qualRc,
			     empty(r.name)   ? String<char>("(empty)") : r.name);
		}
	}
	/// Implementation to be provided by concrete subclasses
	virtual void nextReadImpl(ReadBuf& r, uint32_t& patid) = 0;
	/// Reset state to start over again with the first read
	virtual void reset() { readCnt_ = 0; }
	/**
	 * Whether to reverse reads as they're read in (useful when using
	 * the mirror index)
	 */
	virtual bool reverse() const { return reverse_; }
	/**
	 * Set whether to reverse reads as they're read in (useful when
	 * using the mirror index.
	 */
	virtual void setReverse(bool reverse) {
		reverse_ = reverse;
	}

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

	bool reverse_;         /// reverse patterns before returning them
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
 * Encapsulates a synchronized source of both paired-end patterns and
 * unpaired patterns.
 */
class PairedPatternSource {

public:

	PairedPatternSource(const vector<PatternSource*>& srca,
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
		MUTEX_INIT(lock_);
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
	 * The main member function for dispensing pairs of reads or
	 * singleton reads.  Returns true iff ra and rb contain a new
	 * pair; returns false if ra contains a new unpaired read.
	 */
	bool nextReadPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
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
	volatile uint32_t cur_; // current element in parallel srca_, srcb_ vectors
	vector<PatternSource*> srca_; /// PatternSources for 1st mates and/or unpaired reads
	vector<PatternSource*> srcb_; /// PatternSources for 2nd mates
#ifdef USE_SPINLOCK
	SpinLock spinlock_;
#endif
	MUTEX_T lock_; /// mutex for locking critical regions
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
		buf1_(), buf2_(), patid_(0xffffffff), reverse_(false) { }

	virtual ~PatternSourcePerThread() { }

	/**
	 * Read the next read pair.  Reset reverse_ to indicate that this
	 * object hasn't yet been used to reverse the sequence and quals.
	 */
	virtual void nextReadPair() { reverse_ = false; }

	ReadBuf& bufa()        { return buf1_;         }
	ReadBuf& bufb()        { return buf2_;         }

	uint32_t      patid()  { return patid_;        }
	virtual void  reset()  { patid_ = 0xffffffff;  }
	void    reverseRead()  { buf1_.reverseAll();
	                         buf2_.reverseAll();
	                         reverse_ = !reverse_;     }
	bool          empty()  { return buf1_.empty();     }

	/**
	 * Return true iff the buffers jointly contain a paired-end read.
	 */
	bool paired() {
		bool ret = !buf2_.empty();
	    assert(!ret || !empty());
	    return ret;
	}

	/**
	 * Return true iff the reads in the buffers bufa and bufb are
	 * reversed from their original representation in the input reads.
	 */
	virtual bool reverse() { return reverse_; }

protected:
	ReadBuf  buf1_;    // read buffer for mate a
	ReadBuf  buf2_;    // read buffer for mate b
	uint32_t patid_;   // index of read just read
	bool     reverse_; // whether the read is reversed
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

	/**
	 * Return true iff the reads in the buffers bufa and bufb are
	 * reversed from their original representation in the input reads.
	 */
	virtual bool reverse() { return reverse_ != patsrc_.reverse(); }

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
	TrimmingPatternSource(bool reverse = false,
	                      bool randomizeQuals = false,
	                      bool useSpinlock = true,
	                      const char *dumpfile = NULL,
	                      int trim3 = 0,
	                      int trim5 = 0) :
		PatternSource(reverse, randomizeQuals, useSpinlock, dumpfile),
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
		PatternSource(false, false, useSpinlock, dumpfile),
		numReads_(numReads),
		length_(length),
		seed_(seed),
		rand_(seed),
		reverse_(false)
	{
		if(length_ > 1024) {
			cerr << "Read length for RandomPatternSource may not exceed 1024; got " << length_ << endl;
			exit(1);
		}
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
		fillRandomRead(r, ra, length_, patid, reverse_);
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
		rand_.init(seed_);
	}
	virtual bool reverse() const { return reverse_; }
	virtual void setReverse(bool reverse) {
		reverse_ = reverse;
	}
private:
	uint32_t     numReads_; /// number of reads to dish out
	int          length_;   /// length of reads
	uint32_t     seed_;     /// seed for pseudo-randoms
	RandomSource rand_;     /// pseudo-random generator
	bool         reverse_;  /// whether to reverse reads
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
	                             int thread,
	                             bool reverse) :
		PatternSourcePerThread(),
		numreads_(numreads),
		length_(length),
		numthreads_(numthreads),
		thread_(thread),
		reverse_(reverse),
		rand_(thread_)
	{
		patid_ = thread_;
		if(length_ > 1024) {
			cerr << "Read length for RandomPatternSourcePerThread may not exceed 1024; got " << length_ << endl;
			exit(1);
		}
	}

	virtual void nextReadPair() {
		PatternSourcePerThread::nextReadPair();
		if(patid_ >= numreads_) {
			buf1_.clearAll();
			buf2_.clearAll();
			return;
		}
		RandomPatternSource::fillRandomRead(
			buf1_, rand_.nextU32(), length_, patid_, reverse_);
		RandomPatternSource::fillRandomRead(
			buf2_, rand_.nextU32(), length_, patid_, reverse_);
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
	bool         reverse_;
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
            int thread,
            bool reverse) :
            numreads_(numreads),
            length_(length),
            numthreads_(numthreads),
            thread_(thread),
            reverse_(reverse) { }

	/**
	 * Create a new heap-allocated WrappedPatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new RandomPatternSourcePerThread(
			numreads_, length_, numthreads_, thread_, reverse_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * WrappedPatternSourcePerThreads.
	 */
	virtual std::vector<PatternSourcePerThread*>* create(uint32_t n) const {
		std::vector<PatternSourcePerThread*>* v = new std::vector<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new RandomPatternSourcePerThread(
				numreads_, length_, numthreads_, thread_, reverse_));
			assert(v->back() != NULL);
		}
		return v;
	}

private:
	uint32_t numreads_;
    int length_;
    int numthreads_;
    int thread_;
    bool reverse_;
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
	                    bool reverse = false,
	                    bool randomizeQuals = false,
	                    bool useSpinlock = true,
	                    const char *dumpfile = NULL,
	                    int trim3 = 0,
	                    int trim5 = 0,
		                int npolicy = NS_TO_NS,
		                uint32_t skip = 0,
	                    uint32_t seed = 0) :
		TrimmingPatternSource(false, randomizeQuals, useSpinlock, dumpfile, trim3, trim5),
		reverse_(reverse), cur_(skip), skip_(skip),
		v_(), vrev_(), vrc_(), vrcrev_(), quals_(), qualsrev_(), rand_(seed)
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
			// Count Ns and possibly reject
			for(size_t j = 0; j < s.length(); j++) {
				if(s[j] == 'N' || s[j] == 'n') {
					if(npolicy == NS_TO_NS) {
						// Leave s[j] == 'N'
					} else if(npolicy == NS_TO_RANDS) {
						s[j] = "ACGT"[rand_.nextU32() & 3];
					} else {
						assert_eq(NS_TO_AS, npolicy);
						s[j] = 'A';
					}
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
			{
				vrev_.push_back(s);
				::reverseInPlace(vrev_.back());
				qualsrev_.push_back(vq);
				::reverseInPlace(qualsrev_.back());
			}
			vrc_.push_back(reverseComplement(String<Dna5>(s)));
			{
				vrcrev_.push_back(reverseComplement(String<Dna5>(s)));
				::reverseInPlace(vrcrev_.back());
			}
			ostringstream os;
			os << (names_.size());
			names_.push_back(os.str());
		}
		assert_eq(v_.size(), vrev_.size());
		assert_eq(v_.size(), quals_.size());
		assert_eq(v_.size(), qualsrev_.size());
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
			clear(r.patFw);
			clear(r.patRc);
			clear(r.qualFw);
			clear(r.qualRc);
			clear(r.name);
			return;
		}
		// Copy v_*, quals_* strings into the respective Strings
		if(!reverse_) {
			// not reversed
			r.patFw  = v_[cur_];
			r.patRc  = vrc_[cur_];
			r.qualFw = quals_[cur_];
			r.qualRc = qualsrev_[cur_];
		} else {
			// reversed
			r.patFw  = vrev_[cur_];
			r.patRc  = vrcrev_[cur_];
			r.qualFw = qualsrev_[cur_];
			r.qualRc = quals_[cur_];
		}
		ostringstream os;
		os << cur_;
		r.name = os.str();
		cur_++;
		readCnt_++;
		patid = readCnt_;
		unlock();
	}
	virtual void reset() {
		TrimmingPatternSource::reset();
		cur_ = skip_;
	}
	virtual bool reverse() const { return reverse_; }
	virtual void setReverse(bool reverse) {
		reverse_ = reverse;
	}
private:
	bool   reverse_;
	size_t cur_;
	uint32_t skip_;
	vector<String<Dna5> > v_;        /// forward sequences
	vector<String<Dna5> > vrev_;     /// reversed forward sequences
	vector<String<Dna5> > vrc_;      /// rev-comp sequences
	vector<String<Dna5> > vrcrev_;   /// reversed rev-comp sequences
	vector<String<char> > quals_;    /// quality values parallel to v_
	vector<String<char> > qualsrev_; /// quality values parallel to vrev_
	vector<String<char> > names_;    /// names
	RandomSource rand_;
};

/**
 * Supports reversing all strings as they're read in.  Also supports
 * returning reverse complements interspersed with forward versions of
 * patterns.
 */
class BufferedFilePatternSource : public TrimmingPatternSource {
public:
	BufferedFilePatternSource(const vector<string>& infiles,
	                          bool reverse = false,
	                          bool randomizeQuals = false,
	                          bool useSpinlock = true,
	                          bool __forgiveInput = false,
	                          const char *dumpfile = NULL,
	                          int trim3 = 0,
	                          int trim5 = 0,
	                          uint32_t skip = 0) :
		TrimmingPatternSource(reverse, randomizeQuals, useSpinlock, dumpfile, trim3, trim5),
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
	virtual void read(ReadBuf& r, uint32_t& patid) = 0; // length of name
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() = 0;
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
	                   bool reverse = false,
	                   bool randomizeQuals = false,
	                   bool useSpinlock = true,
	                   const char *dumpfile = NULL,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   int npolicy = NS_TO_NS,
	                   bool __forgiveInput = false,
	                   uint32_t skip = 0,
	                   uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, randomizeQuals, useSpinlock,
		                          __forgiveInput, dumpfile,
		                          trim3, trim5, skip),
		first_(true), reverse_(reverse), policy_(npolicy), rand_(seed)
	{ }
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}
	virtual bool reverse() const { return reverse_; }
	virtual void setReverse(bool reverse) {
		reverse_ = reverse;
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
		const int bufSz = ReadBuf::BUF_SIZE;
		int c;
		int dstLen = 0;
		int nameLen = 0;
		// Pick off the first carat
		if(first_) {
			c = filebuf_.get();
			if(c != '>') {
				if(forgiveInput_) {
					c = FastaPatternSource::skipToNextFastaRecord(filebuf_);
					if(c < 0) {
						r.clearAll(); return;
					}
				} else {
					c = getOverNewline(filebuf_); if(c < 0) {
						r.clearAll(); return;
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
				r.clearAll(); return;
			}
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					c = filebuf_.get(); if(c < 0) {
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
		if(!reverse_) {
			while(c != '>' && c != '#') {
				// Note: can't have a comment in the middle of a sequence,
				// though a comment can end a sequence
				if(isalpha(c) && begin++ >= this->trim5_) {
					if(dstLen + 1 > 1024) {
						cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
							 << "reads and re-run Bowtie";
						exit(1);
					}
					if(c == 'N' || c == 'n') {
						if(policy_ == NS_TO_NS) {
							// Leave c = 'N'
						} else if(policy_ == NS_TO_RANDS) {
							c = "ACGT"[rand_.nextU32() & 3];
						} else {
							assert_eq(NS_TO_AS, policy_);
							c = 'A';
						}
					}
					r.patBufFw [dstLen] = charToDna5[c];
					r.qualBufFw[dstLen] = 'I';
					r.patBufRc [bufSz-dstLen-1] = rcCharToDna5[c];
					r.qualBufRc[bufSz-dstLen-1] = 'I';
					dstLen++;
				}
				if((c = filebuf_.get()) < 0) break;
			}
			dstLen -= this->trim3_;
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
				if(isalpha(c) && begin++ >= this->trim5_) {
					if(dstLen + 1 > 1024) {
						cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
							 << "reads and re-run Bowtie";
						exit(1);
					}
					if(c == 'N' || c == 'n') {
						if(policy_ == NS_TO_NS) {
							// Leave c = 'N'
						} else if(policy_ == NS_TO_RANDS) {
							c = "ACGT"[rand_.nextU32() & 3];
						} else {
							assert_eq(NS_TO_AS, policy_);
							c = 'A';
						}
					}
					r.patBufFw [bufSz-dstLen-1] = charToDna5[c];
					r.qualBufFw[bufSz-dstLen-1] = 'I';
					r.patBufRc [dstLen] = rcCharToDna5[c];
					r.qualBufRc[dstLen] = 'I';
					dstLen++;
				}
				if((c = filebuf_.get()) < 0) break;
			}
			dstLen -= this->trim3_;
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
			itoa10(readCnt_, r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			nameLen = strlen(r.nameBuf);
			_setLength(r.name, nameLen);
		}
		assert_gt(nameLen, 0);
		readCnt_++;
		patid = readCnt_-1;
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
	bool reverse_;
	int policy_;
	RandomSource rand_;
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
	        bool reverse = false,
	        bool useSpinlock = true,
	        const char *dumpfile = NULL,
	        int npolicy = NS_TO_NS,
	        uint32_t skip = 0,
	        uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, false, useSpinlock,
		                          false, dumpfile, 0, 0, skip),
		length_(length), freq_(freq), reverse_(reverse),
		policy_(npolicy), eat_(length_), bufCur_(0), rand_(seed)
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
							r.patBufRc [length_-i-1] = rcCharToDna5[c];
							r.qualBufRc[length_-i-1] = 'I';
						}
						_setBegin (r.patFw,  (Dna5*)r.patBufFw);
						_setLength(r.patFw,  length_);
						_setBegin (r.qualFw, r.qualBufFw);
						_setLength(r.qualFw, length_);
						_setBegin (r.patRc,  (Dna5*)r.patBufRc);
						_setLength(r.patRc,  length_);
						_setBegin (r.qualRc, r.qualBufRc);
						_setLength(r.qualRc, length_);
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
	bool reverse_;      /// reverse reads on read-in
	int policy_;        /// policy for handling Ns

	size_t eat_;        /// number of characters we need to skip before
	                    /// we have flushed all of the ambiguous or
	                    /// non-existent characters out of our read
	                    /// window
	char buf_[1024];    /// read buffer
	size_t bufCur_;     /// buffer cursor; points to where we should
	                    /// insert the next character

	RandomSource rand_;
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
	                   bool reverse = false,
	                   bool randomizeQuals = false,
	                   bool useSpinlock = true,
	                   const char *dumpfile = NULL,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   int npolicy = NS_TO_NS,
	                   bool __forgiveInput = false,
					   bool solexa_quals = false,
					   bool integer_quals = true,
					   uint32_t skip = 0,
	                   uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, randomizeQuals, useSpinlock,
		                          __forgiveInput, dumpfile,
		                          trim3, trim5, skip),
		first_(true), reverse_(reverse),
		solQuals_(solexa_quals),
		intQuals_(integer_quals),
		policy_(npolicy),
		rand_(seed)
	{
		for (int l = 0; l != 128; ++l) {
			table_[l] = (int)(10.0 * log(1.0 + pow(10.0, (l - 64) / 10.0)) / log(10.0) + .499);
			if (table_[l] >= 63) table_[l] = 63;
			if (table_[l] == 0) table_[l] = 1;
		}
	}
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}
	virtual bool reverse() const { return reverse_; }
	virtual void setReverse(bool reverse) {
		reverse_ = reverse;
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
					c = FastqPatternSource::skipToNextFastqRecord(filebuf_, c == '+'); if(c < 0) {
						seqan::clear(r.patFw);
						return;
					}
				} else {
					c = getOverNewline(filebuf_); if(c < 0) {
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
			c = filebuf_.get(); if(c < 0) {
				seqan::clear(r.patFw);
				return;
			}
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					c = filebuf_.get();
					if(c < 0) {
						seqan::clear(r.patFw);
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
		if(!reverse_) {
			while(c != '+') {
				if(isalpha(c)) {
					// If it's past the 5'-end trim point
					if(charsRead >= this->trim5_) {
						if(dstLen + 1 > 1024) {
							cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
								 << "reads and re-run Bowtie";
							exit(1);
						}
						// Add it to the read buffer
						if(c == 'N' || c == 'n') {
							if(policy_ == NS_TO_NS) {
								// Leave c = 'N'
							} else if(policy_ == NS_TO_RANDS) {
								c = "ACGT"[rand_.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, policy_);
								c = 'A';
							}
						}
						r.patBufFw[dstLen] = charToDna5[c];
						r.patBufRc[bufSz-dstLen-1] = rcCharToDna5[c];
						dstLen++;
					}
					charsRead++;
				}
				c = filebuf_.get();
				if(c < 0) {
					// EOF occurred in the middle of a read - abort
					seqan::clear(r.patFw);
					return;
				}
			}
			// Trim from 3' end
			dstLen -= this->trim3_;
			// Set trimmed bounds of buffers
			_setBegin(r.patFw, (Dna5*)r.patBufFw);
			_setLength(r.patFw, dstLen);
			_setBegin(r.patRc, (Dna5*)&r.patBufRc[bufSz-dstLen]);
			_setLength(r.patRc, dstLen);
		} else {
			while(c != '+') {
				if(isalpha(c)) {
					// If it's past the 5'-end trim point
					if(charsRead >= this->trim5_) {
						if(dstLen + 1 > 1024) {
							cerr << "Input file contained a pattern more than 1024 characters long.  Please truncate" << endl
								 << "reads and re-run Bowtie";
							exit(1);
						}
						// Add it to the read buffer
						if(c == 'N' || c == 'n') {
							if(policy_ == NS_TO_NS) {
								// Leave c = 'N'
							} else if(policy_ == NS_TO_RANDS) {
								c = "ACGT"[rand_.nextU32() & 3];
							} else {
								assert_eq(NS_TO_AS, policy_);
								c = 'A';
							}
						}
						r.patBufFw[bufSz-dstLen-1] = charToDna5[c];
						r.patBufRc[dstLen] = rcCharToDna5[c];
						dstLen++;
					}
					charsRead++;
				}
				c = filebuf_.get();
				if(c < 0) {
					// EOF occurred in the middle of a read - abort
					seqan::clear(r.patFw);
					return;
				}
			}
			// Trim from 3' end
			dstLen -= this->trim3_;
			// Set trimmed bounds of buffers
			_setBegin(r.patFw, (Dna5*)&r.patBufFw[bufSz-dstLen]);
			_setLength(r.patFw, dstLen);
			_setBegin(r.patRc, (Dna5*)r.patBufRc);
			_setLength(r.patRc, dstLen);
		}
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
				if(!reverse_) {
					for (unsigned int j = 0; j < s_quals.size(); ++j)
					{
						int iQ = atoi(s_quals[j].c_str());
						int pQ;
						if (solQuals_)
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

						if (qualsRead >= trim5_)
						{
							size_t off = qualsRead - trim5_;
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
						if (solQuals_)
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
						if (qualsRead >= trim5_)
						{
							size_t off = qualsRead - trim5_;
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
			if(!reverse_) {
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
			c = filebuf_.get();
		}
		else
		{
			// Non-integer qualities
			if(!reverse_) {
				while((qualsRead < dstLen + this->trim5_) && c >= 0) {
					c = filebuf_.get();
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
						if (qualsRead >= trim5_) {
							size_t off = qualsRead - trim5_;
							if(off + 1 > 1024) {
								cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
									 << "Please truncate reads and quality values and and re-run Bowtie";
								exit(1);
							}

							if (solQuals_)
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
				assert_eq(qualsRead, dstLen + this->trim5_);
				_setBegin (r.qualFw, (char*)r.qualBufFw);
				_setLength(r.qualFw, dstLen);
				_setBegin (r.qualRc, (char*)&r.qualBufRc[bufSz-dstLen]);
				_setLength(r.qualRc, dstLen);
			} else {
				while((qualsRead < dstLen + this->trim5_) && c >= 0) {
					c = filebuf_.get();
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
						if (qualsRead >= trim5_) {
							size_t off = qualsRead - trim5_;
							if(off + 1 > 1024) {
								cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
									 << "Please truncate reads and quality values and and re-run Bowtie";
								exit(1);
							}

							if (solQuals_)
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
				assert_eq(qualsRead, dstLen + this->trim5_);
				_setBegin (r.qualFw, (char*)&r.qualBufFw[bufSz-dstLen]);
				_setLength(r.qualFw, dstLen);
				_setBegin (r.qualRc, (char*)r.qualBufRc);
				_setLength(r.qualRc, dstLen);
			}
			if(c == '\r' || c == '\n') {
				c = getOverNewline(filebuf_);
			} else {
				c = getToEndOfLine(filebuf_);
			}
		}
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
	bool reverse_;
	bool solQuals_;
	bool intQuals_;
	int policy_;
	int table_[128];
	RandomSource rand_;
};

/**
 * Read a Raw-format file (one sequence per line).
 */
class RawPatternSource : public BufferedFilePatternSource {
public:
	RawPatternSource(const vector<string>& infiles,
	                 bool reverse = false,
	                 bool randomizeQuals = false,
	                 bool useSpinlock = true,
	                 const char *dumpfile = NULL,
	                 int trim3 = 0,
	                 int trim5 = 0,
	                 int npolicy = NS_TO_NS,
	                 uint32_t skip = 0,
	                 uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, randomizeQuals, false, useSpinlock,
		                          dumpfile, trim3, trim5, skip),
		first_(true), reverse_(reverse), policy_(npolicy), rand_(seed)
	{ }
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}
	virtual bool reverse() const { return reverse_; }
	virtual void setReverse(bool reverse) {
		reverse_ = reverse;
	}
protected:
	/// Read another pattern from a Raw input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		const int bufSz = ReadBuf::BUF_SIZE;
		int c;
		int dstLen = 0;
		int nameLen = 0;
		c = getOverNewline(this->filebuf_);
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
		if(!reverse_) {
			while(!isspace(c) && c >= 0) {
				if(isalpha(c) && dstLen >= this->trim5_) {
					size_t len = dstLen - this->trim5_;
					if(len + 1 > 1024) {
						cerr << "Reads file contained a pattern with more than 1024 characters." << endl
							 << "Please truncate reads and and re-run Bowtie";
						exit(1);
					}
					if(c == 'N' || c == 'n') {
						if(policy_ == NS_TO_NS) {
							// Leave c = 'N'
						} else if(policy_ == NS_TO_RANDS) {
							c = "ACGT"[rand_.nextU32() & 3];
						} else {
							assert_eq(NS_TO_AS, policy_);
							c = 'A';
						}
					}
					r.patBufFw [len] = charToDna5[c];
					r.qualBufFw[len] = 'I';
					r.patBufRc [bufSz-len-1] = rcCharToDna5[c];
					r.qualBufRc[bufSz-len-1] = 'I';
					dstLen++;
				} else if(isalpha(c)) dstLen++;
				c = filebuf_.get();
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
			_setBegin (r.patRc,  (Dna5*)&r.patBufRc[bufSz-dstLen]);
			_setLength(r.patRc,  dstLen);
			_setBegin (r.qualRc, &r.qualBufRc[bufSz-dstLen]);
			_setLength(r.qualRc, dstLen);
		} else {
			while(!isspace(c) && c >= 0) {
				if(isalpha(c) && dstLen >= this->trim5_) {
					size_t len = dstLen - this->trim5_;
					if(len + 1 > 1024) {
						cerr << "Reads file contained a pattern with more than 1024 characters.." << endl
							 << "Please truncate reads and and re-run Bowtie";
						exit(1);
					}
					if(c == 'N' || c == 'n') {
						if(policy_ == NS_TO_NS) {
							// Leave c = 'N'
						} else if(policy_ == NS_TO_RANDS) {
							c = "ACGT"[rand_.nextU32() & 3];
						} else {
							assert_eq(NS_TO_AS, policy_);
							c = 'A';
						}
					}
					r.patBufFw [bufSz-len-1] = charToDna5[c];
					r.qualBufFw[bufSz-len-1] = 'I';
					r.patBufRc [len] = rcCharToDna5[c];
					r.qualBufRc[len] = 'I';
					dstLen++;
				} else if(isalpha(c)) dstLen++;
				c = filebuf_.get();
			}
			if(dstLen >= (this->trim3_ + this->trim5_)) {
				dstLen -= (this->trim3_ + this->trim5_);
			} else {
				dstLen = 0;
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
		itoa10(readCnt_, r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		nameLen = strlen(r.nameBuf);
		_setLength(r.name, nameLen);
		readCnt_++;

		if(c == -1) {
			return;
		} else {
			assert(isspace(c));
		}
		patid = readCnt_-1;
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
	bool reverse_;
	int policy_;
	RandomSource rand_;
};

#endif /*PAT_H_*/
