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
#include "hit_set.h"

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
	ReadBuf() { reset(); }

	~ReadBuf() {
		clearAll(); reset();
		// Prevent seqan from trying to free buffers
		_setBegin(patFw, NULL);
		_setBegin(patRc, NULL);
		_setBegin(qual, NULL);
		_setBegin(patFwRev, NULL);
		_setBegin(patRcRev, NULL);
		_setBegin(qualRev, NULL);
		_setBegin(name, NULL);
		for(int j = 0; j < 3; j++) {
			_setBegin(altPatFw[j], NULL);
			_setBegin(altPatFwRev[j], NULL);
			_setBegin(altPatRc[j], NULL);
			_setBegin(altPatRcRev[j], NULL);
			_setBegin(altQual[j], NULL);
			_setBegin(altQualRev[j], NULL);
		}
	}

#define RESET_BUF(str, buf, typ) _setBegin(str, (typ*)buf); _setLength(str, 0); _setCapacity(str, BUF_SIZE);
#define RESET_BUF_LEN(str, buf, len, typ) _setBegin(str, (typ*)buf); _setLength(str, len); _setCapacity(str, BUF_SIZE);

	/// Point all Strings to the beginning of their respective buffers
	/// and set all lengths to 0
	void reset() {
		patid = 0;
		readOrigBufLen = 0;
		alts = 0;
		fuzzy = false;
		RESET_BUF(patFw, patBufFw, Dna5);
		RESET_BUF(patRc, patBufRc, Dna5);
		RESET_BUF(qual, qualBuf, char);
		RESET_BUF(patFwRev, patBufFwRev, Dna5);
		RESET_BUF(patRcRev, patBufRcRev, Dna5);
		RESET_BUF(qualRev, qualBufRev, char);
		RESET_BUF(name, nameBuf, char);
		for(int j = 0; j < 3; j++) {
			RESET_BUF(altPatFw[j], altPatBufFw[j], Dna5);
			RESET_BUF(altPatFwRev[j], altPatBufFwRev[j], Dna5);
			RESET_BUF(altPatRc[j], altPatBufRc[j], Dna5);
			RESET_BUF(altPatRcRev[j], altPatBufRcRev[j], Dna5);
			RESET_BUF(altQual[j], altQualBuf[j], char);
			RESET_BUF(altQualRev[j], altQualBufRev[j], char);
		}
	}

	void clearAll() {
		seqan::clear(patFw);
		seqan::clear(patRc);
		seqan::clear(qual);
		seqan::clear(patFwRev);
		seqan::clear(patRcRev);
		seqan::clear(qualRev);
		seqan::clear(name);
		for(int j = 0; j < 3; j++) {
			seqan::clear(altPatFw[j]);
			seqan::clear(altPatFwRev[j]);
			seqan::clear(altPatRc[j]);
			seqan::clear(altPatRcRev[j]);
			seqan::clear(altQual[j]);
			seqan::clear(altQualRev[j]);
		}
		readOrigBufLen = 0;
	}

	/// Return true iff the read (pair) is empty
	bool empty() const {
		return seqan::empty(patFw);
	}

	/// Return length of the read in the buffer
	uint32_t length() const {
		return seqan::length(patFw);
	}

	/**
	 * Construct patRc in place.
	 */
	void constructRevComps() {
		uint32_t len = length();
		assert_gt(len, 0);
		RESET_BUF_LEN(patRc, patBufRc, len, Dna5);
		for(int j = 0; j < alts; j++) {
			RESET_BUF_LEN(altPatRc[j], altPatBufRc[j], len, Dna5);
		}
		for(uint32_t i = 0; i < len; i++) {
			// Reverse-complement the sequence
			patBufRc[i]  = (patBufFw[len-i-1] == 4) ? 4 : (patBufFw[len-i-1] ^ 3);
			for(int j = 0; j < alts; j++) {
				altPatBufRc[j][i] = (altPatBufFw[j][len-i-1] == 4) ? 4 : (altPatBufFw[j][len-i-1] ^ 3);
			}
		}
	}

	/**
	 * Given patFw, patRc, and qual, construct the *Rev versions in
	 * place.  Assumes constructRevComps() was called previously.
	 */
	void constructReverses() {
		uint32_t len = length();
		assert_gt(len, 0);
		RESET_BUF_LEN(patFwRev, patBufFwRev, len, Dna5);
		RESET_BUF_LEN(patRcRev, patBufRcRev, len, Dna5);
		RESET_BUF_LEN(qualRev, qualBufRev, len, char);
		for(int j = 0; j < alts; j++) {
			RESET_BUF_LEN(altPatFwRev[j], altPatBufFwRev[j], len, Dna5);
			RESET_BUF_LEN(altPatRcRev[j], altPatBufRcRev[j], len, Dna5);
			RESET_BUF_LEN(altQualRev[j], altQualBufRev[j], len, char);
		}
		for(uint32_t i = 0; i < len; i++) {
			patFwRev[i]  = patFw[len-i-1];
			patRcRev[i]  = patRc[len-i-1];
			qualRev[i]   = qual[len-i-1];
			for(int j = 0; j < alts; j++) {
				altPatFwRev[j][i] = altPatFw[j][len-i-1];
				altPatRcRev[j][i] = altPatRc[j][len-i-1];
				altQualRev[j][i]  = altQual[j][len-i-1];
			}
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

	void dump(std::ostream& os) const {
		os << name << " " << patFw << " ";
		// Print out the sequences
		for(int j = 0; j < 3; j++) {
			bool started = false;
			if(seqan::length(altQual[j]) > 0) {
				for(size_t i = 0; i < length(); i++) {
					if(altQual[j][i] != '!') {
						started = true;
					}
					if(started) {
						if(altQual[j][i] == '!') {
							os << '-';
						} else {
							os << altPatFw[j][i];
						}
					}
				}
			}
			cout << " ";
		}
		os << qual << " ";
		// Print out the quality strings
		for(int j = 0; j < 3; j++) {
			bool started = false;
			if(seqan::length(altQual[j]) > 0) {
				for(size_t i = 0; i < length(); i++) {
					if(altQual[j][i] != '!') {
						started = true;
					}
					if(started) {
						os << altQual[j][i];
					}
				}
			}
			if(j == 2) {
				os << endl;
			} else {
				os << " ";
			}
		}
	}

	/**
	 * Write read details to a HitSet object.
	 */
	void toHitSet(HitSet& hs) {
		assert(!empty());
		hs.name = name;
		hs.seq = patFw;
		hs.qual = qual;
	}

	static const int BUF_SIZE = 1024;

	String<Dna5>  patFw;               // forward-strand sequence
	uint8_t       patBufFw[BUF_SIZE];  // forward-strand sequence buffer
	String<Dna5>  patRc;               // reverse-complement sequence
	uint8_t       patBufRc[BUF_SIZE];  // reverse-complement sequence buffer
	String<char>  qual;                // quality values
	char          qualBuf[BUF_SIZE];   // quality value buffer

	String<Dna5>  altPatFw[3];              // forward-strand sequence
	uint8_t       altPatBufFw[3][BUF_SIZE]; // forward-strand sequence buffer
	String<Dna5>  altPatRc[3];              // reverse-complement sequence
	uint8_t       altPatBufRc[3][BUF_SIZE]; // reverse-complement sequence buffer
	String<char>  altQual[3];               // quality values for alternate basecalls
	char          altQualBuf[3][BUF_SIZE];  // quality value buffer for alternate basecalls

	String<Dna5>  patFwRev;               // forward-strand sequence reversed
	uint8_t       patBufFwRev[BUF_SIZE];  // forward-strand sequence buffer reversed
	String<Dna5>  patRcRev;               // reverse-complement sequence reversed
	uint8_t       patBufRcRev[BUF_SIZE];  // reverse-complement sequence buffer reversed
	String<char>  qualRev;                // quality values reversed
	char          qualBufRev[BUF_SIZE];   // quality value buffer reversed

	String<Dna5>  altPatFwRev[3];              // forward-strand sequence reversed
	uint8_t       altPatBufFwRev[3][BUF_SIZE]; // forward-strand sequence buffer reversed
	String<Dna5>  altPatRcRev[3];              // reverse-complement sequence reversed
	uint8_t       altPatBufRcRev[3][BUF_SIZE]; // reverse-complement sequence buffer reversed
	String<char>  altQualRev[3];              // quality values for alternate basecalls reversed
	char          altQualBufRev[3][BUF_SIZE]; // quality value buffer for alternate basecalls reversed

	char          readOrigBuf[FileBuf::LASTN_BUF_SZ];
	size_t        readOrigBufLen;

	String<char>  name;                // read name
	char          nameBuf[BUF_SIZE];   // read name buffer
	uint32_t      patid;               // unique 0-based id based on order in read file(s)
	int           mate;                // 0 = single-end, 1 = mate1, 2 = mate2
	uint32_t      seed;                // random seed
	int           alts;                // number of alternatives
	bool          fuzzy;               // whether to employ fuzziness
	HitSet        hitset;              // holds previously-found hits; for chaining
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
	              const char *dumpfile = NULL,
	              bool verbose = false) :
	    readCnt_(0),
		dumpfile_(dumpfile),
		numWrappers_(0),
		doLocking_(true),
		useSpinlock_(useSpinlock),
		randomizeQuals_(randomizeQuals),
		lock_(),
		verbose_(verbose)
	{
		// Open dumpfile, if specified
		if(dumpfile_ != NULL) {
			out_.open(dumpfile_, ios_base::out);
			if(!out_.good()) {
				cerr << "Could not open pattern dump file \"" << dumpfile_ << "\" for writing" << endl;
				throw 1;
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
			ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name);
			if(!rb.empty()) {
				rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name);
			}
			// Output it, if desired
			if(dumpfile_ != NULL) {
				dumpBuf(ra);
				if(!rb.empty()) {
					dumpBuf(rb);
				}
			}
			if(verbose_) {
				cout << "Parsed mate 1: "; ra.dump(cout);
				cout << "Parsed mate 2: "; rb.dump(cout);
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
			r.seed = genRandSeed(r.patFw, r.qual, r.name);
			// Output it, if desired
			if(dumpfile_ != NULL) {
				dumpBuf(r);
			}
			if(verbose_) {
				cout << "Parsed read: "; r.dump(cout);
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

	/**
	 * Return the number of reads attempted.
	 */
	uint64_t readCnt() const {
		return readCnt_ - 1;
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
				r.qual[i] *= (r.qual[i+1] + 7);
			}
			if(i > 0) {
				r.qual[i] *= (r.qual[i-1] + 11);
			}
			// A user says that g++ complained here about "comparison
			// is always false due to limited range of data type", but
			// I can't see why.  I added the (int) cast to try to avoid
			// the warning.
			if((int)r.qual[i] < 0) r.qual[i] = -(r.qual[i]+1);
			r.qual[i] %= 41;
			assert_leq(r.qual[i], 40);
			r.qual[i] += 33;
		}
	}

	/**
	 * Dump the contents of the ReadBuf to the dump file.
	 */
	void dumpBuf(const ReadBuf& r) {
		assert(dumpfile_ != NULL);
		dump(out_, r.patFw,
		     empty(r.qual) ? String<char>("(empty)") : r.qual,
		     empty(r.name)   ? String<char>("(empty)") : r.name);
		dump(out_, r.patRc,
		     empty(r.qualRev) ? String<char>("(empty)") : r.qualRev,
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
	uint64_t readCnt_;

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
	bool verbose_;
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
	virtual pair<uint64_t,uint64_t> readCnt() const = 0;

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

	/**
	 * Return the number of reads attempted.
	 */
	virtual pair<uint64_t,uint64_t> readCnt() const {
		uint64_t ret = 0llu;
		vector<PatternSource*>::const_iterator it;
		for(it = src_.begin(); it != src_.end(); it++) {
			ret += (*it)->readCnt();
		}
		return make_pair(ret, 0llu);
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

	/**
	 * Return the number of reads attempted.
	 */
	virtual pair<uint64_t,uint64_t> readCnt() const {
		uint64_t rets = 0llu, retp = 0llu;
		for(size_t i = 0; i < srca_.size(); i++) {
			if(srcb_[i] == NULL) {
				rets += srca_[i]->readCnt();
			} else {
				assert_eq(srca_[i]->readCnt(), srcb_[i]->readCnt());
				retp += srca_[i]->readCnt();
			}
		}
		return make_pair(rets, retp);
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

	uint32_t      patid() const { return patid_;        }
	virtual void  reset()       { patid_ = 0xffffffff;  }
	bool          empty() const { return buf1_.empty(); }
	uint32_t length(int mate) const {
		return (mate == 1)? buf1_.length() : buf2_.length();
	}

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
	                      bool verbose = false,
	                      int trim3 = 0,
	                      int trim5 = 0) :
		PatternSource(randomizeQuals, useSpinlock, dumpfile, verbose),
		trim3_(trim3), trim5_(trim5) { }
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
	                    bool verbose = false,
	                    uint32_t seed = 0) :
		PatternSource(false, useSpinlock, dumpfile, verbose),
		numReads_(numReads),
		length_(length),
		seed_(seed)
	{
		if(length_ > 1024) {
			cerr << "Read length for RandomPatternSource may not exceed 1024; got " << length_ << endl;
			throw 1;
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
			r.qualBuf[i]            = c;
		}
		_setBegin (r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, length);
		_setBegin (r.qual, r.qualBuf);
		_setLength(r.qual, length);
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
			throw 1;
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
	                    bool verbose = false,
	                    int trim3 = 0,
	                    int trim5 = 0,
		                uint32_t skip = 0) :
		TrimmingPatternSource(randomizeQuals, useSpinlock, dumpfile,
		                      verbose, trim3, trim5),
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
		r.qual = quals_[cur_];
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
		ra.qual = quals_[cur_];
		cur_++;
		rb.patFw  = v_[cur_];
		rb.qual = quals_[cur_];
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
	vector<String<Dna5> > v_;     /// forward sequences
	vector<String<char> > quals_; /// quality values parallel to v_
	vector<String<char> > names_; /// names
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
	                          bool verbose = false,
	                          int trim3 = 0,
	                          int trim5 = 0,
	                          uint32_t skip = 0) :
		TrimmingPatternSource(randomizeQuals, useSpinlock, dumpfile,
		                      verbose, trim3, trim5),
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
			} else if((in = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
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
		throw 1;
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
	                   bool verbose = false,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool __forgiveInput = false,
	                   uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, useSpinlock,
		                          __forgiveInput, dumpfile, verbose,
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
				throw 1;
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
						 << "reads and re-run Bowtie" << endl;
					throw 1;
				}
				r.patBufFw [dstLen] = charToDna5[c];
				r.qualBuf[dstLen]   = 'I';
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
		_setBegin (r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, dstLen);
		_setBegin (r.qual,  r.qualBuf);
		_setLength(r.qual,  dstLen);

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
		throw 1;
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


/**
 * Tokenize a line of space-separated integer quality values.
 */
static inline bool tokenizeQualLine(FileBuf& filebuf, char *buf, size_t buflen, vector<string>& toks) {
	size_t rd = filebuf.gets(buf, buflen);
	if(rd == 0) return false;
	assert(NULL == strrchr(buf, '\n'));
	tokenize(string(buf), " ", toks);
	return true;
}

extern void wrongQualityFormat(const String<char>& read_name);
extern void tooFewQualities(const String<char>& read_name);
extern void tooManyQualities(const String<char>& read_name);
extern void tooManySeqChars(const String<char>& read_name);

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
	                   bool verbose = false,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool forgiveInput = false,
	                   bool solQuals = false,
	                   bool phred64Quals = false,
	                   bool intQuals = false,
	                   uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, useSpinlock,
		                          forgiveInput, dumpfile, verbose,
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
							 << "reads and re-run Bowtie" << endl;
						throw 1;
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
				vector<string> s_quals;
				if(!tokenizeQualLine(filebuf_, buf, 4096, s_quals)) break;
				for (unsigned int j = 0; j < s_quals.size(); ++j) {
					char c = intToPhred33(atoi(s_quals[j].c_str()), solQuals_);
					assert_geq(c, 33);
					if (qualsRead >= trim5_) {
						size_t off = qualsRead - trim5_;
						if(off >= 1024) tooManyQualities(r.name);
						r.qualBuf[off] = c;
						}
					++qualsRead;
				}
			} // done reading integer quality lines
			if (charsRead > qualsRead) tooFewQualities(r.name);
		} else {
			// Non-integer qualities
			while((qualsRead < dstLen + this->trim5_) && c >= 0) {
				c = filebuf_.get();
				c2 = c;
				if (c == ' ') wrongQualityFormat(r.name);
				if(c < 0) {
					// EOF occurred in the middle of a read - abort
					return -1;
				}
				if(!isspace(c) && c != upto && (upto2 == -1 || c != upto2)) {
					if (qualsRead >= trim5_) {
						size_t off = qualsRead - trim5_;
						if(off >= 1024) tooManyQualities(r.name);
						c = charToPhred33(c, solQuals_, phred64Quals_);
						assert_geq(c, 33);
						r.qualBuf[off] = c;
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
		_setBegin (r.qual, (char*)r.qualBuf);
		_setLength(r.qual, dstLen);
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
			bool verbose = false,
			uint32_t skip = 0,
			uint32_t seed = 0) :
		BufferedFilePatternSource(infiles, false, useSpinlock,
		                          false, dumpfile, verbose, 0, 0, skip),
		length_(length), freq_(freq),
		eat_(length_-1), beginning_(true),
		nameChars_(0), bufCur_(0), subReadCnt_(0llu)
	{
		resetForNextFile();
		assert_lt(length_, (size_t)ReadBuf::BUF_SIZE);
	}

	virtual void reset() {
		BufferedFilePatternSource::reset();
		resetForNextFile();
	}

protected:
	/// Read another pattern from a FASTA input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		while(true) {
			int c = filebuf_.get();
			if(c < 0) {
				seqan::clear(r.patFw);
				return;
			}
			if(c == '>') {
				resetForNextFile();
				c = filebuf_.peek();
				bool sawSpace = false;
				while(c != '\n' && c != '\r') {
					if(!sawSpace) {
						sawSpace = isspace(c);
					}
					if(!sawSpace) {
						nameBuf_[nameChars_++] = c;
					}
					filebuf_.get();
					c = filebuf_.peek();
				}
				while(c == '\n' || c == '\r') {
					filebuf_.get();
					c = filebuf_.peek();
				}
				nameBuf_[nameChars_++] = '_';
			} else {
				int cat = dna4Cat[c];
				if(cat == 2) c = 'N';
				if(cat == 0) {
					// Encountered non-DNA, non-IUPAC char; skip it
					continue;
				} else {
					// DNA char
					buf_[bufCur_++] = c;
					if(bufCur_ == 1024) bufCur_ = 0;
					if(eat_ > 0) {
						eat_--;
						// Try to keep readCnt_ aligned with the offset
						// into the reference; that let's us see where
						// the sampling gaps are by looking at the read
						// name
						if(!beginning_) readCnt_++;
						continue;
					}
					for(size_t i = 0; i < length_; i++) {
						if(length_ - i <= bufCur_) {
							c = buf_[bufCur_ - (length_ - i)];
						} else {
							// Rotate
							c = buf_[bufCur_ - (length_ - i) + 1024];
						}
						r.patBufFw [i] = charToDna5[c];
						r.qualBuf[i] = 'I';
					}
					_setBegin (r.patFw,  (Dna5*)r.patBufFw);
					_setLength(r.patFw,  length_);
					_setBegin (r.qual, r.qualBuf);
					_setLength(r.qual, length_);
					// Set up a default name if one hasn't been set
					for(size_t i = 0; i < nameChars_; i++) {
						r.nameBuf[i] = nameBuf_[i];
					}
					itoa10(readCnt_ - subReadCnt_, &r.nameBuf[nameChars_]);
					_setBegin(r.name, r.nameBuf);
					_setLength(r.name, strlen(r.nameBuf));
					eat_ = freq_-1;
					readCnt_++;
					beginning_ = false;
					patid = readCnt_-1;
					break;
				}
			}
		}
	}
	/// Shouldn't ever be here; it's not sensible to obtain read pairs
	// from a continuous input.
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		cerr << "In FastaContinuousPatternSource.readPair()" << endl;
		throw 1;
	}

	/**
	 * Reset state to be read for the next file.
	 */
	virtual void resetForNextFile() {
		eat_ = length_-1;
		beginning_ = true;
		bufCur_ = nameChars_ = 0;
		subReadCnt_ = readCnt_;
	}

private:
	size_t length_;     /// length of reads to generate
	size_t freq_;       /// frequency to sample reads
	int policy_;        /// policy for handling Ns

	size_t eat_;        /// number of characters we need to skip before
	                    /// we have flushed all of the ambiguous or
	                    /// non-existent characters out of our read
	                    /// window
	bool beginning_;    /// skipping over the first read length?
	char buf_[1024];    /// read buffer
	char nameBuf_[1024];/// read buffer for name of fasta record being
	                    /// split into mers
	size_t nameChars_;  /// number of name characters copied into buf
	size_t bufCur_;     /// buffer cursor; points to where we should
	                    /// insert the next character
	uint64_t subReadCnt_;/// number to subtract from readCnt_ to get
	                    /// the pat id to output (so it resets to 0 for
	                    /// each new sequence)
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
	                   bool verbose = false,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool __forgiveInput = false,
	                   bool solexa_quals = false,
	                   bool phred64Quals = false,
	                   bool integer_quals = false,
	                   bool fuzzy = false,
	                   uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, useSpinlock,
		                          __forgiveInput, dumpfile, verbose,
		                          trim3, trim5, skip),
		first_(true),
		solQuals_(solexa_quals),
		phred64Quals_(phred64Quals),
		intQuals_(integer_quals),
		fuzzy_(fuzzy)
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
		r.fuzzy = fuzzy_;
		r.alts = 0;
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
				throw 1;
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
				throw 1;
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
		uint8_t *sbuf = r.patBufFw;
		int dstLens[] = {0, 0, 0, 0};
		int *dstLenCur = &dstLens[0];
		int trim5 = this->trim5_;
		int altBufIdx = 0;
		while(c != '+') {
			if(fuzzy_ && c == '-') c = 'A';
			if(isalpha(c)) {
				// If it's past the 5'-end trim point
				if(charsRead >= trim5) {
					if((*dstLenCur) >= 1024) tooManySeqChars(r.name);
					sbuf[(*dstLenCur)++] = charToDna5[c];
				}
				charsRead++;
			} else if(fuzzy_ && c == ' ') {
				trim5 = 0; // disable 5' trimming for now
				if(charsRead == 0) {
					c = filebuf_.get();
					continue;
					}
				charsRead = 0;
				if(altBufIdx >= 3) {
					cerr << "At most 3 alternate sequence strings permitted; offending read: " << r.name << endl;
					throw 1;
				}
				// Move on to the next alternate-sequence buffer
				sbuf = r.altPatBufFw[altBufIdx++];
				dstLenCur = &dstLens[altBufIdx];
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
		dstLen = dstLens[0];
		charsRead = dstLen + this->trim5_;
		dstLen -= this->trim3_;
		// Set trimmed bounds of buffers
		_setBegin(r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, dstLen);
		assert_eq('+', c);

		// Chew up the optional name on the '+' line
		peekToEndOfLine(filebuf_);

		// Now read the qualities
		if (intQuals_) {
			assert(!fuzzy_);
			int qualsRead = 0;
			char buf[4096];
			while (qualsRead < charsRead) {
				vector<string> s_quals;
				if(!tokenizeQualLine(filebuf_, buf, 4096, s_quals)) break;
				for (unsigned int j = 0; j < s_quals.size(); ++j) {
					char c = intToPhred33(atoi(s_quals[j].c_str()), solQuals_);
					assert_geq(c, 33);
					if (qualsRead >= trim5_) {
						size_t off = qualsRead - trim5_;
						if(off >= 1024) tooManyQualities(r.name);
						r.qualBuf[off] = c;
					}
					++qualsRead;
				}
			} // done reading integer quality lines
			if (charsRead > qualsRead) tooFewQualities(r.name);
			_setBegin(r.qual, (char*)r.qualBuf);
			_setLength(r.qual, dstLen);
			peekOverNewline(filebuf_);
		} else {
			// Non-integer qualities
			char *qbuf = r.qualBuf;
			altBufIdx = 0;
			trim5 = this->trim5_;
			int qualsRead[4] = {0, 0, 0, 0};
			int *qualsReadCur = &qualsRead[0];
			while((*qualsReadCur) < dstLen + this->trim5_ || fuzzy_) {
				c = filebuf_.get();
				if (!fuzzy_ && c == ' ') {
					wrongQualityFormat(r.name);
				} else if(c == ' ') {
					trim5 = 0; // disable 5' trimming for now
					if((*qualsReadCur) == 0) continue;
					if(altBufIdx >= 3) {
						cerr << "At most 3 alternate quality strings permitted; offending read: " << r.name << endl;
						throw 1;
					}
					qbuf = r.altQualBuf[altBufIdx++];
					qualsReadCur = &qualsRead[altBufIdx];
					continue;
				}
				if(c < 0) {
					// EOF occurred in the middle of a read - abort
					seqan::clear(r.patFw);
					filebuf_.resetLastN();
					return;
				}
				if (c != '\r' && c != '\n') {
					if (*qualsReadCur >= trim5) {
						size_t off = (*qualsReadCur) - trim5;
						if(off >= 1024) tooManyQualities(r.name);
						c = charToPhred33(c, solQuals_, phred64Quals_);
						assert_geq(c, 33);
						qbuf[off] = c;
					}
					(*qualsReadCur)++;
				} else {
					break;
				}
			}
			_setBegin (r.qual, (char*)r.qualBuf);
			_setLength(r.qual, dstLen);

			if(fuzzy_) {
				// Trim from 3' end of alternate basecall and quality strings
				if(this->trim3_ > 0) {
					for(int i = 1; i < 4; i++) {
						assert_eq(qualsRead[i], dstLens[i]);
						qualsRead[i] = dstLens[i] =
							max<int>(0, dstLens[i] - this->trim3_);
					}
				}
				// Shift to RHS, and install in Strings
				assert_eq(0, r.alts);
				for(int i = 1; i < 4; i++) {
					if(qualsRead[i] == 0) continue;
					if(qualsRead[i] > dstLen) {
						// Shift everybody up
						int shiftAmt = qualsRead[i] - dstLen;
						for(int j = 0; j < dstLen; j++) {
							r.altQualBuf[i-1][j]  = r.altQualBuf[i-1][j+shiftAmt];
							r.altPatBufFw[i-1][j] = r.altPatBufFw[i-1][j+shiftAmt];
						}
					} else if (qualsRead[i] < dstLen) {
						// Shift everybody down
						int shiftAmt = dstLen - qualsRead[i];
						for(int j = dstLen-1; j >= shiftAmt; j--) {
							r.altQualBuf[i-1][j]  = r.altQualBuf[i-1][j-shiftAmt];
							r.altPatBufFw[i-1][j] = r.altPatBufFw[i-1][j-shiftAmt];
						}
						// Fill in unset positions
						for(int j = 0; j < shiftAmt; j++) {
							// '!' - indicates no alternate basecall at
							// this position
							r.altQualBuf[i-1][j] = (char)33;
						}
					}
					_setBegin (r.altPatFw[i-1], (Dna5*)r.altPatBufFw[i-1]);
					_setLength(r.altPatFw[i-1], dstLen);
					_setBegin (r.altQual[i-1], (char*)r.altQualBuf[i-1]);
					_setLength(r.altQual[i-1], dstLen);
					r.alts++;
				}
			}

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
		throw 1;
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
	bool fuzzy_;
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
	                 bool verbose = false,
	                 int trim3 = 0,
	                 int trim5 = 0,
	                 uint32_t skip = 0) :
		BufferedFilePatternSource(infiles, randomizeQuals, false, useSpinlock,
		                          dumpfile, verbose, trim3, trim5, skip),
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
				throw 1;
			}
			first_ = false;
		}

		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		while(!isspace(c) && c >= 0) {
			c = filebuf_.get();
			if(isalpha(c) && dstLen >= this->trim5_) {
				size_t len = dstLen - this->trim5_;
				if(len >= 1024) tooManyQualities(String<char>("(no name)"));
				r.patBufFw [len] = charToDna5[c];
				r.qualBuf[len] = 'I';
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
		_setBegin (r.qual, r.qualBuf);
		_setLength(r.qual, dstLen);

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
		throw 1;
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

/**
 * Read a Raw-format file (one sequence per line).  No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class ChainPatternSource : public BufferedFilePatternSource {
public:
	ChainPatternSource(const vector<string>& infiles,
	                   bool useSpinlock,
	                   const char *dumpfile,
	                   bool verbose,
	                   uint32_t skip) :
	BufferedFilePatternSource(
		infiles, false, false, useSpinlock, dumpfile, verbose, 0, 0, skip) { }

protected:

	/// Read another pattern from a Raw input file
	virtual void read(ReadBuf& r, uint32_t& patid) {
		filebuf_.peek();
		if(filebuf_.eof()) {
			filebuf_.resetLastN();
			seqan::clear(r.patFw);
			return;
		}
		do {
			r.hitset.deserialize(filebuf_);
		} while(!r.hitset.initialized() && !filebuf_.eof());
		if(!r.hitset.initialized()) {
			filebuf_.resetLastN();
			seqan::clear(r.patFw);
			return;
		}
		// Now copy the name/sequence/quals into r.name/r.patFw/r.qualFw
		_setBegin(r.name, (char*)r.nameBuf);
		_setCapacity(r.name, seqan::length(r.hitset.name));
		_setLength(r.name, seqan::length(r.hitset.name));
		memcpy(r.nameBuf, seqan::begin(r.hitset.name), seqan::length(r.hitset.name));
		_setBegin (r.patFw, (Dna5*)r.patBufFw);
		_setCapacity(r.patFw, seqan::length(r.hitset.seq));
		_setLength(r.patFw, seqan::length(r.hitset.seq));
		memcpy(r.patBufFw, seqan::begin(r.hitset.seq), seqan::length(r.hitset.seq));
		_setBegin (r.qual, r.qualBuf);
		_setCapacity(r.qual, seqan::length(r.hitset.qual));
		_setLength(r.qual, seqan::length(r.hitset.qual));
		memcpy(r.qualBuf, seqan::begin(r.hitset.qual), seqan::length(r.hitset.qual));

		r.readOrigBufLen = filebuf_.copyLastN(r.readOrigBuf);
		filebuf_.resetLastN();

		readCnt_++;
		patid = readCnt_-1;
	}

	/// Read another read pair
	virtual void readPair(ReadBuf& ra, ReadBuf& rb, uint32_t& patid) {
		// (For now, we shouldn't ever be here)
		cerr << "In ChainPatternSource.readPair()" << endl;
		throw 1;
	}
};

#endif /*PAT_H_*/
