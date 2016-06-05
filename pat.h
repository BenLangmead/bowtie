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
#include "threading.h"
#include "filebuf.h"
#include "qual.h"
#include "hit_set.h"
#include "search_globals.h"

/**
 * Classes and routines for reading reads from various input sources.
 */

using namespace std;
using namespace seqan;

/// Constructs string base-10 representation of integer 'value'
extern char* itoa10(int value, char* result);

typedef uint64_t TReadId;

/**
 * Calculate a per-read random seed based on a combination of
 * the read data (incl. sequence, name, quals) and the global
 * seed in '_randSeed'.
 */
static inline uint32_t genRandSeed(const String<Dna5>& qry,
                                   const String<char>& qual,
                                   const String<char>& name,
                                   uint32_t seed)
{
	// Calculate a per-read random seed based on a combination of
	// the read data (incl. sequence, name, quals) and the global
	// seed
	uint32_t rseed = (seed + 101) * 59 * 61 * 67 * 71 * 73 * 79 * 83;
	size_t qlen = seqan::length(qry);
	// Throw all the characters of the read into the random seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qry[i];
		assert_leq(p, 4);
		size_t off = ((i & 15) << 1);
		rseed ^= (p << off);
	}
	// Throw all the quality values for the read into the random
	// seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qual[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed ^= (p << off);
	}
	// Throw all the characters in the read name into the random
	// seed
	size_t namelen = seqan::length(name);
	for(size_t i = 0; i < namelen; i++) {
		int p = (int)name[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed ^= (p << off);
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
	}

#define RESET_BUF(str, buf, typ) _setBegin(str, (typ*)buf); _setLength(str, 0); _setCapacity(str, BUF_SIZE);
#define RESET_BUF_LEN(str, buf, len, typ) _setBegin(str, (typ*)buf); _setLength(str, len); _setCapacity(str, BUF_SIZE);

	/// Point all Strings to the beginning of their respective buffers
	/// and set all lengths to 0
	void reset() {
		patid = 0;
		readOrigBufLen = 0;
		qualOrigBufLen = 0;
		trimmed5 = trimmed3 = 0;
		color = false;
		primer = '?';
		trimc = '?';
		seed = 0;
		RESET_BUF(patFw, patBufFw, Dna5);
		RESET_BUF(patRc, patBufRc, Dna5);
		RESET_BUF(qual, qualBuf, char);
		RESET_BUF(patFwRev, patBufFwRev, Dna5);
		RESET_BUF(patRcRev, patBufRcRev, Dna5);
		RESET_BUF(qualRev, qualBufRev, char);
		RESET_BUF(name, nameBuf, char);
	}

	void clearAll() {
		seqan::clear(patFw);
		seqan::clear(patRc);
		seqan::clear(qual);
		seqan::clear(patFwRev);
		seqan::clear(patRcRev);
		seqan::clear(qualRev);
		seqan::clear(name);
		trimmed5 = trimmed3 = 0;
		readOrigBufLen = 0;
		qualOrigBufLen = 0;
		color = false;
		primer = '?';
		trimc = '?';
		seed = 0;
	}

	/// Return true iff the read (pair) is empty
	bool empty() const {
		return seqan::empty(patFw);
	}

	/// Return length of the read in the buffer
	uint32_t length() const {
		return (uint32_t)seqan::length(patFw);
	}

	/**
	 * Construct reverse complement of the pattern.  If read is in
	 * colorspace, reverse color string.
	 */
	void constructRevComps() {
		uint32_t len = length();
		assert_gt(len, 0);
		RESET_BUF_LEN(patRc, patBufRc, len, Dna5);
		if(color) {
			for(uint32_t i = 0; i < len; i++) {
				// Reverse the sequence
				patBufRc[i]  = patBufFw[len-i-1];
			}
		} else {
			for(uint32_t i = 0; i < len; i++) {
				// Reverse-complement the sequence
				patBufRc[i]  = (patBufFw[len-i-1] == 4) ? 4 : (patBufFw[len-i-1] ^ 3);
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
		for(uint32_t i = 0; i < len; i++) {
			patFwRev[i]  = patFw[len-i-1];
			patRcRev[i]  = patRc[len-i-1];
			qualRev[i]   = qual[len-i-1];
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

	/**
	 * Dump basic information about this read to the given ostream.
	 */
	void dump(std::ostream& os) const {
		os << name << ' ';
		if(color) {
			for(size_t i = 0; i < seqan::length(patFw); i++) {
				os << "0123."[(int)patFw[i]];
			}
		} else {
			os << patFw;
		}
		os << qual << " ";
	}

	static const int BUF_SIZE = 1024;

	String<Dna5>  patFw;               // forward-strand sequence
	uint8_t       patBufFw[BUF_SIZE];  // forward-strand sequence buffer
	String<Dna5>  patRc;               // reverse-complement sequence
	uint8_t       patBufRc[BUF_SIZE];  // reverse-complement sequence buffer
	String<char>  qual;                // quality values
	char          qualBuf[BUF_SIZE];   // quality value buffer

	String<Dna5>  patFwRev;               // forward-strand sequence reversed
	uint8_t       patBufFwRev[BUF_SIZE];  // forward-strand sequence buffer reversed
	String<Dna5>  patRcRev;               // reverse-complement sequence reversed
	uint8_t       patBufRcRev[BUF_SIZE];  // reverse-complement sequence buffer reversed
	String<char>  qualRev;                // quality values reversed
	char          qualBufRev[BUF_SIZE];   // quality value buffer reversed

	// For remembering the exact input text used to define a read
	char          readOrigBuf[FileBuf::LASTN_BUF_SZ];
	size_t        readOrigBufLen;

	// For when qualities are in a separate file
	char          qualOrigBuf[FileBuf::LASTN_BUF_SZ];
	size_t        qualOrigBufLen;

	String<char>  name;                // read name
	char          nameBuf[BUF_SIZE];   // read name buffer
	uint32_t      patid;               // unique 0-based id based on order in read file(s)
	int           mate;                // 0 = single-end, 1 = mate1, 2 = mate2
	uint32_t      seed;                // random seed
	bool          color;               // whether read is in color space
	char          primer;              // primer base, for csfasta files
	char          trimc;               // trimmed color, for csfasta files
	int           trimmed5;            // amount actually trimmed off 5' end
	int           trimmed3;            // amount actually trimmed off 3' end
	HitSet        hitset;              // holds previously-found hits; for chaining
};

struct BoolTriple {

	BoolTriple() : a(false), b(false), c(false) { }

	BoolTriple(bool a_, bool b_, bool c_) : a(a_), b(b_), c(c_) { }

	bool a;
	bool b;
	bool c;
};

/**
 * All per-thread storage for input read data.
 */
struct PerThreadReadBuf {
	
	PerThreadReadBuf(size_t max_buf) :
		max_buf_(max_buf),
		bufa_(max_buf),
		bufb_(max_buf),
		rdid_()
	{
		bufa_.resize(max_buf);
		bufb_.resize(max_buf);
		reset();
	}
	
	ReadBuf& read_a() { return bufa_[cur_buf_]; }
	ReadBuf& read_b() { return bufb_[cur_buf_]; }
	
	const ReadBuf& read_a() const { return bufa_[cur_buf_]; }
	const ReadBuf& read_b() const { return bufb_[cur_buf_]; }
	
	/**
	 * Return read id for read/pair currently in the buffer.
	 */
	TReadId rdid() const {
		assert_neq(rdid_, std::numeric_limits<TReadId>::max());
		return rdid_ + cur_buf_;
	}
	
	/**
	 * Reset state as though no reads have been read.
	 */
	void reset() {
		cur_buf_ = bufa_.size();
		for(size_t i = 0; i < max_buf_; i++) {
			bufa_[i].reset();
			bufb_[i].reset();
		}
		rdid_ = std::numeric_limits<TReadId>::max();
	}
	
	/**
	 * Advance cursor to next element
	 */
	void next() {
		assert_lt(cur_buf_, bufa_.size());
		cur_buf_++;
	}
	
	/**
	 * Return true when there's nothing left to dish out.
	 */
	bool exhausted() {
		assert_leq(cur_buf_, bufa_.size());
		return cur_buf_ >= bufa_.size()-1;
	}
	
	/**
	 * Just after a new batch has been loaded, use init to
	 * set the cuf_buf_ and rdid_ fields appropriately.
	 */
	void init() {
		cur_buf_ = 0;
	}
	
	/**
	 * Set the read id of the first read in the buffer.
	 */
	void setReadId(TReadId rdid) {
		rdid_ = rdid;
	}
	
	const size_t max_buf_; // max # reads to read into buffer at once
	vector<ReadBuf> bufa_;     // Read buffer for mate as
	vector<ReadBuf> bufb_;     // Read buffer for mate bs
	size_t cur_buf_;       // Read buffer currently active
	TReadId rdid_;         // index of read at offset 0 of bufa_/bufb_
};


/**
 * Encapsulates a synchronized source of patterns; usually a file.
 * Handles dumping patterns to a logfile (useful for debugging).  Also
 * optionally reverses reads and quality strings before returning them,
 * though that is usually more efficiently done by the concrete
 * subclass.  Concrete subclasses should delimit critical sections with
 * ThreadSafe objects.
 */
class PatternSource {
public:
	PatternSource(uint32_t seed,
	              bool randomizeQuals = false,
	              const char *dumpfile = NULL,
	              bool verbose = false) :
		seed_(seed),
		readCnt_(0),
		dumpfile_(dumpfile),
		numWrappers_(0),
		doLocking_(true),
		randomizeQuals_(randomizeQuals),
		mutex_m(),
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
	}

	virtual ~PatternSource() { }

	/**
	 * Call this whenever this PatternSource is wrapped by a new
	 * WrappedPatternSourcePerThread.  This helps us keep track of
	 * whether locks will be contended.
	 */
	//BP: remove?
	void addWrapper() {
		ThreadSafe ts(&mutex_m);
		numWrappers_++;
	}

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats that
	 * can read in a pair of reads in a single transaction with a
	 * single input source.  If paired-end input is given as a pair of
	 * parallel files, this member should throw an error and exit.
	 */
	virtual pair<bool, bool> nextReadPair(ReadBuf& ra, ReadBuf& rb) = 0;

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats
	 * where individual input sources look like single-end-read
	 * sources, e.g., formats where paired-end reads are specified in
	 * parallel read files.
	 */
	virtual bool nextRead(ReadBuf& r) = 0;

	/**
	 * Finishes parsing outside the critical section
	 */
	virtual void finalize(ReadBuf& r) const = 0;
	
	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats
	 * where individual input sources look like single-end-read
	 * sources, e.g., formats where paired-end reads are specified in
	 * parallel read files.
	 */
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		bool batch_a,
		bool lock = true) = 0;
	
	/**
	 * Finishes parsing a given read.  Happens outside the critical section.
	 */
	virtual bool parse(ReadBuf& ra, ReadBuf& rb) const = 0;
	
	/// Reset state to start over again with the first read
	virtual void reset() { readCnt_ = 0; }

	/**
	 * Return the number of reads attempted.
	 */
	//BP: changed, but should I have?
	uint64_t readCount() const { return readCnt_; }

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
		     empty(r.name) ? String<char>("(empty)") : r.name);
		dump(out_, r.patRc,
		     empty(r.qualRev) ? String<char>("(empty)") : r.qualRev,
		     empty(r.name) ? String<char>("(empty)") : r.name);
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

	uint32_t seed_;

	/// The number of reads read by this PatternSource
	volatile uint64_t readCnt_;

	const char *dumpfile_; /// dump patterns to this file before returning them
	ofstream out_;         /// output stream for dumpfile
	int numWrappers_;      /// # threads that own a wrapper for this PatternSource
	bool doLocking_;       /// override whether to lock (true = don't override)
	bool randomizeQuals_;  /// true -> mess up qualities in a random way
	MUTEX_T mutex_m; /// mutex for locking critical regions
	bool verbose_;
};

/**
 * Abstract parent class for synhconized sources of paired-end reads
 * (and possibly also single-end reads).
 */
class PatternComposer {
public:
	PatternComposer(uint32_t seed) {
		seed_ = seed;
	}
	virtual ~PatternComposer() { }

	virtual void reset() = 0;
	
	//virtual pair<bool, bool> nextReadPair(ReadBuf& ra, ReadBuf& rb) = 0;
	
	/**
	 * Member function override by concrete, format-specific classes.
	 */
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf& pt) = 0;
	
	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(ReadBuf& ra, ReadBuf& rb) = 0;

	virtual void free_pmembers( const vector<PatternSource*> &elist) {
    		for (size_t i = 0; i < elist.size(); i++) {
        		if (elist[i] != NULL)
            			delete elist[i];
		}
	}

protected:

	MUTEX_T mutex_m; /// mutex for locking critical regions
	uint32_t seed_;
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class SoloPatternComposer : public PatternComposer {

public:

	SoloPatternComposer(const vector<PatternSource*>& src,
	                        uint32_t seed) :
		PatternComposer(seed), cur_(0), src_(src)
	{
	    for(size_t i = 0; i < src_.size(); i++) {
	    	assert(src_[i] != NULL);
	    }
	}

	virtual ~SoloPatternComposer() { free_pmembers(src_); }

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
	pair<bool, int> nextBatch(PerThreadReadBuf& pt) {
		uint32_t cur = cur_;
		while(cur < src_.size()) {
			// Patterns from srca_[cur_] are unpaired
			pair<bool, int> res;
			do {
				res = src_[cur]->nextBatch(
					pt,
					true,  // batch A (or pairs)
					true); // grab lock below
			} while(!res.first && res.second == 0);
			if(res.second == 0) {
				ThreadSafe ts(&mutex_m);
				if(cur + 1 > cur_) {
					cur_++;
				}
				cur = cur_;
				continue; // on to next pair of PatternSources
			}
			return res;
		}
		assert_leq(cur, src_.size());
		return make_pair(true, 0);
	}
	
	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(ReadBuf& ra, ReadBuf& rb) {
		return src_[0]->parse(ra, rb);
	}

protected:

	volatile uint32_t cur_; // current element in parallel srca_, srcb_ vectors
	vector<PatternSource*> src_; /// PatternSources for paired-end reads
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class DualPatternComposer : public PatternComposer {

public:

	DualPatternComposer(const vector<PatternSource*>& srca,
	                        const vector<PatternSource*>& srcb,
	                        uint32_t seed) :
		PatternComposer(seed), cur_(0), srca_(srca), srcb_(srcb)
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

	virtual ~DualPatternComposer() {
		free_pmembers(srca_);
		free_pmembers(srcb_);
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
	pair<bool, int> nextBatch(PerThreadReadBuf& pt) {
		// 'cur' indexes the current pair of PatternSources
		uint32_t cur = cur_;
		while(cur < srca_.size()) {
			if(srcb_[cur] == NULL) {
				pair<bool, int> res = srca_[cur]->nextBatch(
					pt,
					true,  // batch A (or pairs)
					true); // grab lock below
				bool done = res.first;
				if(!done && res.second == 0) {
					ThreadSafe ts(&mutex_m);
					if(cur + 1 > cur_) cur_++;
					cur = cur_; // Move on to next PatternSource
					continue; // on to next pair of PatternSources
				}
				return make_pair(done, res.second);
			} else {
				pair<bool, int> resa, resb;
				// Lock to ensure that this thread gets parallel reads
				// in the two mate files
				{
					ThreadSafe ts(&mutex_m);
					resa = srca_[cur]->nextBatch(
						pt,
						true,   // batch A
						false); // don't grab lock below
					resb = srcb_[cur]->nextBatch(
						pt,
						false,  // batch B
						false); // don't grab lock below
					assert_eq(srca_[cur]->readCount(),
						  srcb_[cur]->readCount());
				}
				if(resa.second < resb.second) {
					cerr << "Error, fewer reads in file specified with -1 "
					     << "than in file specified with -2" << endl;
					throw 1;
				} else if(resa.second == 0 && resb.second == 0) {
					ThreadSafe ts(&mutex_m);
					if(cur + 1 > cur_) {
						cur_++;
					}
					cur = cur_; // Move on to next PatternSource
					continue; // on to next pair of PatternSources
				} else if(resb.second < resa.second) {
					cerr << "Error, fewer reads in file specified with -2 "
					     << "than in file specified with -1" << endl;
					throw 1;
				}
				assert_eq(resa.first, resb.first);
				assert_eq(resa.second, resb.second);
				return make_pair(resa.first, resa.second);
			}
		}
		assert_leq(cur, srca_.size());
		return make_pair(true, 0);
	}
	
	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(ReadBuf& ra, ReadBuf& rb) {
		return srca_[0]->parse(ra, rb);
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

	PatternSourcePerThread(PatternComposer& patsrc, uint32_t max_buf, uint32_t seed) : 
						patsrc_(patsrc), buf_(max_buf),
      						last_batch_(false), last_batch_size_(0), seed_(seed)	{ }
	
	uint32_t patid() const {
		return buf_.read_a().patid;
	}
	
	virtual void reset() {
		buf_.reset();
	}
	
	bool empty() const {
		return buf_.read_a().empty();
	}
	
	uint32_t length(int mate) const {
		return (mate == 1)? buf_.read_a().length() : buf_.read_a().length();
	}

	/**
	 * Return true iff the buffers jointly contain a paired-end read.
	 */
	bool paired() {
		bool ret = !buf_.read_b().empty();
		assert(!ret || !empty());
		return ret;
	}

	/**
	 * Get the next paired or unpaired read from the wrapped
	 * PatternComposer.  Returns a pair of bools; first indicates
	 * whether we were successful, second indicates whether we're
	 * done.
	 */
	pair<bool, bool> nextReadPair() {
		// Prepare batch
		if(buf_.exhausted()) {
			pair<bool, int> res = nextBatch();
			if(res.first && res.second == 0) {
				return make_pair(false, true);
			}
			last_batch_ = res.first;
			last_batch_size_ = res.second;
			assert_eq(0, buf_.cur_buf_);
		} else {
			buf_.next(); // advance cursor
			assert_gt(buf_.cur_buf_, 0);
		}
		// Parse read/pair
		assert(strlen(buf_.read_a().readOrigBuf) != 0);
		assert(buf_.read_a().empty());
		if(!parse(buf_.read_a(), buf_.read_b())) {
			return make_pair(false, false);
		}
		// Finalize read/pair
		if(strlen(buf_.read_b().readOrigBuf) != 0) {
			finalizePair(buf_.read_a(), buf_.read_b());
		} else {
			finalize(buf_.read_a());
		}
		bool this_is_last = buf_.cur_buf_ == last_batch_size_-1;
		return make_pair(true, this_is_last ? last_batch_ : false);
	}
	virtual void finalize(ReadBuf& ra) {
		ra.mate = 1;
		ra.patid = buf_.rdid();
		ra.constructRevComps();
		ra.constructReverses();
		ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
		//ra.fixMateName(1);
	}


	virtual void finalizePair(ReadBuf& ra, ReadBuf& rb) {
		ra.patid = rb.patid = buf_.rdid();
		
		ra.mate = 1;
		ra.constructRevComps();
		ra.constructReverses();
		ra.fixMateName(1);
		ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
		
		rb.mate = 2;
		rb.constructRevComps();
		rb.constructReverses();
		rb.fixMateName(2);
		rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name, seed_);
	}

	ReadBuf& bufa() { return buf_.read_a(); }
	ReadBuf& bufb() { return buf_.read_b(); }
	
	const ReadBuf& bufa() const { return buf_.read_a(); }
	const ReadBuf& bufb() const { return buf_.read_b(); }

	
private:

	/**
	 * When we've finished fully parsing and dishing out reads in
	 * the current batch, we go get the next one by calling into
	 * the composition layer.
	 */
	std::pair<bool, int> nextBatch() {
		buf_.reset();
		std::pair<bool, int> res = patsrc_.nextBatch(buf_);
		buf_.init();
		return res;
	}
	
	/**
	 * Call into composition layer (which in turn calls into
	 * format layer) to parse the read.
	 */
	bool parse(ReadBuf& ra, ReadBuf& rb) {
		return patsrc_.parse(ra, rb);
	}

	PatternComposer& patsrc_; // pattern composer
	PerThreadReadBuf buf_;    // read data buffer
	bool last_batch_;         // true if this is final batch
	int last_batch_size_;     // # reads read in previous batch
	uint32_t seed_;
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class PatternSourcePerThreadFactory {
public:
	PatternSourcePerThreadFactory(
		PatternComposer& patsrc, uint32_t max_buf, uint32_t seed): 
		patsrc_(patsrc), max_buf_(max_buf), seed_(seed) {}

	virtual PatternSourcePerThread* create() const {
		return new PatternSourcePerThread(patsrc_, max_buf_, seed_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * WrappedPatternSourcePerThreads.
	 */
	virtual std::vector<PatternSourcePerThread*>* create(uint32_t n) const {
		std::vector<PatternSourcePerThread*>* v = new std::vector<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new PatternSourcePerThread(patsrc_, max_buf_, seed_));
			assert(v->back() != NULL);
		}
		return v;
	}

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
private:
	/// Container for obtaining paired reads from PatternSources
	PatternComposer& patsrc_;
	/// Maximum size of batch to read in
	uint32_t max_buf_;
	uint32_t seed_;
};

/**
 * Encapsualtes a source of patterns where each raw pattern is trimmed
 * by some user-defined amount on the 3' and 5' ends.  Doesn't
 * implement the actual trimming - that's up to the concrete
 * descendants.
 */
class TrimmingPatternSource : public PatternSource {
public:
	TrimmingPatternSource(uint32_t seed,
	                      bool randomizeQuals = false,
	                      const char *dumpfile = NULL,
	                      bool verbose = false,
	                      int trim3 = 0,
	                      int trim5 = 0) :
		PatternSource(seed, randomizeQuals, dumpfile, verbose),
		trim3_(trim3), trim5_(trim5) { }
protected:
	int trim3_;
	int trim5_;
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

extern void wrongQualityFormat(const String<char>& read_name);
extern void tooFewQualities(const String<char>& read_name);
extern void tooManyQualities(const String<char>& read_name);
extern void tooManySeqChars(const String<char>& read_name);

/**
 * Encapsulates a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public TrimmingPatternSource {
public:
	VectorPatternSource(uint32_t seed,
	                    const vector<string>& v,
	                    bool color,
	                    bool randomizeQuals = false,
	                    const char *dumpfile = NULL,
	                    bool verbose = false,
	                    int trim3 = 0,
	                    int trim5 = 0,
		                uint32_t skip = 0) :
		TrimmingPatternSource(seed, randomizeQuals,
		                      dumpfile, verbose, trim3, trim5),
		color_(color), cur_(skip), skip_(skip), paired_(false), v_(),
		quals_()
	{
		for(size_t i = 0; i < v.size(); i++) {
			vector<string> ss;
			tokenize(v[i], ":", ss, 2);
			assert_gt(ss.size(), 0);
			assert_leq(ss.size(), 2);
			// Initialize s
			string s = ss[0];
			int mytrim5 = this->trim5_;
			if(color_ && s.length() > 1) {
				// This may be a primer character.  If so, keep it in the
				// 'primer' field of the read buf and parse the rest of the
				// read without it.
				int c = toupper(s[0]);
				if(asc2dnacat[c] > 0) {
					// First char is a DNA char
					int c2 = toupper(s[1]);
					// Second char is a color char
					if(asc2colcat[c2] > 0) {
						mytrim5 += 2; // trim primer and first color
					}
				}
			}
			if(color_) {
				// Convert '0'-'3' to 'A'-'T'
				for(size_t i = 0; i < s.length(); i++) {
					if(s[i] >= '0' && s[i] <= '4') {
						s[i] = "ACGTN"[(int)s[i] - '0'];
					}
					if(s[i] == '.') s[i] = 'N';
				}
			}
			if(s.length() <= (size_t)(trim3_ + mytrim5)) {
				// Entire read is trimmed away
				continue;
			} else {
				// Trim on 5' (high-quality) end
				if(mytrim5 > 0) {
					s.erase(0, mytrim5);
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
			if(vq.length() > (size_t)(trim3_ + mytrim5)) {
				// Trim on 5' (high-quality) end
				if(mytrim5 > 0) {
					vq.erase(0, mytrim5);
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
			trimmed3_.push_back(trim3_);
			trimmed5_.push_back(mytrim5);
			ostringstream os;
			os << (names_.size());
			names_.push_back(os.str());
		}
		assert_eq(v_.size(), quals_.size());
	}
	
	virtual ~VectorPatternSource() { }
	
	virtual bool nextRead(ReadBuf& r) {
		// Let Strings begin at the beginning of the respective bufs
		r.reset();
		ThreadSafe ts(&mutex_m);
		if(cur_ >= v_.size()) {
			// Clear all the Strings, as a signal to the caller that
			// we're out of reads
			r.clearAll();
			assert(r.empty());
			return false;
		}
		// Copy v_*, quals_* strings into the respective Strings
		r.color = color_;
		r.patFw  = v_[cur_];
		r.qual = quals_[cur_];
		r.trimmed3 = trimmed3_[cur_];
		r.trimmed5 = trimmed5_[cur_];
		ostringstream os;
		os << cur_;
		r.name = os.str();
		cur_++;
		readCnt_++;
		r.patid = (uint32_t)readCnt_;
		return true;
	}
	
	/**
	 * This is unused, but implementation is given for completeness.
	 */
	virtual pair<bool, bool> nextReadPair(ReadBuf& ra, ReadBuf& rb) {
		// Let Strings begin at the beginning of the respective bufs
		throw 1;
		return make_pair(false, false);
//		ra.reset();
//		rb.reset();
//		if(!paired_) {
//			paired_ = true;
//			cur_ <<= 1;
//		}
//		ThreadSafe ts(&mutex_m);
//		if(cur_ >= v_.size()-1) {
//			// Clear all the Strings, as a signal to the caller that
//			// we're out of reads
//			ra.clearAll();
//			rb.clearAll();
//			assert(ra.empty());
//			assert(rb.empty());
//			return;
//		}
//		// Copy v_*, quals_* strings into the respective Strings
//		ra.patFw  = v_[cur_];
//		ra.qual = quals_[cur_];
//		ra.trimmed3 = trimmed3_[cur_];
//		ra.trimmed5 = trimmed5_[cur_];
//		cur_++;
//		rb.patFw  = v_[cur_];
//		rb.qual = quals_[cur_];
//		rb.trimmed3 = trimmed3_[cur_];
//		rb.trimmed5 = trimmed5_[cur_];
//		ostringstream os;
//		os << readCnt_;
//		ra.name = os.str();
//		rb.name = os.str();
//		ra.color = rb.color = color_;
//		cur_++;
//		readCnt_++;
//		paired = true;
//		ra.patid = rb.patid = (uint32_t)readCnt_;
	}
	
	virtual void reset() {
		TrimmingPatternSource::reset();
		cur_ = skip_;
		paired_ = false;
	}
	
	virtual void finalize(ReadBuf& r) const { }
	
	/**
	 * Read next batch.  However, batch concept is not very applicable for this
	 * PatternSource where all the info has already been parsed into the fields
	 * in the contsructor.  This essentially modifies the pt as though we read
	 * in some number of patterns.
	 */
	pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		bool batch_a,
		bool lock)
	{
		bool success = true;
		int nread = 0;
		pt.reset();
		ThreadSafe ts(&mutex_m, lock);
		pt.setReadId(readCnt_);
#if 0
		// TODO: set nread to min of pt.size() and total - cur_
		// TODO: implement something like following function
		pt.install_dummies(nread);
#endif
		readCnt_ += nread;
		return make_pair(success, nread);
	}
	
	/**
	 * Finishes parsing outside the critical section
	 */
	virtual bool parse(ReadBuf& ra, ReadBuf& rb) const {
		cerr << "In VectorPatternSource.parse()" << endl;
		throw 1;
		return false;
	}

private:
	bool color_;
	size_t cur_;
	uint32_t skip_;
	bool paired_;
	vector<String<Dna5> > v_;     /// forward sequences
	vector<String<char> > quals_; /// quality values parallel to v_
	vector<String<char> > names_; /// names
	vector<int> trimmed3_;        /// # bases trimmed from 3' end
	vector<int> trimmed5_;        /// # bases trimmed from 5' end
};

/**
 *
 */
class BufferedFilePatternSource : public TrimmingPatternSource {
public:
	BufferedFilePatternSource(uint32_t seed,
	                          const vector<string>& infiles,
	                          const vector<string>* qinfiles,
	                          bool randomizeQuals = false,
	                          const char *dumpfile = NULL,
	                          bool verbose = false,
	                          int trim3 = 0,
	                          int trim5 = 0,
	                          uint32_t skip = 0) :
		TrimmingPatternSource(seed, randomizeQuals,
		                      dumpfile, verbose, trim3, trim5),
		infiles_(infiles),
		filecur_(0),
		fb_(),
		qfb_(),
		skip_(skip),
		first_(true)
	{
		qinfiles_.clear();
		if(qinfiles != NULL) qinfiles_ = *qinfiles;
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size(), false);
		if(qinfiles_.size() > 0 &&
		   qinfiles_.size() != infiles_.size())
		{
			cerr << "Error: Different numbers of input FASTA/quality files ("
			     << infiles_.size() << "/" << qinfiles_.size() << ")" << endl;
			throw 1;
		}
		assert(!fb_.isOpen());
		assert(!qfb_.isOpen());
		open(); // open first file in the list
		filecur_++;
	}

	virtual ~BufferedFilePatternSource() {
		if(fb_.isOpen()) fb_.close();
		if(qfb_.isOpen()) {
			assert_gt(qinfiles_.size(), 0);
			qfb_.close();
		}
	}

	/**
	 * Fill ReadBuf with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 */
	virtual bool nextRead(ReadBuf& r) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		ThreadSafe ts(&mutex_m);
		pair<bool, bool> flags = make_pair(false, false);
		while(!flags.first && !flags.second) {
			flags = readLight(r);
		}
		if(first_ && !flags.first) {
			cerr << "Warning: Could not find any reads in \""
			     << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(!flags.first && filecur_ < infiles_.size()) {
			assert(flags.second);
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			flags.first = flags.second = false;
			while(!flags.first && !flags.second) {
				flags = readLight(r);
			}
			assert_geq(r.patid, skip_);
			if(!flags.first) {
				// No reads could be extracted from this _infile
				cerr << "Warning: Could not find any reads in \""
				     << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		return flags.first;
	}
	
	/**
	 *
	 */
	virtual pair<bool, bool> nextReadPair(ReadBuf& ra, ReadBuf& rb) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		ThreadSafe ts(&mutex_m);
		bool success, eof, paired;
		do {
			BoolTriple result = readPairLight(ra, rb);
			success = result.a;
			eof = result.b;
			paired = result.c;
		} while(!success && !eof);
		if(first_ && !success) {
			cerr << "Warning: Could not find any read pairs in \""
			     << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(!success && filecur_ < infiles_.size()) {
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			do {
				BoolTriple result = readPairLight(ra, rb);
				success = result.a;
				eof = result.b;
				paired = result.c;
			} while(!success && !eof);
			if(!success && eof) {
				cerr << "Warning: Could not find any reads in \""
				     << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		return make_pair(success, paired);
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
	
	/**
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 *
	 * What does the return value signal?
	 * In the event that we read more data, it should at least signal how many
	 * reads were read, and whether we're totally done.  It's debatable whether
	 * it should convey anything about the individual read files, like whether
	 * we finished one of them.
	 */
	pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		bool batch_a,
		bool lock)
	{
		bool done = false;
		int nread = 0;
		
		// synchronization at this level because both reading and manipulation of
		// current file pointer have to be protected
		ThreadSafe ts(&mutex_m, lock);
		pt.setReadId(readCnt_);
		while(true) { // loop that moves on to next file when needed
			do {
				pair<bool, int> ret = nextBatchFromFile(pt, batch_a);
				done = ret.first;
				nread = ret.second;
			} while(!done && nread == 0); // not sure why this would happen
			if(done && filecur_ < infiles_.size()) { // finished with this file
				open();
				resetForNextFile(); // reset state to handle a fresh file
				filecur_++;
				if(nread == 0) {
					continue;
				}
			}
			break;
		}
		assert_geq(nread, 0);
		readCnt_ += nread;
		return make_pair(done, nread);
	}
	
protected:

	/// Read another pattern from the input file; this is overridden
	/// to deal with specific file formats
	virtual pair<bool, bool> readLight(ReadBuf& r) = 0;
	
	/**
	 * Light-parse a batch of unpaired reads from current file into the given
	 * buffer.  Called from CFilePatternSource.nextBatch().
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a) = 0;
	
	/// Read another pattern pair from the input file; this is
	/// overridden to deal with specific file formats
	virtual BoolTriple readPairLight(
		ReadBuf& ra,
		ReadBuf& rb) = 0;
	
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() { }
	
	void open() {
		if(fb_.isOpen()) fb_.close();
		if(qfb_.isOpen()) qfb_.close();
		while(filecur_ < infiles_.size()) {
			// Open read
			FILE *in;
			if(infiles_[filecur_] == "-") {
				in = stdin;
			} else if((in = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
				if(!errs_[filecur_]) {
					cerr << "Warning: Could not open read file \"" << infiles_[filecur_] << "\" for reading; skipping..." << endl;
					errs_[filecur_] = true;
				}
				filecur_++;
				continue;
			}
			fb_.newFile(in);
			// Open quality
			if(!qinfiles_.empty()) {
				FILE *in;
				if(qinfiles_[filecur_] == "-") {
					in = stdin;
				} else if((in = fopen(qinfiles_[filecur_].c_str(), "rb")) == NULL) {
					if(!errs_[filecur_]) {
						cerr << "Warning: Could not open quality file \"" << qinfiles_[filecur_] << "\" for reading; skipping..." << endl;
						errs_[filecur_] = true;
					}
					filecur_++;
					continue;
				}
				qfb_.newFile(in);
			}
			return;
		}
		throw 1;
	}
	
	vector<string> infiles_; /// filenames for read files
	vector<string> qinfiles_; /// filenames for quality files
	vector<bool> errs_; /// whether we've already printed an error for each file
	size_t filecur_;   /// index into infiles_ of next file to read
	FileBuf fb_;  /// read file currently being read from
	FileBuf qfb_; /// quality file currently being read from
	uint32_t skip_;     /// number of reads to skip
	bool first_;
};

/**
 *
 */
class CFilePatternSource : public TrimmingPatternSource {
public:
	CFilePatternSource(
	    uint32_t seed,
	    const vector<string>& infiles,
	    const vector<string>* qinfiles,
	    bool randomizeQuals = false,
	    const char *dumpfile = NULL,
	    bool verbose = false,
	    int trim3 = 0,
	    int trim5 = 0,
	    uint32_t skip = 0) :
		TrimmingPatternSource(seed, randomizeQuals,
		                      dumpfile, verbose, trim3, trim5),
		infiles_(infiles),
		filecur_(0),
		fp_(NULL),
		qfp_(NULL),
		is_open_(false),
		skip_(skip),
		first_(true)
	{
		qinfiles_.clear();
		if(qinfiles != NULL) qinfiles_ = *qinfiles;
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size(), false);
		if(qinfiles_.size() > 0 &&
		   qinfiles_.size() != infiles_.size())
		{
			cerr << "Error: Different numbers of input FASTA/quality files ("
			     << infiles_.size() << "/" << qinfiles_.size() << ")" << endl;
			throw 1;
		}
		open(); // open first file in the list
		filecur_++;
	}

	virtual ~CFilePatternSource() {
		if(is_open_) {
			assert(fp_ != NULL);
			fclose(fp_);
			fp_ = NULL;
			if(qfp_ != NULL) {
				fclose(qfp_);
				qfp_ = NULL;
			}
		}
	}

	/**
	 * Fill ReadBuf with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 */
	virtual bool nextRead(ReadBuf& r) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		ThreadSafe ts(&mutex_m);
		pair<bool, bool> flags = make_pair(false, false);
		while(!flags.first && !flags.second) {
			 flags = readLight(r);
		}
		if(first_ && !flags.first) {
			cerr << "Warning: Could not find any reads in \""
			     << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(!flags.first && filecur_ < infiles_.size()) {
			assert(flags.second);
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			flags.first = flags.second = false;
			while(!flags.first && !flags.second) {
				 flags = readLight(r);
			}
			assert_geq(r.patid, skip_);
			if(!flags.first) {
				cerr << "Warning: Could not find any reads in \""
				     << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		return flags.first;
	}
	
	/**
	 *
	 */
	virtual pair<bool, bool> nextReadPair(ReadBuf& ra, ReadBuf& rb) {
		// We are entering a critical region, because we're
		// manipulating our file handle and filecur_ state
		ThreadSafe ts(&mutex_m);
		bool success, eof, paired;
		do {
			BoolTriple result = readPairLight(ra, rb);
			success = result.a;
			eof = result.b;
			paired = result.c;
		} while(!success && !eof);
		if(first_ && !success) {
			cerr << "Warning: Could not find any read pairs in \""
			     << infiles_[0] << "\"" << endl;
		}
		first_ = false;
		while(!success && filecur_ < infiles_.size()) {
			// Open next file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			do {
				BoolTriple result = readPairLight(ra, rb);
				success = result.a;
				eof = result.b;
				paired = result.c;
			} while(!success && !eof);
			if(!success && eof) {
				cerr << "Warning: Could not find any reads in \""
				     << infiles_[filecur_] << "\"" << endl;
			}
			filecur_++;
		}
		return make_pair(success, paired);
	}

	/**
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 *
	 * Returns pair<bool, int> where bool indicates whether we're
	 * completely done, and int indicates how many reads were read.
	 */
	pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		bool batch_a,
		bool lock)
	{
		bool done = false;
		int nread = 0;
		
		// synchronization at this level because both reading and manipulation of
		// current file pointer have to be protected
		ThreadSafe ts(&mutex_m, lock);
		pt.setReadId(readCnt_);
		while(true) { // loop that moves on to next file when needed
			do {
				pair<bool, int> ret = nextBatchFromFile(pt, batch_a);
				done = ret.first;
				nread = ret.second;
			} while(!done && nread == 0); // not sure why this would happen
			if(done && filecur_ < infiles_.size()) { // finished with this file
				open();
				resetForNextFile(); // reset state to handle a fresh file
				filecur_++;
				if(nread == 0) {
					continue;
				}
			}
			break;
		}
		assert_geq(nread, 0);
		readCnt_ += nread;
		return make_pair(done, nread);
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

	/**
	 * Light-parse a batch of unpaired reads from current file into the given
	 * buffer.  Called from CFilePatternSource.nextBatch().
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a) = 0;
	
	/// Read another pattern from the input file; this is overridden
	/// to deal with specific file formats
	virtual pair<bool, bool> readLight(ReadBuf& r) = 0;
	
	
	/// Read another pattern pair from the input file; this is
	/// overridden to deal with specific file formats
	virtual BoolTriple readPairLight(ReadBuf& ra, ReadBuf& rb) = 0;
	
	/// Reset state to handle a fresh file
	virtual void resetForNextFile() { }
	
	void open() {
		if(is_open_) {
			is_open_ = false;
			fclose(fp_);
			fp_ = NULL;
			if(qfp_ != NULL) {
				fclose(qfp_);
				qfp_ = NULL;
			}
		}
		while(filecur_ < infiles_.size()) {
			// Open read
			if(infiles_[filecur_] == "-") {
				fp_ = stdin;
			} else if((fp_ = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
				if(!errs_[filecur_]) {
					cerr << "Warning: Could not open read file \"" << infiles_[filecur_] << "\" for reading; skipping..." << endl;
					errs_[filecur_] = true;
				}
				filecur_++;
				continue;
			}
			is_open_ = true;
			setvbuf(fp_, buf_, _IOFBF, 64*1024);
			if(!qinfiles_.empty()) {
				if(qinfiles_[filecur_] == "-") {
					qfp_ = stdin;
				} else if((qfp_ = fopen(qinfiles_[filecur_].c_str(), "rb")) == NULL) {
					if(!errs_[filecur_]) {
						cerr << "Warning: Could not open quality file \"" << qinfiles_[filecur_] << "\" for reading; skipping..." << endl;
						errs_[filecur_] = true;
					}
					filecur_++;
					continue;
				}
				assert(qfp_ != NULL);
				setvbuf(qfp_, qbuf_, _IOFBF, 64*1024);
			}
			return;
		}
		throw 1;
	}
	
	vector<string> infiles_; /// filenames for read files
	vector<string> qinfiles_; /// filenames for quality files
	vector<bool> errs_; /// whether we've already printed an error for each file
	size_t filecur_;   /// index into infiles_ of next file to read
	FILE *fp_; /// read file currently being read from
	FILE *qfp_; /// quality file currently being read from
	bool is_open_; /// whether fp_ is currently open
	uint32_t skip_;     /// number of reads to skip
	bool first_;
	char buf_[64*1024]; /// file buffer for sequences
	char qbuf_[64*1024]; /// file buffer for qualities
};

/**
 * Parse a single quality string from fb and store qualities in r.
 * Assume the next character obtained via fb.get() is the first
 * character of the quality string.  When returning, the next
 * character returned by fb.peek() or fb.get() should be the first
 * character of the following line.
 */
int parseQuals(ReadBuf& r,
               FileBuf& fb,
               int readLen,
               int trim3,
               int trim5,
               bool intQuals,
               bool phred64,
               bool solexa64);

/**
 * Synchronized concrete pattern source for a list of FASTA or CSFASTA
 * (if color = true) files.
 */
class FastaPatternSource : public BufferedFilePatternSource {
public:
	FastaPatternSource(uint32_t seed,
	                   const vector<string>& infiles,
	                   const vector<string>* qinfiles,
	                   bool color,
	                   bool randomizeQuals,
	                   const char *dumpfile = NULL,
	                   bool verbose = false,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool solexa64 = false,
	                   bool phred64 = false,
	                   bool intQuals = false,
	                   uint32_t skip = 0) :
		BufferedFilePatternSource(seed, infiles, qinfiles, randomizeQuals,
		                          dumpfile, verbose, trim3,
		                          trim5, skip),
		first_(true), color_(color), solexa64_(solexa64),
		phred64_(phred64), intQuals_(intQuals)
	{ }
	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}


protected:

	/// Read another pattern from a FASTA input file
	pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a)
	{
		throw 1;
		return make_pair(false, 0);
	}

	/// Read another pattern from a FASTA input file
	bool parse(ReadBuf& r, ReadBuf& rb) const {
		r.reset();
		cerr << "In FastaPatternSource.parse()" << endl;
		throw 1;
		return false;
	}

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

	/// Called when we have to bail without having parsed a read.
	void bail(ReadBuf& r) {
		r.clearAll();
		fb_.resetLastN();
		qfb_.resetLastN();
	}

	/// Read another pattern from a FASTA input file
	virtual pair<bool, bool> readLight(ReadBuf& r) {
		int c, qc = 0;
		int dstLen = 0;
		int nameLen = 0;
		bool doquals = qinfiles_.size() > 0;
		assert(!doquals || qfb_.isOpen());
		assert(fb_.isOpen());
		r.color = color_;

		// Skip over between-read comments.  Note that SOLiD's comments use #s
		c = fb_.get();
		if(c < 0) {
			bail(r);
			return make_pair(false, true);
		}
		while(c == '#' || c == ';') {
			c = fb_.peekUptoNewline();
			fb_.resetLastN();
			c = fb_.get();
		}
		assert_eq(1, fb_.lastNCur());
		if(doquals) {
			qc = qfb_.get();
			if(qc < 0) {
				bail(r);
				return make_pair(false, true);
			}
			while(qc == '#' || qc == ';') {
				qc = qfb_.peekUptoNewline();
				qfb_.resetLastN();
				qc = qfb_.get();
			}
			assert_eq(1, qfb_.lastNCur());
		}

		// Pick off the first carat
		if(first_) {
			if(c != '>') {
				cerr << "Error: reads file does not look like a FASTA file" << endl;
				throw 1;
			}
			if(doquals && qc != '>') {
				cerr << "Error: quality file does not look like a FASTA quality file" << endl;
				throw 1;
			}
			first_ = false;
		}
		assert_eq('>', c);
		if(doquals) assert_eq('>', qc);
		c = fb_.get(); // get next char after '>'
		if(doquals) qc = qfb_.get();

		// Read to the end of the id line, sticking everything after the '>'
		// into *name
		bool warning = false;
		while(true) {
			if(c < 0 || qc < 0) {
				bail(r);
				return make_pair(false, true);
			}
			if(c == '\n' || c == '\r') {
				// Break at end of line, after consuming all \r's, \n's
				while(c == '\n' || c == '\r') {
					c = fb_.get();
					if(doquals) qc = qfb_.get();
					if(c < 0 || qc < 0) { bail(r); return make_pair(false, true); }
				}
				break;
			}
			if(doquals && c != qc) {
				cerr << "Warning: one or more mismatched read names between FASTA and quality files" << endl;
				warning = true;
			}
			r.nameBuf[nameLen++] = c;
			c = fb_.get();
			if(doquals) qc = qfb_.get();
		}
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, nameLen);
		if(warning) {
			cerr << "         Offending read name: \"" << r.name << "\"" << endl;
		}

		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		int begin = 0;
		int mytrim5 = this->trim5_;
		if(color_) {
			// This is the primer character, keep it in the
			// 'primer' field of the read buf and keep parsing
			c = toupper(c);
			if(asc2dnacat[c] > 0) {
				// First char is a DNA char
				int c2 = toupper(fb_.peek());
				if(asc2colcat[c2] > 0) {
					// Second char is a color char
					r.primer = c;
					r.trimc = c2;
					mytrim5 += 2;
				}
			}
			if(c < 0) {
				bail(r);
				return make_pair(false, true);
			}
		}
		while(c != '>' && c >= 0) {
			if(color_) {
				if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
				if(c == '.') c = 'N';
			}
			if(asc2dnacat[c] > 0 && begin++ >= mytrim5) {
				if(dstLen + 1 > 1024) tooManySeqChars(r.name);
				r.patBufFw[dstLen] = charToDna5[c];
				if(!doquals) r.qualBuf[dstLen]  = 'I';
				dstLen++;
			}
			if(fb_.peek() == '>') {
				break;
			}
			c = fb_.get();
		}
		dstLen -= this->trim3_;
		if(dstLen < 0) dstLen = 0;
		_setBegin (r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, dstLen);
		r.trimmed3 = this->trim3_;
		r.trimmed5 = mytrim5;
		if(doquals) {
			if(dstLen > 0) {
				parseQuals(r, qfb_, dstLen + r.trimmed3 + r.trimmed5,
						   r.trimmed3, r.trimmed5, intQuals_, phred64_,
						   solexa64_);
			} else {
				// Bail
				qfb_.peekUptoNewline();
				qfb_.get();
				qfb_.resetLastN();
				fb_.resetLastN();
				// Count the read
				readCnt_++;
				r.patid = (uint32_t)readCnt_-1;
				return make_pair(true, false);
			}
		}
		_setBegin (r.qual,  r.qualBuf);
		_setLength(r.qual,  dstLen);
		// Set up a default name if one hasn't been set
		if(nameLen == 0) {
			itoa10((int)readCnt_, r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			nameLen = (int)strlen(r.nameBuf);
			_setLength(r.name, nameLen);
		}
		assert_gt(nameLen, 0);
		readCnt_++;
		r.patid = (uint32_t)(readCnt_-1);
		r.readOrigBufLen = fb_.copyLastN(r.readOrigBuf);
		fb_.resetLastN();
		if(doquals) {
			r.qualOrigBufLen = qfb_.copyLastN(r.qualOrigBuf);
			qfb_.resetLastN();
			if(false) {
				cout << "Name: " << r.name << endl
					 << " Seq: " << r.patFw << " (" << seqan::length(r.patFw) << ")" << endl
					 << "Qual: " << r.qual  << " (" << seqan::length(r.qual) << ")" << endl
					 << "Orig seq:" << endl;
				for(size_t i = 0; i < r.readOrigBufLen; i++) cout << r.readOrigBuf[i];
				cout << "Orig qual:" << endl;
				for(size_t i = 0; i < r.qualOrigBufLen; i++) cout << r.qualOrigBuf[i];
				cout << endl;
			}
		}
		return make_pair(true, false);
	}
	
	/// Read another pair of patterns from a FASTA input file
	virtual BoolTriple readPairLight(ReadBuf& ra, ReadBuf& rb) {
		// (For now, we shouldn't ever be here)
		cerr << "In FastaPatternSource.readPair()" << endl;
		throw 1;
		return BoolTriple();
	}

	/**
	 * Finishes parsing outside the critical section
	 */
	virtual void finalize(ReadBuf& r) const { }
	
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
	bool color_;
	bool solexa64_;
	bool phred64_;
	bool intQuals_;
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

/**
 * Synchronized concrete pattern source for a list of files with tab-
 * delimited name, seq, qual fields (or, for paired-end reads,
 * basename, seq1, qual1, seq2, qual2).
 */
class TabbedPatternSource : public BufferedFilePatternSource {
public:
	TabbedPatternSource(uint32_t seed,
	                    const vector<string>& infiles,
	                    bool color,
	                    bool randomizeQuals = false,
	                    const char *dumpfile = NULL,
	                    bool verbose = false,
	                    int trim3 = 0,
	                    int trim5 = 0,
	                    bool solQuals = false,
	                    bool phred64Quals = false,
	                    bool intQuals = false,
	                    uint32_t skip = 0) :
		BufferedFilePatternSource(seed, infiles, NULL, randomizeQuals,
		                          dumpfile, verbose,
		                          trim3, trim5, skip),
		color_(color),
		solQuals_(solQuals),
		phred64Quals_(phred64Quals),
		intQuals_(intQuals)
	{ }

protected:
	
	/**
	 * "Light" parser.  This is inside the critical section, so the key is to do
	 * just enough parsing so that another function downstream (finalize()) can do
	 * the rest of the parsing.  Really this function's only job is to stick every
	 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
	 * then parses the contents of r.readOrigBuf later.
	 */
	pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a)
	{
		throw 1;
		return make_pair(false, 0);
	}

	/// Read another pattern from a FASTA input file
	bool parse(ReadBuf& r, ReadBuf& rb) const {
		r.reset();
		cerr << "In TabbedPatternSource.parse()" << endl;
		throw 1;
		return false;
	}

	/// Read another pattern from a FASTA input file
	virtual pair<bool, bool> readLight(ReadBuf& r) {
		r.color = color_;
		int trim5 = this->trim5_;
		// fb_ is about to dish out the first character of the
		// name field
		if(parseName(r, NULL, '\t') == -1) {
			peekOverNewline(fb_); // skip rest of line
			r.clearAll();
			return make_pair(false, true);
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// sequence field
		int charsRead = 0;
		int dstLen = parseSeq(r, charsRead, trim5, '\t');
		assert_neq('\t', fb_.peek());
		if(dstLen <= 0) {
			peekOverNewline(fb_); // skip rest of line
			r.clearAll();
			return make_pair(false, true);
		}

		// fb_ is about to dish out the first character of the
		// quality-string field
		char ct = 0;
		if(parseQuals(r, charsRead, dstLen, trim5, ct, '\n') <= 0) {
			peekOverNewline(fb_); // skip rest of line
			r.clearAll();
			return make_pair(false, true);
		}
		r.trimmed3 = this->trim3_;
		r.trimmed5 = trim5;
		assert_eq(ct, '\n');
		assert_neq('\n', fb_.peek());
		r.readOrigBufLen = fb_.copyLastN(r.readOrigBuf);
		fb_.resetLastN();
		// The last character read in parseQuals should have been a
		// '\n'

		readCnt_++;
		r.patid = (uint32_t)(readCnt_-1);
		return make_pair(true, false);
	}

	/// Read another pair of patterns from a FASTA input file
	virtual BoolTriple readPairLight(ReadBuf& ra, ReadBuf& rb) {
		// fb_ is about to dish out the first character of the
		// name field
		int mytrim5_1 = this->trim5_;
		if(parseName(ra, &rb, '\t') == -1) {
			peekOverNewline(fb_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			fb_.resetLastN();
			// !success, eof, !paired
			return BoolTriple(false, true, false);
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// sequence field for the first mate
		int charsRead1 = 0;
		int dstLen1 = parseSeq(ra, charsRead1, mytrim5_1, '\t');
		if(dstLen1 <= -1) {
			peekOverNewline(fb_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			fb_.resetLastN();
			// !success, eof, !paired
			return BoolTriple(false, true, false);
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// quality-string field
		char ct = 0;
		if(parseQuals(ra, charsRead1, dstLen1, mytrim5_1, ct, '\t', '\n') <= 0) {
			peekOverNewline(fb_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			fb_.resetLastN();
			// !success, eof, !paired
			return BoolTriple(false, true, false);
		}
		ra.trimmed3 = this->trim3_;
		ra.trimmed5 = mytrim5_1;
		assert(ct == '\t' || ct == '\n');
		if(ct == '\n') {
			// Unpaired record; return.
			rb.clearAll();
			peekOverNewline(fb_);
			ra.readOrigBufLen = fb_.copyLastN(ra.readOrigBuf);
			fb_.resetLastN();
			readCnt_++;
			ra.patid = (uint32_t)(readCnt_-1);
			// success, !eof, !paired
			return BoolTriple(true, false, false);
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// sequence field for the second mate
		int charsRead2 = 0;
		int mytrim5_2 = this->trim5_;
		int dstLen2 = parseSeq(rb, charsRead2, mytrim5_2, '\t');
		if(dstLen2 <= 0) {
			peekOverNewline(fb_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			fb_.resetLastN();
			// !success, !eof, !paired
			return BoolTriple(false, false, false);
		}
		assert_neq('\t', fb_.peek());

		// fb_ is about to dish out the first character of the
		// quality-string field
		if(parseQuals(rb, charsRead2, dstLen2, mytrim5_2, ct, '\n') <= 0) {
			peekOverNewline(fb_); // skip rest of line
			ra.clearAll();
			rb.clearAll();
			fb_.resetLastN();
			// !success, !eof, !paired
			return BoolTriple(false, false, false);
		}
		assert_eq('\n', ct);
		if(fb_.peek() == '\n') {
			assert(false);
		}
		peekOverNewline(fb_);
		ra.readOrigBufLen = fb_.copyLastN(ra.readOrigBuf);
		fb_.resetLastN();

		rb.trimmed3 = this->trim3_;
		rb.trimmed5 = mytrim5_2;

		// The last character read in parseQuals should have been a
		// '\n'

		readCnt_++;
		ra.patid = rb.patid = (uint32_t)(readCnt_-1);
		// success, !eof, paired
		return BoolTriple(true, false, true);
	}

	/**
	 * Finishes parsing outside the critical section
	 */
	virtual void finalize(ReadBuf& r) const { }
	
	
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
	 * Parse a name from fb_ and store in r.  Assume that the next
	 * character obtained via fb_.get() is the first character of
	 * the sequence and the string stops at the next char upto (could
	 * be tab, newline, etc.).
	 */
	int parseName(ReadBuf& r, ReadBuf* r2, char upto = '\t') {
		// Read the name out of the first field
		int c = 0;
		int nameLen = 0;
		while(true) {
			if((c = fb_.get()) < 0) {
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
			itoa10((int)readCnt_, r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			nameLen = (int)strlen(r.nameBuf);
			_setLength(r.name, nameLen);
			if(r2 != NULL) {
				itoa10((int)readCnt_, (*r2).nameBuf);
				_setBegin((*r2).name, (*r2).nameBuf);
				_setLength((*r2).name, nameLen);
			}
		}
		assert_gt(nameLen, 0);
		return nameLen;
	}

	/**
	 * Parse a single sequence from fb_ and store in r.  Assume
	 * that the next character obtained via fb_.get() is the first
	 * character of the sequence and the sequence stops at the next
	 * char upto (could be tab, newline, etc.).
	 */
	int parseSeq(ReadBuf& r, int& charsRead, int& trim5, char upto = '\t') {
		int begin = 0;
		int dstLen = 0;
		int c = fb_.get();
		assert(c != upto);
		r.color = color_;
		if(color_) {
			// This may be a primer character.  If so, keep it in the
			// 'primer' field of the read buf and parse the rest of the
			// read without it.
			c = toupper(c);
			if(asc2dnacat[c] > 0) {
				// First char is a DNA char
				int c2 = toupper(fb_.peek());
				// Second char is a color char
				if(asc2colcat[c2] > 0) {
					r.primer = c;
					r.trimc = c2;
					trim5 += 2; // trim primer and first color
				}
			}
			if(c < 0) { return -1; }
		}
		while(c != upto) {
			if(color_) {
				if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
			}
			if(c == '.') c = 'N';
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(begin++ >= trim5) {
					assert_neq(0, dna4Cat[c]);
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
			if((c = fb_.get()) < 0) {
				return -1;
			}
		}
		dstLen -= this->trim3_;
		_setBegin (r.patFw,  (Dna5*)r.patBufFw);
		_setLength(r.patFw,  dstLen);
		return dstLen;
	}

	/**
	 * Parse a single quality string from fb_ and store in r.
	 * Assume that the next character obtained via fb_.get() is
	 * the first character of the quality string and the string stops
	 * at the next char upto (could be tab, newline, etc.).
	 */
	int parseQuals(ReadBuf& r, int charsRead, int dstLen, int trim5,
	               char& c2, char upto = '\t', char upto2 = -1)
	{
		int qualsRead = 0;
		int c = 0;
		if (intQuals_) {
			char buf[4096];
			while (qualsRead < charsRead) {
				vector<string> s_quals;
				if(!tokenizeQualLine(fb_, buf, 4096, s_quals)) break;
				for (unsigned int j = 0; j < s_quals.size(); ++j) {
					char c = intToPhred33(atoi(s_quals[j].c_str()), solQuals_);
					assert_geq(c, 33);
					if (qualsRead >= trim5) {
						size_t off = qualsRead - trim5;
						if(off >= 1024) tooManyQualities(r.name);
						r.qualBuf[off] = c;
					}
					++qualsRead;
				}
			} // done reading integer quality lines
			if (charsRead > qualsRead) tooFewQualities(r.name);
		} else {
			// Non-integer qualities
			while((qualsRead < dstLen + trim5) && c >= 0) {
				c = fb_.get();
				c2 = c;
				if (c == ' ') wrongQualityFormat(r.name);
				if(c < 0) {
					// EOF occurred in the middle of a read - abort
					return -1;
				}
				if(!isspace(c) && c != upto && (upto2 == -1 || c != upto2)) {
					if (qualsRead >= trim5) {
						size_t off = qualsRead - trim5;
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
			if(qualsRead != dstLen + trim5) {
				assert(false);
			}
		}
		_setBegin (r.qual, (char*)r.qualBuf);
		_setLength(r.qual, dstLen);
		while(c != upto && (upto2 == -1 || c != upto2)) {
			c = fb_.get();
			c2 = c;
		}
		return qualsRead;
	}

	bool color_;
	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
};

/**
 * Synchronized concrete pattern source for a list of FASTA files where
 * reads need to be extracted from long continuous sequences.
 */
class FastaContinuousPatternSource : public BufferedFilePatternSource {
public:
	FastaContinuousPatternSource(
			uint32_t seed,
			const vector<string>& infiles,
			size_t length,
			size_t freq,
			const char *dumpfile = NULL,
			bool verbose = false,
			uint32_t skip = 0) :
		BufferedFilePatternSource(seed, infiles, NULL, false,
		                          dumpfile, verbose, 0, 0, skip),
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
	
	/**
	 * "Light" parser.  This is inside the critical section, so the key is to do
	 * just enough parsing so that another function downstream (finalize()) can do
	 * the rest of the parsing.  Really this function's only job is to stick every
	 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
	 * then parses the contents of r.readOrigBuf later.
	 */
	pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a)
	{
		throw 1;
		return make_pair(false, 0);
	}

	/// Read another pattern from a FASTA input file
	bool parse(ReadBuf& r, ReadBuf& rb) const {
		assert(r.empty());
		assert(rb.empty());
		cerr << "In FastaContinuousPatternSource.parse()" << endl;
		throw 1;
	}

	/// Read another pattern from a FASTA input file
	virtual pair<bool, bool> readLight(ReadBuf& r) {
		while(true) {
			int c = fb_.get();
			if(c < 0) {
				seqan::clear(r.patFw);
				return make_pair(false, true);
			}
			if(c == '>') {
				resetForNextFile();
				c = fb_.peek();
				bool sawSpace = false;
				while(c != '\n' && c != '\r') {
					if(!sawSpace) {
						sawSpace = isspace(c);
					}
					if(!sawSpace) {
						nameBuf_[nameChars_++] = c;
					}
					fb_.get();
					c = fb_.peek();
				}
				while(c == '\n' || c == '\r') {
					fb_.get();
					c = fb_.peek();
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
					itoa10((int)(readCnt_ - subReadCnt_), &r.nameBuf[nameChars_]);
					_setBegin(r.name, r.nameBuf);
					_setLength(r.name, strlen(r.nameBuf));
					eat_ = freq_-1;
					readCnt_++;
					beginning_ = false;
					r.patid = (uint32_t)(readCnt_-1);
					break;
				}
			}
		}
		return make_pair(true, false);
	}
	/// Shouldn't ever be here; it's not sensible to obtain read pairs
	// from a continuous input.
	virtual BoolTriple readPairLight(ReadBuf& ra, ReadBuf& rb) {
		cerr << "In FastaContinuousPatternSource.readPair()" << endl;
		throw 1;
		return BoolTriple();
	}

	/**
	 * Finishes parsing outside the critical section
	 */
	virtual void finalize(ReadBuf& r) const { }
	

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
class FastqPatternSource : public CFilePatternSource {
public:
	FastqPatternSource(uint32_t seed,
	                   const vector<string>& infiles,
	                   bool color,
	                   bool randomizeQuals = false,
	                   const char *dumpfile = NULL,
	                   bool verbose = false,
	                   int trim3 = 0,
	                   int trim5 = 0,
	                   bool solexa_quals = false,
	                   bool phred64Quals = false,
	                   bool integer_quals = false,
	                   uint32_t skip = 0) :
		CFilePatternSource(seed, infiles, NULL, randomizeQuals,
		                          dumpfile, verbose,
		                          trim3, trim5, skip),
		first_(true),
		solQuals_(solexa_quals),
		phred64Quals_(phred64Quals),
		intQuals_(integer_quals),
		color_(color)
	{ }
	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

protected:
	
	/**
	 * "Light" parser.  This is inside the critical section, so the key is to do
	 * just enough parsing so that another function downstream (finalize()) can do
	 * the rest of the parsing.  Really this function's only job is to stick every
	 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
	 * then parses the contents of r.readOrigBuf later.
	 */
	pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a)
	{
		int c = 0;
		vector<ReadBuf>& readBuf = batch_a ? pt.bufa_ : pt.bufb_;
		size_t& len = readBuf[0].readOrigBufLen;
		len = 0;
		if(first_) {
			c = getc_unlocked(fp_);
			while(c == '\r' || c == '\n') {
				c = getc_unlocked(fp_);
			}
			if(c != '@') {
				cerr << "Error: reads file does not look like a FASTQ file" << endl;
				throw 1;
			}
			readBuf[0].readOrigBuf[len++] = c;
			assert_eq('@', c);
			first_ = false;
		}
		// Note: to reduce the number of times we have to enter the critical
		// section (each entrance has some assocaited overhead), we could populate
		// the buffer with several reads worth of data here, instead of just one.
		bool done = false, aborted = false;
		size_t readi = 0;
		// Read until we run out of input or until we've filled the buffer
		for(; readi < pt.max_buf_ && !done; readi++) {
			char* buf = readBuf[readi].readOrigBuf;
			assert(readi == 0 || strlen(buf) == 0);
			len = readBuf[readi].readOrigBufLen;
			len = 0;
			int newlines = 4;
			while(newlines) {
				c = getc_unlocked(fp_);
				done = c < 0;
				if(c == '\n' || (done && newlines == 1)) {
					newlines--;
					c = '\n';
				} else if(done) {
					//return make_pair(false, true);
					aborted = true; //Unexpected EOF
					break;
				}
				buf[len++] = c;
			}
			readCnt_++;
			readBuf[readi].patid = readCnt_-1;
		}
		if(aborted) {
			readi--;
			readCnt_--;
		}
		return make_pair(done, readi);
		//return make_pair(false, 0);
	}

	/// Read another pattern from a FASTQ input file
	bool parse(ReadBuf& r, ReadBuf& rb) const {
		assert(strlen(r.readOrigBuf) != 0);
		assert(r.empty());
		int c;
		size_t name_len = 0;
		size_t cur = 1;
		r.color = color_;
		r.primer = -1;

		// Parse read name
		while(true) {
			assert(cur < r.readOrigBufLen);
			c = r.readOrigBuf[cur++];
			if(c == '\n' || c == '\r') {
				do {
					c = r.readOrigBuf[cur++];
				} while(c == '\n' || c == '\r');
				break;
			}
			r.nameBuf[name_len++] = c;
		}
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, name_len);

		if(color_) {
			// May be a primer character.  If so, keep it 'primer' field
			// of read buf and parse the rest of the read without it.
	//		c = toupper(c);
	//		if(asc2dnacat[c] > 0) {
	//			// First char is a DNA char
	//			int c2 = toupper(fb_.peek());
	//			// Second char is a color char
	//			if(asc2colcat[c2] > 0) {
	//				r.primer = c;
	//				r.trimc = c2;
	//				mytrim5 += 2; // trim primer and first color
	//			}
	//		}
	//		if(c < 0) { bail(r); return; }
			throw 1;
		}
		
		// Parse sequence
		int nchar = 0;
		uint8_t *seqbuf = r.patBufFw;
		while(c != '+') {
			if(c == '.') {
				c = 'N';
			}
			if(color_ && c >= '0' && c <= '4') {
				c = "ACGTN"[(int)c - '0'];
			}
			if(isalpha(c)) {
				// If it's past the 5'-end trim point
				if(nchar++ >= trim5_) {
					*seqbuf++ = charToDna5[c];
				}
			}
			assert(cur < r.readOrigBufLen);
			c = r.readOrigBuf[cur++];
		}
		int seq_len = (int)(seqbuf - r.patBufFw);
		r.trimmed5 = (int)(nchar - seq_len);
		r.trimmed3 = min(seq_len, trim3_);
		seq_len = max(seq_len - trim3_, 0);
		_setBegin(r.patFw, (Dna5*)r.patBufFw);
		_setLength(r.patFw, seq_len);
		
		assert_eq('+', c);
		do {
			assert(cur < r.readOrigBufLen);
			c = r.readOrigBuf[cur++];
		} while(c != '\n' && c != '\r');
		do {
			assert(cur < r.readOrigBufLen);
			c = r.readOrigBuf[cur++];
		} while(c == '\n' || c == '\r');
		
		// Now we're on the next non-blank line after the + line
		if(seq_len == 0) {
			return true; // done parsing empty read
		}

		int nqual = 0;
		char *qualbuf = r.qualBuf;
		if (intQuals_) {
			// TODO: must implement this for compatibility with other Bowtie
			throw 1; // not yet implemented
		} else {
			while(c != '\r' && c != '\n') {
				c = charToPhred33(c, solQuals_, phred64Quals_);
				if(c == ' ') {
					wrongQualityFormat(r.name);
					return false;
				}
				if(nqual++ >= r.trimmed5) {
					*qualbuf++ = c;
				}
				c = r.readOrigBuf[cur++];
			}
			int qual_len = (int)(qualbuf - r.qualBuf);
			qual_len = max(0, qual_len - r.trimmed3);
			if(qual_len < seq_len) {
				tooFewQualities(r.name);
				return false;
			} else if(qual_len > seq_len+1) {
				tooManyQualities(r.name);
				return false;
			}
			_setBegin(r.qual, (char*)r.qualBuf);
			_setLength(r.qual, seq_len);
		}
		// Set up a default name if one hasn't been set
		if(name_len == 0) {
			itoa10(static_cast<int>(readCnt_), r.nameBuf);
			_setBegin(r.name, r.nameBuf);
			name_len = (int)strlen(r.nameBuf);
			_setLength(r.name, name_len);
		}
		if(strlen(rb.readOrigBuf) != 0 && empty(rb.patFw)) {
			return parse(rb, r);
		}
		return true;
	}

	pair<bool, bool> readLight(ReadBuf& r);

	void finalize(ReadBuf &r) const;

	/// Read another pattern from a FASTQ input file
//	virtual void readHeavy(ReadBuf& r, uint32_t& patid) {
//		const int bufSz = ReadBuf::BUF_SIZE;
//		while(true) {
//			int c;
//			int dstLen = 0;
//			int nameLen = 0;
//			r.color = color_;
//			r.primer = -1;
//			// Pick off the first at
//			if(first_) {
//				c = fb_.get();
//				if(c != '@') {
//					c = getOverNewline(fb_);
//					if(c < 0) { bail(r); return; }
//				}
//				if(c != '@') {
//					cerr << "Error: reads file does not look like a FASTQ file" << endl;
//					throw 1;
//				}
//				assert_eq('@', c);
//				first_ = false;
//			}
//
//			// Read to the end of the id line, sticking everything after the '@'
//			// into *name
//			while(true) {
//				c = fb_.get();
//				if(c < 0) { bail(r); return; }
//				if(c == '\n' || c == '\r') {
//					// Break at end of line, after consuming all \r's, \n's
//					while(c == '\n' || c == '\r') {
//						c = fb_.get();
//						if(c < 0) { bail(r); return; }
//					}
//					break;
//				}
//				r.nameBuf[nameLen++] = c;
//				if(nameLen > bufSz-2) {
//					// Too many chars in read name; print friendly error message
//					_setBegin(r.name, r.nameBuf);
//					_setLength(r.name, nameLen);
//					cerr << "FASTQ read name is too long; read names must be " << (bufSz-2) << " characters or fewer." << endl;
//					cerr << "Beginning of bad read name: " << r.name << endl;
//					throw 1;
//				}
//			}
//			_setBegin(r.name, r.nameBuf);
//			assert_leq(nameLen, bufSz-2);
//			_setLength(r.name, nameLen);
//			// c now holds the first character on the line after the
//			// @name line
//
//			// fb_ now points just past the first character of a
//			// sequence line, and c holds the first character
//			int charsRead = 0;
//			uint8_t *sbuf = r.patBufFw;
//			int dstLens[] = {0, 0, 0, 0};
//			int *dstLenCur = &dstLens[0];
//			int mytrim5 = this->trim5_;
//			if(color_) {
//				// This may be a primer character.  If so, keep it in the
//				// 'primer' field of the read buf and parse the rest of the
//				// read without it.
//				c = toupper(c);
//				if(asc2dnacat[c] > 0) {
//					// First char is a DNA char
//					int c2 = toupper(fb_.peek());
//					// Second char is a color char
//					if(asc2colcat[c2] > 0) {
//						r.primer = c;
//						r.trimc = c2;
//						mytrim5 += 2; // trim primer and first color
//					}
//				}
//				if(c < 0) { bail(r); return; }
//			}
//			int trim5 = mytrim5;
//			if(c == '+') {
//				// Read had length 0; print warning (if not quiet) and quit
//				if(!quiet) {
//					cerr << "Warning: Skipping read (" << r.name << ") because it had length 0" << endl;
//				}
//				peekToEndOfLine(fb_);
//				fb_.get();
//				continue;
//			}
//			while(c != '+') {
//				// Convert color numbers to letters if necessary
//				if(color_) {
//					if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
//				}
//				if(c == '.') c = 'N';
//				if(isalpha(c)) {
//					// If it's past the 5'-end trim point
//					assert_in(toupper(c), "ACGTN");
//					if(charsRead >= trim5) {
//						if((*dstLenCur) >= 1024) tooManySeqChars(r.name);
//						sbuf[(*dstLenCur)++] = charToDna5[c];
//					}
//					charsRead++;
//				}
//				c = fb_.get();
//				if(c < 0) { bail(r); return; }
//			}
//			// Trim from 3' end
//			dstLen = dstLens[0];
//			charsRead = dstLen + mytrim5;
//			dstLen -= this->trim3_;
//			// Set trimmed bounds of buffers
//			_setBegin(r.patFw, (Dna5*)r.patBufFw);
//			_setLength(r.patFw, dstLen);
//			assert_eq('+', c);
//
//			// Chew up the optional name on the '+' line
//			peekToEndOfLine(fb_);
//
//			// Now read the qualities
//			if (intQuals_) {
//				int qualsRead = 0;
//				char buf[4096];
//				if(color_ && r.primer != -1) mytrim5--;
//				while (qualsRead < charsRead) {
//					vector<string> s_quals;
//					if(!tokenizeQualLine(fb_, buf, 4096, s_quals)) break;
//					for (unsigned int j = 0; j < s_quals.size(); ++j) {
//						char c = intToPhred33(atoi(s_quals[j].c_str()), solQuals_);
//						assert_geq(c, 33);
//						if (qualsRead >= mytrim5) {
//							size_t off = qualsRead - mytrim5;
//							if(off >= 1024) tooManyQualities(r.name);
//							r.qualBuf[off] = c;
//						}
//						++qualsRead;
//					}
//				} // done reading integer quality lines
//				if(color_ && r.primer != -1) mytrim5++;
//				if(qualsRead < charsRead-mytrim5) {
//					tooFewQualities(r.name);
//				} else if(qualsRead > charsRead-mytrim5+1) {
//					tooManyQualities(r.name);
//				}
//				if(qualsRead == charsRead-mytrim5+1 && color_ && r.primer != -1) {
//					for(int i = 0; i < qualsRead-1; i++) {
//						r.qualBuf[i] = r.qualBuf[i+1];
//					}
//				}
//				_setBegin(r.qual, (char*)r.qualBuf);
//				_setLength(r.qual, dstLen);
//				peekOverNewline(fb_);
//			} else {
//				// Non-integer qualities
//				char *qbuf = r.qualBuf;
//				trim5 = mytrim5;
//				if(color_ && r.primer != -1) trim5--;
//				int itrim5 = trim5;
//				int qualsRead[4] = {0, 0, 0, 0};
//				int *qualsReadCur = &qualsRead[0];
//				while(true) {
//					c = fb_.get();
//					if(c == ' ') {
//						wrongQualityFormat(r.name);
//					}
//					if(c < 0) { bail(r); return; }
//					if (c != '\r' && c != '\n') {
//						if (*qualsReadCur >= trim5) {
//							size_t off = (*qualsReadCur) - trim5;
//							if(off >= 1024) tooManyQualities(r.name);
//							c = charToPhred33(c, solQuals_, phred64Quals_);
//							assert_geq(c, 33);
//							qbuf[off] = c;
//						}
//						(*qualsReadCur)++;
//					} else {
//						break;
//					}
//				}
//				qualsRead[0] -= this->trim3_;
//				int qRead = 0;
//				if (qualsRead[0] > itrim5)
//					qRead = (int)(qualsRead[0] - itrim5);
//				if(qRead < dstLen) {
//					tooFewQualities(r.name);
//				} else if(qRead > dstLen+1) {
//					tooManyQualities(r.name);
//				}
//				if(qRead == dstLen+1 && color_ && r.primer != -1) {
//					for(int i = 0; i < dstLen; i++) {
//						r.qualBuf[i] = r.qualBuf[i+1];
//					}
//				}
//				_setBegin (r.qual, (char*)r.qualBuf);
//				_setLength(r.qual, dstLen);
//
//				if(c == '\r' || c == '\n') {
//					c = peekOverNewline(fb_);
//				} else {
//					c = peekToEndOfLine(fb_);
//				}
//			}
//			r.readOrigBufLen = fb_.copyLastN(r.readOrigBuf);
//			fb_.resetLastN();
//
//			c = fb_.get();
//			assert(c == -1 || c == '@');
//
//			// Set up a default name if one hasn't been set
//			if(nameLen == 0) {
//				itoa10((int)readCnt_, r.nameBuf);
//				_setBegin(r.name, r.nameBuf);
//				nameLen = (int)strlen(r.nameBuf);
//				_setLength(r.name, nameLen);
//			}
//			r.trimmed3 = this->trim3_;
//			r.trimmed5 = mytrim5;
//			assert_gt(nameLen, 0);
//			readCnt_++;
//			patid = (uint32_t)(readCnt_-1);
//			return;
//		}
//	}
	/// Read another read pair from a FASTQ input file
	virtual BoolTriple readPairLight(ReadBuf& ra, ReadBuf& rb) {
		// (For now, we shouldn't ever be here)
		cerr << "In FastqPatternSource.readPair()" << endl;
		throw 1;
		return BoolTriple();
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

	/**
	 * Do things we need to do if we have to bail in the middle of a
	 * read, usually because we reached the end of the input without
	 * finishing.
	 */
	void bail(ReadBuf& r) {
		seqan::clear(r.patFw);
		r.readOrigBufLen = 0;
	}

	bool first_;
	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	bool color_;
};

/**
 * Read a Raw-format file (one sequence per line).  No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class RawPatternSource : public BufferedFilePatternSource {

public:

	RawPatternSource(uint32_t seed,
	                 const vector<string>& infiles,
	                 bool color,
	                 bool randomizeQuals = false,
	                 const char *dumpfile = NULL,
	                 bool verbose = false,
	                 int trim3 = 0,
	                 int trim5 = 0,
	                 uint32_t skip = 0) :
		BufferedFilePatternSource(seed, infiles, NULL, randomizeQuals, 
		                          dumpfile, verbose, trim3, trim5, skip),
		first_(true), color_(color)
	{ }

	virtual void reset() {
		first_ = true;
		BufferedFilePatternSource::reset();
	}

protected:
	
	/**
	 * "Light" parser.  This is inside the critical section, so the key is to do
	 * just enough parsing so that another function downstream (finalize()) can do
	 * the rest of the parsing.  Really this function's only job is to stick every
	 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
	 * then parses the contents of r.readOrigBuf later.
	 */
	pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a)
	{
		throw 1;
		return make_pair(false, 0);
	}

	/// Read another pattern from a FASTA input file
	bool parse(ReadBuf& r, ReadBuf& rb) const {
		cerr << "In RawPatternSource.parse()" << endl;
		throw 1;
		return false;
	}

	/// Read another pattern from a Raw input file
	virtual pair<bool, bool> readLight(ReadBuf& r) {
		int c;
		int dstLen = 0;
		int nameLen = 0;
		c = getOverNewline(this->fb_);
		if(c < 0) { bail(r); return make_pair(false, false); }
		assert(!isspace(c));
		r.color = color_;
		int mytrim5 = this->trim5_;
		if(first_) {
			// Check that the first character is sane for a raw file
			int cc = c;
			if(color_) {
				if(cc >= '0' && cc <= '4') cc = "ACGTN"[(int)cc - '0'];
				if(cc == '.') cc = 'N';
			}
			if(dna4Cat[cc] == 0) {
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
		if(color_) {
			// This may be a primer character.  If so, keep it in the
			// 'primer' field of the read buf and parse the rest of the
			// read without it.
			c = toupper(c);
			if(asc2dnacat[c] > 0) {
				// First char is a DNA char
				int c2 = toupper(fb_.peek());
				// Second char is a color char
				if(asc2colcat[c2] > 0) {
					r.primer = c;
					r.trimc = c2;
					mytrim5 += 2; // trim primer and first color
				}
			}
			if(c < 0) { bail(r); return make_pair(false, false); }
		}
		// _in now points just past the first character of a sequence
		// line, and c holds the first character
		while(!isspace(c) && c >= 0) {
			if(color_) {
				if(c >= '0' && c <= '4') c = "ACGTN"[(int)c - '0'];
			}
			if(c == '.') c = 'N';
			if(isalpha(c) && dstLen >= mytrim5) {
				size_t len = dstLen - mytrim5;
				if(len >= 1024) tooManyQualities(String<char>("(no name)"));
				r.patBufFw [len] = charToDna5[c];
				r.qualBuf[len] = 'I';
				dstLen++;
			} else if(isalpha(c)) dstLen++;
			if(isspace(fb_.peek())) break;
			c = fb_.get();
		}
		if(dstLen >= (this->trim3_ + mytrim5)) {
			dstLen -= (this->trim3_ + mytrim5);
		} else {
			dstLen = 0;
		}
		_setBegin (r.patFw,  (Dna5*)r.patBufFw);
		_setLength(r.patFw,  dstLen);
		_setBegin (r.qual, r.qualBuf);
		_setLength(r.qual, dstLen);

		c = peekToEndOfLine(fb_);
		r.trimmed3 = this->trim3_;
		r.trimmed5 = mytrim5;
		r.readOrigBufLen = fb_.copyLastN(r.readOrigBuf);
		fb_.resetLastN();

		// Set up name
		itoa10((int)readCnt_, r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		nameLen = (int)strlen(r.nameBuf);
		_setLength(r.name, nameLen);
		readCnt_++;

		r.patid = (uint32_t)(readCnt_-1);
		return make_pair(true, false);
	}
	
	/// Read another read pair from a FASTQ input file
	virtual BoolTriple readPairLight(ReadBuf& ra, ReadBuf& rb) {
		// (For now, we shouldn't ever be here)
		cerr << "In RawPatternSource.readPair()" << endl;
		throw 1;
		return BoolTriple();
	}

	/**
	 * Finishes parsing outside the critical section
	 */
	virtual void finalize(ReadBuf& r) const { }
	
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
	/**
	 * Do things we need to do if we have to bail in the middle of a
	 * read, usually because we reached the end of the input without
	 * finishing.
	 */
	void bail(ReadBuf& r) {
		seqan::clear(r.patFw);
		fb_.resetLastN();
	}
	
	bool first_;
	bool color_;
};

#endif /*PAT_H_*/
