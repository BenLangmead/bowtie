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

#ifdef _WIN32
#define getc_unlocked _fgetc_nolock
#endif

/**
 * Classes and routines for reading reads from various input sources.
 */

using namespace std;
using namespace seqan;

/**
 * C++ version char* style "itoa":
 */
template<typename T>
char* itoa10(const T& value, char* result) {
	// Check that base is valid
	char* out = result;
	T quotient = value;
	if(std::numeric_limits<T>::is_signed) {
		if(quotient <= 0) quotient = -quotient;
	}
	// Now write each digit from most to least significant
	do {
		*out = "0123456789"[quotient % 10];
		++out;
		quotient /= 10;
	} while (quotient > 0);
	// Only apply negative sign for base 10
	if(std::numeric_limits<T>::is_signed) {
		// Avoid compiler warning in cases where T is unsigned
		if (value <= 0 && value != 0) *out++ = '-';
	}
	reverse( result, out );
	*out = 0; // terminator
	return out;
}

typedef uint64_t TReadId;

/**
 * Parameters affecting how reads and read in.
 * Note: Bowtie 2 uses this but Bowtie doesn't yet.
 */
struct PatternParams {
	
	PatternParams() { }

	PatternParams(
		int format_,
		bool fileParallel_,
		uint32_t seed_,
		size_t max_buf_,
		bool solexa64_,
		bool phred64_,
		bool intQuals_,
		int sampleLen_,
		int sampleFreq_,
		size_t skip_,
		int nthreads_,
		bool fixName_) :
		format(format_),
		fileParallel(fileParallel_),
		seed(seed_),
		max_buf(max_buf_),
		solexa64(solexa64_),
		phred64(phred64_),
		intQuals(intQuals_),
		sampleLen(sampleLen_),
		sampleFreq(sampleFreq_),
		skip(skip_),
		nthreads(nthreads_),
		fixName(fixName_) { }

	int format;           // file format
	bool fileParallel;    // true -> wrap files with separate PatternComposers
	uint32_t seed;        // pseudo-random seed
	size_t max_buf;       // number of reads to buffer in one read
	bool solexa64;        // true -> qualities are on solexa64 scale
	bool phred64;         // true -> qualities are on phred64 scale
	bool intQuals;        // true -> qualities are space-separated numbers
	int sampleLen;        // length of sampled reads for FastaContinuous...
	int sampleFreq;       // frequency of sampled reads for FastaContinuous...
	size_t skip;          // skip the first 'skip' patterns
	int nthreads;         // number of threads for locking
	bool fixName;         //
};

/**
 * A buffer for keeping all relevant information about a single read.
 * Each search thread has one.
 */
struct Read {
	Read() { reset(); }

	~Read() {
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
		parsed = false;
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
		parsed = false;
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
	bool          parsed;              // whether read has been fully parsed
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
	
	Read& read_a() { return bufa_[cur_buf_]; }
	Read& read_b() { return bufb_[cur_buf_]; }
	
	const Read& read_a() const { return bufa_[cur_buf_]; }
	const Read& read_b() const { return bufb_[cur_buf_]; }
	
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
	vector<Read> bufa_; // Read buffer for mate as
	vector<Read> bufb_; // Read buffer for mate bs
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
	PatternSource(
		const char *dumpfile = NULL) :
		readCnt_(0),
		dumpfile_(dumpfile),
		mutex()
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
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const = 0;
	
	/// Reset state to start over again with the first read
	virtual void reset() { readCnt_ = 0; }

	/**
	 * Return the number of reads attempted.
	 */
	TReadId readCount() const { return readCnt_; }

protected:

	/**
	 * Dump the contents of the ReadBuf to the dump file.
	 */
	void dumpBuf(const Read& r) {
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

	/// The number of reads read by this PatternSource
	volatile uint64_t readCnt_;

	const char *dumpfile_; /// dump patterns to this file before returning them
	ofstream out_;         /// output stream for dumpfile

	/// Lock enforcing mutual exclusion for (a) file I/O, (b) writing fields
	/// of this or another other shared object.
	MUTEX_T mutex;
};

/**
 * Encapsualtes a source of patterns where each raw pattern is trimmed
 * by some user-defined amount on the 3' and 5' ends.  Doesn't
 * implement the actual trimming - that's up to the concrete
 * descendants.
 */
class TrimmingPatternSource : public PatternSource {
public:
	TrimmingPatternSource(const char *dumpfile = NULL,
	                      int trim3 = 0,
	                      int trim5 = 0) :
		PatternSource(dumpfile),
		trim3_(trim3), trim5_(trim5) { }
protected:
	int trim3_;
	int trim5_;
};

extern void wrongQualityFormat(const String<char>& read_name);
extern void tooFewQualities(const String<char>& read_name);
extern void tooManyQualities(const String<char>& read_name);
extern void tooManySeqChars(const String<char>& read_name);

/**
 * Encapsulates a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public TrimmingPatternSource {
public:
	VectorPatternSource(
		const vector<string>& v,
	    bool color,
	    const char *dumpfile = NULL,
	    int trim3 = 0,
	    int trim5 = 0);
	
	virtual ~VectorPatternSource() { }
	
	/**
	 * Read next batch.  However, batch concept is not very applicable for this
	 * PatternSource where all the info has already been parsed into the fields
	 * in the contsructor.  This essentially modifies the pt as though we read
	 * in some number of patterns.
	 */
	virtual pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		bool batch_a,
		bool lock = true);

	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 */
	virtual void reset() {
		TrimmingPatternSource::reset();
		cur_ = 0;
		paired_ = false;
	}
	
	/**
	 * Finishes parsing outside the critical section
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

private:

	bool color_;                      // colorspace?
	size_t cur_;                      // index for first read of next batch
	bool paired_;                     // whether reads are paired
	std::vector<std::string> tokbuf_; // buffer for storing parsed tokens
	std::vector<std::string> bufs_;   // per-read buffers
	char nametmp_[20];                // temp buffer for constructing name
};

/**
 * Parent class for PatternSources that read from a file.
 * Uses unlocked C I/O, on the assumption that all reading
 * from the file will take place in an otherwise-protected
 * critical section.
 */
class CFilePatternSource : public TrimmingPatternSource {
public:
	CFilePatternSource(
	    const vector<string>& infiles,
	    const vector<string>* qinfiles,
	    const char *dumpfile = NULL,
	    int trim3 = 0,
	    int trim5 = 0) :
		TrimmingPatternSource(dumpfile, trim3, trim5),
		infiles_(infiles),
		filecur_(0),
		fp_(NULL),
		qfp_(NULL),
		is_open_(false),
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
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.  This function gets called by
	 * all the search threads, so we must handle synchronization.
	 *
	 * Returns pair<bool, int> where bool indicates whether we're
	 * completely done, and int indicates how many reads were read.
	 */
	virtual pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		bool batch_a,
		bool lock);
	
	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 * Should only be called by the master thread.
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

	/**
	 * Reset state to handle a fresh file
	 */
	virtual void resetForNextFile() { }
	
	/**
	 * Open the next file in the list of input files.
	 */
	void open();
	
	vector<string> infiles_; /// filenames for read files
	vector<string> qinfiles_; /// filenames for quality files
	vector<bool> errs_; /// whether we've already printed an error for each file
	size_t filecur_;   /// index into infiles_ of next file to read
	FILE *fp_; /// read file currently being read from
	FILE *qfp_; /// quality file currently being read from
	bool is_open_; /// whether fp_ is currently open
	bool first_;
	char buf_[64*1024]; /// file buffer for sequences
	char qbuf_[64*1024]; /// file buffer for qualities
};

/**
 * Synchronized concrete pattern source for a list of FASTA or CSFASTA
 * (if color = true) files.
 */
class FastaPatternSource : public CFilePatternSource {

public:

	FastaPatternSource(
		const vector<string>& infiles,
	    const vector<string>* qinfiles,
	    bool color,
	    const char *dumpfile = NULL,
	    int trim3 = 0,
	    int trim5 = 0,
	    bool solexa64 = false,
	    bool phred64 = false,
	    bool intQuals = false) :
		CFilePatternSource(
			infiles,
			qinfiles,
		    dumpfile,
			trim3,
		    trim5),
		first_(true),
		color_(color),
		solexa64_(solexa64),
		phred64_(phred64),
		intQuals_(intQuals) { }
	
	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 * Should only be called by the master thread.
	 */
	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize FASTA parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a FASTA batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a);

	/**
	 * Reset state to handle a fresh file
	 */
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
class TabbedPatternSource : public CFilePatternSource {
public:
	TabbedPatternSource(
		const vector<string>& infiles,
		bool secondName,  // whether it's --12/--tab5 or --tab6
	    bool color,
	    const char *dumpfile = NULL,
	    int trim3 = 0,
	    int trim5 = 0,
	    bool solQuals = false,
	    bool phred64Quals = false,
	    bool intQuals = false) :
		CFilePatternSource(
			infiles,
			NULL,
		    dumpfile,
		    trim3,
			trim5),
		color_(color),
		solQuals_(solQuals),
		phred64Quals_(phred64Quals),
		intQuals_(intQuals) { }

	/**
	 * Finalize tabbed parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:
	
	/**
	 * Light-parse a batch of tabbed-format reads into given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a);
	
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
	
protected:

	bool color_;        // colorspace reads?
	bool solQuals_;     // base qualities are log odds
	bool phred64Quals_; // base qualities are on -64 scale
	bool intQuals_;     // base qualities are space-separated strings
	bool secondName_;   // true if --tab6, false if --tab5
};

/**
 * Synchronized concrete pattern source for a list of FASTA files where
 * reads need to be extracted from long continuous sequences.
 */
class FastaContinuousPatternSource : public CFilePatternSource {
public:
	FastaContinuousPatternSource(
			const vector<string>& infiles,
			size_t length,
			size_t freq,
			const char *dumpfile = NULL) :
		CFilePatternSource(
			infiles,
			NULL,
		    dumpfile,
			0,
			0),
		length_(length),
		freq_(freq),
		eat_(length_-1),
		beginning_(true),
		bufCur_(0),
		subReadCnt_(0llu)
	{
		assert_gt(freq_, 0);
		resetForNextFile();
		assert_lt(length_, (size_t)Read::BUF_SIZE);
	}

	virtual void reset() {
		CFilePatternSource::reset();
		resetForNextFile();
	}

	/**
	 * Finalize FASTA parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:
	
	/**
	 * Light-parse a batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a);

	/**
	 * Reset state to be read for the next file.
	 */
	virtual void resetForNextFile() {
		eat_ = length_-1;
		//name_prefix_buf_.clear();
		beginning_ = true;
		bufCur_ = 0;
		subReadCnt_ = readCnt_;
	}

private:
	const size_t length_; /// length of reads to generate
	const size_t freq_;   /// frequency to sample reads
	size_t eat_;        /// number of characters we need to skip before
	                    /// we have flushed all of the ambiguous or
	                    /// non-existent characters out of our read
	                    /// window
	bool beginning_;    /// skipping over the first read length?
	char buf_[1024];    /// read buffer
	char name_prefix_buf_[1024]; /// FASTA sequence name buffer
	char name_int_buf_[20]; /// for composing offsets for names
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
	FastqPatternSource(
		const vector<string>& infiles,
	    bool color,
	    const char *dumpfile = NULL,
	    int trim3 = 0,
	    int trim5 = 0,
	    bool solexa_quals = false,
	    bool phred64Quals = false,
	    bool integer_quals = false,
	    uint32_t skip = 0) :
		CFilePatternSource(
			infiles,
			NULL,
		    dumpfile,
		    trim3,
			trim5),
		first_(true),
		solQuals_(solexa_quals),
		phred64Quals_(phred64Quals),
		intQuals_(integer_quals),
		color_(color) { }
	
	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize FASTQ parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:
	
	/**
	 * "Light" parser.  This is inside the critical section, so the key is to do
	 * just enough parsing so that another function downstream (finalize()) can do
	 * the rest of the parsing.  Really this function's only job is to stick every
	 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
	 * then parses the contents of r.readOrigBuf later.
	 */
	virtual pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a);

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
	bool color_;
};

/**
 * Read a Raw-format file (one sequence per line).  No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class RawPatternSource : public CFilePatternSource {

public:

	RawPatternSource(
		const vector<string>& infiles,
	    bool color,
	    const char *dumpfile = NULL,
		int trim3 = 0,
	    int trim5 = 0) :
		CFilePatternSource(
			infiles,
			NULL,
		    dumpfile,
			trim3,
			trim5),
		first_(true),
		color_(color) { }

	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize raw parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:
	
	/**
	 * Light-parse a batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a);

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
	bool color_;
};

/**
 * Abstract parent class for synhconized sources of paired-end reads
 * (and possibly also single-end reads).
 */
class PatternComposer {
public:
	PatternComposer() { }
	
	virtual ~PatternComposer() { }

	virtual void reset() = 0;
	
	/**
	 * Member function override by concrete, format-specific classes.
	 */
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf& pt) = 0;
	
	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) = 0;

	virtual void free_pmembers( const vector<PatternSource*> &elist) {
    		for (size_t i = 0; i < elist.size(); i++) {
        		if (elist[i] != NULL)
            			delete elist[i];
		}
	}

protected:

	MUTEX_T mutex_m; /// mutex for locking critical regions
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class SoloPatternComposer : public PatternComposer {

public:

	SoloPatternComposer(const vector<PatternSource*>& src) :
		PatternComposer(),
		cur_(0),
		src_(src)
	{
	    for(size_t i = 0; i < src_.size(); i++) {
	    	assert(src_[i] != NULL);
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
	pair<bool, int> nextBatch(PerThreadReadBuf& pt);
	
	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) {
		return src_[0]->parse(ra, rb, rdid);
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
	                        const vector<PatternSource*>& srcb) :
		PatternComposer(),
		cur_(0),
		srca_(srca),
		srcb_(srcb)
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
	pair<bool, int> nextBatch(PerThreadReadBuf& pt);
	
	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) {
		return srca_[0]->parse(ra, rb, rdid);
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

	PatternSourcePerThread(
		PatternComposer& composer,
		uint32_t max_buf,
		uint32_t skip,
		uint32_t seed) :
		composer_(composer),
		buf_(max_buf),
      	last_batch_(false),
		last_batch_size_(0),
		skip_(skip),
		seed_(seed) { }

	/**
	 * Get the next paired or unpaired read from the wrapped
	 * PatternComposer.  Returns a pair of bools; first indicates
	 * whether we were successful, second indicates whether we're
	 * done.
	 */
	std::pair<bool, bool> nextReadPair();

	Read& bufa() { return buf_.read_a(); }
	Read& bufb() { return buf_.read_b(); }
	
	const Read& bufa() const { return buf_.read_a(); }
	const Read& bufb() const { return buf_.read_b(); }
	
	TReadId rdid() const { return buf_.rdid(); }
	
	/**
	 * Return true iff the read currently in the buffer is a
	 * paired-end read.
	 */
	bool paired() const {
		// can also do buf_.read_b().mate > 0, but the mate
		// field isn't set until finalize is called, whereas
		// parsed is set by the time parse() is finished.
		return buf_.read_b().parsed;
	}

private:

	/**
	 * When we've finished fully parsing and dishing out reads in
	 * the current batch, we go get the next one by calling into
	 * the composition layer.
	 */
	std::pair<bool, int> nextBatch() {
		buf_.reset();
		std::pair<bool, int> res = composer_.nextBatch(buf_);
		buf_.init();
		return res;
	}

	/**
	 * Once name/sequence/qualities have been parsed for an
	 * unpaired read, set all the other key fields of the Read
	 * struct.
	 */
	void finalize(Read& ra);

	/**
	 * Once name/sequence/qualities have been parsed for a
	 * paired-end read, set all the other key fields of the Read
	 * structs.
	 */
	void finalizePair(Read& ra, Read& rb);
	
	/**
	 * Call into composition layer (which in turn calls into
	 * format layer) to parse the read.
	 */
	bool parse(Read& ra, Read& rb) {
		return composer_.parse(ra, rb, buf_.rdid());
	}

	PatternComposer& composer_; // pattern composer
	PerThreadReadBuf buf_;    // read data buffer
	bool last_batch_;         // true if this is final batch
	int last_batch_size_;     // # reads read in previous batch
	uint32_t skip_;           // skip reads with rdids less than this
	uint32_t seed_;           // pseudo-random seed based on read content
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class PatternSourcePerThreadFactory {
public:
	PatternSourcePerThreadFactory(
		PatternComposer& composer,
		uint32_t max_buf,
		uint32_t skip,
		uint32_t seed):
		composer_(composer),
		max_buf_(max_buf),
		skip_(skip),
		seed_(seed) {}

	/**
	 * Create a new heap-allocated PatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new PatternSourcePerThread(composer_, max_buf_, skip_, seed_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * PatternSourcePerThreads.
	 */
	virtual std::vector<PatternSourcePerThread*>* create(uint32_t n) const {
		std::vector<PatternSourcePerThread*>* v = new std::vector<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new PatternSourcePerThread(composer_, max_buf_, skip_, seed_));
			assert(v->back() != NULL);
		}
		return v;
	}

	/// Free memory associated with a pattern source
	virtual void destroy(PatternSourcePerThread* composer) const {
		assert(composer != NULL);
		// Free the PatternSourcePerThread
		delete composer;
	}

	/// Free memory associated with a pattern source list
	virtual void destroy(std::vector<PatternSourcePerThread*>* composers) const {
		assert(composers != NULL);
		// Free all of the PatternSourcePerThreads
		for(size_t i = 0; i < composers->size(); i++) {
			if((*composers)[i] != NULL) {
				delete (*composers)[i];
				(*composers)[i] = NULL;
			}
		}
		// Free the vector
		delete composers;
	}
	
private:
	/// Container for obtaining paired reads from PatternSources
	PatternComposer& composer_;
	/// Maximum size of batch to read in
	uint32_t max_buf_;
	uint32_t skip_;
	uint32_t seed_;
};

#endif /*PAT_H_*/
