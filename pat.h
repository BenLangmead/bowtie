#ifndef PAT_H_
#define PAT_H_

#include <iostream>
#include <cassert>
#include <cmath>
#include <zlib.h>
#include <sys/stat.h>
#include <stdexcept>
#include <string>
#include <cstring>
#include <ctype.h>
#include <fstream>

#include "alphabet.h"
#include "assert_helpers.h"
#include "ds.h"
#include "ds.h"
#include "filebuf.h"
#include "qual.h"
#include "random_source.h"
#include "read.h"
#include "search_globals.h"
#include "sstring.h"
#include "threading.h"
#include "tokenize.h"
#include "util.h"

#ifdef _WIN32
#define getc_unlocked _fgetc_nolock
#endif

/**
 * Classes and routines for reading reads from various input sources.
 */

using namespace std;

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
		assert_neq(rdid_, std::numeric_limits<TReadId>::max());
	}

	const size_t max_buf_; // max # reads to read into buffer at once
	EList<Read> bufa_; // Read buffer for mate as
	EList<Read> bufb_; // Read buffer for mate bs
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
	PatternSource() :
		readCnt_(0),
		mutex()
	{ }

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
	 * Default format for dumping a read to an output stream.  Concrete
	 * subclasses might want to do something fancier.
	 */
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << name << ": " << seq << " " << qual << endl;
	}

	/// The number of reads read by this PatternSource
	volatile uint64_t readCnt_;

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
	TrimmingPatternSource(int trim3 = 0,
	                      int trim5 = 0) :
		PatternSource(),
		trim3_(trim3), trim5_(trim5) { }
protected:
	int trim3_;
	int trim5_;
};

extern void wrongQualityFormat(const BTString& read_name);
extern void tooFewQualities(const BTString& read_name);
extern void tooManyQualities(const BTString& read_name);
extern void tooManySeqChars(const BTString& read_name);

/**
 * Encapsulates a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public TrimmingPatternSource {
public:
	VectorPatternSource(
		const EList<string>& v,
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

	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		bool batch_a);

	size_t cur_;                      // index for first read of next batch
	bool paired_;                     // whether reads are paired
	EList<std::string> tokbuf_; // buffer for storing parsed tokens
	EList<std::string> bufs_;   // per-read buffers
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
		const EList<string>& infiles,
		const EList<string>* qinfiles,
		int trim3 = 0,
	    int trim5 = 0) :
		TrimmingPatternSource( trim3, trim5),
		infiles_(infiles),
		filecur_(0),
		fp_(NULL),
		qfp_(NULL),
		zfp_(NULL),
		is_open_(false),
		first_(true)
	{
		qinfiles_.clear();
		if(qinfiles != NULL) qinfiles_ = *qinfiles;
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size());
		errs_.fill(0, infiles_.size(), false);
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
			if (compressed_) {
				gzclose(zfp_);
				zfp_ = NULL;
			}
			else if (fp_ != stdin) {
				fclose(fp_);
				fp_ = NULL;
			}
			if(qfp_ != NULL && qfp_ != stdin) {
				fclose(qfp_);
				qfp_ = NULL;
			}
			assert(zfp_ == NULL);
			assert(fp_ == NULL || fp_ == stdin);
			assert(qfp_ == NULL || qfp_ == stdin);
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
		bool batch_a,
		size_t read_idx) = 0;

	/**
	 * Reset state to handle a fresh file
	 */
	virtual void resetForNextFile() { }

	/**
	 * Open the next file in the list of input files.
	 */
	void open();

	int getc_wrapper() {
		return compressed_ ? gzgetc(zfp_) : getc_unlocked(fp_);
	}

	int ungetc_wrapper(int c) {
		return compressed_ ? gzungetc(c, zfp_) : ungetc(c, fp_);
	}

	bool is_gzipped_file(const std::string& filename) {
		struct stat s;
		if (stat(filename.c_str(), &s) != 0) {
			perror("stat");
		}
		else {
			if (S_ISFIFO(s.st_mode))
				return true;
		}
		size_t pos = filename.find_last_of(".");
		std::string ext = (pos == std::string::npos) ? "" : filename.substr(pos + 1);
		if (ext == "" || ext == "gz" || ext == "Z") {
			return true;
		}
		return false;
	}

	EList<string> infiles_; /// filenames for read files
	EList<string> qinfiles_; /// filenames for quality files
	EList<bool> errs_; /// whether we've already printed an error for each file
	size_t filecur_;   /// index into infiles_ of next file to read
	FILE *fp_; /// read file currently being read from
	FILE *qfp_; /// quality file currently being read from
    gzFile zfp_;
	bool is_open_; /// whether fp_ is currently open
	bool first_;
	char buf_[64*1024]; /// file buffer for sequences
	char qbuf_[64*1024]; /// file buffer for qualities
    bool compressed_;

private:

	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		bool batch_a);

};

/**
 * Synchronized concrete pattern source for a list of FASTA
 */
class FastaPatternSource : public CFilePatternSource {

public:

	FastaPatternSource(
		const EList<string>& infiles,
		const EList<string>* qinfiles,
		int trim3 = 0,
		int trim5 = 0) :
		CFilePatternSource(
			infiles,
			qinfiles,
			trim3,
			trim5),
		first_(true) { }

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
		bool batch_a,
		size_t read_idx);

	/**
	 * Reset state to handle a fresh file
	 */
	virtual void resetForNextFile() {
		first_ = true;
	}

	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << ">" << name << endl << seq << endl;
	}

private:

	bool first_;
};


/**
 * Tokenize a line of space-separated integer quality values.
 */
static inline bool tokenizeQualLine(FileBuf& filebuf, char *buf, size_t buflen, EList<string>& toks) {
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
		const EList<string>& infiles,
		bool secondName,  // whether it's --12/--tab5 or --tab6
		int trim3 = 0,
		int trim5 = 0,
		bool solQuals = false,
		bool phred64Quals = false,
		bool intQuals = false) :
		CFilePatternSource(
			infiles,
			NULL,
			trim3,
			trim5),
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
		bool batch_a,
		size_t read_idx);

	/**
	 * Dump a FASTQ-style record for the read.
	 */
	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << "@" << name << endl << seq << endl
		    << "+" << endl << qual << endl;
	}

protected:

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
			const EList<string>& infiles,
			size_t length,
			size_t freq) :
		CFilePatternSource(
			infiles,
			NULL,
			0,
			0),
		length_(length),
		freq_(freq),
		eat_(length_-1),
		beginning_(true),
		bufCur_(0),
		cur_(0llu),
		last_(0llu)
		{
			assert_gt(freq_, 0);
			resetForNextFile();
			assert_lt(length_, 1024);
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
		bool batch_a,
		size_t read_idx);

	/**
	 * Reset state to be read for the next file.
	 */
	virtual void resetForNextFile() {
		eat_ = length_-1;
		name_prefix_buf_.clear();
		beginning_ = true;
		bufCur_ = 0;
		last_ = cur_;
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
	std::string name_prefix_buf_; /// FASTA sequence name buffer
	char name_int_buf_[20]; /// for composing offsets for names
	size_t bufCur_;     /// buffer cursor; points to where we should
	                    /// insert the next character
	uint64_t cur_;
	uint64_t last_;/// number to subtract from readCnt_ to get
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
		const EList<string>& infiles,
		int trim3 = 0,
		int trim5 = 0,
		bool solexa_quals = false,
		bool phred64Quals = false,
		bool integer_quals = false,
		bool interleaved = false,
		uint32_t skip = 0) :
		CFilePatternSource(
			infiles,
			NULL,
			trim3,
			trim5),
		first_(true),
		solQuals_(solexa_quals),
		phred64Quals_(phred64Quals),
		intQuals_(integer_quals),
		interleaved_(interleaved) { }

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
		bool batch_a,
		size_t read_idx);

	virtual void resetForNextFile() {
		first_ = true;
	}

	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
	}

private:

	bool first_;
	bool solQuals_;
	bool phred64Quals_;
	bool intQuals_;
	bool interleaved_;
};

/**
 * Read a Raw-format file (one sequence per line).  No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class RawPatternSource : public CFilePatternSource {

public:

	RawPatternSource(
		const EList<string>& infiles,
		int trim3 = 0,
		int trim5 = 0) :
		CFilePatternSource(
			infiles,
			NULL,
			trim3,
			trim5),
		first_(true) { }

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
		bool batch_a,
		size_t read_idx);

	virtual void resetForNextFile() {
		first_ = true;
	}

	virtual void dump(ostream& out,
	                  const BTDnaString& seq,
	                  const BTString& qual,
	                  const BTString& name)
	{
		out << seq << endl;
	}


private:

	bool first_;
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

	virtual void free_pmembers( const EList<PatternSource*> &elist) {
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

	SoloPatternComposer(const EList<PatternSource*>& src) :
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
	EList<PatternSource*> src_; /// PatternSources for paired-end reads
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class DualPatternComposer : public PatternComposer {

public:

	DualPatternComposer(const EList<PatternSource*>& srca,
	                        const EList<PatternSource*>& srcb) :
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
	EList<PatternSource*> srca_; /// PatternSources for 1st mates and/or unpaired reads
	EList<PatternSource*> srcb_; /// PatternSources for 2nd mates
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
		seed_(seed),
		batch_id_(0) { }

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

	size_t batch_id() const { return batch_id_; }

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
		batch_id_ = (size_t)(buf_.rdid()/buf_.max_buf_);
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
	size_t last_batch_size_;  // # reads read in previous batch
	uint32_t skip_;           // skip reads with rdids less than this
	uint32_t seed_;           // pseudo-random seed based on read content
	size_t batch_id_;	  // identify batches of reads for reordering
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
	virtual EList<PatternSourcePerThread*>* create(uint32_t n) const {
		EList<PatternSourcePerThread*>* v = new EList<PatternSourcePerThread*>;
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
	virtual void destroy(EList<PatternSourcePerThread*>* composers) const {
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

	virtual ~PatternSourcePerThreadFactory() {}

private:
	/// Container for obtaining paired reads from PatternSources
	PatternComposer& composer_;
	/// Maximum size of batch to read in
	uint32_t max_buf_;
	uint32_t skip_;
	uint32_t seed_;
};

#endif /*PAT_H_*/
