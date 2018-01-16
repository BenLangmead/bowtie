#ifndef PAT_H_
#define PAT_H_

#include <cassert>
#include <cmath>
#include <zlib.h>
#include <sys/stat.h>
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
#include "read.h"
#include "search_globals.h"

#ifdef _WIN32
#define getc_unlocked _fgetc_nolock
#endif

/**
 * Classes and routines for reading reads from various input sources.
 */

using namespace std;
using namespace seqan;

typedef uint64_t TReadId;

/**
 * Parameters affecting how reads and read in.
 * Note: Bowtie 2 uses this but Bowtie doesn't yet.
 */
struct PatternParams {
	
	PatternParams() { }

	PatternParams(
		int format_,
		bool color_,
		bool fileParallel_,
		uint32_t seed_,
		size_t max_buf_,
		size_t buffer_sz_,
		bool solexa64_,
		bool phred64_,
		bool intQuals_,
		int trim5_,
		int trim3_,
		int sampleLen_,
		int sampleFreq_,
		size_t skip_,
		int nthreads_,
		int block_bytes_,
		int reads_per_block_,
		bool fixName_) :
		format(format_),
		color(color_),
		fileParallel(fileParallel_),
		seed(seed_),
		max_buf(max_buf_),
		buffer_sz(buffer_sz_),
		solexa64(solexa64_),
		phred64(phred64_),
		intQuals(intQuals_),
		trim5(trim5_),
		trim3(trim3_),
		sampleLen(sampleLen_),
		sampleFreq(sampleFreq_),
		skip(skip_),
		nthreads(nthreads_),
		block_bytes(block_bytes_),
		reads_per_block(reads_per_block_),
		fixName(fixName_) { }

	int format;           // file format
	bool color;           // colorspace?
	bool fileParallel;    // true -> wrap files with separate PatternComposers
	uint32_t seed;        // pseudo-random seed
	size_t max_buf;       // number of reads to buffer in one read
	size_t buffer_sz;     // input buffer size
	bool solexa64;        // true -> qualities are on solexa64 scale
	bool phred64;         // true -> qualities are on phred64 scale
	bool intQuals;        // true -> qualities are space-separated numbers
	int trim5;            // amt to hard clip from 5' end
	int trim3;            // amt to hard clip from 3' end
	int sampleLen;        // length of sampled reads for FastaContinuous...
	int sampleFreq;       // frequency of sampled reads for FastaContinuous...
	size_t skip;          // skip the first 'skip' patterns
	int nthreads;         // number of threads for locking
	int block_bytes;      // # bytes in one input block, 0 if we're not using blocked input
	int reads_per_block;  // # reads per input block, 0 if we're not using blockeds input
	bool fixName;         // whether to fix mate names
};

/**
 * All per-thread storage for input read data.
 */
struct PerThreadReadBuf {
	
	PerThreadReadBuf(size_t max_buf) :
		max_buf_(max_buf),
		bufa_(new Read[max_buf]),
		bufb_(new Read[max_buf]),
		cur_buf_(max_buf),
		rdid_(std::numeric_limits<TReadId>::max())
	{ }
	
	~PerThreadReadBuf() {
		delete[] bufa_;
		delete[] bufb_;
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
		cur_buf_ = max_buf_;
		for(size_t i = 0; i < max_buf_; i++) {
			bufa_[i].reset();
			assert(bufa_[i].empty());
			bufb_[i].reset();
			assert(bufb_[i].empty());
		}
		rdid_ = std::numeric_limits<TReadId>::max();
	}
	
	/**
	 * Advance cursor to next element
	 */
	void next() {
		assert_lt(cur_buf_, max_buf_);
		cur_buf_++;
	}
	
	/**
	 * Return true when there's nothing left to dish out.
	 */
	bool exhausted() {
		assert_leq(cur_buf_, max_buf_);
		return cur_buf_ >= max_buf_ - 1;
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
	Read *bufa_;           // Read buffer for mate as
	Read *bufb_;           // Read buffer for mate bs
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
		const PatternParams& pp,
		const char *dumpfile = NULL) :
		pp_(pp),
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
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const = 0;
	
	/// Reset state to start over again with the first read
	virtual void reset() { readCnt_ = 0; }

	/**
	 * Return the number of reads attempted.
	 */
	TReadId readCount() const { return readCnt_; }

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const = 0;

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
	
	// Parsing parameters
	const PatternParams& pp_;

	/// The number of reads read by this PatternSource
	volatile uint64_t readCnt_;

	const char *dumpfile_; /// dump patterns to this file before returning them
	ofstream out_;         /// output stream for dumpfile

	/// Lock enforcing mutual exclusion for (a) file I/O, (b) writing fields
	/// of this or another other shared object.
	MUTEX_T mutex;
};

/**
 * Encapsulates a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public PatternSource {

public:

	VectorPatternSource(
		const vector<string>& v,
		const PatternParams& pp,
		const char *dumpfile = NULL);
	
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
		PatternSource::reset();
		cur_ = 0;
		paired_ = false;
	}
	
	/**
	 * Finishes parsing outside the critical section
	 */
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const;
	
	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const { return false; }

private:

	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		bool batch_a);

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
class CFilePatternSource : public PatternSource {
public:
	CFilePatternSource(
		const vector<string>& infiles,
		const vector<string>* qinfiles,
		const PatternParams& pp,
		const char *dumpfile = NULL) :
		PatternSource(pp, dumpfile),
		infiles_(infiles),
		filecur_(0),
		fp_(NULL),
		qfp_(NULL),
		zfp_(NULL),
		is_open_(false),
		first_(true),
		buf_(NULL),
		qbuf_(NULL),
		compressed_(false),
		buffer_sz_(pp.buffer_sz)
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
		buf_ = new char[buffer_sz_];
		qbuf_ = new char[buffer_sz_];
		open(); // open first file in the list
		filecur_++;
	}

	/**
	 * Close open file.
	 */
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
		if(buf_ != NULL) {
			delete[] buf_;
			buf_ = NULL;
		}
		if(qbuf_ != NULL) {
			delete[] qbuf_;
			qbuf_ = NULL;
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
		PatternSource::reset();
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
	
	vector<string> infiles_; /// filenames for read files
	vector<string> qinfiles_; /// filenames for quality files
	vector<bool> errs_; /// whether we've already printed an error for each file
	size_t filecur_;   /// index into infiles_ of next file to read
	FILE *fp_; /// read file currently being read from
	FILE *qfp_; /// quality file currently being read from
    gzFile zfp_;
	bool is_open_; /// whether fp_ is currently open
	bool first_;
	char *buf_; /// file buffer for sequences
	char *qbuf_; /// file buffer for qualities
	size_t buffer_sz_; // buffer size for use w/ setvbuf/gzbuffer
    bool compressed_;

private:

	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		bool batch_a);

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
		const PatternParams& pp,
		const char *dumpfile = NULL) :
		CFilePatternSource(
			infiles,
			qinfiles,
			pp,
			dumpfile),
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
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const;

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const {
		return false;  // but wouldn't be too much work to support it
	}

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

	bool first_;
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
		const PatternParams& pp,
		bool secondName, // whether it's --12/--tab5 or --tab6
		const char *dumpfile = NULL) :
		CFilePatternSource(
			infiles,
			NULL,
			pp,
			dumpfile) { }

	/**
	 * Finalize tabbed parsing outside critical section.
	 */
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const;

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const {
		return false;  // but wouldn't be too much work to support it
	}

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
		const PatternParams& pp,
		const char *dumpfile = NULL) :
		CFilePatternSource(
			infiles,
			NULL,
			pp,
			dumpfile),
		length_(pp.sampleLen),
		freq_(pp.sampleFreq),
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
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const;

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const { return false; }

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
		const PatternParams& pp,
		const char *dumpfile = NULL) :
		CFilePatternSource(
			infiles,
			NULL,
			pp,
			dumpfile),
		first_(true) { }
	
	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize FASTQ parsing outside critical section.
	 */
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const;

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const { return true; }

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
		const PatternParams& pp,
		const char *dumpfile = NULL) :
		CFilePatternSource(
			infiles,
			NULL,
			pp,
			dumpfile),
		first_(true) { }

	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize raw parsing outside critical section.
	 */
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const;

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const {
		return false;  // but wouldn't be too much work to support it
	}

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
	
	bool first_;
};

/**
 * Abstract parent class for synhconized sources of paired-end reads
 * (and possibly also single-end reads).
 */
class PatternComposer {

public:

	PatternComposer(const PatternParams& p) : pp_(p), mutex_() { }
	
	virtual ~PatternComposer() { }

	virtual void reset() = 0;
	
	/**
	 * Member function override by concrete, format-specific classes.
	 */
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf& pt) = 0;
	
	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const = 0;
	
	virtual void free_pmembers(const vector<PatternSource*> &elist) {
		for (size_t i = 0; i < elist.size(); i++) {
			if (elist[i] != NULL) {
				delete elist[i];
			}
		}
	}

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const = 0;

protected:

	const PatternParams& pp_;
	MUTEX_T mutex_;
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class SoloPatternComposer : public PatternComposer {

public:

	SoloPatternComposer(
		const vector<PatternSource*>& src,
		const PatternParams& p) :
		PatternComposer(p),
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
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const
	{
		return src_[0]->parse(ra, rb, cura, curb, rdid);
	}

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const {
		return src_[0]->supportsBlocks();
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

	DualPatternComposer(
		const vector<PatternSource*>& srca,
		const vector<PatternSource*>& srcb,
		const PatternParams& p) :
		PatternComposer(p), cur_(0), srca_(srca), srcb_(srcb)
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
	virtual bool parse(
		Read& ra, Read& rb,
		ParsingCursor& cura, ParsingCursor& curb,
		TReadId rdid) const
	{
		return srca_[0]->parse(ra, rb, cura, curb, rdid);
	}

	/**
	 * Returns true iff file format can be parsed in fixed-size blocks.
	 */
	virtual bool supportsBlocks() const {
		return srca_[0]->supportsBlocks();
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
		const PatternParams& pp) :
		composer_(composer),
		pp_(pp),
		blockReads_(pp.block_bytes > 0 && composer.supportsBlocks()),
		buf_(blockReads_ ? pp.reads_per_block : pp.max_buf),
		last_batch_(false),
		last_batch_size_(0) { }


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
		cura_.buf = &(buf_.bufa_[0].readOrigBuf);
		curb_.buf = &(buf_.bufb_[0].readOrigBuf);
		cura_.off = curb_.off = 0;
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
		// advance cursors to next read in case of non-blocked input
		if(!blockReads_) {
			cura_.buf = &ra.readOrigBuf;
			curb_.buf = &rb.readOrigBuf;
			cura_.off = curb_.off = 0;
		}
		assert(!cura_.buf->empty());
		return composer_.parse(ra, rb, cura_, curb_, buf_.rdid());
	}

	PatternComposer& composer_; // pattern composer
	const PatternParams& pp_;   // parsing parameters
	bool blockReads_;           // read input in blocks?
	PerThreadReadBuf buf_;      // read data buffer
	ParsingCursor cura_, curb_; // parsing cursors
	bool last_batch_;           // true if this is final batch
	int last_batch_size_;       // # reads read in previous batch
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class PatternSourcePerThreadFactory {
public:
	PatternSourcePerThreadFactory(
		PatternComposer& composer,
		const PatternParams& pp):
		composer_(composer),
		pp_(pp) { }

	/**
	 * Create a new heap-allocated PatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new PatternSourcePerThread(composer_, pp_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * PatternSourcePerThreads.
	 */
	virtual std::vector<PatternSourcePerThread*>* create(uint32_t n) const {
		std::vector<PatternSourcePerThread*>* v = new std::vector<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new PatternSourcePerThread(composer_, pp_));
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
	// Parsing parameters
	const PatternParams& pp_;
};

#endif /*PAT_H_*/
