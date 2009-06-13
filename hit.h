#ifndef HIT_H_
#define HIT_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "spinlock.h"
#include "threading.h"
#include "bitset.h"
#include "tokenize.h"
#include "pat.h"
#include "formats.h"
#include "filebuf.h"

/**
 * Classes for dealing with reporting alignments.
 */

using namespace std;
using namespace seqan;

/// Constants for the various output modes
enum output_types {
	FULL = 1,
	CONCISE,
	BINARY,
	NONE
};

/// Names of the various output modes
static const std::string output_type_names[] = {
	"Invalid!",
	"Full",
	"Concise",
	"Binary",
	"None"
};

typedef pair<uint32_t,uint32_t> U32Pair;

// Support reads of up to 1024 characters for now
static const int max_read_bp = 1024;

/**
 * Encapsulates a hit, including a text-id/text-offset pair, a pattern
 * id, and a boolean indicating whether it matched as its forward or
 * reverse-complement version.
 */
class Hit {
public:
	Hit() :
		h(make_pair(0, 0)),
		patId(0),
		patName(String<char>("")),
		patSeq(String<char>("")),
		quals(String<char>("")),
		mms(),
		refcs(),
		oms(0),
		fw(true),
		mate(0) { }

	Hit(const Hit& other) {
		this->operator=(other);
		assert(seqan::begin(patName) != seqan::begin(other.patName));
		assert(seqan::begin(patSeq)  != seqan::begin(other.patSeq));
		assert(seqan::begin(quals)   != seqan::begin(other.quals));
	}

	Hit(U32Pair _h,
		uint32_t _patId,
		const String<char>& _patName,
		const String<Dna5>& _patSeq,
		const String<char>& _quals,
		bool _fw,
		const FixedBitset<max_read_bp>& _mms,
		const vector<char>& _refcs,
		uint32_t _oms = 0,
		uint8_t _mate = 0) :
		h(_h),
		patId(_patId),
		oms(_oms),
		fw(_fw),
		mate(_mate)
	{
		patName = _patName;
		patSeq  = _patSeq;
		quals   = _quals;
		mms     = _mms;
		refcs   = _refcs;
		// Enforce the constraints imposed by the binary output format
		if(seqan::length(patName) > 0xffff) {
			cerr << "Error: One or more read names are 2^16 characters or longer.  Please truncate" << endl
			     << "read names and re-run bowtie." << endl;
			exit(1);
		}
		if(mms.count() > 0xff) {
			cerr << "Error: The alignment contains 256 or more mismatches.  bowtie cannot handle" << endl
			     << "alignments with this many alignments.  Please provide smaller reads or consider" << endl
			     << "using a different tool." << endl;
			exit(1);
		}
		if(seqan::length(quals) > 0xffff) {
			cerr << "Error: One or more quality strings are 2^16 characters or longer.  Please" << endl
			     << "truncate reads and re-run bowtie." << endl;
			exit(1);
		}
		if(seqan::length(patSeq) > 0xffff) {
			cerr << "Error: One or more read sequences are 2^16 characters or longer.  Please" << endl
			     << "truncate reads and re-run bowtie." << endl;
			exit(1);
		}
	}

	U32Pair             h;       /// reference index & offset
	uint32_t            patId;   /// read index
	String<char>        patName; /// read name
	String<Dna5>        patSeq;  /// read sequence
	String<char>        quals;   /// read qualities
	FixedBitset<max_read_bp> mms;/// mismatch mask
	vector<char>        refcs;   /// reference characters for mms
	uint32_t            oms;     /// # of other possible mappings; 0 -> this is unique
	bool                fw;      /// orientation of read in alignment
	uint8_t             mate;    /// matedness; 0 = not a mate
	                             ///            1 = upstream mate
	                             ///            2 = downstream mate

	size_t length() const { return seqan::length(patSeq); }

	Hit& operator = (const Hit &other) {
	    this->h       = other.h;
	    this->patId   = other.patId;
		this->patName = other.patName;
		this->patSeq  = other.patSeq;
		this->quals   = other.quals;
	    this->mms     = other.mms;
	    this->refcs   = other.refcs;
	    this->oms     = other.oms;
		this->fw      = other.fw;
		this->mate    = other.mate;
	    return *this;
	}
};

/// Sort by text-id then by text-offset
bool operator< (const Hit& a, const Hit& b);

/**
 * Table for holding recalibration counts, along the lines of the table
 * presented in the SOAPsnp paper in Genome Res.  Each element maps a
 * read allele (o), quality score (q), cycle (c), and reference allele
 * (H), to a count of how many times that combination is observed in a
 * reported alignment.
 *
 * RecalTable is not synchronized, so it's assumed that threads are
 * incrementing the counters from within critical sections.
 */
class RecalTable {
public:
	RecalTable(int maxCycle,
	           int maxQual,
	           int qualShift) : maxCycle_(maxCycle),
	                            maxQual_(maxQual),
	                            qualShift_(qualShift),
	                            shift1_(6 - qualShift_),
	                            shift2_(shift1_ + 2),
	                            shift3_(shift2_ + 2),
	                            ents_(NULL), len_(0)
	{
		if(maxCycle == 0) {
			cerr << "Warning: maximum cycle for recalibration table is 0" << endl;
		} else if(maxQual >> qualShift == 0) {
			cerr << "Warning: maximum quality value " << maxQual << ", when shifted, is 0" << endl;
		} else if(qualShift > 5) {
			cerr << "Warning: quality shift value " << qualShift << " exceeds ceiling of 5" << endl;
		} else {
			try {
				len_ = maxCycle_ * 4 /* subj alleles*/ * 4 /* ref alleles */ * 64 /* quals */;
				ents_ = new uint32_t[len_];
				if(ents_ == NULL) {
					throw std::bad_alloc();
				}
				memset(ents_, 0, len_ << 2);
			} catch(std::bad_alloc& e) {
				cerr << "Error allocating recalibration table with " << len_ << " entries" << endl;
				exit(1);
			}
		}
	}

	~RecalTable() {
		if(ents_ != NULL) delete[] ents_;
	}

	/**
	 * Factor a new alignment into the recalibration table.
	 */
	void commitHit(const Hit& h) {
		// Iterate through the pattern from 5' to 3', calculate the
		// shifted quality value, obtain the reference character, and
		// increment the appropriate counter
		for(int i = 0; i < (int)h.length(); i++) {
			int ii = i;
			if(!h.fw) {
				ii = h.length() - ii - 1;
			}
			int qc = (int)h.patSeq[ii];
			int rc = qc;
			if(h.mms.test(i)) {
				rc = charToDna5[(int)h.refcs[i]];
				assert_neq(rc, qc);
			}
			int q = (int)h.quals[ii]-33;
			assert_lt(q, 64);
			q >>= qualShift_;
			ents_[calcIdx(i, qc, rc, q)]++;
		}
	}

	/**
	 * Print the contents of the recalibration table.
	 */
	void print (std::ostream& out) const {
		if(ents_ == NULL) return;
		const int lim = maxCycle_;
		for(int i = 0; i < lim; i++) {
			out << "t" << i << "\t";
			// Iterate over subject alleles
			for(int j = 0; j < 4; j++) {
				// Iterate over reference alleles
				for(int k = 0; k < 4; k++) {
					// Iterate over qualities
					int lim2 = maxQual_ >> qualShift_;
					for(int l = 0; l < lim2; l++) {
						out << ents_[calcIdx(i, j, k, l)] << '\t';
					}
				}
			}
			out << endl;
		}
	}

protected:

	/**
	 * Calculate index into the ents_ array given cycle, subject
	 * allele, reference allele, and (shifted) quality.
	 */
	int calcIdx(int cyc, int sa, int ra, int q) const {
		int ret = q | (ra << shift1_) | (sa << shift2_) | (cyc << shift3_);
		assert_lt(ret, len_);
		return ret;
	}

	const int maxCycle_;
	const int maxQual_;
	const int qualShift_;
	const int shift1_;
	const int shift2_;
	const int shift3_;
	uint32_t *ents_;
	int len_;
};

/**
 * Encapsulates an object that accepts hits, optionally retains them in
 * a vector, and does something else with them according to
 * descendent's implementation of pure virtual member reportHitImpl().
 */
class HitSink {
public:
	HitSink(OutFileBuf* out,
			const std::string& dumpAlignFaBasename,
			const std::string& dumpAlignFqBasename,
			const std::string& dumpAlignBasename,
			const std::string& dumpUnalFaBasename,
			const std::string& dumpUnalFqBasename,
			const std::string& dumpUnalBasename,
			const std::string& dumpMaxedFaBasename,
			const std::string& dumpMaxedFqBasename,
			const std::string& dumpMaxedBasename,
			bool onePairFile,
			RecalTable *table,
	        vector<string>* refnames = NULL) :
		_outs(),
		_deleteOuts(false),
		recalTable_(table),
		_refnames(refnames),
		_numWrappers(0),
		_locks(),
		dumpAlFaBase_(dumpAlignFaBasename),
		dumpAlFqBase_(dumpAlignFqBasename),
		dumpAlBase_(dumpAlignBasename),
		dumpUnalFaBase_(dumpUnalFaBasename),
		dumpUnalFqBase_(dumpUnalFqBasename),
		dumpUnalBase_(dumpUnalBasename),
		dumpMaxFaBase_(dumpMaxedFaBasename),
		dumpMaxFqBase_(dumpMaxedFqBasename),
		dumpMaxBase_(dumpMaxedBasename),
		onePairFile_(onePairFile),
		first_(true),
		numReported_(0llu),
		numReportedPaired_(0llu),
		quiet_(false),
		ssmode_(ios_base::out)
	{
		_outs.push_back(out);
		_locks.resize(1);
		MUTEX_INIT(_locks[0]);
		MUTEX_INIT(_mainlock);
		initDumps();
	}

	/**
	 * Open a number of output streams; usually one per reference
	 * sequence.  For now, we give then names refXXXXX.map where XXXXX
	 * is the 0-padded reference index.  Someday we may want to include
	 * the name of the reference sequence in the filename somehow.
	 */
	HitSink(size_t numOuts,
			const std::string& dumpAlignFaBasename,
			const std::string& dumpAlignFqBasename,
			const std::string& dumpAlignBasename,
			const std::string& dumpUnalFaBasename,
			const std::string& dumpUnalFqBasename,
			const std::string& dumpUnalBasename,
			const std::string& dumpMaxedFaBasename,
			const std::string& dumpMaxedFqBasename,
			const std::string& dumpMaxedBasename,
			bool onePairFile,
			RecalTable *table,
	        vector<string>* refnames = NULL) :
		_outs(),
		_deleteOuts(true),
		recalTable_(table),
		_refnames(refnames),
		_locks(),
		dumpAlFaBase_(dumpAlignFaBasename),
		dumpAlFqBase_(dumpAlignFqBasename),
		dumpAlBase_(dumpAlignBasename),
		dumpUnalFaBase_(dumpUnalFaBasename),
		dumpUnalFqBase_(dumpUnalFqBasename),
		dumpUnalBase_(dumpUnalBasename),
		dumpMaxFaBase_(dumpMaxedFaBasename),
		dumpMaxFqBase_(dumpMaxedFqBasename),
		dumpMaxBase_(dumpMaxedBasename),
		onePairFile_(onePairFile),
		quiet_(false),
		ssmode_(ios_base::out)
	{
		// Open all files for writing and initialize all locks
		for(size_t i = 0; i < numOuts; i++) {
			_outs.push_back(NULL); // we open output streams lazily
			_locks.resize(i+1);
			MUTEX_INIT(_locks[i]);
		}
		MUTEX_INIT(_mainlock);
   		initDumps();
	}

	/**
	 *
	 */
	virtual ~HitSink() {
		closeOuts();
		if(_deleteOuts) {
			// Delete all non-NULL output streams
			for(size_t i = 0; i < _outs.size(); i++) {
				if(_outs[i] != NULL) {
					delete _outs[i];
				}
			}
		}
		destroyDumps();
	}

	/**
	 * Call this whenever this HitSink is wrapped by a new
	 * HitSinkPerThread.  This helps us keep track of whether the main
	 * lock or any of the per-stream locks will be contended.
	 */
	void addWrapper() {
		_numWrappers++;
	}

	/**
	 * Called by concrete subclasses to figure out which elements of
	 * the _outs/_locks array to use when outputting the alignment.
	 */
	size_t refIdxToStreamIdx(size_t refIdx) {
		if(refIdx >= _outs.size()) return 0;
		return refIdx;
	}

	/// Implementation of hit-report
	virtual void reportHit(const Hit& h) = 0;

	/**
	 * Append a single hit to the given output stream.
	 */
	virtual void append(ostream& o, const Hit& h) = 0;

	/**
	 * Report a batch of hits.
	 */
	virtual void reportHits(vector<Hit>& hs) {
		size_t hssz = hs.size();
		if(hssz == 0) return;
		bool paired = hs[0].mate > 0;
		// Sort reads so that those against the same reference sequence
		// are consecutive.
		if(_outs.size() > 1 && hssz > 2) {
			sort(hs.begin(), hs.end());
		}
		char buf[4096];
		for(size_t i = 0; i < hssz; i++) {
			const Hit& h = hs[i];
			bool diff = false;
			if(i > 0) {
				diff = (refIdxToStreamIdx(h.h.first) != refIdxToStreamIdx(hs[i-1].h.first));
				if(diff) unlock(hs[i-1].h.first);
			}
			ostringstream ss(ssmode_);
			ss.rdbuf()->pubsetbuf(buf, 4096);
			append(ss, h);
			if(i == 0 || diff) {
				lock(h.h.first);
			}
			out(h.h.first).writeChars(buf, ss.tellp());
		}
		unlock(hs[hssz-1].h.first);
		mainlock();
		commitHits(hs);
		first_ = false;
		if(paired) numReportedPaired_ += hssz;
		else       numReported_ += hssz;
		mainunlock();
	}

	void commitHit(const Hit& hit) {
		if(recalTable_ != NULL) {
			recalTable_->commitHit(hit);
		}
	}

	void commitHits(const std::vector<Hit>& hits) {
		if(recalTable_ != NULL) {
			const size_t sz = hits.size();
			for(size_t i = 0; i < sz; i++) {
				commitHit(hits[i]);
			}
		}
	}

	/**
	 * Called when all alignments are complete.  It is assumed that no
	 * synchronization is necessary.
	 */
	void finish() {
		// Close output streams
		closeOuts();
		if(!quiet_) {
			// Print information about how many unpaired and/or paired
			// reads were aligned.
			if(first_) {
				assert_eq(0llu, numReported_);
				cerr << "No results" << endl;
			}
			else if(numReportedPaired_ > 0 && numReported_ == 0) {
				cerr << "Reported " << (numReportedPaired_ >> 1)
					 << " paired-end alignments to " << _outs.size()
					 << " output stream(s)" << endl;
			}
			else if(numReported_ > 0 && numReportedPaired_ == 0) {
				cerr << "Reported " << numReported_
					 << " alignments to " << _outs.size()
					 << " output stream(s)" << endl;
			}
			else {
				assert_gt(numReported_, 0);
				assert_gt(numReportedPaired_, 0);
				cerr << "Reported " << (numReportedPaired_ >> 1)
					 << " paired-end alignments and " << numReported_
					 << " singleton alignments to " << _outs.size()
					 << " output stream(s)" << endl;
			}
		}
		// Print the recalibration table.
		if(recalTable_ != NULL) {
			recalTable_->print(cout);
		}
	}

	/// Returns the alignment output stream; if the stream needs to be
	/// created, create it
	OutFileBuf& out(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
		if(_outs[strIdx] == NULL) {
			assert(_deleteOuts);
			ostringstream oss;
			oss << "ref";
			if     (strIdx < 10)    oss << "0000";
			else if(strIdx < 100)   oss << "000";
			else if(strIdx < 1000)  oss << "00";
			else if(strIdx < 10000) oss << "0";
			oss << strIdx << ".map";
			_outs[strIdx] = new OutFileBuf(oss.str().c_str());
		}
		assert(_outs[strIdx] != NULL);
		return *(_outs[strIdx]);
	}

	/**
	 * Lock the monolithic lock for this HitSink.  This is useful when,
	 * for example, outputting a read to an unaligned-read file.
	 */
	void mainlock() {
		MUTEX_LOCK(_mainlock);
	}

	/**
	 * Unlock the monolithic lock for this HitSink.  This is useful
	 * when, for example, outputting a read to an unaligned-read file.
	 */
	void mainunlock() {
		MUTEX_UNLOCK(_mainlock);
	}

	/**
	 * Return true iff this HitSink dumps aligned reads to an output
	 * stream (i.e., iff --alfa or --alfq are specified).
	 */
	bool dumpsAlignedReads() {
		return dumpAlignFlag_;
	}

	/**
	 * Return true iff this HitSink dumps unaligned reads to an output
	 * stream (i.e., iff --unfa or --unfq are specified).
	 */
	bool dumpsUnalignedReads() {
		return dumpUnalignFlag_;
	}

	/**
	 * Return true iff this HitSink dumps maxed-out reads to an output
	 * stream (i.e., iff --maxfa or --maxfq are specified).
	 */
	bool dumpsMaxedReads() {
		return dumpMaxedFlag_ || dumpUnalignFlag_;
	}

	/**
	 * Return true iff this HitSink dumps either unaligned or maxed-
	 * out reads to an output stream (i.e., iff --unfa, --maxfa,
	 * --unfq, or --maxfq are specified).
	 */
	bool dumpsReads() {
		return dumpAlignFlag_ || dumpUnalignFlag_ || dumpMaxedFlag_;
	}

	/**
	 * Dump an aligned read to all of the appropriate output streams.
	 * Be careful to synchronize correctly - there may be multiple
	 * simultaneous writers.
	 */
	void dumpAlign(PatternSourcePerThread& p) {
		if(!dumpAlignFlag_) return;
		if(!p.paired() || onePairFile_) {
			// Dump unpaired read to the aligned-read FASTA file
			if(!dumpAlFaBase_.empty()) {
				MUTEX_LOCK(dumpAlignFaLock_);
				if(dumpAlFa_ == NULL) {
					dumpAlFa_ = openOf(dumpAlFaBase_, 0, ".fa");
					assert(dumpAlFa_ != NULL);
				}
				printFastaRecord(*dumpAlFa_, p.bufa().name, p.bufa().patFw);
				MUTEX_UNLOCK(dumpAlignFaLock_);
			}
			// Dump unpaired read to the aligned-read FASTQ file
			if(!dumpAlFqBase_.empty()) {
				MUTEX_LOCK(dumpAlignFqLock_);
				if(dumpAlFq_ == NULL) {
					dumpAlFq_ = openOf(dumpAlFqBase_, 0, ".fq");
					assert(dumpAlFq_ != NULL);
				}
				printFastqRecord(*dumpAlFq_, p.bufa().name, p.bufa().patFw, p.bufa().qualFw);
				MUTEX_UNLOCK(dumpAlignFqLock_);
			}
			// Dump unpaired read to an aligned-read file of the same format
			if(!dumpAlBase_.empty()) {
				MUTEX_LOCK(dumpAlignLock_);
				if(dumpAl_ == NULL) {
					dumpAl_ = openOf(dumpAlBase_, 0, "");
					assert(dumpAl_ != NULL);
				}
				dumpAl_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				MUTEX_UNLOCK(dumpAlignLock_);
			}
		} else {
			// Dump paired-end read to the aligned-read FASTA file
			if(!dumpAlFaBase_.empty()) {
				MUTEX_LOCK(dumpAlignFaLockPE_);
				if(dumpAlFa_1_ == NULL) {
					assert(dumpAlFa_2_ == NULL);
					dumpAlFa_1_ = openOf(dumpAlFaBase_, 1, ".fa");
					dumpAlFa_2_ = openOf(dumpAlFaBase_, 2, ".fa");
					assert(dumpAlFa_1_ != NULL && dumpAlFa_2_ != NULL);
				}
				printFastaRecord(*dumpAlFa_1_, p.bufa().name, p.bufa().patFw);
				printFastaRecord(*dumpAlFa_2_, p.bufb().name, p.bufb().patFw);
				MUTEX_UNLOCK(dumpAlignFaLockPE_);
			}
			// Dump paired-end read to the aligned-read FASTQ file
			if(!dumpAlFqBase_.empty()) {
				MUTEX_LOCK(dumpAlignFqLockPE_);
				if(dumpAlFq_1_ == NULL) {
					assert(dumpAlFq_2_ == NULL);
					dumpAlFq_1_ = openOf(dumpAlFqBase_, 1, ".fq");
					dumpAlFq_2_ = openOf(dumpAlFqBase_, 2, ".fq");
					assert(dumpAlFq_1_ != NULL && dumpAlFq_2_ != NULL);
				}
				printFastqRecord(*dumpAlFq_1_, p.bufa().name, p.bufa().patFw, p.bufa().qualFw);
				printFastqRecord(*dumpAlFq_2_, p.bufb().name, p.bufb().patFw, p.bufb().qualFw);
				MUTEX_UNLOCK(dumpAlignFqLockPE_);
			}
			// Dump paired-end read to an aligned-read file (or pair of
			// files) of the same format
			if(!dumpAlBase_.empty()) {
				MUTEX_LOCK(dumpAlignLockPE_);
				if(dumpAl_1_ == NULL) {
					dumpAl_1_ = openOf(dumpAlBase_, 1, "");
		    		dumpAl_2_ = openOf(dumpAlBase_, 2, "");
					assert(dumpAl_1_ != NULL);
					assert(dumpAl_2_ != NULL);
				}
				dumpAl_1_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				dumpAl_2_->write(p.bufb().readOrigBuf, p.bufb().readOrigBufLen);
				MUTEX_UNLOCK(dumpAlignLockPE_);
			}
		}
	}

	/**
	 * Dump an unaligned read to all of the appropriate output streams.
	 * Be careful to synchronize correctly - there may be multiple
	 * simultaneous writers.
	 */
	void dumpUnal(PatternSourcePerThread& p) {
		if(!dumpUnalignFlag_) return;
		if(!p.paired() || onePairFile_) {
			if(!dumpUnalFaBase_.empty()) {
				MUTEX_LOCK(dumpUnalFaLock_);
				if(dumpUnalFa_ == NULL) {
					dumpUnalFa_ = openOf(dumpUnalFaBase_, 0, ".fa");
					assert(dumpUnalFa_ != NULL);
				}
				printFastaRecord(*dumpUnalFa_, p.bufa().name, p.bufa().patFw);
				MUTEX_UNLOCK(dumpUnalFaLock_);
			}
			if(!dumpUnalFqBase_.empty()) {
				MUTEX_LOCK(dumpUnalFqLock_);
				if(dumpUnalFq_ == NULL) {
					dumpUnalFq_ = openOf(dumpUnalFqBase_, 0, ".fq");
					assert(dumpUnalFq_ != NULL);
				}
				printFastqRecord(*dumpUnalFq_, p.bufa().name, p.bufa().patFw, p.bufa().qualFw);
				MUTEX_UNLOCK(dumpUnalFqLock_);
			}
			// Dump unpaired read to an unaligned-read file of the same format
			if(!dumpUnalBase_.empty()) {
				MUTEX_LOCK(dumpUnalLock_);
				if(dumpUnal_ == NULL) {
					dumpUnal_ = openOf(dumpUnalBase_, 0, "");
					assert(dumpUnal_ != NULL);
				}
				dumpUnal_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				MUTEX_UNLOCK(dumpUnalLock_);
			}
		} else {
			if(!dumpUnalFaBase_.empty()) {
				MUTEX_LOCK(dumpUnalFaLockPE_);
				if(dumpUnalFa_1_ == NULL) {
					assert(dumpUnalFa_2_ == NULL);
					dumpUnalFa_1_ = openOf(dumpUnalFaBase_, 1, ".fa");
					dumpUnalFa_2_ = openOf(dumpUnalFaBase_, 2, ".fa");
					assert(dumpUnalFa_1_ != NULL && dumpUnalFa_2_ != NULL);
				}
				printFastaRecord(*dumpUnalFa_1_, p.bufa().name, p.bufa().patFw);
				printFastaRecord(*dumpUnalFa_2_, p.bufb().name, p.bufb().patFw);
				MUTEX_UNLOCK(dumpUnalFaLockPE_);
			}
			if(!dumpUnalFqBase_.empty()) {
				MUTEX_LOCK(dumpUnalFqLockPE_);
				if(dumpUnalFq_1_ == NULL) {
					assert(dumpUnalFq_2_ == NULL);
					dumpUnalFq_1_ = openOf(dumpUnalFqBase_, 1, ".fq");
					dumpUnalFq_2_ = openOf(dumpUnalFqBase_, 2, ".fq");
					assert(dumpUnalFq_1_ != NULL && dumpUnalFq_2_ != NULL);
				}
				printFastqRecord(*dumpUnalFq_1_, p.bufa().name, p.bufa().patFw, p.bufa().qualFw);
				printFastqRecord(*dumpUnalFq_2_, p.bufb().name, p.bufb().patFw, p.bufb().qualFw);
				MUTEX_UNLOCK(dumpUnalFqLockPE_);
			}
			// Dump paired-end read to an unaligned-read file (or pair
			// of files) of the same format
			if(!dumpUnalBase_.empty()) {
				MUTEX_LOCK(dumpUnalLockPE_);
				if(dumpUnal_1_ == NULL) {
					dumpUnal_1_ = openOf(dumpUnalBase_, 1, "");
		    		dumpUnal_2_ = openOf(dumpUnalBase_, 2, "");
					assert(dumpUnal_1_ != NULL);
					assert(dumpUnal_2_ != NULL);
				}
				dumpUnal_1_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				dumpUnal_2_->write(p.bufb().readOrigBuf, p.bufb().readOrigBufLen);
				MUTEX_UNLOCK(dumpUnalLockPE_);
			}
		}
	}

	/**
	 * Dump a maxed-out read to all of the appropriate output streams.
	 * Be careful to synchronize correctly - there may be multiple
	 * simultaneous writers.
	 */
	void dumpMaxed(PatternSourcePerThread& p) {
		if(!dumpMaxedFlag_) {
			if(dumpUnalignFlag_) dumpUnal(p);
			return;
		}
		if(!p.paired() || onePairFile_) {
			if(!dumpMaxFaBase_.empty()) {
				MUTEX_LOCK(dumpMaxFaLock_);
				if(dumpMaxFa_ == NULL) {
					dumpMaxFa_ = openOf(dumpMaxFaBase_, 0, ".fa");
		    		assert(dumpMaxFa_ != NULL);
				}
				printFastaRecord(*dumpMaxFa_, p.bufa().name, p.bufa().patFw);
				MUTEX_UNLOCK(dumpMaxFaLock_);
			}
			if(!dumpMaxFqBase_.empty()) {
				MUTEX_LOCK(dumpMaxFqLock_);
				if(dumpMaxFq_ == NULL) {
		    		dumpMaxFq_ = openOf(dumpMaxFqBase_, 0, ".fq");
		    		assert(dumpMaxFq_ != NULL);
				}
				printFastqRecord(*dumpMaxFq_, p.bufa().name, p.bufa().patFw, p.bufa().qualFw);
				MUTEX_UNLOCK(dumpMaxFqLock_);
			}
			// Dump unpaired read to an maxed-out-read file of the same format
			if(!dumpMaxBase_.empty()) {
				MUTEX_LOCK(dumpMaxLock_);
				if(dumpMax_ == NULL) {
					dumpMax_ = openOf(dumpMaxBase_, 0, "");
					assert(dumpMax_ != NULL);
				}
				dumpMax_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				MUTEX_UNLOCK(dumpMaxLock_);
			}
		} else {
			if(!dumpMaxFaBase_.empty()) {
				MUTEX_LOCK(dumpMaxFaLockPE_);
				if(dumpMaxFa_1_ == NULL) {
					assert(dumpMaxFa_2_ == NULL);
		    		dumpMaxFa_1_ = openOf(dumpMaxFaBase_, 1, ".fa");
		    		dumpMaxFa_2_ = openOf(dumpMaxFaBase_, 2, ".fa");
		    		assert(dumpMaxFa_1_ != NULL && dumpMaxFa_2_ != NULL);
				}
				printFastaRecord(*dumpMaxFa_1_, p.bufa().name, p.bufa().patFw);
				printFastaRecord(*dumpMaxFa_2_, p.bufb().name, p.bufb().patFw);
				MUTEX_UNLOCK(dumpMaxFaLockPE_);
			}
			if(!dumpMaxFqBase_.empty()) {
				MUTEX_LOCK(dumpMaxFqLockPE_);
				if(dumpMaxFq_1_ == NULL) {
					assert(dumpMaxFq_2_ == NULL);
		    		dumpMaxFq_1_ = openOf(dumpMaxFqBase_, 1, ".fq");
		    		dumpMaxFq_2_ = openOf(dumpMaxFqBase_, 2, ".fq");
		    		assert(dumpMaxFq_1_ != NULL && dumpMaxFq_2_ != NULL);
				}
				printFastqRecord(*dumpMaxFq_1_, p.bufa().name, p.bufa().patFw, p.bufa().qualFw);
				printFastqRecord(*dumpMaxFq_2_, p.bufb().name, p.bufb().patFw, p.bufb().qualFw);
				MUTEX_UNLOCK(dumpMaxFqLockPE_);
			}
			// Dump paired-end read to a maxed-out-read file (or pair
			// of files) of the same format
			if(!dumpMaxBase_.empty()) {
				MUTEX_LOCK(dumpMaxLockPE_);
				if(dumpMax_1_ == NULL) {
					dumpMax_1_ = openOf(dumpMaxBase_, 1, "");
		    		dumpMax_2_ = openOf(dumpMaxBase_, 2, "");
					assert(dumpMax_1_ != NULL);
					assert(dumpMax_2_ != NULL);
				}
				dumpMax_1_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				dumpMax_2_->write(p.bufb().readOrigBuf, p.bufb().readOrigBufLen);
				MUTEX_UNLOCK(dumpMaxLockPE_);
			}
		}
	}

protected:

	/**
	 * Close (and flush) all OutFileBufs.
	 */
	void closeOuts() {
		// Flush and close all non-NULL output streams
		for(size_t i = 0; i < _outs.size(); i++) {
			if(_outs[i] != NULL && !_outs[i]->closed()) {
				_outs[i]->close();
			}
		}
	}

	/**
	 * Lock the output buffer for the output stream for reference with
	 * index 'refIdx'.  By default, hits for all references are
	 * directed to the same output stream, but if --refout is
	 * specified, each reference has its own reference stream.
	 */
	void lock(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
		MUTEX_LOCK(_locks[strIdx]);
	}

	/**
	 * Lock the output buffer for the output stream for reference with
	 * index 'refIdx'.  By default, hits for all references are
	 * directed to the same output stream, but if --refout is
	 * specified, each reference has its own reference stream.
	 */
	void unlock(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
		MUTEX_UNLOCK(_locks[strIdx]);
	}

	vector<OutFileBuf*> _outs;        /// the alignment output stream(s)
	bool                _deleteOuts;  /// Whether to delete elements of _outs upon exit
	RecalTable         *recalTable_;  /// recalibration table
	vector<string>*     _refnames;    /// map from reference indexes to names
	int                 _numWrappers; /// # threads owning a wrapper for this HitSink
	vector<MUTEX_T>     _locks;       /// pthreads mutexes for per-file critical sections
	MUTEX_T             _mainlock;    /// pthreads mutexes for fields of this object

	// Output filenames for dumping
	std::string dumpAlFaBase_;
	std::string dumpAlFqBase_;
	std::string dumpAlBase_;
	std::string dumpUnalFaBase_;
	std::string dumpUnalFqBase_;
	std::string dumpUnalBase_;
	std::string dumpMaxFaBase_;
	std::string dumpMaxFqBase_;
	std::string dumpMaxBase_;

	bool onePairFile_;

	// Output streams for dumping
    std::ofstream *dumpAlFa_;   // for single-ended reads
    std::ofstream *dumpAlFa_1_; // for first mates
    std::ofstream *dumpAlFa_2_; // for second mates
    std::ofstream *dumpAlFq_;   // for single-ended reads
    std::ofstream *dumpAlFq_1_; // for first mates
    std::ofstream *dumpAlFq_2_; // for second mates
    std::ofstream *dumpAl_;   // for single-ended reads
    std::ofstream *dumpAl_1_; // for first mates
    std::ofstream *dumpAl_2_; // for second mates
    std::ofstream *dumpUnalFa_;   // for single-ended reads
    std::ofstream *dumpUnalFa_1_; // for first mates
    std::ofstream *dumpUnalFa_2_; // for second mates
    std::ofstream *dumpUnalFq_;   // for single-ended reads
    std::ofstream *dumpUnalFq_1_; // for first mates
    std::ofstream *dumpUnalFq_2_; // for second mates
    std::ofstream *dumpUnal_;   // for single-ended reads
    std::ofstream *dumpUnal_1_; // for first mates
    std::ofstream *dumpUnal_2_; // for second mates
    std::ofstream *dumpMaxFa_;     // for single-ended reads
    std::ofstream *dumpMaxFa_1_;   // for first mates
    std::ofstream *dumpMaxFa_2_;   // for second mates
    std::ofstream *dumpMaxFq_;     // for single-ended reads
    std::ofstream *dumpMaxFq_1_;   // for first mates
    std::ofstream *dumpMaxFq_2_;   // for second mates
    std::ofstream *dumpMax_;     // for single-ended reads
    std::ofstream *dumpMax_1_;   // for first mates
    std::ofstream *dumpMax_2_;   // for second mates

    /**
     * Open an ofstream with given name; output error message and quit
     * if it fails.
     */
    std::ofstream* openOf(const std::string& name,
                          int mateType,
                          const std::string& suffix)
	{
    	std::string s = name;
		size_t dotoff = name.find_last_of(".");
    	if(mateType == 1) {
    		if(dotoff == string::npos) {
    			s += "_1"; s += suffix;
    		} else {
    			s = name.substr(0, dotoff) + "_1" + s.substr(dotoff);
    		}
    	} else if(mateType == 2) {
    		if(dotoff == string::npos) {
    			s += "_2"; s += suffix;
    		} else {
    			s = name.substr(0, dotoff) + "_2" + s.substr(dotoff);
    		}
    	} else if(mateType != 0) {
    		cerr << "Bad mate type " << mateType << endl; exit(1);
    	}
    	std::ofstream* tmp = new ofstream(s.c_str(), ios::out);
    	if(tmp->fail()) {
    		if(mateType == 0) {
    			cerr << "Could not open single-ended aligned/unaligned-read file for writing: " << name << endl;
    		} else {
    			cerr << "Could not open paired-end aligned/unaligned-read file for writing: " << name << endl;
    		}
    		exit(1);
    	}
    	return tmp;
    }

    /**
     * Initialize all the locks for dumping.
     */
    void initDumps() {
        dumpAlFa_     = dumpAlFa_1_   = dumpAlFa_2_   = NULL;
        dumpAlFq_     = dumpAlFq_1_   = dumpAlFq_2_   = NULL;
        dumpAl_       = dumpAl_1_     = dumpAl_2_     = NULL;
        dumpUnalFa_   = dumpUnalFa_1_ = dumpUnalFa_2_ = NULL;
        dumpUnalFq_   = dumpUnalFq_1_ = dumpUnalFq_2_ = NULL;
        dumpUnal_     = dumpUnal_1_   = dumpUnal_2_   = NULL;
        dumpMaxFa_    = dumpMaxFa_1_  = dumpMaxFa_2_  = NULL;
        dumpMaxFq_    = dumpMaxFq_1_  = dumpMaxFq_2_  = NULL;
        dumpMax_      = dumpMax_1_    = dumpMax_2_    = NULL;
    	dumpAlignFlag_   = !dumpAlFaBase_.empty() ||
    	                   !dumpAlFqBase_.empty() ||
    	                   !dumpAlBase_.empty();
    	dumpUnalignFlag_ = !dumpUnalFaBase_.empty() ||
    	                   !dumpUnalFqBase_.empty() ||
    	                   !dumpUnalBase_.empty();
    	dumpMaxedFlag_   = !dumpMaxFaBase_.empty() ||
    	                   !dumpMaxFqBase_.empty() ||
    	                   !dumpMaxBase_.empty();
   		MUTEX_INIT(dumpAlignFaLock_);
   		MUTEX_INIT(dumpAlignFaLockPE_);
   		MUTEX_INIT(dumpAlignFqLock_);
   		MUTEX_INIT(dumpAlignFqLockPE_);
   		MUTEX_INIT(dumpAlignLock_);
   		MUTEX_INIT(dumpAlignLockPE_);
   		MUTEX_INIT(dumpUnalFaLock_);
   		MUTEX_INIT(dumpUnalFaLockPE_);
   		MUTEX_INIT(dumpUnalFqLock_);
   		MUTEX_INIT(dumpUnalFqLockPE_);
   		MUTEX_INIT(dumpUnalLock_);
   		MUTEX_INIT(dumpUnalLockPE_);
   		MUTEX_INIT(dumpMaxFaLock_);
   		MUTEX_INIT(dumpMaxFaLockPE_);
   		MUTEX_INIT(dumpMaxFqLock_);
   		MUTEX_INIT(dumpMaxFqLockPE_);
   		MUTEX_INIT(dumpMaxLock_);
   		MUTEX_INIT(dumpMaxLockPE_);
    }

    void destroyDumps() {
    	if(dumpAlFa_     != NULL) { dumpAlFa_->close();     delete dumpAlFa_; }
    	if(dumpAlFa_1_   != NULL) { dumpAlFa_1_->close();   delete dumpAlFa_1_; }
    	if(dumpAlFa_2_   != NULL) { dumpAlFa_2_->close();   delete dumpAlFa_2_; }
    	if(dumpAlFq_     != NULL) { dumpAlFq_->close();     delete dumpAlFq_; }
    	if(dumpAlFq_1_   != NULL) { dumpAlFq_1_->close();   delete dumpAlFq_1_; }
    	if(dumpAlFq_2_   != NULL) { dumpAlFq_2_->close();   delete dumpAlFq_2_; }
    	if(dumpAl_       != NULL) { dumpAl_->close();       delete dumpAl_; }
    	if(dumpAl_1_     != NULL) { dumpAl_1_->close();     delete dumpAl_1_; }
    	if(dumpAl_2_     != NULL) { dumpAl_2_->close();     delete dumpAl_2_; }

    	if(dumpUnalFa_   != NULL) { dumpUnalFa_->close();   delete dumpUnalFa_; }
    	if(dumpUnalFa_1_ != NULL) { dumpUnalFa_1_->close(); delete dumpUnalFa_1_; }
    	if(dumpUnalFa_2_ != NULL) { dumpUnalFa_2_->close(); delete dumpUnalFa_2_; }
    	if(dumpUnalFq_   != NULL) { dumpUnalFq_->close();   delete dumpUnalFq_; }
    	if(dumpUnalFq_1_ != NULL) { dumpUnalFq_1_->close(); delete dumpUnalFq_1_; }
    	if(dumpUnalFq_2_ != NULL) { dumpUnalFq_2_->close(); delete dumpUnalFq_2_; }
    	if(dumpUnal_     != NULL) { dumpUnal_->close();     delete dumpUnal_; }
    	if(dumpUnal_1_   != NULL) { dumpUnal_1_->close();   delete dumpUnal_1_; }
    	if(dumpUnal_2_   != NULL) { dumpUnal_2_->close();   delete dumpUnal_2_; }

    	if(dumpMaxFa_    != NULL) { dumpMaxFa_->close();    delete dumpMaxFa_; }
    	if(dumpMaxFa_1_  != NULL) { dumpMaxFa_1_->close();  delete dumpMaxFa_1_; }
    	if(dumpMaxFa_2_  != NULL) { dumpMaxFa_2_->close();  delete dumpMaxFa_2_; }
    	if(dumpMaxFq_    != NULL) { dumpMaxFq_->close();    delete dumpMaxFq_; }
    	if(dumpMaxFq_1_  != NULL) { dumpMaxFq_1_->close();  delete dumpMaxFq_1_; }
    	if(dumpMaxFq_2_  != NULL) { dumpMaxFq_2_->close();  delete dumpMaxFq_2_; }
    	if(dumpMax_      != NULL) { dumpMax_->close();      delete dumpMax_; }
    	if(dumpMax_1_    != NULL) { dumpMax_1_->close();    delete dumpMax_1_; }
    	if(dumpMax_2_    != NULL) { dumpMax_2_->close();    delete dumpMax_2_; }
    }

    // Locks for dumping
    MUTEX_T dumpAlignFaLock_;
    MUTEX_T dumpAlignFaLockPE_; // _1 and _2
    MUTEX_T dumpAlignFqLock_;
    MUTEX_T dumpAlignFqLockPE_; // _1 and _2
    MUTEX_T dumpAlignLock_;
    MUTEX_T dumpAlignLockPE_; // _1 and _2
    MUTEX_T dumpUnalFaLock_;
    MUTEX_T dumpUnalFaLockPE_; // _1 and _2
    MUTEX_T dumpUnalFqLock_;
    MUTEX_T dumpUnalFqLockPE_; // _1 and _2
    MUTEX_T dumpUnalLock_;
    MUTEX_T dumpUnalLockPE_; // _1 and _2
    MUTEX_T dumpMaxFaLock_;
    MUTEX_T dumpMaxFaLockPE_;   // _1 and _2
    MUTEX_T dumpMaxFqLock_;
    MUTEX_T dumpMaxFqLockPE_;   // _1 and _2
    MUTEX_T dumpMaxLock_;
    MUTEX_T dumpMaxLockPE_;   // _1 and _2

    // false -> no dumping
    bool dumpAlignFlag_;
    bool dumpUnalignFlag_;
    bool dumpMaxedFlag_;

	volatile bool     first_;       /// true -> first hit hasn't yet been reported
	volatile uint64_t numReported_; /// # single-ended alignments reported
	volatile uint64_t numReportedPaired_; /// # paired-end alignments reported
	bool quiet_;  /// true -> don't print alignment stats at the end
	ios_base::openmode ssmode_;     /// output mode for stringstreams
};

/**
 * A per-thread wrapper for a HitSink.  Incorporates state that a
 * single search thread cares about.
 */
class HitSinkPerThread {
public:
	HitSinkPerThread(HitSink& sink, uint32_t __max = 0xffffffff, bool __keep = false) :
		_sink(sink),
		_bestRemainingStratum(0),
		_numReportableHits(0llu),
		_numValidHits(0llu),
		_keep(__keep),
		_hits(),
		_bufferedHits(),
		_strata(),
		hitsForThisRead_(),
		_max(__max)
	{
		_sink.addWrapper();
	}

	virtual ~HitSinkPerThread() { }

	/// Set whether to retain hits in a vector or not
	void setRetainHits(bool r)  { _keep = r; }
	/// Return whether we're retaining hits or not
	bool retainHits()           { return _keep; }

	/// Return the vector of retained hits
	vector<Hit>& retainedHits()   { return _hits; }
	vector<int>& retainedStrata() { return _strata; }

	/// Finalize current read
	virtual uint32_t finishRead(PatternSourcePerThread& p, bool dump = true) {
		uint32_t ret = finishReadImpl();
		_bestRemainingStratum = 0;
		if(dump && (ret == 0 || (ret > 0 && _bufferedHits.size() == 0))) {
			// Either no reportable hits were found or the number of
			// reportable hits exceeded the -m limit specified by the
			// user
			assert(ret == 0 || ret > _max);
			if(ret > 0) _sink.dumpMaxed(p);
			else        _sink.dumpUnal(p);
		}
		if(_bufferedHits.size() > 0) {
			// Flush buffered hits
			_sink.reportHits(_bufferedHits);
			_sink.dumpAlign(p);
			ret = _bufferedHits.size();
			_bufferedHits.clear();
		} else {
			ret = 0;
		}
		assert_eq(0, _bufferedHits.size());
		return ret;
	}

	virtual uint32_t finishReadImpl() = 0;

	/**
	 * Implementation for hit reporting; update per-thread _hits and
	 * _numReportableHits variables and call the master HitSink to do the actual
	 * reporting
	 */
	virtual void bufferHit(const Hit& h, int stratum) {
		_numReportableHits++;
		if(_keep) {
#ifndef NDEBUG
			// Check whether any of the previous 256 hits are identical
			// to this one.  Really, we should check them all but that
			// would be too slow for long runs.
			size_t len = _hits.size();
			size_t lookBack = 256;
			if(lookBack > len) lookBack = len;
			if(h.mate > 0) lookBack = 0;
			for(size_t i = 0; i < lookBack; i++) {
				if(h.patId    == _hits[len-i-1].patId &&
				   h.h.first  == _hits[len-i-1].h.first &&
				   h.h.second == _hits[len-i-1].h.second &&
				   h.fw       == _hits[len-i-1].fw)
				{
					cerr << "Repeat hit: " << h.patId << (h.fw? "+":"-")
					     << ":<" << h.h.first << "," << h.h.second << ","
					     << h.mms.count() << ">" << endl;
					assert(false); // same!
				}
			}
#endif
			_hits.push_back(h);
			_strata.push_back(stratum);
		}
#ifndef NDEBUG
		// Ensure all buffered hits have the same patid
		for(size_t i = 1; i < _bufferedHits.size(); i++) {
			assert_eq(_bufferedHits[0].patId, _bufferedHits[i].patId);
		}
#endif
		_bufferedHits.push_back(h);
	}

	/**
	 * Concrete subclasses override this to (possibly) report a hit and
	 * return true iff the caller should continue to report more hits.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		_numValidHits++;
		return true;
	}

	/// Return the number of valid hits so far (not necessarily
	/// reported).  It's up to the concrete subclasses
	uint64_t numValidHits()    { return _numValidHits; }

	/// Return the number of hits reported so far
	uint64_t numReportedHits() { return _numReportableHits; }

	/**
	 * Return true if there are no more reportable hits.
	 */
	bool finishedWithStratum(int stratum) {
		bool ret = finishedWithStratumImpl(stratum);
		_bestRemainingStratum = stratum+1;
		return ret;
	}

	/**
	 * Return true if there are no reportable hits with the given cost
	 * (or worse).
	 */
	virtual bool irrelevantCost(uint16_t cost) {
		return false;
	}

	/**
	 * Concrete subclasses override this to determine whether the
	 * search routine should keep searching after having finished
	 * reporting all alignments at the given stratum.
	 */
	virtual bool finishedWithStratumImpl(int stratum) = 0;

	/// Return the maximum number of hits allowed per read
	virtual uint32_t maxHits() = 0;

	/// The mhits maximum
	uint32_t overThresh() { return _max; }

	/// Whether this thread, for this read, knows that we have already
	/// exceeded the mhits maximum
	bool exceededOverThresh() { return hitsForThisRead_ > _max; }

	/// Return whether we span strata
	virtual bool spanStrata() = 0;

	/// Return whether we report only the best possible hits
	virtual bool best() = 0;

	/**
	 * Return true iff the underlying HitSink dumps unaligned or
	 * maxed-out reads.
	 */
	bool dumpsReads() {
		return _sink.dumpsReads();
	}

protected:
	HitSink&    _sink; /// Ultimate destination of reported hits
	/// Least # mismatches in alignments that will be reported in the
	/// future.  Updated by the search routine.
	int         _bestRemainingStratum;
	/// # hits reported to this HitSink so far (not all of which were
	/// necesssary reported to _sink)
	uint64_t    _numReportableHits;
	uint64_t    _numValidHits;
	bool        _keep; /// Whether to retain all reported hits in _hits
	vector<Hit> _hits; /// Repository for retained hits
	/// Buffered hits, to be reported and flushed at end of read-phase
	vector<Hit> _bufferedHits;
	vector<int> _strata; /// Repository for retained strata

	// Following variables are declared in the parent but maintained in
	// the concrete subcalsses
	uint32_t hitsForThisRead_; /// # hits for this read so far
	uint32_t _max;             /// don't report any hits if there were > _max
};

/**
 * Abstract parent factory for HitSinkPerThreads.
 */
class HitSinkPerThreadFactory {
public:
	virtual ~HitSinkPerThreadFactory() { }
	virtual HitSinkPerThread* create() const = 0;
	virtual HitSinkPerThread* createMult(uint32_t m) const = 0;

	/// Free memory associated with a per-thread hit sink
	virtual void destroy(HitSinkPerThread* sink) const {
		assert(sink != NULL);
		// Free the HitSinkPerThread
		delete sink;
	}
};

/**
 * Report first N good alignments encountered; trust search routine
 * to try alignments in something approximating a best-first order.
 * Best used in combination with a stringent alignment policy.
 */
class NGoodHitSinkPerThread : public HitSinkPerThread {

public:
	NGoodHitSinkPerThread(
			HitSink& sink,
	        uint32_t __n,
	        uint32_t __max,
	        bool __keep) :
	        HitSinkPerThread(sink, __max, __keep),
	        _n(__n)
	{
		assert_gt(_n, 0);
	}

	virtual uint32_t maxHits() { return _n; }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return false; // we settle for "good" hits
	}

	/// Finalize current read
	virtual uint32_t finishReadImpl() {
		uint32_t ret = hitsForThisRead_;
		hitsForThisRead_ = 0;
		return ret;
	}

	/**
	 * Report and then return true if we've already reported N good
	 * hits.  Ignore the stratum - it's not relevant for finding "good"
	 * hits.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		hitsForThisRead_++;
		if(hitsForThisRead_ > _max) {
			if(!_bufferedHits.empty()) _bufferedHits.clear();
			return true; // done - report nothing
		}
		if(hitsForThisRead_ <= _n) {
			// Only report hit if we haven't
			bufferHit(h, stratum);
		}
		if(hitsForThisRead_ == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits and max isn't set; stop!
		}
		return false; // not at N or max yet; keep going
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) { return false; }

private:
	uint32_t _n;               /// max # hits to report per read
};

/**
 * Concrete factory for FirstNGoodHitSinkPerThreads.
 */
class NGoodHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NGoodHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t n,
			uint32_t max = 0xffffffff,
			bool keep = false) :
			sink_(sink),
			n_(n),
			max_(max),
			keep_(keep)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NGoodHitSinkPerThread(sink_, n_, max_, keep_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		return new NGoodHitSinkPerThread(sink_, n_ * m, max_ * m, keep_);
	}

private:
	HitSink& sink_;
    uint32_t n_;
    uint32_t max_;
    bool keep_;
};


/**
 * Report the first N best alignments encountered.  Aggregator and
 * search routine collaborate to ensure that there exists no
 * alignment that's better than the ones reported.
 */
class NBestHitSinkPerThread : public HitSinkPerThread {

public:
	NBestHitSinkPerThread(
			HitSink& sink,
	        uint32_t __n,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
	        HitSinkPerThread(sink, __max, __keep),
	        _n(__n)
	{
		assert_gt(_n, 0);
	}

	virtual uint32_t maxHits() { return _n; }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		assert_geq(stratum, _bestRemainingStratum);
		if(stratum == _bestRemainingStratum) {
			// This hit is within th best possible remaining stratum,
			// so it should definitely count
			hitsForThisRead_++;
			if(hitsForThisRead_ > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			if(hitsForThisRead_ <= _n) {
				bufferHit(h, stratum);
			}
			if(hitsForThisRead_ == _n &&
			   (_max == 0xffffffff || _max < _n))
			{
				return true; // already reported N good hits; stop!
			}
		} else {
			// The alignment is not guaranteed to be in the best
			// stratum, but it might be, so save it for possible
			// future reporting
			hitStrata_[stratum].push_back(h);
		}
		return false; // not at N yet; keep going
	}

	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual uint32_t finishReadImpl() {
		if(hitsForThisRead_ == _n || hitsForThisRead_ > _max) {
			uint32_t ret = hitsForThisRead_;
			reset();
			return ret; // already reported N good hits; stop!
		}
		bool done = false;
		for(int j = 0; j < 4; j++) {
			for(size_t i = 0; i < hitStrata_[j].size(); i++) {
				// This hit is within th best possible remaining stratum,
				// so it should definitely count
				hitsForThisRead_++;
				if(hitsForThisRead_ > _max) {
					_bufferedHits.clear();
					done = true;
					break; // done - report nothing
				}
				if(hitsForThisRead_ <= _n) {
					bufferHit(hitStrata_[j][i], j);
				}
				if(hitsForThisRead_ == _n &&
				   (_max == 0xffffffff || _max < _n))
				{
					done = true;
					break; // already reported N good hits; stop!
				}
			}
			if(done) break;
			hitStrata_[j].clear();
		}
		uint32_t ret = hitsForThisRead_;
		reset();
		return ret;
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
#ifndef NDEBUG
		for(int j = 0; j < stratum; j++) {
			assert_eq(0, hitStrata_[j].size());
		}
#endif
		if(hitsForThisRead_ == _n) {
			return true; // already reported N good hits; stop!
		}
		for(int j = stratum; j <= stratum+1; j++) {
			for(size_t i = 0; i < hitStrata_[j].size(); i++) {
				// This hit is within the best possible remaining stratum,
				// so it should definitely count
				hitsForThisRead_++;
				if(hitsForThisRead_ > _max) {
					_bufferedHits.clear();
					return true; // done - report nothing
				}
				if(hitsForThisRead_ <= _n) {
					bufferHit(hitStrata_[j][i], j);
				}
				if(hitsForThisRead_ == _n &&
				   (_max == 0xffffffff || _max < _n))
				{
					hitStrata_[j].clear();
					return true; // already reported N good hits; stop!
				}
			}
			hitStrata_[j].clear();
		}
		return false; // keep going
	}

private:
	void reset() {
		clearAll();
		hitsForThisRead_ = 0;
	}

	void clearAll() {
		for(int j = 0; j < 4; j++) {
			hitStrata_[j].clear();
		}
	}

	uint32_t _n;               /// max # hits to report
	vector<Hit> hitStrata_[4]; /// lower numbered strata are better
};

/**
 * Concrete factory for NBestHitSinkPerThreads.
 */
class NBestHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NBestHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t n,
			uint32_t max = 0xffffffff,
			bool keep = false) :
			sink_(sink),
			n_(n),
			max_(max),
			keep_(keep)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NBestHitSinkPerThread(sink_, n_, max_, keep_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		return new NBestHitSinkPerThread(sink_, n_ * m, max_ * m, keep_);
	}

private:
	HitSink& sink_;
    uint32_t n_;
    uint32_t max_;
    bool keep_;
};

/**
 * Report the first N best alignments encountered in a single
 * alignment stratum (i.e., the best stratum in which there were
 * any hits).  Aggregator and search routine collaborate to ensure
 * that there exists no alignment that's better than the ones
 * reported.
 */
class NBestStratHitSinkPerThread : public HitSinkPerThread {

public:
	NBestStratHitSinkPerThread(
			HitSink& sink,
	        uint32_t __n,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
	        HitSinkPerThread(sink, __max, __keep),
	        _n(__n),
	        bestStratumReported_(999)
	{
		assert_gt(_n, 0);
	}

	virtual uint32_t maxHits() { return _n; }

	virtual bool spanStrata() {
		return false; // we do not span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		assert_geq(stratum, _bestRemainingStratum);
		if(stratum == _bestRemainingStratum &&
		   stratum <= bestStratumReported_)
		{
			// This hit is within th best possible remaining stratum,
			// so it should definitely count
			hitsForThisRead_++;
			if(hitsForThisRead_ > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			if(hitsForThisRead_ <= _n) {
				bufferHit(h, stratum);
			}
			bestStratumReported_ = stratum;
			if(hitsForThisRead_ == _n &&
			   (_max == 0xffffffff || _max < _n))
			{
				return true; // already reported N good hits; stop!
			}
		} else if(stratum <= bestStratumReported_) {
			hitStrata_[stratum].push_back(h);
			bestStratumReported_ = stratum;
		} else {
			// No need to buffer the hit because we're guaranteed not
			// to report it
		}
		return false; // not at N yet; keep going
	}

	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual uint32_t finishReadImpl() {
		if(bestStratumReported_ < 999 &&
		   hitsForThisRead_ < _n &&
		   hitsForThisRead_ <= _max)
		{
			assert_leq(bestStratumReported_, 3);
#ifndef NDEBUG
			// All better strata should be empty
			for(int j = 0; j < bestStratumReported_; j++) {
				assert_eq(0, hitStrata_[j].size());
			}
#endif
			for(size_t i = 0; i < hitStrata_[bestStratumReported_].size(); i++) {
				// This hit is within th best possible remaining stratum,
				// so it should definitely count
				hitsForThisRead_++;
				if(hitsForThisRead_ > _max) {
					_bufferedHits.clear();
					break; // done - report nothing
				}
				if(hitsForThisRead_ <= _n) {
					bufferHit(hitStrata_[bestStratumReported_][i], bestStratumReported_);
				}
				if(hitsForThisRead_ == _n &&
				   (_max == 0xffffffff || _max < _n))
				{
					hitStrata_[bestStratumReported_].clear();
					break; // already reported N good hits; stop!
				}
			}
			hitStrata_[bestStratumReported_].clear();
		} else if(bestStratumReported_ == 999) {
#ifndef NDEBUG
			// All strata should be empty
			for(int j = 0; j < 4; j++) {
				assert_eq(0, hitStrata_[j].size());
			}
#endif
		}
		uint32_t ret = hitsForThisRead_;
		reset();
		return ret;
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
		if(hitsForThisRead_ == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits; stop!
		}
		for(size_t i = 0; i < hitStrata_[stratum].size(); i++) {
			// This hit is within the best possible remaining stratum,
			// so it should definitely count
			hitsForThisRead_++;
			if(hitsForThisRead_ > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			if(hitsForThisRead_ <= _n) {
				bufferHit(hitStrata_[stratum][i], stratum);
			}
			if(hitsForThisRead_ == _n &&
			   (_max == 0xffffffff || _max < _n))
			{
				hitStrata_[stratum].clear();
				return true; // already reported N good hits; stop!
			}
		}
		hitStrata_[stratum].clear();
		// reported at least once in this stratum; do not move to the
		// next stratum
		if(hitsForThisRead_ > 0) return true;
		// didn't report; move to next stratum
		return false;
	}

private:

	void reset() {
		clearAll();
		hitsForThisRead_ = 0;
		bestStratumReported_ = 999;
	}

	void clearAll() {
		for(int j = 0; j < 4; j++) {
			hitStrata_[j].clear();
		}
	}

	uint32_t _n;               /// max # hits to report
	int bestStratumReported_;  /// stratum of best reported hit thus far
	vector<Hit> hitStrata_[4];
};

/**
 * Concrete factory for NBestStratHitSinkPerThread.
 */
class NBestStratHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NBestStratHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t n,
			uint32_t max = 0xffffffff,
			bool keep = false) :
			sink_(sink),
			n_(n),
			max_(max),
			keep_(keep)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NBestStratHitSinkPerThread(sink_, n_, max_, keep_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		return new NBestStratHitSinkPerThread(sink_, n_ * m, max_ * m, keep_);
	}

private:
	HitSink& sink_;
    uint32_t n_;
    uint32_t max_;
    bool keep_;
};

/**
 * Report the first N best alignments encountered in a single
 * alignment stratum assuming that we're receiving the alignments in
 * best-first order.
 */
class NBestFirstStratHitSinkPerThread : public HitSinkPerThread {

public:
	NBestFirstStratHitSinkPerThread(
			HitSink& sink,
	        uint32_t n,
	        uint32_t max,
	        bool keep,
	        uint32_t mult) :
	        HitSinkPerThread(sink, max * mult, keep),
	        n_(n * mult), bestStratum_(999), mult_(mult)
	{
		assert_gt(n_, 0);
	}

	virtual uint32_t maxHits() { return n_; }

	virtual bool spanStrata() {
		return false; // we do not span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		// This hit is within th best possible remaining stratum,
		// so it should definitely count
		hitsForThisRead_++;
		// It doesn't exceed the limit, so buffer it
		if(stratum < bestStratum_) {
			bestStratum_ = stratum;
		}
		if(hitsForThisRead_ > _max) {
			_bufferedHits.clear();
			return true; // done - report nothing
		}
		if(hitsForThisRead_ <= n_) {
			bufferHit(h, stratum);
		}
		if(hitsForThisRead_ == n_ &&
		   (_max == 0xffffffff || _max < n_))
		{
			return true; // already reported N good hits; stop!
		}
		return false; // not at N yet; keep going
	}

	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual uint32_t finishReadImpl() {
		uint32_t ret = hitsForThisRead_;
		hitsForThisRead_ = 0;
		bestStratum_ = 999;
		const size_t sz = _bufferedHits.size();
		for(size_t i = 0; i < sz; i++) {
			// Set 'oms' according to the number of other alignments
			// at this stratum
			_bufferedHits[i].oms = (sz / mult_) - 1;
		}
		return ret;
	}

	/**
	 * If we had any alignments at all and we're now moving on to a new
	 * stratum, then we're done.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
		if(hitsForThisRead_ > 0) {
			return true;
		}
		return false;
	}

	virtual bool irrelevantCost(uint16_t cost) {
		if(hitsForThisRead_) {
			return ((int)cost >> 14) > bestStratum_;
		}
		return false;
	}

private:

	uint32_t n_; /// max # hits to report
	int bestStratum_;
	uint32_t mult_;
};

/**
 * Concrete factory for NBestStratHitSinkPerThread.
 */
class NBestFirstStratHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NBestFirstStratHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t n,
			uint32_t max = 0xffffffff,
			bool keep = false) :
			sink_(sink),
			n_(n),
			max_(max),
			keep_(keep)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NBestFirstStratHitSinkPerThread(sink_, n_, max_, keep_, 1);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		return new NBestFirstStratHitSinkPerThread(sink_, n_, max_, keep_, m);
	}

private:
	HitSink& sink_;
    uint32_t n_;
    uint32_t max_;
    bool keep_;
};

/**
 * Report all valid alignments.
 */
class AllHitSinkPerThread : public HitSinkPerThread {

public:
	AllHitSinkPerThread(
			HitSink& sink,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
		    HitSinkPerThread(sink, __max, __keep) { }

	virtual uint32_t maxHits() { return 0xffffffff; }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and always return true; we're finiding all hits so that
	 * search routine should always continue.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		hitsForThisRead_++;
		if(hitsForThisRead_ > _max) {
			if(!_bufferedHits.empty()) _bufferedHits.clear();
			return true; // done - report nothing
		}
		bufferHit(h, stratum);
		return false; // reporting all; always keep going
	}

	/**
	 * Finalize; do nothing because we haven't buffered anything
	 */
	virtual uint32_t finishReadImpl() {
		uint32_t ret = hitsForThisRead_;
		hitsForThisRead_ = 0;
		return ret;
	}

	/**
	 * Always return false; search routine should not stop.
	 */
	virtual bool finishedWithStratumImpl(int stratum) { return false; }
};

/**
 * Concrete factory for AllHitSinkPerThread.
 */
class AllHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	AllHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t max = 0xffffffff,
			bool keep = false) :
			sink_(sink),
			max_(max),
			keep_(keep)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new AllHitSinkPerThread(sink_, max_, keep_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		return new AllHitSinkPerThread(sink_, max_ * m, keep_);
	}

private:
	HitSink& sink_;
    uint32_t max_;
    bool keep_;
};

/**
 * Report all alignments encountered in the best alignment stratum
 * for which there were any alignments.
 */
class AllStratHitSinkPerThread : public HitSinkPerThread {

public:
	AllStratHitSinkPerThread(
			HitSink& sink,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
	        HitSinkPerThread(sink, __max, __keep),
	        bestStratumReported_(999),
	        reported_(false) { }

	virtual uint32_t maxHits() { return 0xffffffff; }

	virtual bool spanStrata() {
		return false; // we do not span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		assert_geq(stratum, _bestRemainingStratum);
		if(stratum == _bestRemainingStratum &&
		   stratum <= bestStratumReported_)
		{
			// This hit is within the best possible remaining stratum,
			// so it should definitely be reported
			hitsForThisRead_++;
			if(hitsForThisRead_ > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			bufferHit(h, stratum);
			reported_ = true;
			bestStratumReported_ = stratum;
		} else if(stratum <= bestStratumReported_) {
			hitStrata_[stratum].push_back(h);
			bestStratumReported_ = stratum;
		} else {
			// No need to buffer the hit because we're guaranteed not
			// to report it
		}
		return false; // keep going
	}

	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual uint32_t finishReadImpl() {
		if(bestStratumReported_ < 999) {
#ifndef NDEBUG
			// All better strata should be empty
			for(int j = 0; j < bestStratumReported_; j++) {
				assert_eq(0, hitStrata_[j].size());
			}
#endif
			for(size_t i = 0; i < hitStrata_[bestStratumReported_].size(); i++) {
				// This hit is within th best possible remaining stratum,
				// so it should definitely count
				hitsForThisRead_++;
				if(hitsForThisRead_ > _max) {
					_bufferedHits.clear();
					break; // done - report nothing
				}
				bufferHit(hitStrata_[bestStratumReported_][i], bestStratumReported_);
			}
		} else {
#ifndef NDEBUG
			// All strata should be empty
			for(int j = 0; j < 4; j++) {
				assert_eq(0, hitStrata_[j].size());
			}
#endif
		}
		uint32_t ret = hitsForThisRead_;
		reset();
		return ret;
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
#ifndef NDEBUG
		// All better strata should be empty
		for(int j = 0; j < stratum; j++) {
			assert_eq(0, hitStrata_[j].size());
		}
#endif
		for(size_t i = 0; i < hitStrata_[stratum].size(); i++) {
			// This hit is within the best possible remaining stratum,
			// so it should definitely count
			hitsForThisRead_++;
			if(hitsForThisRead_ > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			bufferHit(hitStrata_[stratum][i], stratum);
			reported_ = true;
		}
		hitStrata_[stratum].clear();
		// reported at least once in this stratum; do not move to the
		// next stratum
		if(reported_) return true; // stop
		// didn't report; move to next stratum
		return false;
	}

private:

	void reset() {
		clearAll();
		reported_ = false;
		bestStratumReported_ = 999;
		hitsForThisRead_ = 0;
	}

	void clearAll() {
		for(int j = 0; j < 4; j++) {
			hitStrata_[j].clear();
		}
	}

	int bestStratumReported_;          /// stratum of best reported hit thus far
	bool reported_;
	vector<Hit> hitStrata_[4];
};

/**
 * Concrete factory for AllStratHitSinkPerThread.
 */
class AllStratHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	AllStratHitSinkPerThreadFactory(
			HitSink& sink,
			uint32_t max = 0xffffffff,
			bool keep = false) :
			sink_(sink),
			max_(max),
			keep_(keep)
	{ }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new AllStratHitSinkPerThread(sink_, max_, keep_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		return new AllStratHitSinkPerThread(sink_, max_ * m, keep_);
	}

private:
	HitSink& sink_;
    uint32_t max_;
    bool keep_;
};

#define DECL_HIT_DUMPS \
	const std::string& dumpAlFa, \
	const std::string& dumpAlFq, \
	const std::string& dumpAl, \
	const std::string& dumpUnalFa, \
	const std::string& dumpUnalFq, \
	const std::string& dumpUnal, \
	const std::string& dumpMaxFa, \
	const std::string& dumpMaxFq, \
	const std::string& dumpMax, \
	bool onePairFile, \
	RecalTable *recalTable, \
	std::vector<std::string>* refnames

#define PASS_HIT_DUMPS \
	dumpAlFa, \
	dumpAlFq, \
	dumpAl, \
	dumpUnalFa, \
	dumpUnalFq, \
	dumpUnal, \
	dumpMaxFa, \
	dumpMaxFq, \
	dumpMax, \
	onePairFile, \
	recalTable, \
	refnames

/**
 * Sink that prints lines like this:
 * (pat-id)[-|+]:<hit1-text-id,hit2-text-offset>,<hit2-text-id...
 *
 * Activated with --concise
 */
class ConciseHitSink : public HitSink {
public:
	/**
	 * Construct a single-stream ConciseHitSink (default)
	 */
	ConciseHitSink(OutFileBuf* out,
			       int offBase,
	               DECL_HIT_DUMPS,
	               bool reportOpps = false) :
		HitSink(out, PASS_HIT_DUMPS),
		_reportOpps(reportOpps),
		offBase_(offBase) { }

	/**
	 * Construct a multi-stream ConciseHitSink with one stream per
	 * reference string (see --refout)
	 */
	ConciseHitSink(size_t numOuts,
	               int offBase,
	               DECL_HIT_DUMPS,
	               bool reportOpps = false) :
		HitSink(numOuts, PASS_HIT_DUMPS),
		_reportOpps(reportOpps),
		offBase_(offBase) { }

	/**
	 * Append a verbose, readable hit to the given output stream.
	 */
	static void append(ostream& ss, const Hit& h, int offBase, bool reportOpps) {
		ss << h.patId;
		if(h.mate > 0) {
			assert(h.mate == 1 || h.mate == 2);
			ss << '/' << (int)h.mate;
		}
		ss << (h.fw? "+" : "-") << ":";
    	// .first is text id, .second is offset
		ss << "<" << h.h.first << "," << (h.h.second + offBase) << "," << h.mms.count();
		if(reportOpps) ss << "," << h.oms;
		ss << ">" << endl;
	}

	/**
	 * Append a verbose, readable hit to the given output stream.
	 */
	void append(ostream& ss, const Hit& h) {
		ConciseHitSink::append(ss, h, this->offBase_, this->_reportOpps);
	}

	/**
	 * Report a concise alignment to the appropriate output stream.
	 */
	virtual void reportHit(const Hit& h) {
		ostringstream ss;
		append(ss, h);
		lock(h.h.first);
		out(h.h.first).writeString(ss.str());
		unlock(h.h.first);
		mainlock();
		commitHit(h);
		first_ = false;
		if(h.mate > 0) numReportedPaired_++;
		else           numReported_++;
		mainunlock();
	}

private:
	bool _reportOpps;
	int  offBase_;     /// Add this to reference offsets before outputting.
	                   /// (An easy way to make things 1-based instead of
	                   /// 0-based)
};

/**
 * Sink that prints lines like this:
 * pat-name \t [-|+] \t ref-name \t ref-off \t pat \t qual \t #-alt-hits \t mm-list
 */
class VerboseHitSink : public HitSink {
public:
	/**
	 * Construct a single-stream VerboseHitSink (default)
	 */
	VerboseHitSink(OutFileBuf* out,
	               int offBase,
	               DECL_HIT_DUMPS,
				   int partition = 0) :
	HitSink(out, PASS_HIT_DUMPS),
	_partition(partition),
	offBase_(offBase)
	{ }

	/**
	 * Construct a multi-stream VerboseHitSink with one stream per
	 * reference string (see --refout)
	 */
	VerboseHitSink(size_t numOuts,
	               int offBase,
	               DECL_HIT_DUMPS,
				   int partition = 0) :
	HitSink(numOuts, PASS_HIT_DUMPS),
	_partition(partition),
	offBase_(offBase)
	{ }

	/**
	 * Given a line of output from the VerboseHitSink, parse it into a
	 * Hit object and return the object (by value).  If the reference
	 * and/or read are identified in the output by name, then the patid
	 * and/or h.first fields may be invalid.  TODO: load the id/name
	 * mapping from the .ebwt file so that we can handle either case.
	 */
	static bool readHit(Hit& h,
	                    istream& in,
	                    vector<string>* refnames,
	                    bool verbose = false)
	{
		char buf[4096];
		if(in.eof()) {
			return false;
		}
		if(in.bad()) {
			cerr << "Alignment file set \"bad\" bit" << endl;
			exit(1);
		}
		if(in.fail()) {
			cerr << "A line from the alignment file was longer than 4K" << endl;
			exit(1);
		}
		in.getline(buf, 4096);
		size_t len = in.gcount();
		if(len < 11) {
			// Can't possibly be a well-formed line if it's this short
			return false;
		}

		// Push cur pointer past partition label, if it exists
		char *cur = buf;
		bool sawSpace = true;
		while(*cur != ' ') {
			if(*cur == '\t') {
				sawSpace = false;
				break;
			}
			cur++;
			assert_lt((size_t)(cur - buf), len);
		}
		if(sawSpace) {
			// Keep moving to past the next tab
			while(*cur != '\t') cur++;
			assert_eq('\t', *cur);
			cur++;
		} else {
			cur = buf; // rewind
		}

		// Parse read name and check whether it's totally numeric
		bool readNameIsIdx = true;
		h.patId = 0;
		char *readName = cur;
		while(*cur != '\t') {
			if(*cur < '0' || *cur > '9') {
				readNameIsIdx = false;
			} else if(readNameIsIdx) {
				h.patId *= 10;
				h.patId += (*cur - '0');
			}
			assert_lt((size_t)(cur - readName), len);
			cur++;
		}
		*cur = '\0';
		size_t readNameLen = cur - readName;
		cur++;

		// Copy read name into h.patName
		seqan::resize(h.patName, readNameLen);
		char *patName = (char *)h.patName.data_begin;
		for(size_t i = 0; i < readNameLen; i++) {
			patName[i] = readName[i];
		}
		h.mate = 0;
		if(readNameLen >= 2 && patName[readNameLen-2] == '/') {
			if     (patName[readNameLen-1] == '1') h.mate = 1;
			else if(patName[readNameLen-1] == '2') h.mate = 2;
		}

		// Parse orientation
		assert(*cur == '+' || *cur == '-');
		h.fw = (*cur == '+');
		cur++;
		assert_eq('\t', *cur);
		cur++;

		// Parse reference sequence id
		bool refIsIdx = true;
		uint32_t refIdx = 0;
		char *refName = cur;
		while(*cur != '\t') {
			if(*cur < '0' || *cur > '9') {
				refIsIdx = false;
			} else if(readNameIsIdx) {
				refIdx *= 10;
				refIdx += (*cur - '0');
			}
			assert_lt((size_t)(cur - readName), len);
			cur++;
		}
		*cur = '\0';
		//size_t refNameLen = (size_t)(cur - refName);
		cur++;

		if(!refIsIdx && refnames != NULL) {
			bool found = false;
			for(size_t i = 0; i < refnames->size(); i++) {
				if((*refnames)[i] == refName) {
					found = true;
					refIdx = i;
					break;
				}
			}
			cerr << "Could not find an id to map reference name \"" << refName << "\" to." << endl;
			exit(1);
		}

		// Parse reference sequence offset
		uint32_t refOff = 0;
		ASSERT_ONLY(char *refOffStr = cur);
		while(*cur != '\t') {
			if(*cur < '0' || *cur > '9') {
				assert(false);
			}
			refOff *= 10;
			refOff += (*cur - '0');
			assert_leq((size_t)(cur - refOffStr), len);
			cur++;
		}
		*cur = '\0'; // terminate reference-offset sequence
		//size_t refOffLen = (size_t)(cur - refOffStr);
		cur++;

		// Fill in h.h now that we have refIdx and refOff
		h.h = make_pair<uint32_t>(refIdx, refOff);

		// Parse read sequence
		char *readSeq = cur;
		while(*cur != '\t') {
			cur++;
			assert_leq((size_t)(cur - refOffStr), len);
		}
		*cur = '\0'; // terminate read sequence
		size_t readSeqLen = cur - readSeq;
		assert_gt(readSeqLen, 0);
		cur++;

		// Copy read sequence into h.patSeq
		seqan::resize(h.patSeq, readSeqLen);
		uint8_t *readSeqDest = (uint8_t *)h.patSeq.data_begin;
		for(size_t i = 0; i < readSeqLen; i++) {
			assert_neq(0, dna4Cat[(int)readSeq[i]]);
			readSeqDest[i] = charToDna5[(int)readSeq[i]];
		}

		// Parse read qualities
		char *readQuals = cur;
		while(*cur != '\t') {
			cur++;
			assert_leq((size_t)(cur - refOffStr), len);
		}
		*cur = '\0'; // terminate read sequence
		size_t readQualsLen = cur - readQuals;
		assert_gt(readQualsLen, 0);
		cur++;
		assert_eq(readSeqLen, readQualsLen);

		// Copy quality values into h.quals
		seqan::resize(h.quals, readSeqLen);
		char *readQualsDest = (char *)h.quals.data_begin;
		for(size_t i = 0; i < readQualsLen; i++) {
			assert_geq(readQuals[i], 33);
			readQualsDest[i] = readQuals[i];
		}

		// Parse # other hits at this stratum (an underestimate)
		ASSERT_ONLY(char *omsStr = cur);
		h.oms = 0;
		// (the oms field is always followed by a tab, even when the
		// (following) mismatch field is empty
		while(*cur != '\t') {
			// Must be a number
			assert(*cur >= '0' && *cur <= '9');
			h.oms *= 10;
			h.oms += (*cur - '0');
			cur++;
			assert_leq((size_t)(cur - omsStr), len);
		}
		*cur = '\0'; // terminate read sequence
		ASSERT_ONLY(size_t omsLen = cur - omsStr);
		assert_gt(omsLen, 0);
		cur++;

		// Parse the # other hits at this stratum estimate
		h.refcs.resize(readSeqLen, 0);
		assert_eq(readSeqLen, h.refcs.size());

		// (h.mm is fixed-width so we don't need to resize it)
		if(*cur == '\0') {
			return true;
		}
		assert(*cur >= '0' && *cur <= '9');
		while(true) {
			uint32_t i = 0;
			ASSERT_ONLY(char *offStr = cur);
			while(*cur != ':') {
				// Must be a number
				assert(*cur >= '0' && *cur <= '9');
				i *= 10;
				i += (*cur - '0');
				cur++;
				assert_leq((size_t)(cur - offStr), len);
			}
			assert_lt(i, readSeqLen);
			cur++; // now points to reference base
			assert_eq(1, (int)dna4Cat[(int)(*cur)]);
			ASSERT_ONLY(int rc = *cur);
			h.refcs[i] = *cur;
			h.mms.set(i);
			cur++; assert_eq('>', *cur);
			cur++;
			assert_gt((int)dna4Cat[(int)(*cur)], 0);
			assert_neq((int)dna4Cat[(int)(*cur)], rc);
			if(h.fw) {
				assert_eq("ACGTN"[(int)h.patSeq[i]], *cur);
			} else {
				assert_eq("ACGTN"[(int)h.patSeq[readSeqLen - i - 1]], *cur);
			}
			cur++;
			if(*cur == '\0') {
				break;
			}
			assert_eq(',', *cur);
			cur++;
		}
		return true;
	}

	/**
	 * Append a verbose, readable hit to the given output stream.
	 */
	static void append(ostream& ss,
	                   const Hit& h,
	                   const vector<string>* refnames,
	                   int partition,
	                   int offBase)
	{
		bool spill = false;
		int spillAmt = 0;
		uint32_t pdiv = 0xffffffff;
		uint32_t pmod = 0xffffffff;
		do {
			bool dospill = false;
			if(spill) {
				// The read spilled over a partition boundary and so
				// needs to be printed more than once
				assert(partition > 0);
				spill = false;
				dospill = true;
				spillAmt++;
			}
			assert(!spill);
			if(partition > 0) {
				// Output a partitioning key
				// First component of the key is the reference index
				ss << h.h.first << "\t";
				ostringstream ss2, ss3;
				// Next component of the key is the partition id
				if(!dospill) {
					pdiv = (h.h.second + offBase) / partition;
					pmod = (h.h.second + offBase) % partition;
				}
				assert_neq(0xffffffff, pdiv);
				assert_neq(0xffffffff, pmod);
				if(dospill) assert_gt(spillAmt, 0);
				ss2 << (pdiv + (dospill ? spillAmt : 0));
				if((pmod + h.length()) >= ((uint32_t)partition * (spillAmt + 1))) {
					// Spills into the next partition so we need to
					// output another alignment for that partition
					spill = true;
				}
				// Print partition id with leading 0s so that Hadoop
				// can do lexicographical sort (modern Hadoop versions
				// seen to support numeric)
				string s2 = ss2.str();
				size_t partDigits = 1;
				if(partition >= 10) partDigits++;
				if(partition >= 100) partDigits++;
				if(partition >= 1000) partDigits++;
				if(partition >= 10000) partDigits++;
				if(partition >= 100000) partDigits++;
				for(size_t i = s2.length(); i < (10-partDigits); i++) {
					ss << "0";
				}
				ss << s2.c_str() << "\t";
				// Print offset with leading 0s
				ss3 << (h.h.second + offBase);
				string s3 = ss3.str();
				for(size_t i = s3.length(); i < 9; i++) {
					ss << "0";
				}
				ss << s3;
				ss << "\t" << (h.fw? "+":"-");
				ss << "\t" << h.patSeq;
				ss << "\t" << h.quals;
				ss << "\t" << h.oms;
			} else {
				assert(!dospill);
				ss << h.patName << "\t" << (h.fw? "+":"-") << "\t";
				// .first is text id, .second is offset
				if(refnames != NULL && h.h.first < refnames->size()) {
					ss << (*refnames)[h.h.first];
				} else {
					ss << h.h.first;
				}
				ss << "\t" << (h.h.second + offBase);
				ss << "\t" << h.patSeq;
				ss << "\t" << h.quals;
				ss << "\t" << h.oms;
				ss << "\t";
				// Output mismatch column
				bool firstmiss = true;
				size_t c = 0;
				for (unsigned int i = 0; i < h.mms.size(); ++ i) {
					if (h.mms.test(i)) {
						if (!firstmiss) ss << ",";
						ss << i;
						if(h.refcs.size() > 0) {
							assert_gt(h.refcs.size(), i);
							ASSERT_ONLY(char cc = toupper(h.refcs[i]));
							assert(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T');
							char refChar = toupper(h.refcs[i]);
							char qryChar = (h.fw ? h.patSeq[i] : h.patSeq[length(h.patSeq)-i-1]);
							assert_neq(refChar, qryChar);
							ss << ":" << refChar << ">" << qryChar;
						}
						firstmiss = false;
						c++;
					}
				}
			}
			ss << endl;
		} while(spill);
	}

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(ostream& ss, const Hit& h) {
		VerboseHitSink::append(ss, h, this->_refnames, this->_partition, this->offBase_);
	}


	/**
	 * Report a verbose, human-readable alignment to the appropriate
	 * output stream.
	 */
	virtual void reportHit(const Hit& h) {
		ostringstream ss;
		append(ss, h);
		// Make sure to grab lock before writing to output stream
		lock(h.h.first);
		out(h.h.first).writeString(ss.str());
		unlock(h.h.first);
		mainlock();
		commitHit(h);
		first_ = false;
		if(h.mate > 0) numReportedPaired_++;
		else           numReported_++;
		mainunlock();
	}

private:
	int      _partition;   /// partition size, or 0 if partitioning is disabled
	int      offBase_;     /// Add this to reference offsets before outputting.
	                       /// (An easy way to make things 1-based instead of
	                       /// 0-based)
};

/**
 * Sink for outputting alignments in a binary format.
 */
class BinaryHitSink : public HitSink {
public:

	/**
	 * Construct a single-stream BinaryHitSink (default)
	 */
	BinaryHitSink(OutFileBuf* out,
	              int offBase,
	              DECL_HIT_DUMPS) :
	HitSink(out, PASS_HIT_DUMPS),
	offBase_(offBase)
	{
		ssmode_ |= ios_base::binary;
	}

	/**
	 * Construct a multi-stream BinaryHitSink with one stream per
	 * reference string (see --refout)
	 */
	BinaryHitSink(size_t numOuts,
	              int offBase,
	              DECL_HIT_DUMPS) :
	HitSink(numOuts, PASS_HIT_DUMPS),
	offBase_(offBase)
	{
		ssmode_ |= ios_base::binary;
	}

	/**
	 * Append a binary alignment to the output stream corresponding to
	 * the reference sequence involved.
	 */
	static void append(ostream& o,
					   const Hit& h,
					   const vector<string>* refnames,
					   int offBase)
	{
		// true iff we're going to print the reference sequence name
		bool refName = refnames != NULL &&
		                h.h.first < refnames->size();
		uint16_t pnamelen = (uint16_t)length(h.patName);
		// Write read name
		o.write((const char *)&pnamelen, 2);
		o.write(begin(h.patName), pnamelen);
		// Write fw/refname flags
		assert_lt(h.mate, 3);
		uint8_t flags = (h.fw ? 1 : 0) | (refName? 2 : 0) | (h.mate << 2);
		o.write((const char *)&flags, 1);
		if(refName) {
			// Write reference name as string
			uint16_t rnamelen = (uint16_t)(*refnames)[h.h.first].length();
			o.write((const char *)&rnamelen, 2);
			o.write((*refnames)[h.h.first].c_str(), rnamelen);
		} else {
			// Write reference name as index into global reference name
			// list
			o.write((const char *)&h.h.first, 4);
		}
		// Write reference offset
		uint32_t offset = h.h.second + offBase;
		o.write((const char *)&offset, 4);
		// Write pattern sequence
		uint16_t plen = (uint16_t)length(h.patSeq);
		assert_lt(plen, 1024);
		o.write((const char *)&plen, 2);
		for(size_t i = 0; i < plen; i += 2) {
			uint8_t twoChrs = (uint8_t)h.patSeq[i];
			if(i+1 < plen) {
				twoChrs |= ((uint8_t)h.patSeq[i+1] << 4);
			}
			o.write((const char *)&twoChrs, 1);
		}
		// Write quals sequence
		uint16_t qlen = (uint16_t)length(h.quals);
		assert_eq(plen, qlen);
		o.write(begin(h.quals), qlen);
		// Write oms
		o.write((const char *)&h.oms, 4);
		// Write # mismatches
		uint8_t numMms = h.mms.count();
		o.write((const char *)&numMms, 1);
		// Output mismatches
		size_t c = 0;
		for (uint8_t i = 0; i < h.mms.size(); ++ i) {
			if (h.mms.test(i)) {
				o.write((const char *)&i, 1);
				assert_gt(h.refcs.size(), i);
				assert_eq(1, dna4Cat[(int)h.refcs[i]]);
				uint8_t refChar = charToDna5[(int)h.refcs[i]];
				assert_leq(refChar, 4);
				uint8_t qryChar = (h.fw ? (int)h.patSeq[i] :
				                          (int)h.patSeq[length(h.patSeq)-i-1]);
				assert_leq(refChar, 4);
				assert_neq(refChar, qryChar);
				uint8_t both = refChar | (qryChar << 4);
				o.write((const char *)&both, 1);
				c++;
			}
		}
	}

	/**
	 * Append a binary alignment to the output stream corresponding to
	 * the reference sequence involved.
	 */
	virtual void append(ostream& o, const Hit& h) {
		BinaryHitSink::append(o, h, this->_refnames, this->offBase_);
	}

	/**
	 * Read a binary-encoded hit (written by append() above) from an
	 * input stream.
	 */
	static bool readHit(Hit& h,
	                    istream& in,
	                    vector<string>* refnames,
	                    bool verbose = false)
	{
		if(!in.good()) return false;
		uint16_t pnamelen;
		in.read((char *)&pnamelen, 2);
		if(in.eof()) return false;
		if(verbose) cout << "name(" << pnamelen << ")=\"";
		seqan::resize(h.patName, pnamelen);
		in.read((char*)begin(h.patName), pnamelen);
		if(verbose) cout << h.patName << "\", ";
		// Parse read name
		h.patId = 0;
		for(int i = 0; i < pnamelen; i++) {
			// Skip over mate designation
			if(h.patName[i] == '/' && i == pnamelen-2) break;
			// Check whether we can continue interpreting this as a
			// number
			if(h.patName[i] < '0' || h.patName[i] > '9') {
				h.patId = 0;
				break;
			}
			h.patId *= 10;
			h.patId += (h.patName[i] - '0');
		}
		if(verbose) cout << "patid(" << h.patId << ")=\"";
		// Write fw/refname flags
		uint8_t flags;
		in.read((char *)&flags, 1);
		h.fw         = ((flags & 1) != 0);
		bool refName = ((flags & 2) != 0);
		h.mate       = (flags >> 2) & 3;
		if(verbose) cout << "fw=" << (h.fw ? "true":"false") << ", ";
		if(refName) {
			// Read and ignore reference name
			uint16_t rnamelen;
			char buf[2048];
			in.read((char *)&rnamelen, 2);
			in.read(buf, rnamelen);
			buf[rnamelen] = '\0';
			if(verbose) cout << "refname=\"" << buf << "\", ";
			string name(buf);
			// Add to the refnames vector; this isn't efficient - if
			// this ends up being a problem, we should use a map
			// instead of a vector
			if(refnames != NULL) {
				bool found = true;
				h.h.first = 0;
				for(int i = 0; i < rnamelen; i++) {
					if(buf[i] < '0' || buf[i] > '9') {
						h.h.first = 0;
						found = false;
						break;
					}
					h.h.first *= 10;
					h.h.first += (buf[i] - '0');
				}
				if(!found) {
					for(size_t i = 0; i < refnames->size(); i++) {
						if((*refnames)[i] == name) {
							h.h.first = i;
							found = true;
							break;
						}
					}
					if(!found) {
						h.h.first = refnames->size();
						refnames->push_back(name);
					}
				}
			} else {
				h.h.first = 0;
				for(int i = 0; i < rnamelen; i++) {
					if(buf[i] < '0' || buf[i] > '9') {
						h.h.first = 0;
						break;
					}
					h.h.first *= 10;
					h.h.first += (buf[i] - '0');
				}
			}
			if(verbose) cout << "refidx=\"" << h.h.first << "\", ";
		} else {
			// Read reference id as index into global reference name list
			in.read((char*)&h.h.first, 4);
			if(verbose) cout << "refidx=" << h.h.first << ", ";
		}
		// Read reference offset
		in.read((char*)&h.h.second, 4);
		if(verbose) cout << "refoff=" << h.h.second << ", ";
		// Read pattern length
		uint16_t plen;
		in.read((char*)&plen, 2);
		if(verbose) cout << "plen=" << plen << ", ";
		h.refcs.resize(plen, 0);
		assert_gt(plen, 1);
		assert_lt(plen, 1024);
		// Read pattern sequence
		seqan::resize(h.patSeq, plen);
		for(size_t i = 0; i < plen; i += 2) {
			uint8_t twoChars;
			in.read((char *)&twoChars, 1);
			assert_lt((twoChars & 15), 5);
			assert_lt((twoChars >> 4), 5);
			h.patSeq[i] = (twoChars & 15);
			if(i+1 < plen) {
				h.patSeq[i+1] = (twoChars >> 4);
			}
		}
		if(verbose) cout << "pat=\"" << h.patSeq << "\", ";
		// Read quals sequence
		uint16_t qlen = plen;
		seqan::resize(h.quals, qlen);
		in.read((char*)begin(h.quals), qlen);
		if(verbose) cout << "quals=\"" << h.quals << "\", ";
#ifndef NDEBUG
		for(size_t i = 0; i < length(h.quals); i++) {
			assert_geq(h.quals[i], 33);
			assert_leq(h.quals[i], 126);
		}
#endif
		// Read oms
		in.read((char *)&h.oms, 4);
		if(verbose) cout << "oms=" << h.oms << ", ";
		// Read # mismatches
		uint8_t numMms;
		in.read((char *)&numMms, 1);
		if(verbose) cout << "nummms=" << ((int)numMms) << ", ";
		// Read mismatches
		for(uint8_t i = 0; i < numMms; i++) {
			uint8_t ii;
			// Read the read offset (from 5' end)
			in.read((char*)&ii, 1);
			h.mms.set(ii);
			uint8_t both;
			// Read the reference and query characters involved
			in.read((char*)&both, 1);
			uint8_t refChar = both & 15;
			uint8_t qryChar = both >> 4;
			assert_neq(refChar, qryChar);
			h.refcs[ii] = "ACGTN"[refChar];
			if(verbose) cout << ((int)ii) << ":" << "ACGTN"[refChar] << ">" << "ACGTN"[qryChar] << ", ";
		}
		if(verbose) cout << endl;
		// Done
		return true;
	}

	/**
	 * Report a single hit to the appropriate output stream.
	 */
	virtual void reportHit(const Hit& h) {
		ostringstream ss;
		append(ss, h);
		lock(h.h.first);
		out(h.h.first).writeString(ss.str());
		unlock(h.h.first);
		mainlock();
		commitHit(h);
		first_ = false;
		if(h.mate > 0) numReportedPaired_++;
		else           numReported_++;
		mainunlock();
	}

private:
	int offBase_;          /// Add this to reference offsets before outputting.
                           /// (An easy way to make things 1-based instead of
                           /// 0-based)
};

/**
 * Sink that does nothing.
 */
class StubHitSink : public HitSink {
public:
	StubHitSink() : HitSink(new OutFileBuf(".tmp"), "", "", "", "", "", "", "", "", "", false, NULL) { quiet_ = true; }
	virtual void reportHit(const Hit& h) { }
	virtual void append(ostream& o, const Hit& h) { }
};

#endif /*HIT_H_*/

