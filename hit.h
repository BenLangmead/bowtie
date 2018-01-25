#ifndef HIT_H_
#define HIT_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "threading.h"
#include "bitset.h"
#include "tokenize.h"
#include "pat.h"
#include "formats.h"
#include "filebuf.h"
#include "edit.h"
#include "sstring.h"
#include <algorithm>

/**
 * Classes for dealing with reporting alignments.
 */

using namespace std;
using namespace seqan;

/// Constants for the various output modes
enum output_types {
	OUTPUT_FULL = 1,
	OUTPUT_BINARY,
	OUTPUT_CHAIN,
	OUTPUT_SAM,
	OUTPUT_NONE
};

/// Names of the various output modes
static const std::string output_type_names[] = {
	"Invalid!",
	"Full",
	"Binary",
	"None"
};

typedef pair<TIndexOffU,TIndexOffU> UPair;

/**
 * Encapsulates a hit, including a text-id/text-offset pair, a pattern
 * id, and a boolean indicating whether it matched as its forward or
 * reverse-complement version.
 */
class Hit {
public:
	Hit() : stratum(-1) { }

	UPair             h;       /// reference index & offset
	UPair             mh;      /// reference index & offset for mate
	uint32_t            patId;   /// read index
	String<char>        patName; /// read name
	String<Dna5>        patSeq;  /// read sequence
	String<Dna5>        colSeq;  /// original color sequence, not decoded
	String<char>        quals;   /// read qualities
	String<char>        colQuals;/// original color qualities, not decoded
	FixedBitset<1024>   mms;     /// nucleotide mismatch mask
	FixedBitset<1024>   cmms;    /// color mismatch mask (if relevant)
	vector<char>        refcs;   /// reference characters for mms
	vector<char>        crefcs;  /// reference characters for cmms
	uint32_t            oms;     /// # of other possible mappings; 0 -> this is unique
	bool                fw;      /// orientation of read in alignment
	bool                mfw;     /// orientation of mate in alignment
	uint16_t            mlen;    /// length of mate
	int8_t              stratum; /// stratum of hit (= mismatches in seed)
	uint32_t            cost;    /// total cost, factoring in stratum and quality penalty
	uint8_t             mate;    /// matedness; 0 = not a mate
	                             ///            1 = upstream mate
	                             ///            2 = downstream mate
	bool                color;   /// read is in colorspace?
	char                primer;  /// primer base, for csfasta files
	char                trimc;   /// trimmed color, for csfasta files
	uint32_t            seed;    /// pseudo-random seed for aligned read

	/**
	 * Return true if this Hit is internally consistent.  Otherwise,
	 * throw an assertion.
	 */
	bool repOk() const {
		assert_geq(cost, (uint32_t)(stratum << 14));
		return true;
	}

	size_t length() const { return seqan::length(patSeq); }

	Hit& operator = (const Hit &other) {
		this->h       = other.h;
		this->mh      = other.mh;
		this->patId   = other.patId;
		this->patName = other.patName;
		this->patSeq  = other.patSeq;
		this->colSeq  = other.colSeq;
		this->quals   = other.quals;
		this->colQuals= other.colQuals;
		this->mms     = other.mms;
		this->cmms    = other.cmms;
		this->refcs   = other.refcs;
		this->crefcs  = other.crefcs;
		this->oms     = other.oms;
		this->fw      = other.fw;
		this->mfw     = other.mfw;
		this->mlen    = other.mlen;
		this->stratum = other.stratum;
		this->cost    = other.cost;
		this->mate    = other.mate;
		this->color   = other.color;
		this->cmms    = other.cmms;
		this->seed    = other.seed;
		return *this;
	}
};

/**
 * Compare hits a and b; a < b if its cost is less than B's.  If
 * there's a tie, break on position and orientation.
 */
class HitCostCompare {
public:
	bool operator() (const Hit& a, const Hit& b) {
		if(a.cost < b.cost) return true;
		if(a.cost > b.cost) return false;
		if(a.h < b.h) return true;
		if(a.h > b.h) return false;
		if(a.fw < b.fw) return true;
		if(a.fw > b.fw) return false;
		return false;
	}
};

/// Sort by text-id then by text-offset
bool operator< (const Hit& a, const Hit& b);

/**
 * Encapsulates an object that accepts hits, optionally retains them in
 * a vector, and does something else with them according to
 * descendent's implementation of pure virtual member reportHitImpl().
 */
class HitSink {
public:
	explicit HitSink(
		OutFileBuf& out,
		const std::string& dumpAl,
		const std::string& dumpUnal,
		const std::string& dumpMax,
		bool onePairFile,
		bool sampleMax,
		vector<string>* refnames,
		size_t nthreads,
		size_t perThreadBufSize,
		bool reorder) :
		out_(out),
		_refnames(refnames),
		mutex_(),
		dumpAlBase_(dumpAl),
		dumpUnalBase_(dumpUnal),
		dumpMaxBase_(dumpMax),
		onePairFile_(onePairFile),
		sampleMax_(sampleMax),
		quiet_(false),
		nthreads_((nthreads > 0) ? nthreads : 1),
		ptBufs_(),
		ptCounts_(nthreads_),
		perThreadBufSize_(perThreadBufSize),
		ptNumAligned_(NULL),
		reorder_(reorder)
	{
		size_t nelt = 5 * nthreads_;
		ptNumAligned_ = new uint64_t[nelt];
		std::memset(reinterpret_cast<void*>(const_cast<uint64_t*>(ptNumAligned_)), 0, sizeof(uint64_t) * nelt);
		ptNumReported_ = ptNumAligned_ + nthreads_;
		ptNumReportedPaired_ = ptNumReported_ + nthreads_;
		ptNumUnaligned_ = ptNumReportedPaired_ + nthreads_;
		ptNumMaxed_ = ptNumUnaligned_ + nthreads_;
		ptBufs_.resize(nthreads_);
		ptCounts_.resize(nthreads_, 0);
		initDumps();
	}

	/**
	 * Destroy HitSinkobject;
	 */
	virtual ~HitSink() {
		if(ptNumAligned_ != NULL) {
			delete[] ptNumAligned_;
			ptNumAligned_ = NULL;
		}
		closeOuts();
		destroyDumps();
	}

	/**
	 * Append a single hit to the given output stream.
	 */
	virtual void append(BTString& o, const Hit& h, int mapq, int xms) = 0;

	/**
	 * Add a number of alignments to the tally.  Tally shouldn't
	 * include reads that fail to align either because they had 0
	 * alignments or because of -m.
	 */
	void tallyAlignments(size_t threadId, size_t numAl, bool paired) {
        assert(!paired || (numAl % 2) == 0);
        ptNumAligned_[threadId] ++;
		if(paired) {
			ptNumReportedPaired_[threadId] += numAl;
		} else {
			ptNumReported_[threadId] += numAl;
		}
	}

	/**
	 * Report a batch of hits from a vector, perhaps subsetting it.
	 */
	virtual void reportHits(
		const Hit *hptr,
		vector<Hit> *hsptr,
		size_t start,
		size_t end,
		size_t threadId,
		int mapq,
		int xms,
		bool tally, size_t rdid)
	{
		assert_geq(end, start);
		assert(nthreads_ > 1 || threadId == 0);
		if(end == start) {
			return;
		}
		const Hit& firstHit = (hptr == NULL) ? (*hsptr)[start] : *hptr;
		bool paired = firstHit.mate > 0;
		maybeFlush(threadId);
		BTString& o = ptBufs_[threadId];
		// Per-thread buffering is active
		for(size_t i = start; i < end; i++) {
			const Hit& h = (hptr == NULL) ? (*hsptr)[i] : *hptr;
			assert(h.repOk());
			append(o, h, mapq, xms);
			if(nthreads_ == 1) {
				out_.writeString(o);
				o.clear();
			}
		}
		ptCounts_[threadId]++;
		if(tally) {
			tallyAlignments(threadId, end - start, paired);
		}
	}

	/**
	 * Called when all alignments are complete.  It is assumed that no
	 * synchronization is necessary
	 */
	void finish(bool hadoopOut) {
		// Flush all per-thread buffers
		flushAll();
		
		// Close all output streams
		closeOuts();
		
		// Print information about how many unpaired and/or paired
		// reads were aligned.
		if(!quiet_) {
			uint64_t numReported = 0, numReportedPaired = 0;
			uint64_t numAligned = 0, numUnaligned = 0;
			uint64_t numMaxed = 0;
			for(size_t i = 0; i < nthreads_; i++) {
				numReported += ptNumReported_[i];
				numReportedPaired += ptNumReportedPaired_[i];
				numAligned += ptNumAligned_[i];
				numUnaligned += ptNumUnaligned_[i];
				numMaxed += ptNumMaxed_[i];
			}
			
			uint64_t tot = numAligned + numUnaligned + numMaxed;
			double alPct = 0.0, unalPct = 0.0, maxPct = 0.0;
			if(tot > 0) {
				alPct   = 100.0 * (double)numAligned / (double)tot;
				unalPct = 100.0 * (double)numUnaligned / (double)tot;
				maxPct  = 100.0 * (double)numMaxed / (double)tot;
			}
			cerr << "# reads processed: " << tot << endl;
			cerr << "# reads with at least one reported alignment: "
			     << numAligned << " (" << fixed << setprecision(2)
			     << alPct << "%)" << endl;
			cerr << "# reads that failed to align: "
			     << numUnaligned << " (" << fixed << setprecision(2)
			     << unalPct << "%)" << endl;
			if(numMaxed > 0) {
				if(sampleMax_) {
					cerr << "# reads with alignments sampled due to -M: "
						 << numMaxed << " (" << fixed << setprecision(2)
						 << maxPct << "%)" << endl;
				} else {
					cerr << "# reads with alignments suppressed due to -m: "
						 << numMaxed << " (" << fixed << setprecision(2)
						 << maxPct << "%)" << endl;
				}
			}
			if(numReported == 0 && numReportedPaired == 0) {
				cerr << "No alignments" << endl;
			}
			else if(numReportedPaired > 0 && numReported == 0) {
				cerr << "Reported " << (numReportedPaired >> 1)
					 << " paired-end alignments" << endl;
			}
			else if(numReported > 0 && numReportedPaired == 0) {
				cerr << "Reported " << numReported
					 << " alignments" << endl;
			}
			else {
				assert_gt(numReported + numReportedPaired, 0);
				cerr << "Reported " << (numReportedPaired >> 1)
					 << " paired-end alignments and " << numReported
					 << " singleton alignments" << endl;
			}
			if(hadoopOut) {
				cerr << "reporter:counter:Bowtie,Reads with reported alignments," << numAligned << endl;
				cerr << "reporter:counter:Bowtie,Reads with no alignments," << numUnaligned << endl;
				cerr << "reporter:counter:Bowtie,Reads exceeding -m limit," << numMaxed << endl;
				cerr << "reporter:counter:Bowtie,Unpaired alignments reported," << numReported << endl;
				cerr << "reporter:counter:Bowtie,Paired alignments reported," << numReportedPaired << endl;
			}
		}
	}

	/**
	 * Returns alignment output stream.
	 */
	OutFileBuf& out() { return out_; }

	/**
	 * Return true iff this HitSink dumps aligned reads to an output
	 * stream (i.e., iff --alfa or --alfq are specified).
	 */
	bool dumpsAlignedReads() const { return dumpAlignFlag_; }

	/**
	 * Return true iff this HitSink dumps unaligned reads to an output
	 * stream (i.e., iff --unfa or --unfq are specified).
	 */
	bool dumpsUnalignedReads() const { return dumpUnalignFlag_; }

	/**
	 * Return true iff this HitSink dumps maxed-out reads to an output
	 * stream (i.e., iff --maxfa or --maxfq are specified).
	 */
	bool dumpsMaxedReads() const { return dumpMaxedFlag_ || dumpUnalignFlag_; }

	/**
	 * Return true iff this HitSink dumps either unaligned or maxed-
	 * out reads to an output stream (i.e., iff --unfa, --maxfa,
	 * --unfq, or --maxfq are specified).
	 */
	bool dumpsReads() const {
		return dumpAlignFlag_ || dumpUnalignFlag_ || dumpMaxedFlag_;
	}

	/**
	 * Dump an aligned read to all of the appropriate output streams.
	 * Be careful to synchronize correctly - there may be multiple
	 * simultaneous writers.
	 */
	void dumpAlign(PatternSourcePerThread& p) {
		if(!dumpAlignFlag_) return;
		const bool paired = p.bufa().mate > 0;
		if(!paired || onePairFile_) {
			// Dump unpaired read to an aligned-read file of the same format
			if(!dumpAlBase_.empty()) {
				ThreadSafe _ts(&dumpAlignLock_);
				if(dumpAl_ == NULL) {
					assert(dumpAlQv_ == NULL);
					dumpAl_ = openOf(dumpAlBase_, 0, "");
					assert(dumpAl_ != NULL);
					if(p.bufa().qualOrigBufLen > 0) {
						dumpAlQv_ = openOf(dumpAlBase_ + ".qual", 0, "");
						assert(dumpAlQv_ != NULL);
					}
				}
				dumpAl_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				if(dumpAlQv_ != NULL) {
					dumpAlQv_->write(p.bufa().qualOrigBuf, p.bufa().qualOrigBufLen);
				}
			}
		} else {
			// Dump paired-end read to an aligned-read file (or pair of
			// files) of the same format
			if(!dumpAlBase_.empty()) {
				ThreadSafe _ts(&dumpAlignLockPE_);
				if(dumpAl_1_ == NULL) {
					assert(dumpAlQv_1_ == NULL);
					assert(dumpAlQv_2_ == NULL);
					dumpAl_1_ = openOf(dumpAlBase_, 1, "");
					dumpAl_2_ = openOf(dumpAlBase_, 2, "");
					assert(dumpAl_1_ != NULL);
					assert(dumpAl_2_ != NULL);
					if(p.bufa().qualOrigBufLen > 0) {
						dumpAlQv_1_ = openOf(dumpAlBase_ + ".qual", 1, "");
						dumpAlQv_2_ = openOf(dumpAlBase_ + ".qual", 2, "");
						assert(dumpAlQv_1_ != NULL);
						assert(dumpAlQv_2_ != NULL);
					}
				}
				dumpAl_1_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				dumpAl_2_->write(p.bufb().readOrigBuf, p.bufb().readOrigBufLen);
				if(dumpAlQv_1_ != NULL) {
					dumpAlQv_1_->write(p.bufa().qualOrigBuf, p.bufa().qualOrigBufLen);
					dumpAlQv_2_->write(p.bufb().qualOrigBuf, p.bufb().qualOrigBufLen);
				}
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
		const bool paired = p.bufa().mate > 0;
		if(!paired || onePairFile_) {
			// Dump unpaired read to an unaligned-read file of the same format
			if(!dumpUnalBase_.empty()) {
				ThreadSafe _ts(&dumpUnalLock_);
				if(dumpUnal_ == NULL) {
					assert(dumpUnalQv_ == NULL);
					dumpUnal_ = openOf(dumpUnalBase_, 0, "");
					assert(dumpUnal_ != NULL);
					if(p.bufa().qualOrigBufLen > 0) {
						dumpUnalQv_ = openOf(dumpUnalBase_ + ".qual", 0, "");
						assert(dumpUnalQv_ != NULL);
					}
				}
				dumpUnal_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				if(dumpUnalQv_ != NULL) {
					dumpUnalQv_->write(p.bufa().qualOrigBuf, p.bufa().qualOrigBufLen);
				}
			}
		} else {
			// Dump paired-end read to an unaligned-read file (or pair
			// of files) of the same format
			if(!dumpUnalBase_.empty()) {
				ThreadSafe _ts(&dumpUnalLockPE_);
				if(dumpUnal_1_ == NULL) {
					assert(dumpUnal_1_ == NULL);
					assert(dumpUnal_2_ == NULL);
					dumpUnal_1_ = openOf(dumpUnalBase_, 1, "");
					dumpUnal_2_ = openOf(dumpUnalBase_, 2, "");
					assert(dumpUnal_1_ != NULL);
					assert(dumpUnal_2_ != NULL);
					if(p.bufa().qualOrigBufLen > 0) {
						dumpUnalQv_1_ = openOf(dumpUnalBase_ + ".qual", 1, "");
						dumpUnalQv_2_ = openOf(dumpUnalBase_ + ".qual", 2, "");
					}
				}
				dumpUnal_1_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				dumpUnal_2_->write(p.bufb().readOrigBuf, p.bufb().readOrigBufLen);
				if(dumpUnalQv_1_ != NULL) {
					dumpUnalQv_1_->write(p.bufa().qualOrigBuf, p.bufa().qualOrigBufLen);
					dumpUnalQv_2_->write(p.bufb().qualOrigBuf, p.bufb().qualOrigBufLen);
				}
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
		const bool paired = p.bufa().mate > 0;
		if(!paired || onePairFile_) {
			// Dump unpaired read to an maxed-out-read file of the same format
			if(!dumpMaxBase_.empty()) {
				ThreadSafe _ts(&dumpMaxLock_);
				if(dumpMax_ == NULL) {
					dumpMax_ = openOf(dumpMaxBase_, 0, "");
					assert(dumpMax_ != NULL);
					if(p.bufa().qualOrigBufLen > 0) {
						dumpMaxQv_ = openOf(dumpMaxBase_ + ".qual", 0, "");
					}
				}
				dumpMax_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				if(dumpMaxQv_ != NULL) {
					dumpMaxQv_->write(p.bufa().qualOrigBuf, p.bufa().qualOrigBufLen);
				}
			}
		} else {
			// Dump paired-end read to a maxed-out-read file (or pair
			// of files) of the same format
			if(!dumpMaxBase_.empty()) {
				ThreadSafe _ts(&dumpMaxLockPE_);
				if(dumpMax_1_ == NULL) {
					assert(dumpMaxQv_1_ == NULL);
					assert(dumpMaxQv_2_ == NULL);
					dumpMax_1_ = openOf(dumpMaxBase_, 1, "");
					dumpMax_2_ = openOf(dumpMaxBase_, 2, "");
					assert(dumpMax_1_ != NULL);
					assert(dumpMax_2_ != NULL);
					if(p.bufa().qualOrigBufLen > 0) {
						dumpMaxQv_1_ = openOf(dumpMaxBase_ + ".qual", 1, "");
						dumpMaxQv_2_ = openOf(dumpMaxBase_ + ".qual", 2, "");
					}
				}
				dumpMax_1_->write(p.bufa().readOrigBuf, p.bufa().readOrigBufLen);
				dumpMax_2_->write(p.bufb().readOrigBuf, p.bufb().readOrigBufLen);
				if(dumpMaxQv_1_ != NULL) {
					dumpMaxQv_1_->write(p.bufa().qualOrigBuf, p.bufa().qualOrigBufLen);
					dumpMaxQv_2_->write(p.bufb().qualOrigBuf, p.bufb().qualOrigBufLen);
				}
			}
		}
	}

	/**
	 * Report a maxed-out read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportMaxed(
		vector<Hit>& hs,
		size_t threadId,
		PatternSourcePerThread& p)
	{
		ptNumMaxed_[threadId]++;
	}

	/**
	 * Report an unaligned read.  Typically we do nothing, but we might
	 * want to print a placeholder when output is chained.
	 */
	virtual void reportUnaligned(
		size_t threadId,
		PatternSourcePerThread& p)
	{
		ptNumUnaligned_[threadId]++;
	}

protected:

	/**
	 * Flush thread's output buffer and reset both buffer and count.
	 */
	void flush(size_t threadId, bool finalBatch) {
		{
			ThreadSafe _ts(&mutex_); // flush
			out_.writeString(ptBufs_[threadId]);
		}
		ptCounts_[threadId] = 0;
		ptBufs_[threadId].clear();
	}
	
	/**
	 * Flush all output buffers.
	 */
	void flushAll() {
		for(size_t i = 0; i < nthreads_; i++) {
			flush(i, i == nthreads_ - 1);
		}
	}

	/**
	 * If the thread's output buffer is currently full, flush it and
	 * reset both buffer and count.
	 */
	void maybeFlush(size_t threadId) {
		if(ptCounts_[threadId] >= perThreadBufSize_) {
			flush(threadId, false /* final batch? */);
		}
	}
	
	/**
	 * Close (and flush) all OutFileBufs.
	 */
	void closeOuts() {
		out_.close();
	}

	OutFileBuf&         out_;        /// the alignment output stream(s)
	vector<string>*     _refnames;    /// map from reference indexes to names
	MUTEX_T             mutex_;       /// pthreads mutexes for per-file critical sections
	
	// used for output read buffer	
	size_t nthreads_;
	std::vector<BTString> ptBufs_;
	std::vector<size_t> ptCounts_;
	int perThreadBufSize_;
    bool reorder_;

	// Output filenames for dumping
	std::string dumpAlBase_;
	std::string dumpUnalBase_;
	std::string dumpMaxBase_;

	bool onePairFile_;
	bool sampleMax_;

	// Output streams for dumping sequences
	std::ofstream *dumpAl_;       // for single-ended reads
	std::ofstream *dumpAl_1_;     // for first mates
	std::ofstream *dumpAl_2_;     // for second mates
	std::ofstream *dumpUnal_;     // for single-ended reads
	std::ofstream *dumpUnal_1_;   // for first mates
	std::ofstream *dumpUnal_2_;   // for second mates
	std::ofstream *dumpMax_;      // for single-ended reads
	std::ofstream *dumpMax_1_;    // for first mates
	std::ofstream *dumpMax_2_;    // for second mates

	// Output streams for dumping qualities
	std::ofstream *dumpAlQv_;     // for single-ended reads
	std::ofstream *dumpAlQv_1_;   // for first mates
	std::ofstream *dumpAlQv_2_;   // for second mates
	std::ofstream *dumpUnalQv_;   // for single-ended reads
	std::ofstream *dumpUnalQv_1_; // for first mates
	std::ofstream *dumpUnalQv_2_; // for second mates
	std::ofstream *dumpMaxQv_;    // for single-ended reads
	std::ofstream *dumpMaxQv_1_;  // for first mates
	std::ofstream *dumpMaxQv_2_;  // for second mates

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
			cerr << "Bad mate type " << mateType << endl; throw 1;
		}
		std::ofstream* tmp = new ofstream(s.c_str(), ios::out);
		if(tmp->fail()) {
			if(mateType == 0) {
				cerr << "Could not open single-ended aligned/unaligned-read file for writing: " << name << endl;
			} else {
				cerr << "Could not open paired-end aligned/unaligned-read file for writing: " << name << endl;
			}
			throw 1;
		}
		return tmp;
	}

	/**
	 * Initialize all the locks for dumping.
	 */
	void initDumps() {
		dumpAl_       = dumpAl_1_     = dumpAl_2_     = NULL;
		dumpUnal_     = dumpUnal_1_   = dumpUnal_2_   = NULL;
		dumpMax_      = dumpMax_1_    = dumpMax_2_    = NULL;
		dumpAlQv_     = dumpAlQv_1_   = dumpAlQv_2_   = NULL;
		dumpUnalQv_   = dumpUnalQv_1_ = dumpUnalQv_2_ = NULL;
		dumpMaxQv_    = dumpMaxQv_1_  = dumpMaxQv_2_  = NULL;
		dumpAlignFlag_   = !dumpAlBase_.empty();
		dumpUnalignFlag_ = !dumpUnalBase_.empty();
		dumpMaxedFlag_   = !dumpMaxBase_.empty();
	}

	void destroyDumps() {
		if(dumpAl_       != NULL) { dumpAl_->close();       delete dumpAl_; }
		if(dumpAl_1_     != NULL) { dumpAl_1_->close();     delete dumpAl_1_; }
		if(dumpAl_2_     != NULL) { dumpAl_2_->close();     delete dumpAl_2_; }
		if(dumpUnal_     != NULL) { dumpUnal_->close();     delete dumpUnal_; }
		if(dumpUnal_1_   != NULL) { dumpUnal_1_->close();   delete dumpUnal_1_; }
		if(dumpUnal_2_   != NULL) { dumpUnal_2_->close();   delete dumpUnal_2_; }
		if(dumpMax_      != NULL) { dumpMax_->close();      delete dumpMax_; }
		if(dumpMax_1_    != NULL) { dumpMax_1_->close();    delete dumpMax_1_; }
		if(dumpMax_2_    != NULL) { dumpMax_2_->close();    delete dumpMax_2_; }
		if(dumpAlQv_     != NULL) { dumpAlQv_->close();     delete dumpAlQv_; }
		if(dumpAlQv_1_   != NULL) { dumpAlQv_1_->close();   delete dumpAlQv_1_; }
		if(dumpAlQv_2_   != NULL) { dumpAlQv_2_->close();   delete dumpAlQv_2_; }
		if(dumpUnalQv_   != NULL) { dumpUnalQv_->close();   delete dumpUnalQv_; }
		if(dumpUnalQv_1_ != NULL) { dumpUnalQv_1_->close(); delete dumpUnalQv_1_; }
		if(dumpUnalQv_2_ != NULL) { dumpUnalQv_2_->close(); delete dumpUnalQv_2_; }
		if(dumpMaxQv_    != NULL) { dumpMaxQv_->close();    delete dumpMaxQv_; }
		if(dumpMaxQv_1_  != NULL) { dumpMaxQv_1_->close();  delete dumpMaxQv_1_; }
		if(dumpMaxQv_2_  != NULL) { dumpMaxQv_2_->close();  delete dumpMaxQv_2_; }
	}

	// Locks for dumping
	MUTEX_T dumpAlignLock_;
	MUTEX_T dumpAlignLockPE_; // _1 and _2
	MUTEX_T dumpUnalLock_;
	MUTEX_T dumpUnalLockPE_; // _1 and _2
	MUTEX_T dumpMaxLock_;
	MUTEX_T dumpMaxLockPE_;   // _1 and _2

	// false -> no dumping
	bool dumpAlignFlag_;
	bool dumpUnalignFlag_;
	bool dumpMaxedFlag_;

	volatile bool first_;       /// true -> first hit hasn't yet been reported
	volatile uint64_t *ptNumAligned_;
	volatile uint64_t *ptNumReported_;
	volatile uint64_t *ptNumReportedPaired_;
	volatile uint64_t *ptNumUnaligned_;
	volatile uint64_t *ptNumMaxed_;

	bool quiet_;  /// true -> don't print alignment stats at the end
};

/**
 * A per-thread wrapper for a HitSink.  Incorporates state that a
 * single search thread cares about.
 */
class HitSinkPerThread {
public:
	explicit HitSinkPerThread(
		HitSink& sink,
		uint32_t max,
		uint32_t n,
		int defaultMapq,
		size_t threadId) :
		_sink(sink),
		_bestRemainingStratum(0),
		_numValidHits(0llu),
		_hits(),
		_bufferedHits(),
		hitsForThisRead_(),
		_max(max),
		_n(n),
		defaultMapq_(defaultMapq),
		threadId_(threadId)
	{
		assert_gt(_n, 0);
	}

	virtual ~HitSinkPerThread() { }

	/// Return the vector of retained hits
	vector<Hit>& retainedHits()   { return _hits; }

	/// Finalize current read
	virtual uint32_t finishRead(PatternSourcePerThread& p, bool report, bool dump) {
		uint32_t ret = finishReadImpl();
		_bestRemainingStratum = 0;
		if(!report) {
			_bufferedHits.clear();
			return 0;
		}
		bool maxed = (ret > _max);
		bool unal = (ret == 0);
		if(dump && (unal || maxed)) {
			// Either no reportable hits were found or the number of
			// reportable hits exceeded the -m limit specified by the
			// user
			assert(ret == 0 || ret > _max);
			if(maxed) _sink.dumpMaxed(p);
			else      _sink.dumpUnal(p);
		}
		ret = 0;
		if(maxed) {
			// Report that the read maxed-out; useful for chaining output
			if(dump) _sink.reportMaxed(_bufferedHits, threadId_, p);
			_bufferedHits.clear();
		} else if(unal) {
			// Report that the read failed to align; useful for chaining output
			if(dump) _sink.reportUnaligned(threadId_, p);
		} else {
			// Flush buffered hits
			assert_gt(_bufferedHits.size(), 0);
			if(_bufferedHits.size() > _n) {
				_bufferedHits.resize(_n);
			}
			int mapq = defaultMapq_;
			int xms = (int)(_bufferedHits.size());
			const bool paired = p.bufa().mate > 0;
			if(paired) {
				xms /= 2;
			}
			xms++;
			_sink.reportHits(NULL, &_bufferedHits, 0, _bufferedHits.size(),
			                 threadId_, mapq, xms, true, p.rdid());
			_sink.dumpAlign(p);
			ret = (uint32_t)_bufferedHits.size();
			_bufferedHits.clear();
		}
		assert_eq(0, _bufferedHits.size());
		return ret;
	}

	virtual uint32_t finishReadImpl() = 0;

	/**
	 * Add a hit to the internal buffer.  Not yet reporting the hits.
	 */
	virtual void bufferHit(const Hit& h, int stratum) {
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
		assert(h.repOk());
		_numValidHits++;
		return true;
	}

	/// Return the number of valid hits so far.
	uint64_t numValidHits()    { return _numValidHits; }

	/**
	 * Return true if there are no more reportable hits.
	 */
	bool finishedWithStratum(int stratum) {
		bool ret = finishedWithStratumImpl(stratum);
		_bestRemainingStratum = stratum+1;
		return ret;
	}

	/**
	 * Use the given set of hits as a starting point.  By default, we don't
	 */
	virtual bool setHits(HitSet& hs) {
		if(!hs.empty()) {
			cerr << "Error: default setHits() called with non-empty HitSet" << endl;
			throw 1;
		}
		return false;
	}

	/**
	 * Return true if there are no reportable hits with the given cost
	 * (or worse).
	 */
	virtual bool irrelevantCost(uint16_t cost) const {
		return false;
	}

	/**
	 * Concrete subclasses override this to determine whether the
	 * search routine should keep searching after having finished
	 * reporting all alignments at the given stratum.
	 */
	virtual bool finishedWithStratumImpl(int stratum) = 0;

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
	bool dumpsReads() const {
		return _sink.dumpsReads();
	}

	/**
	 * Return true iff there are currently no buffered hits.
	 */
	bool empty() const {
		return _bufferedHits.empty();
	}

	/**
	 * Return the number of currently buffered hits.
	 */
	bool size() const {
		return _bufferedHits.size();
	}

	/**
	 * Return max # hits to report (*2 in paired-end mode because mates
	 * count separately)
	 */
	virtual uint32_t maxHits() const {
		return _n;
	}

protected:
	HitSink&    _sink; /// Ultimate destination of reported hits
	/// Least # mismatches in alignments that will be reported in the
	/// future.  Updated by the search routine.
	int         _bestRemainingStratum;
	/// # hits reported to this HitSink so far (not all of which were
	/// necesssary reported to _sink)
	uint64_t    _numValidHits;
	vector<Hit> _hits; /// Repository for retained hits
	/// Buffered hits, to be reported and flushed at end of read-phase
	vector<Hit> _bufferedHits;

	// Following variables are declared in the parent but maintained in
	// the concrete subcalsses
	uint32_t hitsForThisRead_; /// # hits for this read so far
	uint32_t _max; /// don't report any hits if there were > _max
	uint32_t _n;   /// report at most _n hits
	int defaultMapq_;
	size_t threadId_;
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
		uint32_t n,
		uint32_t max,
		int defaultMapq,
		size_t threadId) :
		HitSinkPerThread(sink, max, n, defaultMapq, threadId)
	{ }

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
			return true; // done - report nothing
		}
		//if(hitsForThisRead_ <= _n) {
			// Only report hit if we haven't
			bufferHit(h, stratum);
		//}
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
};

/**
 * Concrete factory for FirstNGoodHitSinkPerThreads.
 */
class NGoodHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NGoodHitSinkPerThreadFactory(
		HitSink& sink,
		uint32_t n,
		uint32_t max,
		int defaultMapq,
		size_t threadId) :
		sink_(sink),
		n_(n),
		max_(max),
		defaultMapq_(defaultMapq),
		threadId_(threadId) { }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NGoodHitSinkPerThread(sink_, n_, max_, defaultMapq_, threadId_);
	}
	
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		uint32_t max = max_ * (max_ == 0xffffffff ? 1 : m);
		uint32_t n = n_ * (n_ == 0xffffffff ? 1 : m);
		return new NGoodHitSinkPerThread(sink_, n, max, defaultMapq_, threadId_);
	}

private:

	HitSink& sink_;
	uint32_t n_;
	uint32_t max_;
	int defaultMapq_;
	size_t threadId_;
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
		uint32_t mult,
		int defaultMapq,
		size_t threadId) :
		HitSinkPerThread(sink, max, n, defaultMapq, threadId),
		bestStratum_(999),
		mult_(mult) { }

	/**
	 * false -> we do not allow strata to be spanned
	 */
	virtual bool spanStrata() {
		return false; // we do not span strata
	}

	/**
	 * true -> we report best hits
	 */
	virtual bool best() {
		return true;
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
			return true; // done - report nothing
		}
		//if(hitsForThisRead_ <= _n) {
			bufferHit(h, stratum);
		//}
		if(hitsForThisRead_ == _n &&
		   (_max == 0xffffffff || _max < _n))
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
			_bufferedHits[i].oms = ((uint32_t)sz / mult_) - 1;
		}
		return ret;
	}

	/**
	 * If we had any alignments at all and we're now moving on to a new
	 * stratum, then we're done.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
		return hitsForThisRead_ > 0;
	}

	/**
	 * If there have been any hits reported so far, classify any
	 * subsequent alignments with higher strata as irrelevant.
	 */
	virtual bool irrelevantCost(uint16_t cost) const {
		if(hitsForThisRead_) {
			// irrelevant iff at worse stratum
			return ((int)cost >> 14) > bestStratum_;
		}
		return false;
	}

private:

	int bestStratum_; /// best stratum observed so far
	uint32_t mult_; /// number of batched-up alignments
};

/**
 * Concrete factory for NBestStratHitSinkPerThread.
 */
class NBestFirstStratHitSinkPerThreadFactory : public HitSinkPerThreadFactory {
public:
	NBestFirstStratHitSinkPerThreadFactory(
		HitSink& sink,
		uint32_t n,
		uint32_t max,
		int defaultMapq,
		size_t threadId) :
		sink_(sink),
		n_(n),
		max_(max),
		defaultMapq_(defaultMapq),
		threadId_(threadId) { }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new NBestFirstStratHitSinkPerThread(sink_, n_, max_, 1, defaultMapq_, threadId_);
	}
	
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		uint32_t max = max_ * (max_ == 0xffffffff ? 1 : m);
		uint32_t n = n_ * (n_ == 0xffffffff ? 1 : m);
		return new NBestFirstStratHitSinkPerThread(sink_, n, max, m, defaultMapq_, threadId_);
	}

private:
	HitSink& sink_;
	uint32_t n_;
	uint32_t max_;
	int defaultMapq_;
	size_t threadId_;
};

/**
 * Report all valid alignments.
 */
class AllHitSinkPerThread : public HitSinkPerThread {

public:
	AllHitSinkPerThread(
		HitSink& sink,
		uint32_t max,
		int defaultMapq,
		size_t threadId) :
		HitSinkPerThread(sink, max, 0xffffffff, defaultMapq, threadId) { }

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
		uint32_t max,
		int defaultMapq,
		size_t threadId) :
		sink_(sink),
		max_(max),
		defaultMapq_(defaultMapq),
		threadId_(threadId) { }

	/**
	 * Allocate a new NGoodHitSinkPerThread object on the heap,
	 * using the parameters given in the constructor.
	 */
	virtual HitSinkPerThread* create() const {
		return new AllHitSinkPerThread(sink_, max_, defaultMapq_, threadId_);
	}
	virtual HitSinkPerThread* createMult(uint32_t m) const {
		uint32_t max = max_ * (max_ == 0xffffffff ? 1 : m);
		return new AllHitSinkPerThread(sink_, max, defaultMapq_, threadId_);
	}

private:
	HitSink& sink_;
	uint32_t max_;
	int defaultMapq_;
	size_t threadId_;
};

/**
 * Print the given string.  If ws = true, print only up to and not
 * including the first space or tab.  Useful for printing reference
 * names.
 */
inline void printUptoWs(BTString& o, const std::string& str, bool ws) {
	if(!ws) {
		o << str;
	} else {
		size_t pos = str.find_first_of(" \t");
		if(pos != string::npos) {
			o << str.substr(0, pos);
		} else {
			o << str;
		}
	}
}

/**
 * Sink that prints lines like this:
 * pat-name \t [-|+] \t ref-name \t ref-off \t pat \t qual \t #-alt-hits \t mm-list
 */
class VerboseHitSink : public HitSink {

public:

	/**
	 * Construct a single-stream VerboseHitSink (default)
	 */
	VerboseHitSink(
		OutFileBuf& out,
		int offBase,
		bool colorSeq,
		bool colorQual,
		bool printCost,
		const Bitset& suppressOuts,
		bool fullRef,
		const std::string& dumpAl,
		const std::string& dumpUnal,
		const std::string& dumpMax,
		bool onePairFile,
		bool sampleMax,
		std::vector<std::string>* refnames,
		size_t nthreads,
		size_t perThreadBufSize,
		int partition = 0) :
		HitSink(
			out,
			dumpAl,
			dumpUnal,
			dumpMax,
			onePairFile,
			sampleMax,
			refnames,
			nthreads,
			perThreadBufSize,
			false),
		partition_(partition),
		offBase_(offBase),
		colorSeq_(colorSeq),
		colorQual_(colorQual),
		cost_(printCost),
		suppress_(suppressOuts),
		fullRef_(fullRef)
		{ }

	static void append(
		BTString& o,
		const Hit& h,
		const vector<string>* refnames,
		bool fullRef,
		int partition,
		int offBase,
		bool colorSeq,
		bool colorQual,
		bool cost,
		const Bitset& suppress);

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(BTString& o, const Hit& h, int mapq, int xms) {
		VerboseHitSink::append(o, h, _refnames,
		                       fullRef_, partition_, offBase_,
		                       colorSeq_, colorQual_, cost_,
		                       suppress_);
	}

	/**
	 * See hit.cpp
	 */
	virtual void reportMaxed(
		vector<Hit>& hs,
		size_t threadId,
		PatternSourcePerThread& p);

private:
	int      partition_;   /// partition size, or 0 if partitioning is disabled
	int      offBase_;     /// Add this to reference offsets before outputting.
	                       /// (An easy way to make things 1-based instead of
	                       /// 0-based)
	bool     colorSeq_;    /// true -> print colorspace alignment sequence in colors
	bool     colorQual_;   /// true -> print colorspace quals as originals, not decoded
	bool     cost_;        /// true -> print statum and cost
	Bitset   suppress_;    /// output fields to suppress
	bool fullRef_;         /// print full reference name
};

#endif /*HIT_H_*/
