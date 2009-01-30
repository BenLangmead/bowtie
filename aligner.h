/**
 * aligner.h
 *
 * A generic class providing a stateful way to find alignments.
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <iostream>
#include <stdint.h>
#include "seqan/sequence.h"
#include "assert_helpers.h"
#include "ebwt.h"
#include "pat.h"
#include "range.h"
#include "range_source.h"
#include "range_chaser.h"

/**
 * State machine for carrying out an alignment, which usually consists
 * of a series of phases that conduct different alignments using
 * different backtracking constraints.
 *
 * Each Aligner should have a dedicated PatternSourcePerThread.
 */
class Aligner {
public:
	Aligner(bool rangeMode, uint32_t seed) :
		patsrc_(NULL), bufa_(NULL), bufb_(NULL),
		rangeMode_(rangeMode), seed_(seed)
	{ }

	virtual ~Aligner() { }
	/// Advance the range search by one memory op
	virtual bool advance() = 0;
	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() = 0;

	/// Prepare Aligner for the next read
	virtual void setQuery(PatternSourcePerThread *patsrc) {
		assert(patsrc != NULL);
		patsrc_ = patsrc;
		bufa_ = &patsrc->bufa();
		assert(bufa_ != NULL);
		bufb_ = &patsrc->bufb();
		alen_ = bufa_->length();
		blen_ = (bufb_ != NULL) ? bufb_->length() : 0;
		rand_.init(seed_ + genRandSeed(bufa_->patFw, bufa_->qualFw, bufa_->name));
	}

protected:

	// Current read pair
	PatternSourcePerThread* patsrc_;
	ReadBuf* bufa_;
	uint32_t alen_;
	ReadBuf* bufb_;
	uint32_t blen_;
	// RandomSource for choosing alignments to report from ranges
	bool rangeMode_;
	uint32_t seed_;
	RandomSource rand_;
};

/**
 * Abstract parent factory class for constructing aligners of all kinds.
 */
class AlignerFactory {
public:
	virtual ~AlignerFactory() { }
	virtual Aligner* create() const = 0;

	/**
	 * Allocate a vector of n Aligners; use destroy(std::vector...) to
	 * free the memory.
	 */
	virtual std::vector<Aligner*>* create(uint32_t n) const {
		std::vector<Aligner*>* v = new std::vector<Aligner*>;
		for(uint32_t i = 0; i < n; i++) {
			v->push_back(create());
			assert(v->back() != NULL);
		}
		return v;
	}

	/// Free memory associated with the aligner
	virtual void destroy(Aligner* al) const {
		assert(al != NULL);
		// Free the Aligner
		delete al;
	}

	/// Free memory associated with an aligner list
	virtual void destroy(std::vector<Aligner*>* als) const {
		assert(als != NULL);
		// Free all of the Aligners
		for(size_t i = 0; i < als->size(); i++) {
			if((*als)[i] != NULL) {
				delete (*als)[i];
				(*als)[i] = NULL;
			}
		}
		// Free the vector
		delete als;
	}
};

/**
 * Coordinates multiple aligners of the same type (i.e. either all
 * single-end or all paired-end).
 */
class MultiAligner {
public:
	MultiAligner(
			uint32_t n,
			uint32_t qUpto,
			const AlignerFactory& alignFact,
			const PatternSourcePerThreadFactory& patsrcFact) :
			n_(n), qUpto_(qUpto),
			alignFact_(alignFact), patsrcFact_(patsrcFact),
			aligners_(NULL), patsrcs_(NULL)
	{
		aligners_ = alignFact_.create(n_);
		assert(aligners_ != NULL);
		patsrcs_ = patsrcFact_.create(n_);
		assert(patsrcs_ != NULL);
	}

	/// Free memory associated with the aligners and their pattern sources.
	virtual ~MultiAligner() {
		alignFact_.destroy(aligners_);
		patsrcFact_.destroy(patsrcs_);
	}

	/**
	 * Advance an array of aligners in parallel, using prefetches to
	 * try to hide all the latency.
	 */
	void run() {
		bool done = false;
		while(!done) {
			done = true;
			for(uint32_t i = 0; i < n_; i++) {
				if(!(*aligners_)[i]->done()) {
					// Advance an aligner already in progress
					done = false;
					(*aligners_)[i]->advance();
				} else if(qUpto_ > 0) {
					// Get a new read and initialize an aligner with it
					(*patsrcs_)[i]->nextReadPair();
					if(!(*patsrcs_)[i]->empty()) {
						qUpto_--;
						(*aligners_)[i]->setQuery((*patsrcs_)[i]);
						assert(!(*aligners_)[i]->done());
						done = false;
					} else {
						// No more reads; if done == true, it remains
						// true
					}
				} else {
					// Past read limit; if done == true, it remains true
				}
			}
		}
	}

protected:
	uint32_t n_;     /// Number of aligners
	uint32_t qUpto_; /// Number of reads to align before stopping
	const AlignerFactory&                  alignFact_;
	const PatternSourcePerThreadFactory&   patsrcFact_;
	std::vector<Aligner *>*                aligners_;
	std::vector<PatternSourcePerThread *>* patsrcs_;
};

/**
 * Coordinates multiple single-end and paired-end aligners, routing
 * reads to one or the other type as appropriate.
 */
class MixedMultiAligner {
public:
	MixedMultiAligner(
			uint32_t n,
			uint32_t qUpto,
			const AlignerFactory& alignSEFact,
			const AlignerFactory& alignPEFact,
			const PatternSourcePerThreadFactory& patsrcFact) :
			n_(n), qUpto_(qUpto),
			alignSEFact_(alignSEFact),
			alignPEFact_(alignPEFact),
			patsrcFact_(patsrcFact),
			alignersSE_(NULL),
			alignersPE_(NULL),
			seOrPe_(NULL),
			patsrcs_(NULL)
	{
		// Instantiate all single-end aligners
		alignersSE_ = alignSEFact_.create(n_);
		assert(alignersSE_ != NULL);
		// Instantiate all paired-end aligners
		alignersPE_ = alignPEFact_.create(n_);
		assert(alignersPE_ != NULL);
		// Allocate array of boolean flags indicating whether each of
		// the slots is currently using the single-end or paired-end
		// aligner
		seOrPe_ = new bool[n_];
		// Instantiate all read sources
		patsrcs_ = patsrcFact_.create(n_);
		assert(patsrcs_ != NULL);
	}

	/// Free memory associated with the aligners and their pattern sources.
	virtual ~MixedMultiAligner() {
		alignSEFact_.destroy(alignersSE_);
		alignPEFact_.destroy(alignersPE_);
		patsrcFact_.destroy(patsrcs_);
		delete[] seOrPe_;
	}

	/**
	 * Advance an array of aligners in parallel, using prefetches to
	 * try to hide all the latency.
	 */
	void run(bool verbose = false) {
		bool done = false;
		while(!done) {
			done = true;
			for(uint32_t i = 0; i < n_; i++) {
				Aligner *al = seOrPe_[i] ? (*alignersSE_)[i] :
				                           (*alignersPE_)[i];
				if(!al->done()) {
					// Advance an aligner already in progress; this is
					// the common case
					done = false;
					al->advance();
				} else if(qUpto_ > 0) {
					// Feed a new read to a vacant aligner
					PatternSourcePerThread *ps = (*patsrcs_)[i];
					// Get a new read
					ps->nextReadPair();
					if(!ps->empty()) {
						qUpto_--;
						if(ps->paired()) {
							// Read currently in buffer is paired-end
							if(verbose) cout << "Paired input: " << ps->bufa().patFw << "," << ps->bufb().patFw << endl;
							(*alignersPE_)[i]->setQuery(ps);
							seOrPe_[i] = false; // false -> paired
						} else {
							// Read currently in buffer is single-end
							if(verbose) cout << "Unpaired input: " << ps->bufa().patFw << endl;
							(*alignersSE_)[i]->setQuery(ps);
							seOrPe_[i] = true; // true = unpaired
						}
						done = false;
					} else {
						// No more reads; if done == true, it remains
						// true
					}
				} else {
					// Past read limit; if done == true, it remains true
				}
			}
		}
	}

protected:
	uint32_t n_;     /// Number of aligners
	uint32_t qUpto_; /// Number of reads to align before stopping
	const AlignerFactory&                  alignSEFact_;
	const AlignerFactory&                  alignPEFact_;
	const PatternSourcePerThreadFactory&   patsrcFact_;
	std::vector<Aligner *>*                alignersSE_;
	std::vector<Aligner *>*                alignersPE_;
	bool *                                 seOrPe_;
	std::vector<PatternSourcePerThread *>* patsrcs_;
};

/**
 * An aligner for finding exact matches of unpaired reads.  Always
 * tries the forward-oriented version of the read before the reverse-
 * oriented read.
 */
template<typename TRangeSource, typename TContMan>
class UnpairedAlignerV1 : public Aligner {
	typedef RangeSourceDriver<TRangeSource, TContMan> TDriver;
public:
	UnpairedAlignerV1(
		EbwtSearchParams<String<Dna> >* params,
		TDriver* driverFw, TDriver* driverRc,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		vector<String<Dna5> >& os,
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(rangeMode, seed),
		doneFw_(true), done_(true),
		chaseFw_(false), chaseRc_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		params_(params),
		rchase_(rand_),
		driverFw_(driverFw),
		driverRc_(driverRc)
	{
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
		assert(driverFw_ != NULL);
		assert(driverRc_ != NULL);
	}

	virtual ~UnpairedAlignerV1() {
		delete driverFw_; driverFw_ = NULL;
		delete driverRc_; driverRc_ = NULL;
		delete params_;   params_   = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		driverFw_->setQuery(patsrc);
		driverRc_->setQuery(patsrc);
		doneFw_  = false;
		done_    = false;
		chaseFw_ = false;
		chaseRc_ = false;
	}

	/**
	 * Helper for reporting an alignment.
	 */
	inline bool report(const Range& ra,
	                   uint32_t first,
	                   uint32_t second,
	                   uint32_t tlen,
	                   bool fw)
	{
		bool ebwtFw = ra.ebwt->fw();
		return params_->reportHit(
				fw ? (ebwtFw? bufa_->patFw   : bufa_->patFwRev) :
				     (ebwtFw? bufa_->patRc   : bufa_->patRcRev),
				fw ? (ebwtFw? &bufa_->qualFw : &bufa_->qualFwRev) :
				     (ebwtFw? &bufa_->qualRc : &bufa_->qualRcRev),
				&bufa_->name,
				ebwtFw,
				ra.mms,                   // mismatch positions
				ra.refcs,                 // reference characters for mms
				ra.numMms,                // # mismatches
				make_pair(first, second), // position
				make_pair(ra.top, ra.bot),// arrows
				tlen,                     // textlen
				alen_,                    // qlen
				ra.stratum,               // alignment stratum
				ra.bot - ra.top - 1);     // # other hits
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual bool advance() {
		assert(!done_);
		assert(!chaseFw_ || !chaseRc_);
		if(chaseFw_) {
			assert(!rangeMode_);
			assert(driverFw_->foundRange());
			if(!rchase_.foundOff() && !rchase_.done()) {
				rchase_.advance();
				return false;
			}
			if(rchase_.foundOff()) {
				done_ = report(driverFw_->range(), rchase_.off().first,
				               rchase_.off().second, rchase_.tlen(), true);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				chaseFw_ = false;
				doneFw_ = driverFw_->done();
			}
		} else if(chaseRc_) {
			assert(!rangeMode_);
			assert(driverRc_->foundRange());
			if(!rchase_.foundOff() && !rchase_.done()) {
				rchase_.advance();
				return false;
			}
			if(rchase_.foundOff()) {
				done_ = report(driverRc_->range(), rchase_.off().first,
				               rchase_.off().second, rchase_.tlen(), false);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				chaseRc_ = false;
				done_ = driverRc_->done();
			}
		}
		// Still advancing a
		if(!done_ && !chaseFw_ && !chaseRc_) {
			if(!doneFw_) {
				driverFw_->advance();
				if(driverFw_->foundRange()) {
					const Range& ra = driverFw_->range();
					if(rangeMode_) {
						done_ = report(ra, ra.top, ra.bot, 0, true);
					} else {
						rchase_.setTopBot(ra.top, ra.bot, alen_, driverFw_->curEbwt());
						if(rchase_.foundOff()) {
							done_ = report(ra, rchase_.off().first,
							               rchase_.off().second, rchase_.tlen(),
							               true);
							rchase_.reset();
						}
						if(!rchase_.done()) {
							// Keep chasing this range
							chaseFw_ = true;
						}
					}
				}
				if(!doneFw_ && !chaseFw_) {
					doneFw_ = driverFw_->done();
				}
			} else {
				driverRc_->advance();
				// Advance the RangeSource for the reverse-complement read
				if(driverRc_->foundRange()) {
					const Range& ra = driverRc_->range();
					if(rangeMode_) {
						done_ = report(ra, ra.top, ra.bot, 0, false);
					} else {
						rchase_.setTopBot(ra.top, ra.bot, alen_, driverRc_->curEbwt());
						if(rchase_.foundOff()) {
							done_ = report(ra, rchase_.off().first,
							               rchase_.off().second, rchase_.tlen(),
							               false);
							rchase_.reset();
						}
						if(!rchase_.done()) {
							// Keep chasing this range
							chaseRc_ = true;
						}
					}
				}
				if(!done_ && !chaseRc_) {
					done_ = driverRc_->done();
				}
			}
		}
		if(done_) {
			sinkPt_->finishRead(*patsrc_, true);
		}
		return done_;
	}

	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() {
		return done_;
	}

protected:
	// Progress state
	bool doneFw_;
	bool done_;
	bool chaseFw_;
	bool chaseRc_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> >* params_;

	// State for getting alignments from ranges statefully
	WideRandomScanningRangeChaser<String<Dna> > rchase_;

	// Range-finding state
	TDriver* driverFw_;
	TDriver* driverRc_;
};

/**
 * An aligner for finding exact matches of paired reads.
 */
template<typename TRangeSource, typename TContMan>
class PairedAlignerV1 : public Aligner {

	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef std::vector<U32Pair> U32PairVec;
	typedef std::vector<Range> TRangeVec;
	typedef RangeSourceDriver<TRangeSource, TContMan> TDriver;

public:
	PairedAlignerV1(
		EbwtSearchParams<String<Dna> >* params,
		TDriver* driver1Fw, TDriver* driver1Rc,
		TDriver* driver2Fw, TDriver* driver2Rc,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		bool fw1, bool fw2,
		uint32_t minInsert,
		uint32_t maxInsert,
		uint32_t symCeiling,
		vector<String<Dna5> >& os,
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(rangeMode, seed),
		doneFw_(true), doneFwFirst_(true), done_(true),
		chase1Fw_(false), chase1Rc_(false),
		chase2Fw_(false), chase2Rc_(false),
		delayedChase1Fw_(false), delayedChase1Rc_(false),
		delayedChase2Fw_(false), delayedChase2Rc_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		params_(params),
		minInsert_(minInsert),
		maxInsert_(maxInsert),
		symCeiling_(symCeiling),
		fw1_(fw1), fw2_(fw2),
		rchase_(rand_),
		driver1Fw_(driver1Fw), driver1Rc_(driver1Rc),
		offs1FwSz_(0), offs1RcSz_(0),
		driver2Fw_(driver2Fw), driver2Rc_(driver2Rc),
		offs2FwSz_(0), offs2RcSz_(0),

		chaseL_fw_       (fw1_ ? chase1Fw_        : chase1Rc_),
		chaseR_fw_       (fw2_ ? chase2Fw_        : chase2Rc_),
		delayedchaseL_fw_(fw1_ ? delayedChase1Fw_ : delayedChase1Rc_),
		delayedchaseR_fw_(fw2_ ? delayedChase2Fw_ : delayedChase2Rc_),
		drL_fw_          (fw1_ ? *driver1Fw_      : *driver1Rc_),
		drR_fw_          (fw2_ ? *driver2Fw_      : *driver2Rc_),
		offsL_fw_        (fw1_ ? offs1Fw_         : offs1Rc_),
		offsR_fw_        (fw2_ ? offs2Fw_         : offs2Rc_),
		rangesL_fw_      (fw1_ ? ranges1Fw_       : ranges1Rc_),
		rangesR_fw_      (fw2_ ? ranges2Fw_       : ranges2Rc_),
		offsLsz_fw_      (fw1_ ? offs1FwSz_       : offs1RcSz_),
		offsRsz_fw_      (fw2_ ? offs2FwSz_       : offs2RcSz_),

		chaseL_rc_       (fw2_ ? chase2Rc_        : chase2Fw_),
		chaseR_rc_       (fw1_ ? chase1Rc_        : chase1Fw_),
		delayedchaseL_rc_(fw2_ ? delayedChase2Rc_ : delayedChase2Fw_),
		delayedchaseR_rc_(fw1_ ? delayedChase1Rc_ : delayedChase1Fw_),
		drL_rc_          (fw2_ ? *driver2Rc_      : *driver2Fw_),
		drR_rc_          (fw1_ ? *driver1Rc_      : *driver1Fw_),
		offsL_rc_        (fw2_ ? offs2Rc_         : offs2Fw_),
		offsR_rc_        (fw1_ ? offs1Rc_         : offs1Fw_),
		rangesL_rc_      (fw2_ ? ranges2Rc_       : ranges2Fw_),
		rangesR_rc_      (fw1_ ? ranges1Rc_       : ranges1Fw_),
		offsLsz_rc_      (fw2_ ? offs2RcSz_       : offs2FwSz_),
		offsRsz_rc_      (fw1_ ? offs1RcSz_       : offs1FwSz_),

		chaseL_       (&chaseL_fw_),
		chaseR_       (&chaseR_fw_),
		delayedchaseL_(&delayedchaseL_fw_),
		delayedchaseR_(&delayedchaseR_fw_),
		drL_          (&drL_fw_),
		drR_          (&drR_fw_),
		offsL_        (&offsL_fw_),
		offsR_        (&offsR_fw_),
		rangesL_      (&rangesL_fw_),
		rangesR_      (&rangesR_fw_),
		offsLsz_      (&offsLsz_fw_),
		offsRsz_      (&offsRsz_fw_),
		donePair_     (&doneFw_),
		fwL_(fw1),
		fwR_(fw2)
	{
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
		assert(driver1Fw_ != NULL);
		assert(driver1Rc_ != NULL);
		assert(driver2Fw_ != NULL);
		assert(driver2Rc_ != NULL);
	}

	virtual ~PairedAlignerV1() {
		delete driver1Fw_; driver1Fw_ = NULL;
		delete driver1Rc_; driver1Rc_ = NULL;
		delete driver2Fw_; driver2Fw_ = NULL;
		delete driver2Rc_; driver2Rc_ = NULL;
		delete params_;    params_    = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		assert(!patsrc->bufa().empty());
		Aligner::setQuery(patsrc); // set fields & random seed
		assert(!patsrc->bufb().empty());
		// Give all of the drivers pointers to the relevant read info
		driver1Fw_->setQuery(patsrc, true  /* mate1 */);
		driver1Rc_->setQuery(patsrc, true  /* mate1 */);
		driver2Fw_->setQuery(patsrc, false /* mate2 */);
		driver2Rc_->setQuery(patsrc, false /* mate2 */);
		// Neither orientation is done
		doneFw_   = false;
		doneFwFirst_ = true;
		done_     = false;
		// No ranges are being chased yet
		chase1Fw_ = false;
		chase1Rc_ = false;
		chase2Fw_ = false;
		chase2Rc_ = false;
		delayedChase1Fw_ = false;
		delayedChase1Rc_ = false;
		delayedChase2Fw_ = false;
		delayedChase2Rc_ = false;
		// Clear all intermediate ranges
		offs1Fw_.clear(); offs1Rc_.clear();
		offs2Fw_.clear(); offs2Rc_.clear();
		ranges1Fw_.clear(); ranges1Rc_.clear();
		ranges2Fw_.clear(); ranges2Rc_.clear();
		offs1FwSz_ = offs1RcSz_ = offs2FwSz_ = offs2RcSz_ = 0;
		chaseL_        = &chaseL_fw_;
		chaseR_        = &chaseR_fw_;
		delayedchaseL_ = &delayedchaseL_fw_;
		delayedchaseR_ = &delayedchaseR_fw_;
		drL_           = &drL_fw_;
		drR_           = &drR_fw_;
		offsL_         = &offsL_fw_;
		offsR_         = &offsR_fw_;
		rangesL_       = &rangesL_fw_;
		rangesR_       = &rangesR_fw_;
		offsLsz_       = &offsLsz_fw_;
		offsRsz_       = &offsRsz_fw_;
		donePair_      = &doneFw_;
		fwL_           = fw1_;
		fwR_           = fw2_;
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 *
	 * A call to this function does one of many things:
	 * 1. Advance a RangeSourceDriver and check if it found a new range
	 * 2. Advance a RowChaseDriver and check if it found a reference
	 *    offset for a an alignment in a range
	 */
	virtual bool advance() {
		assert(!done_);
		bool verbose = false;
		if(doneFw_ && doneFwFirst_) {
			chaseL_        = &chaseL_rc_;
			chaseR_        = &chaseR_rc_;
			delayedchaseL_ = &delayedchaseL_rc_;
			delayedchaseR_ = &delayedchaseR_rc_;
			drL_           = &drL_rc_;
			drR_           = &drR_rc_;
			offsL_         = &offsL_rc_;
			offsR_         = &offsR_rc_;
			rangesL_       = &rangesL_rc_;
			rangesR_       = &rangesR_rc_;
			offsLsz_       = &offsLsz_rc_;
			offsRsz_       = &offsRsz_rc_;
			donePair_      = &done_;
			fwL_           = !fw2_;
			fwR_           = !fw1_;
			doneFwFirst_   = false;
		}
		advanceOrientation(!doneFw_, verbose);
		if(done_) {
			sinkPt_->finishRead(*patsrc_, true);
		}
		return done_;
	}

	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() {
		return done_;
	}

protected:

	/**
	 * Helper for reporting a pair of alignments.  As of now, we report
	 * a paired alignment by reporting two consecutive alignments, one
	 * for each mate.
	 */
	bool report(const Range& rL, // range for mate1
	            const Range& rR, // range for mate2
	            uint32_t first,  // ref idx
	            uint32_t upstreamOff,// offset for mate1
	            uint32_t dnstreamOff,// offset for mate2
	            uint32_t tlen,   // length of ref
	            bool fwL,        // whether mate1 is in fw orientation
	            bool fwR,        // whether mate2 is in fw orientation
	            bool pairFw)    // whether the pair is being mapped to fw strand
	{
		assert_lt(upstreamOff, dnstreamOff);
		uint32_t spreadL = rL.bot - rL.top;
		uint32_t spreadR = rR.bot - rR.top;
		uint32_t oms = min(spreadL, spreadR) - 1;
		bool ebwtFwL = rL.ebwt->fw();
		bool ebwtFwR = rR.ebwt->fw();
		ReadBuf* bufL = pairFw ? bufa_ : bufb_;
		ReadBuf* bufR = pairFw ? bufb_ : bufa_;
		uint32_t lenL = pairFw ? alen_ : blen_;
		uint32_t lenR = pairFw ? blen_ : alen_;
		bool ret;
		params_->setFw(fwL);
		// Print upstream mate first
		ret = params_->reportHit(
				fwL ? (ebwtFwL?  bufL->patFw  :  bufL->patFwRev) :
					  (ebwtFwL?  bufL->patRc  :  bufL->patRcRev),
				fwL ? (ebwtFwL? &bufL->qualFw : &bufL->qualFwRev) :
					  (ebwtFwL? &bufL->qualRc : &bufL->qualRcRev),
				&bufL->name,
				ebwtFwL,
				rL.mms,                       // mismatch positions
				rL.refcs,                     // reference characters for mms
				rL.numMms,                    // # mismatches
				make_pair(first, upstreamOff),// position
				make_pair(rL.top, rL.bot),    // arrows
				tlen,                         // textlen
				lenL,                         // qlen
				rL.stratum,                   // alignment stratum
				oms);                         // # other hits
		assert(!ret);
		params_->setFw(fwR);
		ret = params_->reportHit(
				fwR ? (ebwtFwR?  bufR->patFw  :  bufR->patFwRev) :
					  (ebwtFwR?  bufR->patRc  :  bufR->patRcRev),
				fwR ? (ebwtFwR? &bufR->qualFw : &bufR->qualFwRev) :
					  (ebwtFwR? &bufR->qualRc : &bufR->qualRcRev),
				&bufR->name,
				ebwtFwR,
				rR.mms,                       // mismatch positions
				rR.refcs,                     // reference characters for mms
				rR.numMms,                    // # mismatches
				make_pair(first, dnstreamOff),// position
				make_pair(rR.top, rR.bot),    // arrows
				tlen,                         // textlen
				lenR,                         // qlen
				rR.stratum,                   // alignment stratum
				oms);                         // # other hits
		return ret;
	}

	/**
	 * Given a new reference position where one mate aligns, add it
	 * to the list of positions for that mate and reconcile it against
	 * all of the positions for the other mate, reporting paired
	 * alignments wherever the mating constraints are met.
	 *
	 * The parameters ending in 1 pertain to the upstream mate.
	 */
	bool reconcileAndAdd(const U32Pair& h,
	                     bool newFromL,
	                     U32PairVec& offsL,
	                     U32PairVec& offsR,
	                     TRangeVec& rangesL,
	                     TRangeVec& rangesR,
		                 TDriver& drL,
		                 TDriver& drR,
		                 bool fwL,
		                 bool fwR,
		                 bool pairFw,
		                 bool verbose = false)
	{
		assert_eq(offsL.size(), rangesL.size());
		assert_eq(offsR.size(), rangesR.size());
		// For each known hit for the other mate, check if this new
		// alignment can be mated with it.  If so, report the mates.
		size_t offsDstSz = newFromL ? offsR.size() : offsL.size();
		if(verbose) cout << "reconcileAndAdd called" << endl;
		if(offsDstSz > 0) {
			// Start in a random spot in the offsR array and scan
			// linearly
			uint32_t rand = rand_.nextU32() % offsDstSz;
			for(size_t i = 0; i < offsDstSz; i++) {
				rand++;
				if(rand == offsDstSz) rand = 0;
				const U32Pair& h2 = newFromL ? offsR[rand] : offsL[rand];
				if(h.first == h2.first) {
					if(verbose) cout << "Found pair on some reference" << endl;
					// Incoming hit hits same reference as buffered
					uint32_t left  = newFromL ? h.second : h2.second;
					uint32_t right = newFromL ? h2.second : h.second;
					if(right > left) {
						if(verbose) cout << "...and in right orientation" << endl;
						uint32_t gap = right - left;
						if(gap >= minInsert_ && gap <= maxInsert_) {
							if(verbose) cout << "...and with the right sized gap; reporting" << endl;
							// Gap between the two alignments satisfies
							// the paired-end policy, so we can report
							// them
							const Range& rL = newFromL ? drL.range() : rangesL[rand];
							const Range& rR = newFromL ? rangesR[rand] : drR.range();
							if(report(
							       rL, // upstream (left) range
							       rR, // downstream (right) range
								   h.first,     // reference target
								   left, right, // up/downstream offsets
								   // _plen array is same both Fw and Bw
								   drL.curEbwt()->_plen[h.first],
								   fwL, // true -> upstream mate is fw
								   fwR, // true -> downstream mate is fw
								   pairFw)) return true;
						}
					}
				}
			}
		}
		// Push new reference locus and range onto the appropriate
		// parallel mate buffers
		if(newFromL) { offsL.push_back(h); rangesL.push_back(drL.range()); }
		else         { offsR.push_back(h); rangesR.push_back(drR.range()); }
		return false;
	}

	/**
	 * Advance paired-end alignment.
	 */
	void advanceOrientation(bool pairFw, bool verbose = false) {
		assert(!done_);
		assert(!*donePair_);
		assert(!*chaseL_ || !*chaseR_);
		if(*chaseL_) {
			assert(!rangeMode_);
			assert(!*delayedchaseL_);
			assert(drL_->foundRange());
			if(!rchase_.foundOff() && !rchase_.done()) {
				// Keep trying to resolve the reference loci for
				// alignments in this range
				rchase_.advance();
				return;
			} else if(rchase_.foundOff()) {
				// Resolve this against the reference loci
				// determined for the other mate
				done_ = reconcileAndAdd(rchase_.off(), true /* new entry is from 1 */,
				                        *offsL_, *offsR_, *rangesL_, *rangesR_, *drL_, *drR_,
				                        fwL_, fwR_, pairFw, verbose);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				*chaseL_ = false;
				if(verbose) cout << "Done with case for first mate" << endl;
				if(*delayedchaseR_) {
					// Start chasing the delayed range
					if(verbose) cout << "Resuming delayed chase for second mate" << endl;
					assert(drR_->foundRange());
					uint32_t top = drR_->range().top; uint32_t bot = drR_->range().bot;
					rchase_.setTopBot(top, bot, drR_->qlen(), drR_->curEbwt());
					*chaseR_ = true;
					*delayedchaseR_ = false;
				}
			}
		} else if(*chaseR_) {
			assert(!rangeMode_);
			assert(!*delayedchaseR_);
			assert(drR_->foundRange());
			if(!rchase_.foundOff() && !rchase_.done()) {
				// Keep trying to resolve the reference loci for
				// alignments in this range
				rchase_.advance();
				return;
			} else if(rchase_.foundOff()) {
				// Resolve this against the reference loci
				// determined for the other mate
				done_ = reconcileAndAdd(rchase_.off(), false /* new entry is from 2 */,
				                        *offsL_, *offsR_, *rangesL_, *rangesR_, *drL_, *drR_,
				                        fwL_, fwR_, pairFw, verbose);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				*chaseR_ = false;
				if(verbose) cout << "Done with case for second mate" << endl;
				if(*delayedchaseL_) {
					// Start chasing the delayed range
					if(verbose) cout << "Resuming delayed chase for first mate" << endl;
					assert(drL_->foundRange());
					uint32_t top = drL_->range().top; uint32_t bot = drL_->range().bot;
					rchase_.setTopBot(top, bot, drL_->qlen(), drL_->curEbwt());
					*chaseL_ = true;
					*delayedchaseL_ = false;
				}
			}
		}
		if(!done_ && !*donePair_ && !*chaseL_ && !*chaseR_) {
			// Search for more ranges for whichever mate currently has
			// fewer ranges
			if((*offsLsz_ < *offsRsz_ || drR_->done()) && !drL_->done()) {
				// If there are no more ranges for the other mate and
				// there are no candidate alignments either, then we're
				// not going to find a paired alignment in this
				// orientation.
				if(drR_->done() && *offsRsz_ == 0) {
					// Give up on this orientation
					if(verbose) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 1" << endl;
					*donePair_ = true;
					return;
				}
				assert(!*delayedchaseL_);
				drL_->advance();
				if(drL_->foundRange()) {
					// Add the size of this range to the total for this mate
					*offsLsz_ += (drL_->range().bot - drL_->range().top);
					if(*offsLsz_ > symCeiling_ && *offsRsz_ > symCeiling_) {
						// Too many candidates for both mates; abort
						// without any more searching
						done_ = true;
						return;
					}
					if(*offsRsz_ == 0) {
						// Delay chasing this range; we delay to avoid
						// needlessly chasing rows in this range when
						// the other mate doesn't end up aligning
						// anywhere
						if(verbose) cout << "Delaying a chase for first mate" << endl;
						*delayedchaseL_ = true;
					} else {
						// Start chasing this range
						if(verbose) cout << "Chasing a range for first mate" << endl;
						rchase_.setTopBot(drL_->range().top, drL_->range().bot, drL_->qlen(), drL_->curEbwt());
						*chaseL_ = true;
					}
				}
			} else if(!drR_->done()) {
				// If there are no more ranges for the other mate and
				// there are no candidate alignments either, then we're
				// not going to find a paired alignment in this
				// orientation.
				if(drL_->done() && *offsLsz_ == 0) {
					// Give up on this orientation
					if(verbose) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 2" << endl;
					*donePair_ = true;
					return;
				}
				assert(!*delayedchaseR_);
				drR_->advance();
				if(drR_->foundRange()) {
					// Add the size of this range to the total for this mate
					*offsRsz_ += (drR_->range().bot - drR_->range().top);
					if(*offsLsz_ > symCeiling_ && *offsRsz_ > symCeiling_) {
						// Too many candidates for both mates; abort
						// without any more searching
						done_ = true;
						return;
					}
					if(*offsLsz_ == 0) {
						// Delay chasing this range; we delay to avoid
						// needlessly chasing rows in this range when
						// the other mate doesn't end up aligning
						// anywhere
						if(verbose) cout << "Delaying a chase for second mate" << endl;
						*delayedchaseR_ = true;
					} else {
						// Start chasing this range
						if(verbose) cout << "Chasing a range for second mate" << endl;
						rchase_.setTopBot(drR_->range().top, drR_->range().bot, drR_->qlen(), drR_->curEbwt());
						*chaseR_ = true;
					}
				}
			} else {
				// Finished processing ranges for both mates
				assert(drL_->done() && drR_->done());
				*donePair_ = true;
			}
		}
	}

	// Progress state
	bool doneFw_;   // finished with forward orientation of both mates?
	bool doneFwFirst_;
	bool done_;

	bool chase1Fw_;
	bool chase1Rc_;
	bool chase2Fw_;
	bool chase2Rc_;

	bool delayedChase1Fw_;
	bool delayedChase1Rc_;
	bool delayedChase2Fw_;
	bool delayedChase2Rc_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> >* params_;

	// Paired-end boundaries
	const uint32_t minInsert_;
	const uint32_t maxInsert_;

	// If both reads align >= symCeiling times, then immediately give
	// up on reporting a paired alignment
	const uint32_t symCeiling_;

	// Orientation of upstream/downstream mates when aligning to
	// forward strand
	const bool fw1_;
	const bool fw2_;

	// State for getting alignments from ranges statefully
	WideRandomScanningRangeChaser<String<Dna> > rchase_;

	// Range-finding state for first mate
	TDriver*      driver1Fw_;
	TDriver*      driver1Rc_;
	U32PairVec    offs1Fw_;
	TRangeVec     ranges1Fw_;
	uint32_t      offs1FwSz_; // total size of all ranges found in this category
	U32PairVec    offs1Rc_;
	TRangeVec     ranges1Rc_;
	uint32_t      offs1RcSz_; // total size of all ranges found in this category

	// Range-finding state for second mate
	TDriver*      driver2Fw_;
	TDriver*      driver2Rc_;
	U32PairVec    offs2Fw_;
	TRangeVec     ranges2Fw_;
	uint32_t      offs2FwSz_; // total size of all ranges found in this category
	U32PairVec    offs2Rc_;
	TRangeVec     ranges2Rc_;
	uint32_t      offs2RcSz_; // total size of all ranges found in this category

	bool&       chaseL_fw_;
	bool&       chaseR_fw_;
	bool&       delayedchaseL_fw_;
	bool&       delayedchaseR_fw_;
	TDriver&    drL_fw_;
	TDriver&    drR_fw_;
	U32PairVec& offsL_fw_;
	U32PairVec& offsR_fw_;
	TRangeVec&  rangesL_fw_;
	TRangeVec&  rangesR_fw_;
	uint32_t&   offsLsz_fw_;
	uint32_t&   offsRsz_fw_;

	bool&       chaseL_rc_;
	bool&       chaseR_rc_;
	bool&       delayedchaseL_rc_;
	bool&       delayedchaseR_rc_;
	TDriver&    drL_rc_;
	TDriver&    drR_rc_;
	U32PairVec& offsL_rc_;
	U32PairVec& offsR_rc_;
	TRangeVec&  rangesL_rc_;
	TRangeVec&  rangesR_rc_;
	uint32_t&   offsLsz_rc_;
	uint32_t&   offsRsz_rc_;

	bool*       chaseL_;
	bool*       chaseR_;
	bool*       delayedchaseL_;
	bool*       delayedchaseR_;
	TDriver*    drL_;
	TDriver*    drR_;
	U32PairVec* offsL_;
	U32PairVec* offsR_;
	TRangeVec*  rangesL_;
	TRangeVec*  rangesR_;
	uint32_t*   offsLsz_;
	uint32_t*   offsRsz_;
	bool*       donePair_;
	bool        fwL_;
	bool        fwR_;
};

#endif /* ALIGNER_H_ */
