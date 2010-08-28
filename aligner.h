/**
 * aligner.h
 *
 * A generic class providing a stateful way to find alignments.
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <iostream>
#include <set>
#include <stdint.h>
#include "seqan/sequence.h"
#include "assert_helpers.h"
#include "ebwt.h"
#include "pat.h"
#include "range.h"
#include "range_source.h"
#include "range_chaser.h"
#include "ref_aligner.h"
#include "reference.h"
#include "aligner_metrics.h"
#include "search_globals.h"

/**
 * State machine for carrying out an alignment, which usually consists
 * of a series of phases that conduct different alignments using
 * different backtracking constraints.
 *
 * Each Aligner should have a dedicated PatternSourcePerThread.
 */
class Aligner {
public:
	Aligner(bool _done, bool rangeMode) :
		done(_done), patsrc_(NULL), bufa_(NULL), bufb_(NULL),
		rangeMode_(rangeMode)
	{ }

	virtual ~Aligner() { }
	/// Advance the range search by one memory op
	virtual bool advance() = 0;

	/// Prepare Aligner for the next read
	virtual void setQuery(PatternSourcePerThread *patsrc) {
		assert(patsrc != NULL);
		patsrc_ = patsrc;
		bufa_ = &patsrc->bufa();
		assert(bufa_ != NULL);
		bufb_ = &patsrc->bufb();
		alen_ = bufa_->length();
		blen_ = (bufb_ != NULL) ? bufb_->length() : 0;
		rand_.init(bufa_->seed);
	}

	/**
	 * Set to true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	bool done;

protected:

	// Current read pair
	PatternSourcePerThread* patsrc_;
	ReadBuf* bufa_;
	uint32_t alen_;
	ReadBuf* bufb_;
	uint32_t blen_;
	bool rangeMode_;
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
				if(!(*aligners_)[i]->done) {
					// Advance an aligner already in progress
					done = false;
					(*aligners_)[i]->advance();
				} else {
					// Get a new read and initialize an aligner with it
					(*patsrcs_)[i]->nextReadPair();
					if(!(*patsrcs_)[i]->empty() && (*patsrcs_)[i]->patid() < qUpto_) {
						(*aligners_)[i]->setQuery((*patsrcs_)[i]);
						assert(!(*aligners_)[i]->done);
						done = false;
					} else {
						// No more reads; if done == true, it remains
						// true
					}
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
		for(uint32_t i = 0; i < n_; i++) {
			seOrPe_[i] = true;
		}
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
		bool first = true;
		if(n_ == 1) {
			Aligner *al = seOrPe_[0] ? (*alignersSE_)[0] : (*alignersPE_)[0];
			PatternSourcePerThread *ps = (*patsrcs_)[0];
			while(!done) {
				done = true;
				if(!first && !al->done) {
					// Advance an aligner already in progress; this is
					// the common case
					done = false;
					al->advance();
				} else {
					// Get a new read
					ps->nextReadPair();
					if(ps->patid() < qUpto_ && !ps->empty()) {
						if(ps->paired()) {
							// Read currently in buffer is paired-end
							(*alignersPE_)[0]->setQuery(ps);
							al = (*alignersPE_)[0];
							seOrPe_[0] = false; // false -> paired
						} else {
							// Read currently in buffer is single-end
							(*alignersSE_)[0]->setQuery(ps);
							al = (*alignersSE_)[0];
							seOrPe_[0] = true; // true = unpaired
						}
						done = false;
					} else {
						// No more reads; if done == true, it remains
						// true
					}
				}
				first = false;
			}
		} else {
			while(!done) {
				done = true;
				for(uint32_t i = 0; i < n_; i++) {
					Aligner *al = seOrPe_[i] ? (*alignersSE_)[i] :
											   (*alignersPE_)[i];
					if(!first && !al->done) {
						// Advance an aligner already in progress; this is
						// the common case
						done = false;
						al->advance();
					} else {
						// Feed a new read to a vacant aligner
						PatternSourcePerThread *ps = (*patsrcs_)[i];
						// Get a new read
						ps->nextReadPair();
						if(ps->patid() < qUpto_ && !ps->empty()) {
							if(ps->paired()) {
								// Read currently in buffer is paired-end
								(*alignersPE_)[i]->setQuery(ps);
								seOrPe_[i] = false; // false -> paired
							} else {
								// Read currently in buffer is single-end
								(*alignersSE_)[i]->setQuery(ps);
								seOrPe_[i] = true; // true = unpaired
							}
							done = false;
						} else {
							// No more reads; if done == true, it remains
							// true
						}
					}
				}
				first = false;
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
template<typename TRangeSource>
class UnpairedAlignerV2 : public Aligner {
	typedef RangeSourceDriver<TRangeSource> TDriver;
public:
	UnpairedAlignerV2(
		EbwtSearchParams<String<Dna> >* params,
		TDriver* driver,
		RangeChaser<String<Dna> >* rchase,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		vector<String<Dna5> >& os,
		const BitPairReference* refs,
		bool rangeMode,
		bool verbose,
		bool quiet,
		int maxBts,
		ChunkPool *pool,
		int *btCnt = NULL,
		AlignerMetrics *metrics = NULL) :
		Aligner(true, rangeMode),
		refs_(refs),
		doneFirst_(true),
		firstIsFw_(true),
		chase_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		params_(params),
		rchase_(rchase),
		driver_(driver),
		verbose_(verbose),
		quiet_(quiet),
		maxBts_(maxBts),
		pool_(pool),
		btCnt_(btCnt),
		metrics_(metrics)
	{
		assert(pool_   != NULL);
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
		assert(driver_ != NULL);
	}

	virtual ~UnpairedAlignerV2() {
		delete driver_;  driver_  = NULL;
		delete params_;  params_  = NULL;
		delete rchase_;  rchase_  = NULL;
		delete[] btCnt_; btCnt_   = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		if(metrics_ != NULL) {
			metrics_->nextRead(patsrc->bufa().patFw);
		}
		pool_->reset(&patsrc->bufa().name, patsrc->patid());
		if(patsrc->bufa().length() < 4) {
			if(!quiet_) {
				cerr << "Warning: Skipping read " << patsrc->bufa().name
				     << " because it is less than 4 characters long" << endl;
			}
			this->done = true;
			sinkPt_->finishRead(*patsrc_, true, true);
			return;
		}
		driver_->setQuery(patsrc, NULL);
		this->done = driver_->done;
		doneFirst_ = false;
		// Reset #-backtrack countdown
		if(btCnt_ != NULL) *btCnt_ = maxBts_;
		if(sinkPt_->setHits(patsrc->bufa().hitset)) {
			this->done = true;
			sinkPt_->finishRead(*patsrc_, true, true);
		}
		// Grab a bit from the pseudo-random seed to determine whether
		// to start with forward or reverse complement
		firstIsFw_ = ((patsrc->bufa().seed & 0x10) == 0);
		chase_ = false;
	}

	/**
	 * Helper for reporting an alignment.
	 */
	inline bool report(const Range& ra,
	                   uint32_t first,
	                   uint32_t second,
	                   uint32_t tlen)
	{
		bool ebwtFw = ra.ebwt->fw();
		params_->setFw(ra.fw);
		assert_eq(bufa_->color, color);
		return params_->reportHit(
				ra.fw ? (ebwtFw? bufa_->patFw    : bufa_->patFwRev) :
				        (ebwtFw? bufa_->patRc    : bufa_->patRcRev),
				ra.fw ? (ebwtFw? &bufa_->qual    : &bufa_->qualRev) :
				        (ebwtFw? &bufa_->qualRev : &bufa_->qual),
				&bufa_->name,
				bufa_->color,
				colorExEnds,
				snpPhred,
				refs_,
				ra.ebwt->rmap(),
				ebwtFw,
				ra.mms,                   // mismatch positions
				ra.refcs,                 // reference characters for mms
				ra.numMms,                // # mismatches
				make_pair(first, second), // position
				make_pair(0, 0),          // (bogus) mate position
				true,                     // (bogus) mate orientation
				0,                        // (bogus) mate length
				make_pair(ra.top, ra.bot),// arrows
				tlen,                     // textlen
				alen_,                    // qlen
				ra.stratum,               // alignment stratum
				ra.cost,                  // cost, including qual penalty
				ra.bot - ra.top - 1,      // # other hits
				patsrc_->patid(),         // pattern id
				bufa_->seed,              // pseudo-random seed
				0);                       // mate (0 = unpaired)
	}

	/**
	 * Advance the aligner.  Return true iff we're
	 * done with this read.
	 */
	virtual bool advance() {
		assert(!this->done);
		if(chase_) {
			assert(!rangeMode_);
			assert(driver_->foundRange);
			assert(!sinkPt_->irrelevantCost(driver_->range().cost));
			if(!rchase_->foundOff() && !rchase_->done) {
				rchase_->advance();
				return false;
			}
			if(rchase_->foundOff()) {
				this->done = report(driver_->range(), rchase_->off().first,
				                    rchase_->off().second, rchase_->tlen());
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				chase_ = false;
				driver_->foundRange = false;
				this->done = driver_->done;
			}
		}
		// Still advancing a
		if(!this->done && !chase_) {
			assert(!driver_->done || driver_->foundRange);
			if(driver_->foundRange) {
				const Range& ra = driver_->range();
				assert(!sinkPt_->irrelevantCost(ra.cost));
				assert(ra.repOk());
				if(rangeMode_) {
					this->done = report(ra, ra.top, ra.bot, 0);
					driver_->foundRange = false;
				} else {
					rchase_->setTopBot(ra.top, ra.bot, alen_, rand_, ra.ebwt);
					if(rchase_->foundOff()) {
						this->done = report(
								ra, rchase_->off().first,
								rchase_->off().second, rchase_->tlen());
						rchase_->reset();
					}
					if(!rchase_->done && !sinkPt_->irrelevantCost(ra.cost)) {
						// Keep chasing this range
						chase_ = true;
					} else {
						driver_->foundRange = false;
					}
				}
			} else {
				this->done = sinkPt_->irrelevantCost(driver_->minCost);
				if(!this->done) {
					driver_->advance(ADV_COST_CHANGES);
				} else {
					// No longer necessarily true with chain input
					//assert(!sinkPt_->spanStrata());
				}
			}
			if(driver_->done && !driver_->foundRange && !chase_) {
				this->done = true;
			}
		}
		if(this->done) {
			sinkPt_->finishRead(*patsrc_, true, true);
		}
		return this->done;
	}

protected:

	// Reference sequences (needed for colorspace decoding)
	const BitPairReference* refs_;

	// Progress state
	bool doneFirst_;
	bool firstIsFw_;
	bool chase_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> >* params_;

	// State for getting alignments from ranges statefully
	RangeChaser<String<Dna> >* rchase_;

	// Range-finding state
	TDriver* driver_;

	bool verbose_; // be talkative
	bool quiet_; // don't print informational/warning info

	const int maxBts_;
	ChunkPool *pool_;
	int *btCnt_;
	AlignerMetrics *metrics_;
};

/**
 * An aligner for finding paired alignments while operating entirely
 * within the Burrows-Wheeler domain.
 */
template<typename TRangeSource>
class PairedBWAlignerV1 : public Aligner {

	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef std::vector<U32Pair> U32PairVec;
	typedef std::vector<Range> TRangeVec;
	typedef RangeSourceDriver<TRangeSource> TDriver;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	PairedBWAlignerV1(
		EbwtSearchParams<String<Dna> >* params,
		TDriver* driver1Fw, TDriver* driver1Rc,
		TDriver* driver2Fw, TDriver* driver2Rc,
		RefAligner<String<Dna5> >* refAligner,
		RangeChaser<String<Dna> >* rchase,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		bool fw1, bool fw2,
		uint32_t minInsert,
		uint32_t maxInsert,
		bool dontReconcile,
		uint32_t symCeiling,
		uint32_t mixedThresh,
		uint32_t mixedAttemptLim,
		const BitPairReference* refs,
		bool rangeMode,
		bool verbose,
		bool quiet,
		int maxBts,
		ChunkPool *pool,
		int *btCnt) :
		Aligner(true, rangeMode),
		refs_(refs),
		patsrc_(NULL), qlen1_(0), qlen2_(0), doneFw_(true),
		doneFwFirst_(true),
		chase1Fw_(false), chase1Rc_(false),
		chase2Fw_(false), chase2Rc_(false),
		delayedChase1Fw_(false), delayedChase1Rc_(false),
		delayedChase2Fw_(false), delayedChase2Rc_(false),
		refAligner_(refAligner),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		params_(params),
		minInsert_(minInsert),
		maxInsert_(maxInsert),
		dontReconcile_(dontReconcile),
		symCeiling_(symCeiling),
		mixedThresh_(mixedThresh),
		mixedAttemptLim_(mixedAttemptLim),
		mixedAttempts_(0),
		fw1_(fw1), fw2_(fw2),
		rchase_(rchase),
		verbose_(verbose),
		quiet_(quiet),
		maxBts_(maxBts),
		pool_(pool),
		btCnt_(btCnt),
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
		offsLarr_fw_     (fw1_ ? offs1FwArr_      : offs1RcArr_),
		offsRarr_fw_     (fw2_ ? offs2FwArr_      : offs2RcArr_),
		rangesLarr_fw_   (fw1_ ? ranges1FwArr_    : ranges1RcArr_),
		rangesRarr_fw_   (fw2_ ? ranges2FwArr_    : ranges2RcArr_),
		offsLsz_fw_      (fw1_ ? offs1FwSz_       : offs1RcSz_),
		offsRsz_fw_      (fw2_ ? offs2FwSz_       : offs2RcSz_),

		chaseL_rc_       (fw2_ ? chase2Rc_        : chase2Fw_),
		chaseR_rc_       (fw1_ ? chase1Rc_        : chase1Fw_),
		delayedchaseL_rc_(fw2_ ? delayedChase2Rc_ : delayedChase2Fw_),
		delayedchaseR_rc_(fw1_ ? delayedChase1Rc_ : delayedChase1Fw_),
		drL_rc_          (fw2_ ? *driver2Rc_      : *driver2Fw_),
		drR_rc_          (fw1_ ? *driver1Rc_      : *driver1Fw_),
		offsLarr_rc_     (fw2_ ? offs2RcArr_      : offs2FwArr_),
		offsRarr_rc_     (fw1_ ? offs1RcArr_      : offs1FwArr_),
		rangesLarr_rc_   (fw2_ ? ranges2RcArr_    : ranges2FwArr_),
		rangesRarr_rc_   (fw1_ ? ranges1RcArr_    : ranges1FwArr_),
		offsLsz_rc_      (fw2_ ? offs2RcSz_       : offs2FwSz_),
		offsRsz_rc_      (fw1_ ? offs1RcSz_       : offs1FwSz_),

		chaseL_       (&chaseL_fw_),
		chaseR_       (&chaseR_fw_),
		delayedchaseL_(&delayedchaseL_fw_),
		delayedchaseR_(&delayedchaseR_fw_),
		drL_          (&drL_fw_),
		drR_          (&drR_fw_),
		offsLarr_     (offsLarr_fw_),
		offsRarr_     (offsRarr_fw_),
		rangesLarr_   (rangesLarr_fw_),
		rangesRarr_   (rangesRarr_fw_),
		offsLsz_      (&offsLsz_fw_),
		offsRsz_      (&offsRsz_fw_),
		donePair_     (&doneFw_),
		fwL_(fw1),
		fwR_(fw2),
		verbose2_(false)
	{
		assert(pool_      != NULL);
		assert(sinkPt_    != NULL);
		assert(params_    != NULL);
		assert(driver1Fw_ != NULL);
		assert(driver1Rc_ != NULL);
		assert(driver2Fw_ != NULL);
		assert(driver2Rc_ != NULL);
	}

	virtual ~PairedBWAlignerV1() {
		delete driver1Fw_; driver1Fw_ = NULL;
		delete driver1Rc_; driver1Rc_ = NULL;
		delete driver2Fw_; driver2Fw_ = NULL;
		delete driver2Rc_; driver2Rc_ = NULL;
		delete params_;    params_    = NULL;
		delete rchase_;    rchase_    = NULL;
		delete[] btCnt_;   btCnt_     = NULL;
		delete refAligner_; refAligner_ = NULL;
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
		patsrc_ = patsrc;
		pool_->reset(&patsrc->bufa().name, patsrc->patid());
		if(patsrc->bufa().length() < 4 || patsrc->bufb().length() < 4) {
			if(!quiet_) {
				cerr << "Warning: Skipping pair " << patsrc->bufa().name
					 << " because a mate is less than 4 characters long" << endl;
			}
			this->done = true;
			sinkPt_->finishRead(*patsrc_, true, true);
			return;
		}
		driver1Fw_->setQuery(patsrc, NULL);
		driver1Rc_->setQuery(patsrc, NULL);
		driver2Fw_->setQuery(patsrc, NULL);
		driver2Rc_->setQuery(patsrc, NULL);
		qlen1_ = patsrc_->bufa().length();
		qlen2_ = patsrc_->bufb().length();
		if(btCnt_ != NULL) (*btCnt_) = maxBts_;
		// Neither orientation is done
		doneFw_   = false;
		doneFwFirst_ = true;
		this->done   = false;
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
		for(size_t i = 0; i < 32; i++) {
			offs1FwArr_[i].clear();   offs1RcArr_[i].clear();
			offs2FwArr_[i].clear();   offs2RcArr_[i].clear();
			ranges1FwArr_[i].clear(); ranges1RcArr_[i].clear();
			ranges2FwArr_[i].clear(); ranges2RcArr_[i].clear();
		}
		offs1FwSz_ = offs1RcSz_ = offs2FwSz_ = offs2RcSz_ = 0;
		chaseL_        = &chaseL_fw_;
		chaseR_        = &chaseR_fw_;
		delayedchaseL_ = &delayedchaseL_fw_;
		delayedchaseR_ = &delayedchaseR_fw_;
		drL_           = &drL_fw_;
		drR_           = &drR_fw_;
		offsLarr_      = offsLarr_fw_;
		offsRarr_      = offsRarr_fw_;
		rangesLarr_    = rangesLarr_fw_;
		rangesRarr_    = rangesRarr_fw_;
		offsLsz_       = &offsLsz_fw_;
		offsRsz_       = &offsRsz_fw_;
		donePair_      = &doneFw_;
		fwL_           = fw1_;
		fwR_           = fw2_;
		mixedAttempts_ = 0;
		pairs_fw_.clear();
		pairs_rc_.clear();
#ifndef NDEBUG
		allTopsL_fw_.clear();
		allTopsR_fw_.clear();
		allTopsL_rc_.clear();
		allTopsR_rc_.clear();
#endif
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
		assert(!this->done);
		if(doneFw_ && doneFwFirst_) {
			if(verbose2_) cout << "--" << endl;
			chaseL_        = &chaseL_rc_;
			chaseR_        = &chaseR_rc_;
			delayedchaseL_ = &delayedchaseL_rc_;
			delayedchaseR_ = &delayedchaseR_rc_;
			drL_           = &drL_rc_;
			drR_           = &drR_rc_;
			offsLarr_      = offsLarr_rc_;
			offsRarr_      = offsRarr_rc_;
			rangesLarr_    = rangesLarr_rc_;
			rangesRarr_    = rangesRarr_rc_;
			offsLsz_       = &offsLsz_rc_;
			offsRsz_       = &offsRsz_rc_;
			donePair_      = &this->done;
			fwL_           = !fw2_;
			fwR_           = !fw1_;
			doneFwFirst_   = false;
			mixedAttempts_ = 0;
		}
		bool chasing = *chaseL_ || *chaseR_;
		if(chasing && !rchase_->foundOff() && !rchase_->done) {
			rchase_->advance();
			return false;
		}
		advanceOrientation(!doneFw_);
		if(this->done) {
			if(verbose2_) cout << "----" << endl;
			sinkPt_->finishRead(*patsrc_, true, true);
		}
		return this->done;
	}

protected:

	/**
	 * Helper for reporting a pair of alignments.  As of now, we report
	 * a paired alignment by reporting two consecutive alignments, one
	 * for each mate.
	 */
	bool report(const Range& rL, // range for upstream mate
	            const Range& rR, // range for downstream mate
	            uint32_t first,  // ref idx
	            uint32_t upstreamOff, // offset for upstream mate
	            uint32_t dnstreamOff, // offset for downstream mate
	            uint32_t tlen, // length of ref
	            bool pairFw,   // whether the pair is being mapped to fw strand
	            bool ebwtFwL,
	            bool ebwtFwR,
	            const ReferenceMap* rmap)
	{
		assert_lt(upstreamOff, dnstreamOff);
		uint32_t spreadL = rL.bot - rL.top;
		uint32_t spreadR = rR.bot - rR.top;
		uint32_t oms = min(spreadL, spreadR) - 1;
		ReadBuf* bufL = pairFw ? bufa_ : bufb_;
		ReadBuf* bufR = pairFw ? bufb_ : bufa_;
		uint32_t lenL = pairFw ? alen_ : blen_;
		uint32_t lenR = pairFw ? blen_ : alen_;
		bool ret;
		assert(!params_->sink().exceededOverThresh());
		params_->setFw(rL.fw);
		assert_eq(bufL->color, color);
		// Print upstream mate first
		ret = params_->reportHit(
				rL.fw ? (ebwtFwL?  bufL->patFw  :  bufL->patFwRev) :
					    (ebwtFwL?  bufL->patRc  :  bufL->patRcRev),
				rL.fw ? (ebwtFwL? &bufL->qual    : &bufL->qualRev) :
				        (ebwtFwL? &bufL->qualRev : &bufL->qual),
				&bufL->name,
				bufL->color,
				colorExEnds,
				snpPhred,
				refs_,
				rmap,
				ebwtFwL,
				rL.mms,                       // mismatch positions
				rL.refcs,                     // reference characters for mms
				rL.numMms,                    // # mismatches
				make_pair(first, upstreamOff),// position
				make_pair(first, dnstreamOff),// mate position
				rR.fw,                        // mate orientation
				lenR,                         // mate length
				make_pair(rL.top, rL.bot),    // arrows
				tlen,                         // textlen
				lenL,                         // qlen
				rL.stratum,                   // alignment stratum
				rL.cost,                      // cost, including quality penalty
				oms,                          // # other hits
				bufL->patid,
				bufL->seed,
				pairFw ? 1 : 2);
		if(ret) {
			return true; // can happen when -m is set
		}
		params_->setFw(rR.fw);
		assert_eq(bufR->color, color);
		ret = params_->reportHit(
				rR.fw ? (ebwtFwR?  bufR->patFw  :  bufR->patFwRev) :
					    (ebwtFwR?  bufR->patRc  :  bufR->patRcRev),
				rR.fw ? (ebwtFwR? &bufR->qual    : &bufR->qualRev) :
				        (ebwtFwR? &bufR->qualRev : &bufR->qual),
				&bufR->name,
				bufR->color,
				colorExEnds,
				snpPhred,
				refs_,
				rmap,
				ebwtFwR,
				rR.mms,                       // mismatch positions
				rR.refcs,                     // reference characters for mms
				rR.numMms,                    // # mismatches
				make_pair(first, dnstreamOff),// position
				make_pair(first, upstreamOff),// mate position
				rL.fw,                        // mate orientation
				lenL,                         // mate length
				make_pair(rR.top, rR.bot),    // arrows
				tlen,                         // textlen
				lenR,                         // qlen
				rR.stratum,                   // alignment stratum
				rR.cost,                      // cost, including quality penalty
				oms,                          // # other hits
				bufR->patid,
				bufR->seed,
				pairFw ? 2 : 1);
		return ret;
	}

	bool report(const Range& rL, // range for upstream mate
	            const Range& rR, // range for downstream mate
	            uint32_t first,  // ref idx
	            uint32_t upstreamOff, // offset for upstream mate
	            uint32_t dnstreamOff, // offset for downstream mate
	            uint32_t tlen, // length of ref
	            bool pairFw,   // whether the pair is being mapped to fw strand
	            const ReferenceMap* rmap)
	{
		return report(rL, rR, first, upstreamOff,
		              dnstreamOff, tlen,
		              pairFw, rL.ebwt->fw(), rR.ebwt->fw(), rmap);
	}

	/**
	 * Given a vector of reference positions where one of the two mates
	 * (the "anchor" mate) has aligned, look directly at the reference
	 * sequence for instances where the other mate (the "outstanding"
	 * mate) aligns such that mating constraint is satisfied.
	 *
	 * This function picks up to 'pick' anchors at random from the
	 * 'offs' array.  It returns the number that it actually picked.
	 */
	bool resolveOutstandingInRef(const bool off1,
	                             const U32Pair& off,
	                             const uint32_t tlen,
	                             const Range& range)
	{
		assert(refs_->loaded());
		assert_lt(off.first, refs_->numRefs());
		// If matchRight is true, then we're trying to align the other
		// mate to the right of the already-aligned mate.  Otherwise,
		// to the left.
		bool matchRight = (off1 ? !doneFw_ : doneFw_);
		// Sequence and quals for mate to be matched
		bool fw = off1 ? fw2_ : fw1_; // whether outstanding mate is fw/rc
		if(doneFw_) fw = !fw;
		// 'seq' gets sequence of outstanding mate w/r/t the forward
		// reference strand
		const String<Dna5>& seq  = fw ? (off1 ? patsrc_->bufb().patFw   :
		                                        patsrc_->bufa().patFw)  :
		                                (off1 ? patsrc_->bufb().patRc   :
		                                        patsrc_->bufa().patRc);
		// 'seq' gets qualities of outstanding mate w/r/t the forward
		// reference strand
		const String<char>& qual = fw ? (off1 ? patsrc_->bufb().qual  :
		                                        patsrc_->bufa().qual) :
		                                (off1 ? patsrc_->bufb().qualRev  :
		                                        patsrc_->bufa().qualRev);
		uint32_t qlen = seqan::length(seq);  // length of outstanding mate
		uint32_t alen = (off1 ? patsrc_->bufa().length() :
		                        patsrc_->bufb().length());
		int minins = minInsert_;
		int maxins = maxInsert_;
		if(fw1_) {
			minins = max<int>(0, minins - patsrc_->bufa().trimmed5);
			maxins = max<int>(0, maxins - patsrc_->bufa().trimmed5);
		} else {
			minins = max<int>(0, minins - patsrc_->bufa().trimmed3);
			maxins = max<int>(0, maxins - patsrc_->bufa().trimmed3);
		}
		if(fw2_) {
			minins = max<int>(0, minins - patsrc_->bufb().trimmed3);
			maxins = max<int>(0, maxins - patsrc_->bufb().trimmed3);
		} else {
			minins = max<int>(0, minins - patsrc_->bufb().trimmed5);
			maxins = max<int>(0, maxins - patsrc_->bufb().trimmed5);
		}
		assert_geq(minins, 0);
		assert_geq(maxins, 0);
		// Don't even try if either of the mates is longer than the
		// maximum insert size.
		if((uint32_t)maxins <= max(qlen, alen)) {
			return false;
		}
		const uint32_t tidx = off.first;
		const uint32_t toff = off.second;
		// Set begin/end to be a range of all reference
		// positions that are legally permitted to be involved in
		// the alignment of the outstanding mate.  It's up to the
		// callee to worry about how to scan these positions.
		uint32_t begin, end;
		assert_geq(maxins, minins);
		uint32_t insDiff = maxins - minins;
		if(matchRight) {
			end = toff + maxins;
			begin = toff + 1;
			if(qlen < alen) begin += alen-qlen;
			if(end > insDiff + qlen) {
				begin = max<uint32_t>(begin, end - insDiff - qlen);
			}
			end = min<uint32_t>(refs_->approxLen(tidx), end);
			begin = min<uint32_t>(refs_->approxLen(tidx), begin);
		} else {
			if(toff + alen < (uint32_t)maxins) {
				begin = 0;
			} else {
				begin = toff + alen - maxins;
			}
			uint32_t mi = min<uint32_t>(alen, qlen);
			end = toff + mi - 1;
			end = min<uint32_t>(end, toff + alen - minins + qlen - 1);
			if(toff + alen + qlen < (uint32_t)minins + 1) end = 0;
		}
		// Check if there's not enough space in the range to fit an
		// alignment for the outstanding mate.
		if(end - begin < qlen) return false;
		std::vector<Range> ranges;
		std::vector<uint32_t> offs;
		refAligner_->find(1, tidx, refs_, seq, qual, begin, end, ranges,
		                  offs, doneFw_ ? &pairs_rc_ : &pairs_fw_,
		                  toff, fw);
		assert_eq(ranges.size(), offs.size());
		for(size_t i = 0; i < ranges.size(); i++) {
			Range& r = ranges[i];
			r.fw = fw;
			r.cost |= (r.stratum << 14);
			r.mate1 = !off1;
			const uint32_t result = offs[i];
			// Just copy the known range's top and bot for now
			r.top = range.top;
			r.bot = range.bot;
			bool ebwtLFw = matchRight ? range.ebwt->fw() : true;
			bool ebwtRFw = matchRight ? true : range.ebwt->fw();
			if(report(
				matchRight ? range : r, // range for upstream mate
				matchRight ? r : range, // range for downstream mate
				tidx,                   // ref idx
				matchRight ? toff : result, // upstream offset
				matchRight ? result : toff, // downstream offset
				tlen,       // length of ref
				!doneFw_,   // whether the pair is being mapped to fw strand
				ebwtLFw,
				ebwtRFw,
				range.ebwt->rmap())) return true;
		}
		return false;
	}

	/**
	 * Advance paired-end alignment.
	 */
	void advanceOrientation(bool pairFw) {
		assert(!this->done);
		assert(!*donePair_);
		assert(!*chaseL_ || !*chaseR_);
		if(*chaseL_) {
			assert(!rangeMode_);
			assert(!*delayedchaseL_);
			assert(drL_->foundRange);
			assert(rchase_->foundOff() || rchase_->done);
			if(rchase_->foundOff()) {
				// Resolve this against the reference loci
				// determined for the other mate
				const bool overThresh = (*offsLsz_ + *offsRsz_) > mixedThresh_;
				if(!this->done && (overThresh || dontReconcile_)) {
					// Because the total size of both ranges exceeds
					// our threshold, we're now operating in "mixed
					// mode"
					const Range& r = drL_->range();
					assert(r.repOk());
					if(verbose_) cout << "Making an attempt to find the outstanding mate" << endl;
					this->done = resolveOutstandingInRef(
							pairFw, rchase_->off(),
					        r.ebwt->_plen[rchase_->off().first], r);
					if(++mixedAttempts_ > mixedAttemptLim_) {
						// Give up on this pair
						*donePair_ = true;
						return;
					}
				}
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				*chaseL_ = false;
				drL_->foundRange = false;
				if(verbose_) cout << "Done with chase for first mate" << endl;
				if(*delayedchaseR_) {
					// Start chasing the delayed range
					if(verbose_) cout << "Resuming delayed chase for second mate" << endl;
					assert(drR_->foundRange);
					const Range& r = drR_->range();
					assert(r.repOk());
					uint32_t top = r.top;
					uint32_t bot = r.bot;
					uint32_t qlen = doneFw_? qlen1_ : qlen2_;
					rchase_->setTopBot(top, bot, qlen, rand_, r.ebwt);
					*chaseR_ = true;
					*delayedchaseR_ = false;
				}
			}
		} else if(*chaseR_) {
			assert(!rangeMode_);
			assert(!*delayedchaseR_);
			assert(drR_->foundRange);
			assert(rchase_->foundOff() || rchase_->done);
			if(rchase_->foundOff()) {
				// Resolve this against the reference loci
				// determined for the other mate
				const bool overThresh = (*offsLsz_ + *offsRsz_) > mixedThresh_;
				if(!this->done && (overThresh || dontReconcile_)) {
					// Because the total size of both ranges exceeds
					// our threshold, we're now operating in "mixed
					// mode"
					const Range& r = drR_->range();
					if(verbose_) cout << "Making an attempt to find the outstanding mate" << endl;
					this->done = resolveOutstandingInRef(
							!pairFw, rchase_->off(),
					        r.ebwt->_plen[rchase_->off().first], r);
					if(++mixedAttempts_ > mixedAttemptLim_) {
						// Give up on this pair
						*donePair_ = true;
						return;
					}
				}
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				*chaseR_ = false;
				drR_->foundRange = false;
				if(verbose_) cout << "Done with chase for second mate" << endl;
				if(*delayedchaseL_) {
					// Start chasing the delayed range
					if(verbose_) cout << "Resuming delayed chase for first mate" << endl;
					assert(drL_->foundRange);
					const Range& r = drL_->range();
					assert(r.repOk());
					uint32_t top = r.top;
					uint32_t bot = r.bot;
					uint32_t qlen = doneFw_? qlen2_ : qlen1_;
					rchase_->setTopBot(top, bot, qlen, rand_, r.ebwt);
					*chaseL_ = true;
					*delayedchaseL_ = false;
				}
			}
		}
		if(!this->done && !*donePair_ && !*chaseL_ && !*chaseR_) {
			// Search for more ranges for whichever mate currently has
			// fewer candidate alignments
			if((*offsLsz_ < *offsRsz_ || drR_->done) && !drL_->done) {
				// If there are no more ranges for the other mate and
				// there are no candidate alignments either, then we're
				// not going to find a paired alignment in this
				// orientation.
				if(drR_->done && *offsRsz_ == 0) {
					// Give up on this orientation
					if(verbose_) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 1" << endl;
					*donePair_ = true;
					if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << endl;
					return;
				}
				assert(!*delayedchaseL_);
				if(!drL_->foundRange) drL_->advance(ADV_FOUND_RANGE);
				if(drL_->foundRange) {
#ifndef NDEBUG
					{
						std::set<int64_t>& s = (pairFw ? allTopsL_fw_ : allTopsL_rc_);
						int64_t t = drL_->range().top + 1; // add 1 to avoid 0
						if(!drL_->range().ebwt->fw()) t = -t; // invert for bw index
						assert(s.find(t) == s.end());
						s.insert(t);
					}
#endif
					// Add the size of this range to the total for this mate
					*offsLsz_ += (drL_->range().bot - drL_->range().top);
					if(*offsRsz_ == 0 && (!dontReconcile_ || *offsLsz_ > 3)) {
						// Delay chasing this range; we delay to avoid
						// needlessly chasing rows in this range when
						// the other mate doesn't end up aligning
						// anywhere
						if(verbose_) cout << "Delaying a chase for first mate" << endl;
						*delayedchaseL_ = true;
					} else {
						// Start chasing this range
						if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << " " << drL_->range().top << endl;
						if(verbose_) cout << "Chasing a range for first mate" << endl;
						if(*offsLsz_ > symCeiling_ && *offsRsz_ > symCeiling_) {
							// Too many candidates for both mates; abort
							// without any more searching
							*donePair_ = true;
							return;
						}
						// If this is the first range for both mates,
						// choose the smaller range to chase down first
						if(*delayedchaseR_ && (*offsRsz_ < *offsLsz_)) {
							assert(drR_->foundRange);
							*delayedchaseR_ = false;
							*delayedchaseL_ = true;
							*chaseR_ = true;
							const Range& r = drR_->range();
							assert(r.repOk());
							uint32_t qlen = doneFw_? qlen1_ : qlen2_;
							rchase_->setTopBot(r.top, r.bot, qlen, rand_, r.ebwt);
						} else {
							// Use Burrows-Wheeler for this pair (as
							// usual)
							*chaseL_ = true;
							const Range& r = drL_->range();
							uint32_t qlen = doneFw_? qlen2_ : qlen1_;
							rchase_->setTopBot(r.top, r.bot, qlen, rand_, r.ebwt);
						}
					}
				}
			} else if(!drR_->done) {
				// If there are no more ranges for the other mate and
				// there are no candidate alignments either, then we're
				// not going to find a paired alignment in this
				// orientation.
				if(drL_->done && *offsLsz_ == 0) {
					// Give up on this orientation
					if(verbose_) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 2" << endl;
					if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << endl;
					*donePair_ = true;
					return;
				}
				assert(!*delayedchaseR_);
				if(!drR_->foundRange) drR_->advance(ADV_FOUND_RANGE);
				if(drR_->foundRange) {
#ifndef NDEBUG
					{
						std::set<int64_t>& s = (pairFw ? allTopsR_fw_ : allTopsR_rc_);
						int64_t t = drR_->range().top + 1; // add 1 to avoid 0
						if(!drR_->range().ebwt->fw()) t = -t; // invert for bw index
						assert(s.find(t) == s.end());
						s.insert(t);
					}
#endif
					// Add the size of this range to the total for this mate
					*offsRsz_ += (drR_->range().bot - drR_->range().top);
					if(*offsLsz_ == 0 && (!dontReconcile_ || *offsRsz_ > 3)) {
						// Delay chasing this range; we delay to avoid
						// needlessly chasing rows in this range when
						// the other mate doesn't end up aligning
						// anywhere
						if(verbose_) cout << "Delaying a chase for second mate" << endl;
						*delayedchaseR_ = true;
					} else {
						// Start chasing this range
						if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << " " << drR_->range().top << endl;
						if(verbose_) cout << "Chasing a range for second mate" << endl;
						if(*offsLsz_ > symCeiling_ && *offsRsz_ > symCeiling_) {
							// Too many candidates for both mates; abort
							// without any more searching
							*donePair_ = true;
							return;
						}
						// If this is the first range for both mates,
						// choose the smaller range to chase down first
						if(*delayedchaseL_ && *offsLsz_ < *offsRsz_) {
							assert(drL_->foundRange);
							*delayedchaseL_ = false;
							*delayedchaseR_ = true;
							*chaseL_ = true;
							const Range& r = drL_->range();
							assert(r.repOk());
							uint32_t qlen = doneFw_? qlen2_ : qlen1_;
							rchase_->setTopBot(r.top, r.bot, qlen, rand_, r.ebwt);
						} else {
							// Use Burrows-Wheeler for this pair (as
							// usual)
							*chaseR_ = true;
							const Range& r = drR_->range();
							assert(r.repOk());
							uint32_t qlen = doneFw_? qlen1_ : qlen2_;
							rchase_->setTopBot(r.top, r.bot, qlen, rand_, r.ebwt);
						}
					}
				}
			} else {
				// Finished processing ranges for both mates
				assert(drL_->done && drR_->done);
				*donePair_ = true;
			}
		}
	}

	const BitPairReference* refs_;

	PatternSourcePerThread *patsrc_;
	uint32_t qlen1_;
	uint32_t qlen2_;

	// Progress state
	bool doneFw_;   // finished with forward orientation of both mates?
	bool doneFwFirst_;

	bool chase1Fw_;
	bool chase1Rc_;
	bool chase2Fw_;
	bool chase2Rc_;

	bool delayedChase1Fw_;
	bool delayedChase1Rc_;
	bool delayedChase2Fw_;
	bool delayedChase2Rc_;

	// For searching for outstanding mates
	RefAligner<String<Dna5> >* refAligner_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> >* params_;

	// Paired-end boundaries
	const uint32_t minInsert_;
	const uint32_t maxInsert_;

	// Don't attempt pairwise all-versus-all style of mate
	// reconciliation; just rely on mixed mode
	const bool dontReconcile_;

	// If both mates in a given orientation align >= symCeiling times,
	// then immediately give up
	const uint32_t symCeiling_;

	// If the total number of alignments for both mates in a given
	// orientation exceeds mixedThresh, then switch to mixed mode
	const uint32_t mixedThresh_;
	const uint32_t mixedAttemptLim_;
	uint32_t mixedAttempts_;

	// Orientation of upstream/downstream mates when aligning to
	// forward strand
	const bool fw1_;
	const bool fw2_;

	// State for getting alignments from ranges statefully
	RangeChaser<String<Dna> >* rchase_;

	// true -> be talkative
	bool verbose_;
	// true -> suppress warnings
	bool quiet_;

	int maxBts_;
	ChunkPool *pool_;
	int *btCnt_;

	// Range-finding state for first mate
	TDriver*      driver1Fw_;
	TDriver*      driver1Rc_;
	U32PairVec    offs1FwArr_[32];
	TRangeVec     ranges1FwArr_[32];
	uint32_t      offs1FwSz_; // total size of all ranges found in this category
	U32PairVec    offs1RcArr_[32];
	TRangeVec     ranges1RcArr_[32];
	uint32_t      offs1RcSz_; // total size of all ranges found in this category

	// Range-finding state for second mate
	TDriver*      driver2Fw_;
	TDriver*      driver2Rc_;
	U32PairVec    offs2FwArr_[32];
	TRangeVec     ranges2FwArr_[32];
	uint32_t      offs2FwSz_; // total size of all ranges found in this category
	U32PairVec    offs2RcArr_[32];
	TRangeVec     ranges2RcArr_[32];
	uint32_t      offs2RcSz_; // total size of all ranges found in this category

	bool&       chaseL_fw_;
	bool&       chaseR_fw_;
	bool&       delayedchaseL_fw_;
	bool&       delayedchaseR_fw_;
	TDriver&    drL_fw_;
	TDriver&    drR_fw_;
	U32PairVec* offsLarr_fw_;
	U32PairVec* offsRarr_fw_;
	TRangeVec*  rangesLarr_fw_;
	TRangeVec*  rangesRarr_fw_;
	uint32_t&   offsLsz_fw_;
	uint32_t&   offsRsz_fw_;

	bool&       chaseL_rc_;
	bool&       chaseR_rc_;
	bool&       delayedchaseL_rc_;
	bool&       delayedchaseR_rc_;
	TDriver&    drL_rc_;
	TDriver&    drR_rc_;
	U32PairVec* offsLarr_rc_;
	U32PairVec* offsRarr_rc_;
	TRangeVec*  rangesLarr_rc_;
	TRangeVec*  rangesRarr_rc_;
	uint32_t&   offsLsz_rc_;
	uint32_t&   offsRsz_rc_;

	bool*       chaseL_;
	bool*       chaseR_;
	bool*       delayedchaseL_;
	bool*       delayedchaseR_;
	TDriver*    drL_;
	TDriver*    drR_;
	U32PairVec* offsLarr_;
	U32PairVec* offsRarr_;
	TRangeVec*  rangesLarr_;
	TRangeVec*  rangesRarr_;
	uint32_t*   offsLsz_;
	uint32_t*   offsRsz_;
	bool*       donePair_;
	bool        fwL_;
	bool        fwR_;

	/// For keeping track of paired alignments that have already been
	/// found for the forward and reverse-comp pair orientations
	TSetPairs   pairs_fw_;
	TSetPairs   pairs_rc_;

#ifndef NDEBUG
	std::set<int64_t> allTopsL_fw_;
	std::set<int64_t> allTopsR_fw_;
	std::set<int64_t> allTopsL_rc_;
	std::set<int64_t> allTopsR_rc_;
#endif

	bool verbose2_;
};

/**
 * Helper struct that holds a Range together with the coordinates where it al
 */
struct RangeWithCoords {
	Range r;
	U32Pair h;
};

/**
 * An aligner for finding paired alignments while operating entirely
 * within the Burrows-Wheeler domain.
 */
template<typename TRangeSource>
class PairedBWAlignerV2 : public Aligner {

	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef std::vector<U32Pair> U32PairVec;
	typedef std::vector<Range> TRangeVec;
	typedef RangeSourceDriver<TRangeSource> TDriver;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	PairedBWAlignerV2(
		EbwtSearchParams<String<Dna> >* params,
		EbwtSearchParams<String<Dna> >* paramsSe1,
		EbwtSearchParams<String<Dna> >* paramsSe2,
		TDriver* driver,
		RefAligner<String<Dna5> >* refAligner,
		RangeChaser<String<Dna> >* rchase,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		HitSinkPerThread* sinkPtSe1,
		HitSinkPerThread* sinkPtSe2,
		bool fw1, bool fw2,
		uint32_t minInsert,
		uint32_t maxInsert,
		uint32_t mixedAttemptLim,
		const BitPairReference* refs,
		bool rangeMode,
		bool verbose,
		bool quiet,
		int maxBts,
		ChunkPool *pool,
		int *btCnt) :
		Aligner(true, rangeMode),
		refs_(refs),
		patsrc_(NULL),
		qlen1_(0), qlen2_(0),
		chase_(false),
		donePe_(false),
		doneSe1_(false),
		doneSe2_(false),
		refAligner_(refAligner),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		sinkPtSe1_(sinkPtSe1),
		sinkPtSe2_(sinkPtSe2),
		params_(params),
		paramsSe1_(paramsSe1),
		paramsSe2_(paramsSe2),
		minInsert_(minInsert),
		maxInsert_(maxInsert),
		mixedAttemptLim_(mixedAttemptLim),
		mixedAttempts_(0),
		fw1_(fw1), fw2_(fw2),
		rchase_(rchase),
		driver_(driver),
		pool_(pool),
		verbose_(verbose),
		quiet_(quiet),
		maxBts_(maxBts),
		btCnt_(btCnt)
	{
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
		assert(driver_ != NULL);
	}

	virtual ~PairedBWAlignerV2() {
		delete driver_; driver_ = NULL;
		delete params_; params_ = NULL;
		if(paramsSe1_ != NULL) {
			delete paramsSe1_; paramsSe1_ = NULL;
			delete paramsSe2_; paramsSe2_ = NULL;
		}
		delete rchase_; rchase_ = NULL;
		delete[] btCnt_; btCnt_ = NULL;
		delete refAligner_; refAligner_ = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
		if(sinkPtSe1_ != NULL) {
			sinkPtFactory_.destroy(sinkPtSe1_); sinkPtSe1_ = NULL;
			sinkPtFactory_.destroy(sinkPtSe2_); sinkPtSe2_ = NULL;
		}
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		assert(!patsrc->bufa().empty());
		Aligner::setQuery(patsrc); // set fields & random seed
		assert(!patsrc->bufb().empty());
		// Give all of the drivers pointers to the relevant read info
		patsrc_ = patsrc;
		pool_->reset(&patsrc->bufa().name, patsrc->patid());
		if(patsrc->bufa().length() < 4 || patsrc->bufb().length() < 4) {
			if(!quiet_) {
				cerr << "Warning: Skipping pair " << patsrc->bufa().name
				     << " because a mate is less than 4 characters long" << endl;
			}
			this->done = true;
			sinkPt_->finishRead(*patsrc_, true, true);
			return;
		}
		driver_->setQuery(patsrc, NULL);
		qlen1_ = patsrc_->bufa().length();
		qlen2_ = patsrc_->bufb().length();
		if(btCnt_ != NULL) (*btCnt_) = maxBts_;
		mixedAttempts_ = 0;
		// Neither orientation is done
		this->done = false;
		// No ranges are being chased yet
		chase_ = false;
		donePe_ = doneSe1_ = doneSe2_ = false;
		pairs_fw_.clear();
		pairs_rc_.clear();
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
		assert(!this->done);
		if(chase_) {
			assert(!rangeMode_); // chasing ranges
			if(!rchase_->foundOff() && !rchase_->done) {
				rchase_->advance();
				return false;
			}
			assert(rchase_->foundOff() || rchase_->done);
			if(rchase_->foundOff()) {
				const Range& r = driver_->range();
				assert(r.repOk());
				resolveOutstanding(
					rchase_->off(),
					r.ebwt->_plen[rchase_->off().first], r);
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				chase_ = false;
				this->done = driver_->done;
			}
		}

		if(!this->done && !chase_) {
			// Search for more ranges for whichever mate currently has
			// fewer candidate alignments
			if(!driver_->done) {
				if(!this->done) {
					//
					// Check whether any of the PE/SE possibilities
					// have become impossible due to the minCost
					//
					if(!donePe_) {
						assert(!this->done);
						donePe_ = sinkPt_->irrelevantCost(driver_->minCost);
						if(donePe_ && (!sinkPt_->empty() || sinkPtSe1_ == NULL)) {
							// Paired-end alignment(s) were found, no
							// more will be found, and no unpaired
							// alignments are requested, so stop
							this->done = true;
						}
						if(donePe_ && sinkPtSe1_ != NULL) {
							// Note: removeMate affects minCost
							if(doneSe1_) driver_->removeMate(1);
							if(doneSe2_) driver_->removeMate(2);
						}
					}
					if(!this->done && sinkPtSe1_ != NULL) {
						if(!doneSe1_) {
							doneSe1_ = sinkPtSe1_->irrelevantCost(driver_->minCost);
							if(doneSe1_ && donePe_) driver_->removeMate(1);
						}
						if(!doneSe2_) {
							doneSe2_ = sinkPtSe2_->irrelevantCost(driver_->minCost);
							if(doneSe2_ && donePe_) driver_->removeMate(2);
						}
						// Do Se1 again, because removing Se2 may have
						// nudged minCost over the threshold
						if(!doneSe1_) {
							doneSe1_ = sinkPtSe1_->irrelevantCost(driver_->minCost);
							if(doneSe1_ && donePe_) driver_->removeMate(1);
						}
						if(doneSe1_ && doneSe2_) assert(donePe_);
						this->done = donePe_ && doneSe1_ && doneSe2_;
					}

					if(!this->done) {
						if(sinkPtSe1_ != NULL) {
							assert(doneSe1_ || !sinkPtSe1_->irrelevantCost(driver_->minCost));
							assert(doneSe2_ || !sinkPtSe2_->irrelevantCost(driver_->minCost));
						}
						assert(donePe_ || !sinkPt_->irrelevantCost(driver_->minCost));
						driver_->advance(ADV_COST_CHANGES);
					}
				}
				if(driver_->foundRange) {
					// Use Burrows-Wheeler for this pair (as usual)
					chase_ = true;
					driver_->foundRange = false;
					const Range& r = driver_->range();
					assert(r.repOk());
					rchase_->setTopBot(r.top, r.bot,
					                   r.mate1 ? qlen1_ : qlen2_,
					                   rand_, r.ebwt);
				}
			} else {
				this->done = true;
			}
		}

		if(this->done) {
			bool reportedPe = (sinkPt_->finishRead(*patsrc_, true, true) > 0);
			if(sinkPtSe1_ != NULL) {
				sinkPtSe1_->finishRead(*patsrc_, !reportedPe, false);
				sinkPtSe2_->finishRead(*patsrc_, !reportedPe, false);
			}
		}
		return this->done;
	}

protected:

	/**
	 * Helper for reporting a pair of alignments.  As of now, we report
	 * a paired alignment by reporting two consecutive alignments, one
	 * for each mate.
	 */
	bool report(const Range& rL, // range for upstream mate
	            const Range& rR, // range for downstream mate
	            uint32_t first,  // ref idx
	            uint32_t upstreamOff, // offset for upstream mate
	            uint32_t dnstreamOff, // offset for downstream mate
	            uint32_t tlen, // length of ref
	            bool pairFw,   // whether the pair is being mapped to fw strand
	            bool ebwtFwL,
	            bool ebwtFwR,
	            const ReferenceMap *rmap)
	{
		assert_lt(upstreamOff, dnstreamOff);
		uint32_t spreadL = rL.bot - rL.top;
		uint32_t spreadR = rR.bot - rR.top;
		uint32_t oms = min(spreadL, spreadR) - 1;
		ReadBuf* bufL = pairFw ? bufa_ : bufb_;
		ReadBuf* bufR = pairFw ? bufb_ : bufa_;
		uint32_t lenL = pairFw ? alen_ : blen_;
		uint32_t lenR = pairFw ? blen_ : alen_;
		bool ret;
		assert(!params_->sink().exceededOverThresh());
		params_->setFw(rL.fw);
		assert_eq(bufL->color, color);
		// Print upstream mate first
		ret = params_->reportHit(
				rL.fw ? (ebwtFwL?  bufL->patFw  :  bufL->patFwRev) :
				        (ebwtFwL?  bufL->patRc  :  bufL->patRcRev),
				rL.fw ? (ebwtFwL? &bufL->qual    : &bufL->qualRev) :
				        (ebwtFwL? &bufL->qualRev : &bufL->qual),
				&bufL->name,
				bufL->color,
				colorExEnds,
				snpPhred,
				refs_,
				rmap,
				ebwtFwL,
				rL.mms,                       // mismatch positions
				rL.refcs,                     // reference characters for mms
				rL.numMms,                    // # mismatches
				make_pair(first, upstreamOff),// position
				make_pair(first, dnstreamOff),// mate position
				rR.fw,                        // mate orientation
				lenR,                         // mate length
				make_pair(rL.top, rL.bot),    // arrows
				tlen,                         // textlen
				lenL,                         // qlen
				rL.stratum,                   // alignment stratum
				rL.cost,                      // cost, including quality penalty
				oms,                          // # other hits
				bufL->patid,
				bufL->seed,
				pairFw ? 1 : 2);
		if(ret) {
			return true; // can happen when -m is set
		}
		params_->setFw(rR.fw);
		assert_eq(bufR->color, color);
		ret = params_->reportHit(
				rR.fw ? (ebwtFwR?  bufR->patFw  :  bufR->patFwRev) :
				        (ebwtFwR?  bufR->patRc  :  bufR->patRcRev),
				rR.fw ? (ebwtFwR? &bufR->qual    : &bufR->qualRev) :
				        (ebwtFwR? &bufR->qualRev : &bufR->qual),
				&bufR->name,
				bufR->color,
				colorExEnds,
				snpPhred,
				refs_,
				rmap,
				ebwtFwR,
				rR.mms,                       // mismatch positions
				rR.refcs,                     // reference characters for mms
				rR.numMms,                    // # mismatches
				make_pair(first, dnstreamOff),// position
				make_pair(first, upstreamOff),// mate position
				rL.fw,                        // mate orientation
				lenL,                         // mate length
				make_pair(rR.top, rR.bot),    // arrows
				tlen,                         // textlen
				lenR,                         // qlen
				rR.stratum,                   // alignment stratum
				rR.cost,                      // cost, including quality penalty
				oms,                          // # other hits
				bufR->patid,
				bufR->seed,
				pairFw ? 2 : 1);
		return ret;
	}

	/**
	 * Helper for reporting a pair of alignments.  As of now, we report
	 * a paired alignment by reporting two consecutive alignments, one
	 * for each mate.
	 */
	void reportSe(const Range& r, U32Pair h, uint32_t tlen) {
		EbwtSearchParams<String<Dna> >*params = (r.mate1 ? paramsSe1_ : paramsSe2_);
		assert(!(r.mate1 ? doneSe1_ : doneSe2_));
		params->setFw(r.fw);
		ReadBuf* buf = r.mate1 ? bufa_ : bufb_;
		bool ebwtFw = r.ebwt->fw();
		uint32_t len = r.mate1 ? alen_ : blen_;
		assert_eq(buf->color, color);
		// Print upstream mate first
		if(params->reportHit(
			r.fw ? (ebwtFw?  buf->patFw   :  buf->patFwRev) :
			       (ebwtFw?  buf->patRc   :  buf->patRcRev),
			r.fw ? (ebwtFw? &buf->qual    : &buf->qualRev) :
			       (ebwtFw? &buf->qualRev : &buf->qual),
			&buf->name,
			buf->color,
			colorExEnds,
			snpPhred,
			refs_,
			r.ebwt->rmap(),
			ebwtFw,
			r.mms,                   // mismatch positions
			r.refcs,                 // reference characters for mms
			r.numMms,                // # mismatches
			h,                       // position
			make_pair(0, 0),         // (bogus) mate coords
			true,                    // (bogus) mate orientation
			0,                       // (bogus) mate length
			make_pair(r.top, r.bot), // arrows
			tlen,                    // textlen
			len,                     // qlen
			r.stratum,               // alignment stratum
			r.cost,                  // cost, including quality penalty
			r.bot - r.top - 1,       // # other hits
			buf->patid,
			buf->seed,
			0))
		{
			if(r.mate1) doneSe1_ = true;
			else        doneSe2_ = true;
			if(donePe_) driver_->removeMate(r.mate1 ? 1 : 2);
		}
	}

	void resolveOutstanding(const U32Pair& off,
	                        const uint32_t tlen,
	                        const Range& range)
	{
		assert(!this->done);
		if(!donePe_) {
			bool ret = resolveOutstandingInRef(off, tlen, range);
			if(++mixedAttempts_ > mixedAttemptLim_ || ret) {
				// Give up on this pair
				donePe_ = true;
				if(sinkPtSe1_ != NULL) {
					if(doneSe1_) driver_->removeMate(1);
					if(doneSe2_) driver_->removeMate(2);
				}
			}
			this->done = (donePe_ && (!sinkPt_->empty() || sinkPtSe1_ == NULL || (doneSe1_ && doneSe2_)));
		}
		if(!this->done && sinkPtSe1_ != NULL) {
			bool doneSe = (range.mate1 ? doneSe1_ : doneSe2_);
			if(!doneSe) {
				// Hold onto this single-end alignment in case we don't
				// find any paired alignments
				reportSe(range, off, tlen);
			}
			this->done = doneSe1_ && doneSe2_ && donePe_;
		}
	}

	/**
	 * Given a vector of reference positions where one of the two mates
	 * (the "anchor" mate) has aligned, look directly at the reference
	 * sequence for instances where the other mate (the "outstanding"
	 * mate) aligns such that mating constraint is satisfied.
	 *
	 * This function picks up to 'pick' anchors at random from the
	 * 'offs' array.  It returns the number that it actually picked.
	 */
	bool resolveOutstandingInRef(const U32Pair& off,
	                             const uint32_t tlen,
	                             const Range& range)
	{
		assert(!donePe_);
		assert(refs_->loaded());
		assert_lt(off.first, refs_->numRefs());
		// pairFw = true if the anchor indicates that the pair will
		// align in its forward orientation (i.e. with mate1 to the
		// left of mate2)
		bool pairFw = (range.mate1)? (range.fw == fw1_) : (range.fw == fw2_);
		// matchRight = true, if the opposite mate will be to the right
		// of the anchor mate
		bool matchRight = (pairFw ? range.mate1 : !range.mate1);
		// fw = orientation of the opposite mate
		bool fw = range.mate1 ? fw2_ : fw1_; // whether outstanding mate is fw/rc
		if(!pairFw) fw = !fw;
		// 'seq' = sequence for opposite mate
		const String<Dna5>& seq  =
			fw ? (range.mate1 ? patsrc_->bufb().patFw   :
		                        patsrc_->bufa().patFw)  :
		         (range.mate1 ? patsrc_->bufb().patRc   :
		                        patsrc_->bufa().patRc);
		// 'qual' = qualities for opposite mate
		const String<char>& qual =
			fw ? (range.mate1 ? patsrc_->bufb().qual  :
			                    patsrc_->bufa().qual) :
			     (range.mate1 ? patsrc_->bufb().qualRev :
			                    patsrc_->bufa().qualRev);
		uint32_t qlen = seqan::length(seq);  // length of outstanding mate
		uint32_t alen = (range.mate1 ? patsrc_->bufa().length() :
		                               patsrc_->bufb().length());
		int minins = minInsert_;
		int maxins = maxInsert_;
		if(fw1_) {
			minins = max<int>(0, minins - patsrc_->bufa().trimmed5);
			maxins = max<int>(0, maxins - patsrc_->bufa().trimmed5);
		} else {
			minins = max<int>(0, minins - patsrc_->bufa().trimmed3);
			maxins = max<int>(0, maxins - patsrc_->bufa().trimmed3);
		}
		if(fw2_) {
			minins = max<int>(0, minins - patsrc_->bufb().trimmed3);
			maxins = max<int>(0, maxins - patsrc_->bufb().trimmed3);
		} else {
			minins = max<int>(0, minins - patsrc_->bufb().trimmed5);
			maxins = max<int>(0, maxins - patsrc_->bufb().trimmed5);
		}
		assert_geq(minins, 0);
		assert_geq(maxins, 0);
		// Don't even try if either of the mates is longer than the
		// maximum insert size.
		if((uint32_t)maxins <= max(qlen, alen)) {
			return false;
		}
		const uint32_t tidx = off.first;  // text id where anchor mate hit
		const uint32_t toff = off.second; // offset where anchor mate hit
		// Set begin/end to the range of reference positions where
		// outstanding mate may align while fulfilling insert-length
		// constraints.
		uint32_t begin, end;
		assert_geq(maxins, minins);
		uint32_t insDiff = maxins - minins;
		if(matchRight) {
			end = toff + maxins;
			begin = toff + 1;
			if(qlen < alen) begin += alen-qlen;
			if(end > insDiff + qlen) {
				begin = max<uint32_t>(begin, end - insDiff - qlen);
			}
			end = min<uint32_t>(refs_->approxLen(tidx), end);
			begin = min<uint32_t>(refs_->approxLen(tidx), begin);
		} else {
			if(toff + alen < (uint32_t)maxins) {
				begin = 0;
			} else {
				begin = toff + alen - maxins;
			}
			uint32_t mi = min<uint32_t>(alen, qlen);
			end = toff + mi - 1;
			end = min<uint32_t>(end, toff + alen - minins + qlen - 1);
			if(toff + alen + qlen < (uint32_t)(minins + 1)) end = 0;
		}
		// Check if there's not enough space in the range to fit an
		// alignment for the outstanding mate.
		if(end - begin < qlen) return false;
		std::vector<Range> ranges;
		std::vector<uint32_t> offs;
		refAligner_->find(1, tidx, refs_, seq, qual, begin, end, ranges,
		                  offs, pairFw ? &pairs_fw_ : &pairs_rc_,
		                  toff, fw);
		assert_eq(ranges.size(), offs.size());
		for(size_t i = 0; i < ranges.size(); i++) {
			Range& r = ranges[i];
			r.fw = fw;
			r.cost |= (r.stratum << 14);
			r.mate1 = !range.mate1;
			const uint32_t result = offs[i];
			// Just copy the known range's top and bot for now
			r.top = range.top;
			r.bot = range.bot;
			bool ebwtLFw = matchRight ? range.ebwt->fw() : true;
			bool ebwtRFw = matchRight ? true : range.ebwt->fw();
			if(report(
				matchRight ? range : r, // range for upstream mate
				matchRight ? r : range, // range for downstream mate
				tidx,                   // ref idx
				matchRight ? toff : result, // upstream offset
				matchRight ? result : toff, // downstream offset
				tlen,       // length of ref
				pairFw,     // whether the pair is being mapped to fw strand
				ebwtLFw,
				ebwtRFw,
				range.ebwt->rmap())) return true;
		}
		return false;
	}

	const BitPairReference* refs_;

	PatternSourcePerThread *patsrc_;
	uint32_t qlen1_, qlen2_;
	bool chase_;

	// true -> we're no longer shooting for paired-end alignments;
	// just collecting single-end ones
	bool donePe_, doneSe1_, doneSe2_;

	// For searching for outstanding mates
	RefAligner<String<Dna5> >* refAligner_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;
	HitSinkPerThread* sinkPtSe1_, * sinkPtSe2_;

	// State for alignment
	EbwtSearchParams<String<Dna> >* params_;
	// for single-end:
	EbwtSearchParams<String<Dna> >* paramsSe1_, * paramsSe2_;

	// Paired-end boundaries
	const uint32_t minInsert_;
	const uint32_t maxInsert_;

	const uint32_t mixedAttemptLim_;
	uint32_t mixedAttempts_;

	// Orientation of upstream/downstream mates when aligning to
	// forward strand
	const bool fw1_, fw2_;

	// State for getting alignments from ranges statefully
	RangeChaser<String<Dna> >* rchase_;

	// Range-finding state for first mate
	TDriver* driver_;

	// Pool for distributing chunks of best-first path descriptor memory
	ChunkPool *pool_;

	bool verbose_;
	bool quiet_;

	int maxBts_; // maximum allowed # backtracks
	int *btCnt_; // current backtrack count

	/// For keeping track of paired alignments that have already been
	/// found for the forward and reverse-comp pair orientations
	TSetPairs pairs_fw_, pairs_rc_;
};

#endif /* ALIGNER_H_ */
