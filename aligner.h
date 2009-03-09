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

/**
 * State machine for carrying out an alignment, which usually consists
 * of a series of phases that conduct different alignments using
 * different backtracking constraints.
 *
 * Each Aligner should have a dedicated PatternSourcePerThread.
 */
class Aligner {
public:
	Aligner(bool _done, bool rangeMode, uint32_t seed) :
		done(_done), patsrc_(NULL), bufa_(NULL), bufb_(NULL),
		rangeMode_(rangeMode), seed_(seed)
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
		qseed_ = seed_ + genRandSeed(bufa_->patFw, bufa_->qualFw, bufa_->name);
		rand_.init(qseed_);
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
	// RandomSource for choosing alignments to report from ranges
	bool rangeMode_;
	uint32_t seed_;
	uint32_t qseed_; // query-specific seed
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
				if(!al->done) {
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
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(true, rangeMode, seed),
		doneFirst_(true),
		firstIsFw_(true),
		chase_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		params_(params),
		rchase_(rchase),
		driver_(driver)
	{
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
		assert(driver_ != NULL);
	}

	virtual ~UnpairedAlignerV2() {
		delete driver_; driver_ = NULL;
		delete params_; params_ = NULL;
		delete rchase_; rchase_ = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		driver_->setQuery(patsrc);
		rchase_->initRand(qseed_);
		this->done = driver_->done;
		doneFirst_ = false;
		firstIsFw_ = ((qseed_ & 0x10) == 0);
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
		return params_->reportHit(
				ra.fw ? (ebwtFw? bufa_->patFw   : bufa_->patFwRev) :
				        (ebwtFw? bufa_->patRc   : bufa_->patRcRev),
				ra.fw ? (ebwtFw? &bufa_->qualFw : &bufa_->qualFwRev) :
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
				ra.bot - ra.top - 1,      // # other hits
				patsrc_->patid());
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual bool advance() {
		assert(!this->done);
		if(chase_) {
			assert(!rangeMode_);
			assert(driver_->foundRange);
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
				this->done = driver_->done;
			}
		}
		// Still advancing a
		if(!this->done && !chase_) {
			assert(!driver_->done);
			driver_->advance();
			if(driver_->foundRange) {
				const Range& ra = driver_->range();
				if(rangeMode_) {
					this->done = report(ra, ra.top, ra.bot, 0);
				} else {
					rchase_->setTopBot(ra.top, ra.bot, alen_, ra.ebwt);
					if(rchase_->foundOff()) {
						this->done = report(
								ra, rchase_->off().first,
								rchase_->off().second, rchase_->tlen());
						rchase_->reset();
					}
					if(!rchase_->done) {
						// Keep chasing this range
						chase_ = true;
					}
				}
			}
			if(driver_->done && !chase_) {
				this->done = true;
			}
		}
		if(this->done) {
			sinkPt_->finishRead(*patsrc_, true);
		}
		return this->done;
	}

protected:
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
};

/**
 * An aligner for finding exact matches of unpaired reads.  Always
 * tries the forward-oriented version of the read before the reverse-
 * oriented read.
 */
template<typename TRangeSource>
class UnpairedAlignerV1 : public Aligner {
	typedef RangeSourceDriver<TRangeSource> TDriver;
public:
	UnpairedAlignerV1(
		EbwtSearchParams<String<Dna> >* params,
		TDriver* driverFw, TDriver* driverRc,
		RangeChaser<String<Dna> >* rchase,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		HitSinkPerThread* sinkPt,
		vector<String<Dna5> >& os,
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(true, rangeMode, seed),
		doneFirst_(true),
		firstIsFw_(true),
		chaseFw_(false), chaseRc_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPt),
		params_(params),
		rchase_(rchase),
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
		delete rchase_;   rchase_   = NULL;
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		driverFw_->setQuery(patsrc);
		driverRc_->setQuery(patsrc);
		rchase_->initRand(qseed_);
		this->done = false;
		doneFirst_ = false;
		firstIsFw_ = ((qseed_ & 0x10) == 0);
		chaseFw_ = false;
		chaseRc_ = false;
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
		return params_->reportHit(
				ra.fw ? (ebwtFw? bufa_->patFw   : bufa_->patFwRev) :
				        (ebwtFw? bufa_->patRc   : bufa_->patRcRev),
				ra.fw ? (ebwtFw? &bufa_->qualFw : &bufa_->qualFwRev) :
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
				ra.bot - ra.top - 1,      // # other hits
				patsrc_->patid());
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual bool advance() {
		assert(!this->done);
		assert(!chaseFw_ || !chaseRc_);
		if(chaseFw_) {
			assert(!rangeMode_);
			assert(driverFw_->foundRange);
			if(!rchase_->foundOff() && !rchase_->done) {
				params_->setFw(true);
				rchase_->advance();
				return false;
			}
			if(rchase_->foundOff()) {
				this->done = report(driverFw_->range(), rchase_->off().first,
				                    rchase_->off().second, rchase_->tlen());
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				chaseFw_ = false;
				if(doneFirst_) this->done = driverFw_->done;
				else           doneFirst_ = driverFw_->done;
			}
		} else if(chaseRc_) {
			assert(!rangeMode_);
			assert(driverRc_->foundRange);
			if(!rchase_->foundOff() && !rchase_->done) {
				params_->setFw(false);
				rchase_->advance();
				return false;
			}
			if(rchase_->foundOff()) {
				this->done = report(driverRc_->range(), rchase_->off().first,
				                    rchase_->off().second, rchase_->tlen());
				rchase_->reset();
			} else {
				assert(rchase_->done);
				// Forget this range; keep looking for ranges
				chaseRc_ = false;
				if(doneFirst_) this->done = driverRc_->done;
				else           doneFirst_ = driverRc_->done;
			}
		}
		// Still advancing a
		if(!this->done && !chaseFw_ && !chaseRc_) {
			bool fw = (doneFirst_ != firstIsFw_);
			if(fw) {
				params_->setFw(true);
				driverFw_->advance();
				if(driverFw_->foundRange) {
					const Range& ra = driverFw_->range();
					assert(ra.fw == fw);
					if(rangeMode_) {
						this->done = report(ra, ra.top, ra.bot, 0);
					} else {
						rchase_->setTopBot(ra.top, ra.bot, alen_, ra.ebwt);
						if(rchase_->foundOff()) {
							this->done = report(ra, rchase_->off().first,
							               rchase_->off().second, rchase_->tlen());
							rchase_->reset();
						}
						if(!rchase_->done) {
							// Keep chasing this range
							chaseFw_ = true;
						}
					}
				}
				if(driverFw_->done && !chaseFw_) {
					if(doneFirst_ ) this->done = true;
					else            doneFirst_ = true;
				}
			} else {
				params_->setFw(false);
				driverRc_->advance();
				// Advance the RangeSource for the reverse-complement read
				if(driverRc_->foundRange) {
					const Range& ra = driverRc_->range();
					if(rangeMode_) {
						this->done = report(ra, ra.top, ra.bot, 0);
					} else {
						rchase_->setTopBot(ra.top, ra.bot, alen_, ra.ebwt);
						if(rchase_->foundOff()) {
							this->done = report(ra, rchase_->off().first,
							               rchase_->off().second, rchase_->tlen());
							rchase_->reset();
						}
						if(!rchase_->done) {
							// Keep chasing this range
							chaseRc_ = true;
						}
					}
				}
				if(driverRc_->done && !chaseRc_) {
					if(doneFirst_ ) this->done = true;
					else            doneFirst_ = true;
				}
			}
		}
		if(this->done) {
			sinkPt_->finishRead(*patsrc_, true);
		}
		return this->done;
	}

protected:
	// Progress state
	bool doneFirst_;
	bool firstIsFw_;
	bool chaseFw_;
	bool chaseRc_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> >* params_;

	// State for getting alignments from ranges statefully
	RangeChaser<String<Dna> >* rchase_;

	// Range-finding state
	TDriver* driverFw_;
	TDriver* driverRc_;
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
		uint32_t seed) :
		Aligner(true, rangeMode, seed),
		refs_(refs), patsrc_(NULL), doneFw_(true),
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
		assert(sinkPt_ != NULL);
		assert(params_ != NULL);
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
		driver1Fw_->setQuery(patsrc, true  /* mate1 */);
		driver1Rc_->setQuery(patsrc, true  /* mate1 */);
		driver2Fw_->setQuery(patsrc, false /* mate2 */);
		driver2Rc_->setQuery(patsrc, false /* mate2 */);
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
		rchase_->initRand(qseed_);
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
		bool verbose = false;
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
		advanceOrientation(!doneFw_, verbose);
		if(this->done) {
			if(verbose2_) cout << "----" << endl;
			sinkPt_->finishRead(*patsrc_, true);
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
	            bool ebwtFwR)
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
		// Print upstream mate first
		ret = params_->reportHit(
				rL.fw ? (ebwtFwL?  bufL->patFw  :  bufL->patFwRev) :
					    (ebwtFwL?  bufL->patRc  :  bufL->patRcRev),
				rL.fw ? (ebwtFwL? &bufL->qualFw : &bufL->qualFwRev) :
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
				oms,                          // # other hits
				bufL->patid,
				pairFw ? 1 : 2);
		if(ret) {
			return true; // can happen when -m is set
		}
		params_->setFw(rR.fw);
		ret = params_->reportHit(
				rR.fw ? (ebwtFwR?  bufR->patFw  :  bufR->patFwRev) :
					    (ebwtFwR?  bufR->patRc  :  bufR->patRcRev),
				rR.fw ? (ebwtFwR? &bufR->qualFw : &bufR->qualFwRev) :
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
				oms,                          // # other hits
				bufR->patid,
				pairFw ? 2 : 1);
		return ret;
	}

	bool report(const Range& rL, // range for upstream mate
	            const Range& rR, // range for downstream mate
	            uint32_t first,  // ref idx
	            uint32_t upstreamOff, // offset for upstream mate
	            uint32_t dnstreamOff, // offset for downstream mate
	            uint32_t tlen, // length of ref
	            bool pairFw)   // whether the pair is being mapped to fw strand
	{
		return report(rL, rR, first, upstreamOff, dnstreamOff, tlen,
		              pairFw, rL.ebwt->fw(), rR.ebwt->fw());
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
	                     U32PairVec* offsLarr,
	                     U32PairVec* offsRarr,
	                     TRangeVec* rangesLarr,
	                     TRangeVec* rangesRarr,
		                 TDriver& drL,
		                 TDriver& drR,
		                 bool pairFw,
		                 bool verbose = false)
	{
		U32PairVec& offsL   = offsLarr  [h.first & 31];
		U32PairVec& offsR   = offsRarr  [h.first & 31];
		TRangeVec&  rangesL = rangesLarr[h.first & 31];
		TRangeVec&  rangesR = rangesRarr[h.first & 31];
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
					// Incoming hit hits same reference as buffered
					uint32_t left  = newFromL ? h.second : h2.second;
					uint32_t right = newFromL ? h2.second : h.second;
					if(right > left) {
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
								   rL.ebwt->_plen[h.first],
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
	 * Given a vector of reference positions where one of the two mates
	 * (the "anchor" mate) has aligned, look directly at the reference
	 * sequence for instances where the other mate (the "outstanding"
	 * mate) aligns such that mating constraint is satisfied.
	 *
	 * This function picks up to 'pick' anchors at random from the
	 * 'offs' array.  It returns the number that it actually picked.
	 */
	bool resolveOutstandingInRef(const bool offs1,
	                             const U32Pair& off,
	                             const uint32_t tlen,
	                             const Range& range)
	{
		assert(refs_->loaded());
		assert_lt(off.first, refs_->numRefs());
		// If matchRight is true, then we're trying to align the other
		// mate to the right of the already-aligned mate.  Otherwise,
		// to the left.
		bool matchRight = (offs1 ? !doneFw_ : doneFw_);
		// Sequence and quals for mate to be matched
		bool fw = offs1 ? fw2_ : fw1_; // whether outstanding mate is fw/rc
		if(doneFw_) fw = !fw;
		// 'seq' gets sequence of outstanding mate w/r/t the forward
		// reference strand
		const String<Dna5>& seq  = fw ? (offs1 ? patsrc_->bufb().patFw   :
		                                         patsrc_->bufa().patFw)  :
		                                (offs1 ? patsrc_->bufb().patRc   :
		                                         patsrc_->bufa().patRc);
		// 'seq' gets qualities of outstanding mate w/r/t the forward
		// reference strand
		const String<char>& qual = fw ? (offs1 ? patsrc_->bufb().qualFw  :
		                                         patsrc_->bufa().qualFw) :
		                                (offs1 ? patsrc_->bufb().qualRc  :
		                                         patsrc_->bufa().qualRc);
		uint32_t qlen = seqan::length(seq);  // length of outstanding mate
		uint32_t alen = (offs1 ? patsrc_->bufa().length() :
		                         patsrc_->bufb().length());
		// Don't even try if either of the mates is longer than the
		// maximum insert size; this seems to be compatible with what
		// Maq does.
		if(maxInsert_ <= max(qlen, alen)) {
			return false;
		}
		const uint32_t tidx = off.first;
		const uint32_t toff = off.second;
		// Set begin/end to be a range of all reference
		// positions that are legally permitted to be involved in
		// the alignment of the outstanding mate.  It's up to the
		// callee to worry about how to scan these positions.
		uint32_t begin, end;
		if(matchRight) {
			begin = toff + 1;
			end = toff + maxInsert_;
			end = min<uint32_t>(refs_->approxLen(tidx), end);
		} else {
			if(toff + alen < maxInsert_) {
				begin = 0;
			} else {
				begin = toff + alen - maxInsert_;
			}
			end = toff + min<uint32_t>(alen, qlen) - 1;
		}
		// Check if there's not enough space in the range to fit an
		// alignment for the outstanding mate.
		if(end - begin < qlen) return false;
		std::vector<Range> ranges;
		std::vector<uint32_t> offs;
		refAligner_->find(1, tidx, refs_, seq, qual, begin, end, ranges,
		                  offs, doneFw_ ? &pairs_rc_ : &pairs_fw_,
		                  toff, !matchRight);
		assert_eq(ranges.size(), offs.size());
		for(size_t i = 0; i < ranges.size(); i++) {
			Range& r = ranges[i];
			r.fw = fw;
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
				    ebwtRFw)) return true;
		}
		return false;
	}

	/**
	 * Advance paired-end alignment.
	 */
	void advanceOrientation(bool pairFw, bool verbose = false) {
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
				if(!dontReconcile_ && !overThresh) {
					this->done = reconcileAndAdd(
							rchase_->off(), true /* new entry is from 1 */,
							offsLarr_, offsRarr_, rangesLarr_, rangesRarr_,
							*drL_, *drR_, pairFw, verbose);
				}
				if(!this->done && (overThresh || dontReconcile_)) {
					// Because the total size of both ranges exceeds
					// our threshold, we're now operating in "mixed
					// mode"
					const Range& r = drL_->range();
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
				if(verbose) cout << "Done with case for first mate" << endl;
				if(*delayedchaseR_) {
					// Start chasing the delayed range
					if(verbose) cout << "Resuming delayed chase for second mate" << endl;
					assert(drR_->foundRange);
					const Range& r = drR_->range();
					uint32_t top = r.top;
					uint32_t bot = r.bot;
					rchase_->setTopBot(top, bot, drR_->qlen(), r.ebwt);
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
				if(!dontReconcile_ && !overThresh) {
					this->done = reconcileAndAdd(
							rchase_->off(), false /* new entry is from 2 */,
							offsLarr_, offsRarr_, rangesLarr_, rangesRarr_,
							*drL_, *drR_, pairFw, verbose);
				}
				if(!this->done && (overThresh || dontReconcile_)) {
					// Because the total size of both ranges exceeds
					// our threshold, we're now operating in "mixed
					// mode"
					const Range& r = drR_->range();
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
				if(verbose) cout << "Done with case for second mate" << endl;
				if(*delayedchaseL_) {
					// Start chasing the delayed range
					if(verbose) cout << "Resuming delayed chase for first mate" << endl;
					assert(drL_->foundRange);
					const Range& r = drL_->range();
					uint32_t top = r.top;
					uint32_t bot = r.bot;
					rchase_->setTopBot(top, bot, drL_->qlen(), r.ebwt);
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
					if(verbose) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 1" << endl;
					*donePair_ = true;
					if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << endl;
					return;
				}
				assert(!*delayedchaseL_);
				drL_->advance();
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
						if(verbose) cout << "Delaying a chase for first mate" << endl;
						*delayedchaseL_ = true;
					} else {
						// Start chasing this range
						if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << " " << drL_->range().top << endl;
						if(verbose) cout << "Chasing a range for first mate" << endl;
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
							rchase_->setTopBot(r.top, r.bot, drR_->qlen(), r.ebwt);
						} else {
							// Use Burrows-Wheeler for this pair (as
							// usual)
							*chaseL_ = true;
							const Range& r = drL_->range();
							rchase_->setTopBot(r.top, r.bot, drL_->qlen(), r.ebwt);
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
					if(verbose) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 2" << endl;
					if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << endl;
					*donePair_ = true;
					return;
				}
				assert(!*delayedchaseR_);
				drR_->advance();
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
						if(verbose) cout << "Delaying a chase for second mate" << endl;
						*delayedchaseR_ = true;
					} else {
						// Start chasing this range
						if(verbose2_) cout << *offsLsz_ << " " << *offsRsz_ << " " << drR_->range().top << endl;
						if(verbose) cout << "Chasing a range for second mate" << endl;
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
							rchase_->setTopBot(r.top, r.bot, drL_->qlen(), r.ebwt);
						} else {
							// Use Burrows-Wheeler for this pair (as
							// usual)
							*chaseR_ = true;
							const Range& r = drR_->range();
							rchase_->setTopBot(r.top, r.bot, drR_->qlen(), r.ebwt);
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

#endif /* ALIGNER_H_ */
