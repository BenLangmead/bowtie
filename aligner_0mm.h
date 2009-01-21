/*
 * aligner_0mm.h
 */

#ifndef ALIGNER_0MM_H_
#define ALIGNER_0MM_H_

#include "aligner.h"
#include "hit.h"
#include "row_chaser.h"

/**
 * An aligner for finding exact matches of paired reads.
 */
class PairedExactAlignerV1 : public Aligner {
	typedef std::pair<uint32_t,uint32_t> U32Pair;
public:
	PairedExactAlignerV1(
		const Ebwt<String<Dna> >& ebwt,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		vector<String<Dna5> >& os,
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(&ebwt, NULL, rangeMode, seed),
		doneFw_(true), done_(true), firstFw_(true), firstRc_(true),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPtFactory.create()),
		params_(*sinkPt_, os, true, true, true, rangeMode),
		r1Fw_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c1Fw_(),
		r1Rc_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c1Rc_(),
		driver1Fw_(ebwt, params_, r1Fw_, c1Fw_, true,  sink, sinkPt_, os, verbose, seed),
		driver1Rc_(ebwt, params_, r1Rc_, c1Rc_, false, sink, sinkPt_, os, verbose, seed),
		r2Fw_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c2Fw_(),
		r2Rc_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c2Rc_(),
		driver2Fw_(ebwt, params_, r2Fw_, c2Fw_, true,  sink, sinkPt_, os, verbose, seed),
		driver2Rc_(ebwt, params_, r2Rc_, c2Rc_, false, sink, sinkPt_, os, verbose, seed),
		fwRs_(), rcRs_()
	{ }

	virtual ~PairedExactAlignerV1() {
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		if(patsrc->bufb().empty()) {
			// Can't handle single-ended alignments; use
			// UnpairedExactAlignerV1 for that
			return;
		}
		driver1Fw_.setQuery(patsrc, true  /* mate1 */);
		driver1Rc_.setQuery(patsrc, true  /* mate1 */);
		driver2Fw_.setQuery(patsrc, false /* mate2 */);
		driver2Rc_.setQuery(patsrc, false /* mate2 */);
		doneFw_  = false;
		done_    = false;
		firstFw_ = true;
		firstRc_ = true;
		// Clear all intermediate ranges
		fwRs_.clear();
		rcRs_.clear();
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
		if(!doneFw_) {
//			if(fwRs_.empty() && !driver1Fw_.done()) {
//				driver1Fw_.advance();
//			} else if
//			if(driver1Fw_.foundRange()) {
//				done_ = reportSingleEndHitFromRange(driver1Fw_.range(), params_, true);
//			}
			if(!done_) doneFw_ = driver1Fw_.done();
		} else {
			driver1Rc_.advance();
			// Advance the RangeSource for the reverse-complement read
			if(driver1Rc_.foundRange()) {
				done_ = reportSingleEndHitFromRange(driver1Rc_.range(), params_, false);
			}
			if(!done_) done_ = driver1Rc_.done();
		}
		if(done_) {
			sinkPt_->finishRead(*patsrc_, NULL, NULL, NULL, NULL);
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
	bool doneFw_;   // finished with forward orientation of both mates?
	bool done_;
	bool firstFw_;
	bool firstRc_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> > params_;

	// First mate
	GreedyDFSRangeSource           r1Fw_;
	GreedyDFSContinuationManager   c1Fw_;
	GreedyDFSRangeSource           r1Rc_;
	GreedyDFSContinuationManager   c1Rc_;
	RangeSourceDriver<GreedyDFSRangeSource, GreedyDFSContinuationManager> driver1Fw_;
	RangeSourceDriver<GreedyDFSRangeSource, GreedyDFSContinuationManager> driver1Rc_;

	// Second mate
	GreedyDFSRangeSource           r2Fw_;
	GreedyDFSContinuationManager   c2Fw_;
	GreedyDFSRangeSource           r2Rc_;
	GreedyDFSContinuationManager   c2Rc_;
	RangeSourceDriver<GreedyDFSRangeSource, GreedyDFSContinuationManager> driver2Fw_;
	RangeSourceDriver<GreedyDFSRangeSource, GreedyDFSContinuationManager> driver2Rc_;

	// Range buffers
	std::vector<U32Pair> fwRs_;
	std::vector<U32Pair> rcRs_;
};

/**
 * An aligner for finding exact matches of unpaired reads.  Always
 * tries the forward-oriented version of the read before the reverse-
 * oriented read.
 */
class UnpairedExactAlignerV1 : public Aligner {
public:
	UnpairedExactAlignerV1(
		const Ebwt<String<Dna> >& ebwt,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		vector<String<Dna5> >& os,
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(&ebwt, NULL, rangeMode, seed),
		doneFw_(true), done_(true), firstFw_(true), firstRc_(true),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPtFactory.create()),
		params_(*sinkPt_, os, true, true, true, rangeMode),
		rFw_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		    false, NULL, NULL, verbose, seed, &os, false, false, false),
		cFw_(),
		rRc_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		    false, NULL, NULL, verbose, seed, &os, false, false, false),
		cRc_(),
		driverFw_(ebwt, params_, rFw_, cFw_, true,  sink, sinkPt_, os, verbose, seed),
		driverRc_(ebwt, params_, rRc_, cRc_, false, sink, sinkPt_, os, verbose, seed)
	{ }

	virtual ~UnpairedExactAlignerV1() {
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		driverFw_.setQuery(patsrc);
		driverRc_.setQuery(patsrc);
		doneFw_  = false;
		done_    = false;
		firstFw_ = true;
		firstRc_ = true;
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual bool advance() {
		assert(!done_);
		if(!doneFw_) {
			driverFw_.advance();
			if(driverFw_.foundRange()) {
				done_ = reportSingleEndHitFromRange(driverFw_.range(), params_, true);
			}
			if(!done_) doneFw_ = driverFw_.done();
		} else {
			driverRc_.advance();
			// Advance the RangeSource for the reverse-complement read
			if(driverRc_.foundRange()) {
				done_ = reportSingleEndHitFromRange(driverRc_.range(), params_, false);
			}
			if(!done_) done_ = driverRc_.done();
		}
		if(done_) {
			sinkPt_->finishRead(*patsrc_, NULL, NULL, NULL, NULL);
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
	bool firstFw_;
	bool firstRc_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> > params_;
	GreedyDFSRangeSource           rFw_;
	GreedyDFSContinuationManager   cFw_;
	GreedyDFSRangeSource           rRc_;
	GreedyDFSContinuationManager   cRc_;
	RangeSourceDriver<GreedyDFSRangeSource, GreedyDFSContinuationManager> driverFw_;
	RangeSourceDriver<GreedyDFSRangeSource, GreedyDFSContinuationManager> driverRc_;
};

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class UnpairedExactAlignerV1Factory : public AlignerFactory {
public:
	UnpairedExactAlignerV1Factory(
			Ebwt<String<Dna> >& ebwt,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			vector<String<Dna5> >& os,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwt_(ebwt),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			os_(os),
			rangeMode_(rangeMode),
			verbose_(verbose),
			seed_(seed)
	{ }

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		return new UnpairedExactAlignerV1(
			ebwt_, sink_, sinkPtFactory_, os_, rangeMode_, verbose_, seed_);
	}

	/**
	 * Create a new vector of new UnpairedExactAlignerV1s.
	 */
	virtual std::vector<Aligner*>* create(uint32_t n) const {
		std::vector<Aligner*>* v = new std::vector<Aligner*>;
		for(uint32_t i = 0; i < n; i++) {
			v->push_back(new UnpairedExactAlignerV1(
				ebwt_, sink_, sinkPtFactory_, os_, rangeMode_, verbose_, seed_));
			assert(v->back() != NULL);
		}
		return v;
	}

private:
	Ebwt<String<Dna> >& ebwt_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	vector<String<Dna5> >& os_;
	bool rangeMode_;
	bool verbose_;
	uint32_t seed_;
};

#endif /* ALIGNER_0MM_H_ */
