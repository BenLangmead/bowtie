/*
 * aligner_0mm.h
 */

#ifndef ALIGNER_0MM_H_
#define ALIGNER_0MM_H_

#include <utility>
#include <vector>
#include "aligner.h"
#include "hit.h"
#include "row_chaser.h"
#include "range_chaser.h"

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class UnpairedExactAlignerV1Factory : public AlignerFactory {
public:
	UnpairedExactAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache* cacheFw,
			RangeCache* cacheBw,
			uint32_t cacheLimit,
			vector<String<Dna5> >& os,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			os_(os),
			rangeMode_(rangeMode),
			verbose_(verbose),
			seed_(seed)
	{
		assert(ebwtFw.isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.create();
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_, true, true, true, rangeMode_);
		GreedyDFSRangeSource *rFw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff, BacktrackLimits(), 0, true,
			false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSource *rRc = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff, BacktrackLimits(), 0, true,
			false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSourceDriver * driverFw = new GreedyDFSRangeSourceDriver(
			*params, rFw, true, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		GreedyDFSRangeSourceDriver * driverRc = new GreedyDFSRangeSourceDriver(
			*params, rRc, false, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		return new UnpairedAlignerV1<GreedyDFSRangeSource>(
			params, driverFw, driverRc, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, rangeMode_, verbose_,
			seed_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	vector<String<Dna5> >& os_;
	bool rangeMode_;
	bool verbose_;
	uint32_t seed_;
};

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class PairedExactAlignerV1Factory : public AlignerFactory {
public:
	PairedExactAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			bool mate1fw,
			bool mate2fw,
			uint32_t peInner,
			uint32_t peOuter,
			bool dontReconcile,
			uint32_t symCeil,
			uint32_t mixedThresh,
			uint32_t mixedAttemptLim,
			RangeCache* cacheFw,
			RangeCache* cacheBw,
			uint32_t cacheLimit,
			BitPairReference* refs,
			vector<String<Dna5> >& os,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			mate1fw_(mate1fw),
			mate2fw_(mate2fw),
			peInner_(peInner),
			peOuter_(peOuter),
			dontReconcile_(dontReconcile),
			symCeil_(symCeil),
			mixedThresh_(mixedThresh),
			mixedAttemptLim_(mixedAttemptLim),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			refs_(refs), os_(os),
			rangeMode_(rangeMode),
			verbose_(verbose),
			seed_(seed)
	{
		assert(ebwtFw.isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2);
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_, true, true, true, rangeMode_);

		GreedyDFSRangeSource *r1Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff, BacktrackLimits(), 0, true,
			false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSource *r1Rc = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff, BacktrackLimits(), 0, true,
			false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSourceDriver * driver1Fw = new GreedyDFSRangeSourceDriver(
			*params, r1Fw, true, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		GreedyDFSRangeSourceDriver * driver1Rc = new GreedyDFSRangeSourceDriver(
			*params, r1Rc, false, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);

		GreedyDFSRangeSource *r2Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff, BacktrackLimits(), 0, true,
			false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSource *r2Rc = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff, BacktrackLimits(), 0, true,
			false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSourceDriver * driver2Fw = new GreedyDFSRangeSourceDriver(
			*params, r2Fw, true, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		GreedyDFSRangeSourceDriver * driver2Rc = new GreedyDFSRangeSourceDriver(
			*params, r2Rc, false, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);

		RefAligner<String<Dna5> >* refAligner = new ExactRefAligner<String<Dna5> >(0);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		return new PairedBWAlignerV1<GreedyDFSRangeSource>(
			params,
			driver1Fw, driver1Rc, driver2Fw, driver2Rc, refAligner,
			rchase, sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
			peInner_, peOuter_, dontReconcile_, symCeil_, mixedThresh_,
			mixedAttemptLim_, refs_, rangeMode_, verbose_, seed_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	const bool mate1fw_;
	const bool mate2fw_;
	const uint32_t peInner_;
	const uint32_t peOuter_;
	const bool dontReconcile_;
	const uint32_t symCeil_;
	const uint32_t mixedThresh_;
	const uint32_t mixedAttemptLim_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	BitPairReference* refs_;
	vector<String<Dna5> >& os_;
	const bool rangeMode_;
	const bool verbose_;
	const uint32_t seed_;
};

#endif /* ALIGNER_0MM_H_ */
