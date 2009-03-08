/*
 * aligner_1mm.h
 */

#ifndef ALIGNER_1MM_H_
#define ALIGNER_1MM_H_

#include <utility>
#include <vector>
#include "aligner.h"
#include "hit.h"
#include "range_source.h"
#include "row_chaser.h"
#include "range_chaser.h"
#include "ref_aligner.h"

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class Unpaired1mmAlignerV1Factory : public AlignerFactory {
	typedef SingleRangeSourceDriver<GreedyDFSRangeSource> TRangeSrcDr;
	typedef ListRangeSourceDriver<GreedyDFSRangeSource> TListRangeSrcDr;
	typedef CostAwareRangeSourceDriver<GreedyDFSRangeSource> TCostAwareRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Unpaired1mmAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache *cacheFw,
			RangeCache *cacheBw,
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
		assert(ebwtBw != NULL);
		assert(ebwtBw->isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {

		HitSinkPerThread* sinkPt = sinkPtFactory_.create();
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_, true, true, true, rangeMode_);

		// Source for ranges from mirror index & forward read
		GreedyDFSRangeSource *rFw_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from forward index & forward read
		GreedyDFSRangeSource *rFw_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		//
		GreedyDFSRangeSourceDriver * drFw_Bw = new GreedyDFSRangeSourceDriver(
			*params, rFw_Bw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		//
		GreedyDFSRangeSourceDriver * drFw_Fw = new GreedyDFSRangeSourceDriver(
			*params, rFw_Fw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec drVec;
		drVec.push_back(drFw_Bw);
		drVec.push_back(drFw_Fw);

		// Source for ranges from forward index & reverse-complement read
		GreedyDFSRangeSource *rRc_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & reverse-complement read
		GreedyDFSRangeSource *rRc_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		//
		GreedyDFSRangeSourceDriver * drRc_Fw = new GreedyDFSRangeSourceDriver(
			*params, rRc_Fw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		//
		GreedyDFSRangeSourceDriver * drRc_Bw = new GreedyDFSRangeSourceDriver(
			*params, rRc_Bw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		drVec.push_back(drRc_Fw); drVec.push_back(drRc_Bw);
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(seed_, drVec);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		// Set up the aligner
		return new UnpairedAlignerV2<GreedyDFSRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, rangeMode_, verbose_, seed_);
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
class Paired1mmAlignerV1Factory : public AlignerFactory {
	typedef SingleRangeSourceDriver<GreedyDFSRangeSource> TRangeSrcDr;
	typedef ListRangeSourceDriver<GreedyDFSRangeSource> TListRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Paired1mmAlignerV1Factory(
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
			RangeCache *cacheFw,
			RangeCache *cacheBw,
			uint32_t cacheLimit,
			BitPairReference* refs,
			vector<String<Dna5> >& os,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
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
		assert(ebwtBw != NULL);
		assert(ebwtFw.isInMemory());
		assert(ebwtBw->isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2);
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_, true, true, true, rangeMode_);

		GreedyDFSRangeSource *r1Fw_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSource *r1Fw_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSourceDriver * dr1Fw_Fw = new GreedyDFSRangeSourceDriver(
			*params, r1Fw_Fw, true, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		GreedyDFSRangeSourceDriver * dr1Fw_Bw = new GreedyDFSRangeSourceDriver(
			*params, r1Fw_Bw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr1FwVec;
		dr1FwVec.push_back(dr1Fw_Fw); dr1FwVec.push_back(dr1Fw_Bw);
		// Overall range source driver for the forward orientation of
		// the first mate
		TListRangeSrcDr* dr1Fw = new TListRangeSrcDr(dr1FwVec);

		GreedyDFSRangeSource *r1Rc_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSource *r1Rc_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSourceDriver * dr1Rc_Fw = new GreedyDFSRangeSourceDriver(
			*params, r1Rc_Fw, false, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		GreedyDFSRangeSourceDriver * dr1Rc_Bw = new GreedyDFSRangeSourceDriver(
			*params, r1Rc_Bw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr1RcVec;
		dr1RcVec.push_back(dr1Rc_Fw); dr1RcVec.push_back(dr1Rc_Bw);
		// Overall range source driver for the reverse-comp orientation
		// of the first mate
		TListRangeSrcDr* dr1Rc = new TListRangeSrcDr(dr1RcVec);

		GreedyDFSRangeSource *r2Fw_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSource *r2Fw_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSourceDriver * dr2Fw_Fw = new GreedyDFSRangeSourceDriver(
			*params, r2Fw_Fw, true, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		GreedyDFSRangeSourceDriver * dr2Fw_Bw = new GreedyDFSRangeSourceDriver(
			*params, r2Fw_Bw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr2FwVec;
		dr2FwVec.push_back(dr2Fw_Fw); dr2FwVec.push_back(dr2Fw_Bw);
		// Overall range source driver for the forward orientation of
		// the first mate
		TListRangeSrcDr* dr2Fw = new TListRangeSrcDr(dr2FwVec);

		GreedyDFSRangeSource *r2Rc_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSource *r2Rc_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		GreedyDFSRangeSourceDriver * dr2Rc_Fw = new GreedyDFSRangeSourceDriver(
			*params, r2Rc_Fw, false, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		GreedyDFSRangeSourceDriver * dr2Rc_Bw = new GreedyDFSRangeSourceDriver(
			*params, r2Rc_Bw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr2RcVec;
		dr2RcVec.push_back(dr2Rc_Fw); dr2RcVec.push_back(dr2Rc_Bw);
		// Overall range source driver for the reverse-comp orientation
		// of the first mate
		TListRangeSrcDr* dr2Rc = new TListRangeSrcDr(dr2RcVec);

		RefAligner<String<Dna5> >* refAligner = new OneMMRefAligner<String<Dna5> >(0);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		return new PairedBWAlignerV1<GreedyDFSRangeSource>(
			params, dr1Fw, dr1Rc, dr2Fw, dr2Rc, refAligner, rchase,
			sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
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

#endif /* ALIGNER_1MM_H_ */
