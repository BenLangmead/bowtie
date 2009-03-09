/*
 * aligner_23mm.h
 */

#ifndef ALIGNER_23MM_H_
#define ALIGNER_23MM_H_

#include <utility>
#include <vector>
#include "aligner.h"
#include "hit.h"
#include "range_source.h"
#include "row_chaser.h"
#include "range_chaser.h"
#include "ref_aligner.h"

/**
 * Concrete factory for constructing unpaired 2- or 3-mismatch aligners.
 */
class Unpaired23mmAlignerV1Factory : public AlignerFactory {
	typedef SingleRangeSourceDriver<GreedyDFSRangeSource> TRangeSrcDr;
	typedef CostAwareRangeSourceDriver<GreedyDFSRangeSource> TCostAwareRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Unpaired23mmAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool two,
			bool doFw,
			bool doRc,
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
			two_(two),
			doFw_(doFw), doRc_(doRc),
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

		// Source for ranges from forward index & forward read
		GreedyDFSRangeSource *rFw_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & forward read
		GreedyDFSRangeSource *rFw_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & forward read
		GreedyDFSRangeSource *rFw_FwHalf = new GreedyDFSRangeSource(
		    &ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, true /*half-and-half*/, false);
		// Driver wrapper for rFw_Fw
		GreedyDFSRangeSourceDriver * drFw_Fw = new GreedyDFSRangeSourceDriver(
			*params, rFw_Fw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rFw_Bw
		GreedyDFSRangeSourceDriver * drFw_Bw = new GreedyDFSRangeSourceDriver(
			*params, rFw_Bw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rFw_Fw
		GreedyDFSRangeSourceDriver * drFw_FwHalf = new GreedyDFSRangeSourceDriver(
			*params, rFw_FwHalf, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			two_ ? PIN_TO_HI_HALF_EDGE : PIN_TO_BEGINNING,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec drVec;
		if(doFw_) {
			drVec.push_back(drFw_Fw);
			drVec.push_back(drFw_Bw);
			drVec.push_back(drFw_FwHalf);
		}

		// Source for ranges from forward index & reverse-complement read
		GreedyDFSRangeSource *rRc_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & reverse-complement read
		GreedyDFSRangeSource *rRc_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from forward index & reverse-complement read
		GreedyDFSRangeSource *rRc_FwHalf = new GreedyDFSRangeSource(
		    &ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, true /*half-and-half*/, false);
		// Driver wrapper for rRc_Fw
		GreedyDFSRangeSourceDriver * drRc_Fw = new GreedyDFSRangeSourceDriver(
			*params, rRc_Fw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rRc_Bw
		GreedyDFSRangeSourceDriver * drRc_Bw = new GreedyDFSRangeSourceDriver(
			*params, rRc_Bw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rRc_Fw
		GreedyDFSRangeSourceDriver * drRc_FwHalf = new GreedyDFSRangeSourceDriver(
			*params, rRc_FwHalf, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			two_ ? PIN_TO_HI_HALF_EDGE : PIN_TO_BEGINNING,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		if(doRc_) {
			drVec.push_back(drRc_Fw);
			drVec.push_back(drRc_Bw);
			drVec.push_back(drRc_FwHalf);
		}
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(seed_, drVec);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		return new UnpairedAlignerV2<GreedyDFSRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, rangeMode_, verbose_, seed_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	bool two_;
	bool doFw_;
	bool doRc_;
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
 * Concrete factory for constructing paired 2- or 3-mismatch aligners.
 */
class Paired23mmAlignerV1Factory : public AlignerFactory {
	typedef SingleRangeSourceDriver<GreedyDFSRangeSource> TRangeSrcDr;
	typedef ListRangeSourceDriver<GreedyDFSRangeSource> TListRangeSrcDr;
	typedef CostAwareRangeSourceDriver<GreedyDFSRangeSource> TCostAwareRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Paired23mmAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool two,
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
			two_(two),
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

		// Source for ranges from forward index & forward read
		GreedyDFSRangeSource *r1Fw_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true,  false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & forward read
		GreedyDFSRangeSource *r1Fw_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & forward read
		GreedyDFSRangeSource *r1Fw_BwHalf = new GreedyDFSRangeSource(
		     ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, true /*half-and-half*/, false);
		// Driver wrapper for rFw_Fw
		GreedyDFSRangeSourceDriver * dr1Fw_Fw = new GreedyDFSRangeSourceDriver(
			*params, r1Fw_Fw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rFw_Bw
		GreedyDFSRangeSourceDriver * dr1Fw_Bw = new GreedyDFSRangeSourceDriver(
			*params, r1Fw_Bw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rFw_Fw
		GreedyDFSRangeSourceDriver * dr1Fw_BwHalf = new GreedyDFSRangeSourceDriver(
			*params, r1Fw_BwHalf, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			two_ ? PIN_TO_HI_HALF_EDGE : PIN_TO_BEGINNING,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr1FwVec;
		dr1FwVec.push_back(dr1Fw_Bw);
		dr1FwVec.push_back(dr1Fw_Fw);
		dr1FwVec.push_back(dr1Fw_BwHalf);
		// Overall range source driver for the forward orientation of
		// the first mate
#if 0
		TListRangeSrcDr* dr1Fw = new TListRangeSrcDr(dr1FwVec);
#else
		TCostAwareRangeSrcDr* dr1Fw = new TCostAwareRangeSrcDr(seed_, dr1FwVec);
#endif

		// Source for ranges from forward index & reverse-complement read
		GreedyDFSRangeSource *r1Rc_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & reverse-complement read
		GreedyDFSRangeSource *r1Rc_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from forward index & reverse-complement read
		GreedyDFSRangeSource *r1Rc_FwHalf = new GreedyDFSRangeSource(
		    &ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, true /*half-and-half*/, false);

		// Driver wrapper for rRc_Fw
		GreedyDFSRangeSourceDriver * dr1Rc_Fw = new GreedyDFSRangeSourceDriver(
			*params, r1Rc_Fw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rRc_Bw
		GreedyDFSRangeSourceDriver * dr1Rc_Bw = new GreedyDFSRangeSourceDriver(
			*params, r1Rc_Bw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rRc_Fw
		GreedyDFSRangeSourceDriver * dr1Rc_FwHalf = new GreedyDFSRangeSourceDriver(
			*params, r1Rc_FwHalf, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			two_ ? PIN_TO_HI_HALF_EDGE : PIN_TO_BEGINNING,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr1RcVec;
		dr1RcVec.push_back(dr1Rc_Fw);
		dr1RcVec.push_back(dr1Rc_Bw);
		dr1RcVec.push_back(dr1Rc_FwHalf);
		// Overall range source driver for the reverse-comp orientation
		// of the first mate
#if 0
		TListRangeSrcDr* dr1Rc = new TListRangeSrcDr(dr1RcVec);
#else
		TCostAwareRangeSrcDr* dr1Rc = new TCostAwareRangeSrcDr(seed_, dr1RcVec);
#endif

		// Source for ranges from forward index & forward read
		GreedyDFSRangeSource *r2Fw_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & forward read
		GreedyDFSRangeSource *r2Fw_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & forward read
		GreedyDFSRangeSource *r2Fw_BwHalf = new GreedyDFSRangeSource(
		     ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, true /*half-and-half*/, false);
		// Driver wrapper for rFw_Fw
		GreedyDFSRangeSourceDriver * dr2Fw_Fw = new GreedyDFSRangeSourceDriver(
			*params, r2Fw_Fw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rFw_Bw
		GreedyDFSRangeSourceDriver * dr2Fw_Bw = new GreedyDFSRangeSourceDriver(
			*params, r2Fw_Bw, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rFw_Fw
		GreedyDFSRangeSourceDriver * dr2Fw_BwHalf = new GreedyDFSRangeSourceDriver(
			*params, r2Fw_BwHalf, true, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			two_ ? PIN_TO_HI_HALF_EDGE : PIN_TO_BEGINNING,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr2FwVec;
		dr2FwVec.push_back(dr2Fw_Bw);
		dr2FwVec.push_back(dr2Fw_Fw);
		dr2FwVec.push_back(dr2Fw_BwHalf);
		// Overall range source driver for the forward orientation of
		// the first mate
#if 0
		TListRangeSrcDr* dr2Fw = new TListRangeSrcDr(dr2FwVec);
#else
		TCostAwareRangeSrcDr* dr2Fw = new TCostAwareRangeSrcDr(seed_, dr2FwVec);
#endif

		// Source for ranges from forward index & reverse-complement read
		GreedyDFSRangeSource *r2Rc_Fw = new GreedyDFSRangeSource(
			&ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			true, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from mirror index & reverse-complement read
		GreedyDFSRangeSource *r2Rc_Bw = new GreedyDFSRangeSource(
			 ebwtBw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, false, false);
		// Source for ranges from forward index & reverse-complement read
		GreedyDFSRangeSource *r2Rc_FwHalf = new GreedyDFSRangeSource(
		    &ebwtFw_, *params, 0xffffffff /*max dist*/, BacktrackLimits(), 0 /*reportPartials*/,
			false, false, NULL, NULL, verbose_, seed_, &os_, false, true /*half-and-half*/, false);

		// Driver wrapper for rRc_Fw
		GreedyDFSRangeSourceDriver * dr2Rc_Fw = new GreedyDFSRangeSourceDriver(
			*params, r2Rc_Fw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rRc_Bw
		GreedyDFSRangeSourceDriver * dr2Rc_Bw = new GreedyDFSRangeSourceDriver(
			*params, r2Rc_Bw, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		// Driver wrapper for rRc_Fw
		GreedyDFSRangeSourceDriver * dr2Rc_FwHalf = new GreedyDFSRangeSourceDriver(
			*params, r2Rc_FwHalf, false, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			two_ ? PIN_TO_HI_HALF_EDGE : PIN_TO_BEGINNING,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, seed_);
		TRangeSrcDrPtrVec dr2RcVec;
		dr2RcVec.push_back(dr2Rc_Fw);
		dr2RcVec.push_back(dr2Rc_Bw);
		dr2RcVec.push_back(dr2Rc_FwHalf);
		// Overall range source driver for the reverse-comp orientation
		// of the first mate
#if 0
		TListRangeSrcDr* dr2Rc = new TListRangeSrcDr(dr2RcVec);
#else
		TCostAwareRangeSrcDr* dr2Rc = new TCostAwareRangeSrcDr(seed_, dr2RcVec);
#endif

		RefAligner<String<Dna5> >* refAligner;
		if(two_) {
			refAligner = new TwoMMRefAligner<String<Dna5> >(0);
		} else {
			refAligner = new ThreeMMRefAligner<String<Dna5> >(0);
		}

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
	bool two_;
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

#endif /* ALIGNER_23MM_H_ */
