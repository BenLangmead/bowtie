/*
 * aligner_seed_mm.h
 */

#ifndef ALIGNER_SEED_MM_H_
#define ALIGNER_SEED_MM_H_

#include <utility>
#include <vector>
#include "aligner.h"
#include "hit.h"
#include "row_chaser.h"
#include "range_chaser.h"

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class UnpairedSeedAlignerV1Factory : public AlignerFactory {

	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;

public:
	UnpairedSeedAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool doFw,
			bool doRc,
			uint32_t seedMms,
			uint32_t seedLen,
			int qualCutoff,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache* cacheFw,
			RangeCache* cacheBw,
			uint32_t cacheLimit,
			vector<String<Dna5> >& os,
			bool strandFix,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			doFw_(doFw), doRc_(doRc),
			seedMms_(seedMms),
			seedLen_(seedLen),
			qualCutoff_(qualCutoff),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			os_(os),
			strandFix_(strandFix),
			maqPenalty_(true),
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

		TRangeSrcDrPtrVec drVec;
		if(seedMms_ == 0) {
			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, true,  qualCutoff_, true, false, seed_, false, maqPenalty_);
			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, qualCutoff_, true, false, seed_, false, maqPenalty_);
			EbwtRangeSourceDriver * driverFw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, true, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				os_, verbose_, seed_, true);
			EbwtRangeSourceDriver * driverRc = new EbwtRangeSourceDriver(
				*params, rRc_Fw, false, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				os_, verbose_, seed_, true);
			if(doFw_) drVec.push_back(driverFw);
			if(doRc_) drVec.push_back(driverRc);

		} else if(seedMms_ == 1) {

			bool fw = true;
			bool mate1 = true;

			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, seed_, false, maqPenalty_);
			EbwtRangeSource *rFw_BwSeed = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, seed_, false, maqPenalty_);
			EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false, verbose_, seed_, false, maqPenalty_);

			EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, seed_, mate1);
			EbwtRangeSourceDriver * drFw_BwSeed = new EbwtRangeSourceDriver(
				*params, rFw_BwSeed, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, seed_, mate1);
			EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen, fw, true, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, seed_, mate1);
			EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, verbose_, mate1);

			fw = false;

			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, true,  verbose_, seed_, false, maqPenalty_);
			EbwtRangeSource *rRc_FwSeed = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, true,  verbose_, seed_, false, maqPenalty_);
			EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
				 ebwtBw_, fw, qualCutoff_, false, verbose_, seed_, false, maqPenalty_);

			EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
				*params, rRc_Fw, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, seed_, mate1);
			EbwtRangeSourceDriver * drRc_FwSeed = new EbwtRangeSourceDriver(
				*params, rRc_FwSeed, fw, false, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, seed_, mate1);
			EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen, fw, true, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, seed_, mate1);
			EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, verbose_, true);

			if(doFw_) {
				drVec.push_back(drFw_Bw);
				drVec.push_back(drFw_Seed);
			}
			if(doRc_) {
				drVec.push_back(drRc_Fw);
				drVec.push_back(drRc_Seed);
			}
		} else {
			cerr << "Unsupported --stateful mode: " << seedMms_ << endl;
		}
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(seed_, strandFix_, drVec);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		return new UnpairedAlignerV2<EbwtRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, rangeMode_, verbose_,
			seed_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	bool doFw_;
	bool doRc_;
	const uint32_t seedMms_;
	const uint32_t seedLen_;
	const int qualCutoff_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	vector<String<Dna5> >& os_;
	bool strandFix_;
	bool maqPenalty_;
	bool rangeMode_;
	bool verbose_;
	uint32_t seed_;
};

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class PairedSeedAlignerV1Factory : public AlignerFactory {
public:
	PairedSeedAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool doFw,
			bool doRc,
			uint32_t seedMms,
			uint32_t seedLen,
			int qualCutoff,
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
			bool strandFix,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			doFw_(doFw),
			doRc_(doRc),
			seedMms_(seedMms),
			seedLen_(seedLen),
			qualCutoff_(qualCutoff),
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
			strandFix_(strandFix),
			rangeMode_(rangeMode),
			verbose_(verbose),
			seed_(seed)
	{
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

		EbwtRangeSource *r1Fw_Bw = new EbwtRangeSource(
			 ebwtBw_, true,  qualCutoff_, true, false, seed_, false, false);
		EbwtRangeSource *r1Rc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, qualCutoff_, true, false, seed_, false, false);

		EbwtRangeSourceDriver * driver1Fw = new EbwtRangeSourceDriver(
			*params, r1Fw_Bw, true, false, sink_, sinkPt,
			seedLen_,   // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			os_, verbose_, seed_, true);
		EbwtRangeSourceDriver * driver1Rc = new EbwtRangeSourceDriver(
			*params, r1Rc_Fw, false, false, sink_, sinkPt,
			seedLen_,   // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			os_, verbose_, seed_, true);

		EbwtRangeSource *r2Fw_Bw = new EbwtRangeSource(
			 ebwtBw_, true,  qualCutoff_, true, false, seed_, false, false);
		EbwtRangeSource *r2Rc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, qualCutoff_, true, false, seed_, false, false);

		EbwtRangeSourceDriver * driver2Fw = new EbwtRangeSourceDriver(
			*params, r2Fw_Bw, true, false, sink_, sinkPt,
			seedLen_,   // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			os_, verbose_, seed_, false);
		EbwtRangeSourceDriver * driver2Rc = new EbwtRangeSourceDriver(
			*params, r2Rc_Fw, false, false, sink_, sinkPt,
			seedLen_,   // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			PIN_TO_SEED_EDGE, // "
			os_, verbose_, seed_, false);

		RefAligner<String<Dna5> >* refAligner = new ExactRefAligner<String<Dna5> >(0);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		return new PairedBWAlignerV1<EbwtRangeSource>(
			params,
			driver1Fw, driver1Rc, driver2Fw, driver2Rc, refAligner,
			rchase, sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
			peInner_, peOuter_, dontReconcile_, symCeil_, mixedThresh_,
			mixedAttemptLim_, refs_, rangeMode_, verbose_, seed_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	const bool doFw_;
	const bool doRc_;
	const uint32_t seedMms_;
	const uint32_t seedLen_;
	const int qualCutoff_;
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
	const bool strandFix_;
	const bool rangeMode_;
	const bool verbose_;
	const uint32_t seed_;
};

#endif /* ALIGNER_SEED_MM_H_ */
