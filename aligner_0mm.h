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

	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;

public:
	UnpairedExactAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool doFw,
			bool doRc,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache* cacheFw,
			RangeCache* cacheBw,
			uint32_t cacheLimit,
			ChunkPool *pool,
			BitPairReference* refs,
			vector<String<Dna5> >& os,
			bool maqPenalty,
			bool qualOrder,
			bool strandFix,
			bool rangeMode,
			bool verbose,
			bool quiet,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			doFw_(doFw), doRc_(doRc),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			pool_(pool),
			refs_(refs),
			os_(os),
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			strandFix_(strandFix),
			rangeMode_(rangeMode),
			verbose_(verbose),
			quiet_(quiet),
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
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_, true, true);

		const int halfAndHalf = 0;
		const bool seeded = false;

		EbwtRangeSource *rFw = new EbwtRangeSource(
			&ebwtFw_, true,  0xffffffff, true, verbose_, quiet_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rRc = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true, verbose_, quiet_, halfAndHalf, seeded, maqPenalty_, qualOrder_);

		EbwtRangeSourceDriver * driverFw = new EbwtRangeSourceDriver(
			*params, rFw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, quiet_, true, pool_, NULL);
		EbwtRangeSourceDriver * driverRc = new EbwtRangeSourceDriver(
			*params, rRc, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (not applicable)
			PIN_TO_LEN, // whole alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, quiet_, true, pool_, NULL);
		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(doFw_) drVec->push_back(driverFw);
		if(doRc_) drVec->push_back(driverRc);
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(strandFix_, drVec, verbose_, quiet_, false);
		delete drVec;

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		return new UnpairedAlignerV2<EbwtRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, refs_,
			rangeMode_, verbose_, quiet_, INT_MAX, pool_, NULL, NULL);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	bool doFw_;
	bool doRc_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	ChunkPool *pool_;
	BitPairReference* refs_;
	vector<String<Dna5> >& os_;
	bool maqPenalty_;
	bool qualOrder_;
	bool strandFix_;
	bool rangeMode_;
	bool verbose_;
	bool quiet_;
	uint32_t seed_;
};

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class PairedExactAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	PairedExactAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool color,
			bool doFw,
			bool doRc,
			bool v1,
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
			ChunkPool *pool,
			BitPairReference* refs,
			vector<String<Dna5> >& os,
			bool reportSe,
			bool maqPenalty,
			bool strandFix,
			bool qualOrder,
			bool rangeMode,
			bool verbose,
			bool quiet,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			color_(color),
			doFw_(doFw),
			doRc_(doRc),
			v1_(v1),
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
			pool_(pool),
			refs_(refs), os_(os),
			reportSe_(reportSe),
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			strandFix_(strandFix),
			rangeMode_(rangeMode),
			verbose_(verbose),
			quiet_(quiet),
			seed_(seed)
	{
		assert(ebwtFw.isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2);
		HitSinkPerThread* sinkPtSe1 = NULL, * sinkPtSe2 = NULL;
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_, true, true);
		EbwtSearchParams<String<Dna> >* paramsSe1 = NULL, * paramsSe2 = NULL;
		if(reportSe_) {
			sinkPtSe1 = sinkPtFactory_.create();
			sinkPtSe2 = sinkPtFactory_.create();
			paramsSe1 =
				new EbwtSearchParams<String<Dna> >(*sinkPtSe1, os_, true, true);
			paramsSe2 =
				new EbwtSearchParams<String<Dna> >(*sinkPtSe2, os_, true, true);
		}

		const int halfAndHalf = 0;
		const bool seeded = false;

		bool do1Fw = true;
		bool do1Rc = true;
		bool do2Fw = true;
		bool do2Rc = true;
		if(!doFw_) {
			if(mate1fw_) do1Fw = false;
			else         do1Rc = false;
			if(mate2fw_) do2Fw = false;
			else         do2Rc = false;
		}
		if(!doRc_) {
			if(mate1fw_) do1Rc = false;
			else         do1Fw = false;
			if(mate2fw_) do2Rc = false;
			else         do2Fw = false;
		}

		EbwtRangeSource *r1Fw = NULL;
		EbwtRangeSource *r1Rc = NULL;
		TRangeSrcDr * driver1Fw = NULL;
		TRangeSrcDr * driver1Rc = NULL;
		EbwtRangeSource *r2Fw = NULL;
		EbwtRangeSource *r2Rc = NULL;
		TRangeSrcDr * driver2Fw = NULL;
		TRangeSrcDr * driver2Rc = NULL;
		if(do1Fw) {
			r1Fw = new EbwtRangeSource(
				&ebwtFw_, true,  0xffffffff, true, verbose_, quiet_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
			driver1Fw = new EbwtRangeSourceDriver(
				*params, r1Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, verbose_, quiet_, true, pool_, NULL);
		}
		if(do2Fw) {
			r2Fw = new EbwtRangeSource(
				&ebwtFw_, true,  0xffffffff, true, verbose_, quiet_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
			driver2Fw = new EbwtRangeSourceDriver(
				*params, r2Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, verbose_, quiet_, false, pool_, NULL);
		}
		if(do1Rc) {
			r1Rc = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true, verbose_, quiet_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
			driver1Rc = new EbwtRangeSourceDriver(
				*params, r1Rc, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, verbose_, quiet_, true, pool_, NULL);
		}
		if(do2Rc) {
			r2Rc = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true, verbose_, quiet_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
			driver2Rc = new EbwtRangeSourceDriver(
				*params, r2Rc, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_LEN, // whole alignment is unrevisitable
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				PIN_TO_LEN, // "
				os_, verbose_, quiet_, false, pool_, NULL);
		}

		RefAligner<String<Dna5> >* refAligner
			= new ExactRefAligner<String<Dna5> >(color_, verbose_, quiet_);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		if(v1_) {
			PairedBWAlignerV1<EbwtRangeSource>* al = new PairedBWAlignerV1<EbwtRangeSource>(
				params,
				driver1Fw == NULL ? (new StubRangeSourceDriver<EbwtRangeSource>()) : driver1Fw,
				driver1Rc == NULL ? (new StubRangeSourceDriver<EbwtRangeSource>()) : driver1Rc,
				driver2Fw == NULL ? (new StubRangeSourceDriver<EbwtRangeSource>()) : driver2Fw,
				driver2Rc == NULL ? (new StubRangeSourceDriver<EbwtRangeSource>()) : driver2Rc,
				refAligner,
				rchase, sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
				peInner_, peOuter_, dontReconcile_, symCeil_, mixedThresh_,
				mixedAttemptLim_, refs_, rangeMode_, verbose_,
				quiet_, INT_MAX, pool_, NULL);
			return al;
		} else {
			TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
			if(driver1Fw != NULL) drVec->push_back(driver1Fw);
			if(driver1Rc != NULL) drVec->push_back(driver1Rc);
			if(driver2Fw != NULL) drVec->push_back(driver2Fw);
			if(driver2Rc != NULL) drVec->push_back(driver2Rc);
			PairedBWAlignerV2<EbwtRangeSource>* al = new PairedBWAlignerV2<EbwtRangeSource>(
				params, paramsSe1, paramsSe2,
				new TCostAwareRangeSrcDr(strandFix_, drVec, verbose_, quiet_, true),
				refAligner,
				rchase, sink_, sinkPtFactory_, sinkPt,
				sinkPtSe1, sinkPtSe2, mate1fw_, mate2fw_,
				peInner_, peOuter_,
				mixedAttemptLim_, refs_, rangeMode_,
				verbose_, quiet_, INT_MAX, pool_, NULL);
			delete drVec;
			return al;
		}
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	bool color_;
	bool doFw_;
	bool doRc_;
	bool v1_;
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
	ChunkPool *pool_;
	BitPairReference* refs_;
	vector<String<Dna5> >& os_;
	const bool reportSe_;
	const bool maqPenalty_;
	const bool qualOrder_;
	const bool strandFix_;
	const bool rangeMode_;
	const bool verbose_;
	const bool quiet_;
	const uint32_t seed_;
};

#endif /* ALIGNER_0MM_H_ */
