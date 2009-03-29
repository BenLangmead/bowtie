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
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef ListRangeSourceDriver<EbwtRangeSource> TListRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Unpaired1mmAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool doFw,
			bool doRc,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			RangeCache *cacheFw,
			RangeCache *cacheBw,
			uint32_t cacheLimit,
			vector<String<Dna5> >& os,
			bool maqPenalty,
			bool qualOrder,
			bool strandFix,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			doFw_(doFw),
			doRc_(doRc),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			os_(os),
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			strandFix_(strandFix),
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
		AllocOnlyPool<Branch> *bpool = new AllocOnlyPool<Branch>((1 << 20), "branch");
		RangeStatePool *rpool = new RangeStatePool(1 * 1024 * 1024);
		AllocOnlyPool<Edit> *epool = new AllocOnlyPool<Edit>((1 << 20), "edit");

		const bool halfAndHalf = false;
		const bool seeded = false;

		EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, true,  false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rFw_Fw = new EbwtRangeSource(
			&ebwtFw_, true, 0xffffffff, false, false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);

		EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
			*params, rFw_Bw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);
		//
		EbwtRangeSourceDriver * drFw_Fw = new EbwtRangeSourceDriver(
			*params, rFw_Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);
		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(doFw_) {
			drVec->push_back(drFw_Bw);
			drVec->push_back(drFw_Fw);
		}

		EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true,  false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rRc_Bw = new EbwtRangeSource(
			 ebwtBw_, false, 0xffffffff, false, false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);

		EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
			*params, rRc_Fw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);
		//
		EbwtRangeSourceDriver * drRc_Bw = new EbwtRangeSourceDriver(
			*params, rRc_Bw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);
		if(doRc_) {
			drVec->push_back(drRc_Fw);
			drVec->push_back(drRc_Bw);
		}
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(seed_, strandFix_, drVec, verbose_);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

		// Set up the aligner
		return new UnpairedAlignerV2<EbwtRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, rangeMode_, verbose_,
			seed_, INT_MAX, bpool, rpool, epool, NULL, NULL);
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
	vector<String<Dna5> >& os_;
	const bool maqPenalty_;
	const bool qualOrder_;
	bool strandFix_;
	bool rangeMode_;
	bool verbose_;
	uint32_t seed_;
};

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class Paired1mmAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef ListRangeSourceDriver<EbwtRangeSource> TListRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
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
			bool maqPenalty,
			bool qualOrder,
			bool strandFix,
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
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			strandFix_(strandFix),
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
		AllocOnlyPool<Branch> *bpool = new AllocOnlyPool<Branch>((1 << 20), "branch");
		RangeStatePool *rpool = new RangeStatePool(1 * 1024 * 1024);
		AllocOnlyPool<Edit> *epool = new AllocOnlyPool<Edit>((1 << 20), "edit");

		const bool halfAndHalf = false;
		const bool seeded = false;

		EbwtRangeSource *r1Fw_Bw = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, true,  false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *r1Fw_Fw = new EbwtRangeSource(
			&ebwtFw_, true, 0xffffffff, false, false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);

		EbwtRangeSourceDriver * dr1Fw_Bw = new EbwtRangeSourceDriver(
			*params, r1Fw_Bw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);
		EbwtRangeSourceDriver * dr1Fw_Fw = new EbwtRangeSourceDriver(
			*params, r1Fw_Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);

#if 1
		TRangeSrcDrPtrVec dr1FwVec;
		dr1FwVec.push_back(dr1Fw_Bw);
		dr1FwVec.push_back(dr1Fw_Fw);
		// Overall range source driver for the forward orientation of
		// the first mate
		TListRangeSrcDr* dr1Fw = new TListRangeSrcDr(dr1FwVec);
#else
		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		drVec->push_back(dr1Fw_Bw);
		drVec->push_back(dr1Fw_Fw);
#endif

		EbwtRangeSource *r1Rc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true,  false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *r1Rc_Bw = new EbwtRangeSource(
			 ebwtBw_, false, 0xffffffff, false, false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);

		EbwtRangeSourceDriver * dr1Rc_Fw = new EbwtRangeSourceDriver(
			*params, r1Rc_Fw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);
		EbwtRangeSourceDriver * dr1Rc_Bw = new EbwtRangeSourceDriver(
			*params, r1Rc_Bw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, true, bpool, rpool, epool, NULL);
#if 1
		TRangeSrcDrPtrVec dr1RcVec;
		dr1RcVec.push_back(dr1Rc_Fw);
		dr1RcVec.push_back(dr1Rc_Bw);
		// Overall range source driver for the reverse-comp orientation
		// of the first mate
		TListRangeSrcDr* dr1Rc = new TListRangeSrcDr(dr1RcVec);
#else
		drVec->push_back(dr1Rc_Fw);
		drVec->push_back(dr1Rc_Bw);
#endif

		EbwtRangeSource *r2Fw_Bw = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, true,  false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *r2Fw_Fw = new EbwtRangeSource(
			&ebwtFw_, true, 0xffffffff, false, false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);

		EbwtRangeSourceDriver * dr2Fw_Bw = new EbwtRangeSourceDriver(
			*params, r2Fw_Bw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, false, bpool, rpool, epool, NULL);
		EbwtRangeSourceDriver * dr2Fw_Fw = new EbwtRangeSourceDriver(
			*params, r2Fw_Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, false, bpool, rpool, epool, NULL);
#if 1
		TRangeSrcDrPtrVec dr2FwVec;
		dr2FwVec.push_back(dr2Fw_Bw);
		dr2FwVec.push_back(dr2Fw_Fw);
		// Overall range source driver for the forward orientation of
		// the first mate
		TListRangeSrcDr* dr2Fw = new TListRangeSrcDr(dr2FwVec);
#else
		drVec->push_back(dr2Fw_Bw);
		drVec->push_back(dr2Fw_Fw);
#endif

		EbwtRangeSource *r2Rc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true,  false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *r2Rc_Bw = new EbwtRangeSource(
			 ebwtBw_, false, 0xffffffff, false, false, seed_, halfAndHalf, seeded, maqPenalty_, qualOrder_);

		EbwtRangeSourceDriver * dr2Rc_Fw = new EbwtRangeSourceDriver(
			*params, r2Rc_Fw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right-hand half alignment is unrevisitable
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, false, bpool, rpool, epool, NULL);
		EbwtRangeSourceDriver * dr2Rc_Bw = new EbwtRangeSourceDriver(
			*params, r2Rc_Bw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_LEN, // allow 1 mismatch in rest of read
			PIN_TO_LEN, // "
			PIN_TO_LEN, // "
			os_, verbose_, seed_, false, bpool, rpool, epool, NULL);
#if 1
		TRangeSrcDrPtrVec dr2RcVec;
		dr2RcVec.push_back(dr2Rc_Fw);
		dr2RcVec.push_back(dr2Rc_Bw);
		// Overall range source driver for the reverse-comp orientation
		// of the first mate
		TListRangeSrcDr* dr2Rc = new TListRangeSrcDr(dr2RcVec);
#else
		drVec->push_back(dr2Rc_Fw);
		drVec->push_back(dr2Rc_Bw);
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(seed_, strandFix_, drVec, verbose_);
#endif

		RefAligner<String<Dna5> >* refAligner = new OneMMRefAligner<String<Dna5> >(0);

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(seed_, cacheLimit_, cacheFw_, cacheBw_);

#if 1
		return new PairedBWAlignerV1<EbwtRangeSource>(
			params, dr1Fw, dr1Rc, dr2Fw, dr2Rc, refAligner, rchase,
			sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
			peInner_, peOuter_, dontReconcile_, symCeil_, mixedThresh_,
			mixedAttemptLim_, refs_, rangeMode_, verbose_, seed_,
			INT_MAX, bpool, rpool, epool, NULL);
#else
		return new PairedBWAlignerV2<EbwtRangeSource>(
			params, dr, refAligner, rchase,
			sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
			peInner_, peOuter_,
			mixedAttemptLim_, refs_, rangeMode_, verbose_, seed_);
#endif
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
	const bool maqPenalty_;
	const bool qualOrder_;
	const bool strandFix_;
	const bool rangeMode_;
	const bool verbose_;
	const uint32_t seed_;
};

#endif /* ALIGNER_1MM_H_ */
