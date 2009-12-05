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
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
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
			two_(two),
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
		assert(ebwtBw != NULL);
		assert(ebwtBw->isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {

		HitSinkPerThread* sinkPt = sinkPtFactory_.create();
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_);

		const bool seeded = false;

		EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, true, verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rFw_Fw = new EbwtRangeSource(
			&ebwtFw_, true, 0xffffffff, false, verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
			 ebwtBw_, true, 0xffffffff, false, verbose_, quiet_, 2,  seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rFw_FwHalf = NULL;
		if(!two_) {
			rFw_FwHalf = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, verbose_, quiet_, 3,  seeded, maqPenalty_, qualOrder_);
		}

		// Driver wrapper for rFw_Bw
		EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
			*params, rFw_Bw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, quiet_, true, pool_, NULL);
		// Driver wrapper for rFw_Fw
		EbwtRangeSourceDriver * drFw_Fw = new EbwtRangeSourceDriver(
			*params, rFw_Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, quiet_, true, pool_, NULL);
		// Driver wrapper for rFw_Fw
		EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
			*params, rFw_BwHalf, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			PIN_TO_HI_HALF_EDGE,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, quiet_, true, pool_, NULL);
		// Driver wrapper for rFw_Fw
		EbwtRangeSourceDriver * drFw_FwHalf = NULL;
		if(!two_) {
			drFw_FwHalf = new EbwtRangeSourceDriver(
				*params, rFw_FwHalf, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
		}

		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(doFw_) {
			drVec->push_back(drFw_Bw);
			drVec->push_back(drFw_Fw);
			drVec->push_back(drFw_BwHalf);
			if(!two_) {
				drVec->push_back(drFw_FwHalf);
			}
		}

		EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, true,  verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rRc_Bw = new EbwtRangeSource(
			 ebwtBw_, false, 0xffffffff, false, verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
			&ebwtFw_, false, 0xffffffff, false, verbose_, quiet_, 2,  seeded, maqPenalty_, qualOrder_);
		EbwtRangeSource *rRc_BwHalf = NULL;
		if(!two_) {
			rRc_BwHalf = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, verbose_, quiet_, 3,  seeded, maqPenalty_, qualOrder_);
		}

		// Driver wrapper for rRc_Fw
		EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
			*params, rRc_Fw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, quiet_, true, pool_, NULL);
		// Driver wrapper for rRc_Bw
		EbwtRangeSourceDriver * drRc_Bw = new EbwtRangeSourceDriver(
			*params, rRc_Bw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			false,      // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
			PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, quiet_, true, pool_, NULL);
		// Driver wrapper for rRc_Fw
		EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
			*params, rRc_FwHalf, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
			0,          // seedLen (0 = whole read is seed)
			true,       // nudgeLeft (true for Fw index, false for Bw)
			PIN_TO_BEGINNING,    // nothing's unrevisitable
			PIN_TO_HI_HALF_EDGE,
			two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
			PIN_TO_LEN,
			os_, verbose_, quiet_, true, pool_, NULL);
		EbwtRangeSourceDriver * drRc_BwHalf = NULL;
		if(!two_) {
			drRc_BwHalf = new EbwtRangeSourceDriver(
				*params, rRc_BwHalf, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
		}
		if(doRc_) {
			drVec->push_back(drRc_Fw);
			drVec->push_back(drRc_Bw);
			drVec->push_back(drRc_FwHalf);
			if(!two_) {
				drVec->push_back(drRc_BwHalf);
			}
		}
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
	bool two_;
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
	const bool maqPenalty_;
	const bool qualOrder_;
	const bool strandFix_;
	const bool rangeMode_;
	const bool verbose_;
	const bool quiet_;
	uint32_t seed_;
};

/**
 * Concrete factory for constructing paired 2- or 3-mismatch aligners.
 */
class Paired23mmAlignerV1Factory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef ListRangeSourceDriver<EbwtRangeSource> TListRangeSrcDr;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
public:
	Paired23mmAlignerV1Factory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool color,
			bool doFw,
			bool doRc,
			bool v1,
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
			ChunkPool *pool,
			BitPairReference* refs,
			vector<String<Dna5> >& os,
			bool reportSe,
			bool maqPenalty,
			bool qualOrder,
			bool strandFix,
			bool rangeMode,
			bool verbose,
			bool quiet,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			color_(color),
			doFw_(doFw),
			doRc_(doRc),
			v1_(v1),
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
		assert(ebwtBw != NULL);
		assert(ebwtFw.isInMemory());
		assert(ebwtBw->isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2);
		HitSinkPerThread* sinkPtSe1 = NULL, * sinkPtSe2 = NULL;
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_);
		EbwtSearchParams<String<Dna> >* paramsSe1 = NULL, * paramsSe2 = NULL;
		if(reportSe_) {
			sinkPtSe1 = sinkPtFactory_.create();
			sinkPtSe2 = sinkPtFactory_.create();
			paramsSe1 =
				new EbwtSearchParams<String<Dna> >(*sinkPtSe1, os_);
			paramsSe2 =
				new EbwtSearchParams<String<Dna> >(*sinkPtSe2, os_);
		}

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

		TRangeSrcDrPtrVec *dr1FwVec = new TRangeSrcDrPtrVec();

		if(do1Fw) {
			EbwtRangeSource *r1Fw_Bw = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, true,  verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r1Fw_Fw = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r1Fw_BwHalf = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, false, verbose_, quiet_, 2, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r1Fw_FwHalf = two_ ? NULL : new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, verbose_, quiet_, 3, seeded, maqPenalty_, qualOrder_);

			// Driver wrapper for rFw_Bw
			EbwtRangeSourceDriver * dr1Fw_Bw = new EbwtRangeSourceDriver(
				*params, r1Fw_Bw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr1Fw_Fw = new EbwtRangeSourceDriver(
				*params, r1Fw_Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr1Fw_BwHalf = new EbwtRangeSourceDriver(
				*params, r1Fw_BwHalf, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
			dr1FwVec->push_back(dr1Fw_Bw);
			dr1FwVec->push_back(dr1Fw_Fw);
			dr1FwVec->push_back(dr1Fw_BwHalf);
			if(!two_) {
				// Driver wrapper for rFw_Fw
				EbwtRangeSourceDriver * dr1Fw_FwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r1Fw_FwHalf, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, verbose_, quiet_, true, pool_, NULL);
				dr1FwVec->push_back(dr1Fw_FwHalf);
			}
		}

		TRangeSrcDrPtrVec *dr1RcVec;
		if(v1_) {
			dr1RcVec = new TRangeSrcDrPtrVec();
		} else {
			dr1RcVec = dr1FwVec;
		}

		if(do1Rc) {
			EbwtRangeSource *r1Rc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true,  verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r1Rc_Bw = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r1Rc_FwHalf = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, false, verbose_, quiet_, 2, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r1Rc_BwHalf = two_ ? NULL : new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, verbose_, quiet_, 3, seeded, maqPenalty_, qualOrder_);

			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr1Rc_Fw = new EbwtRangeSourceDriver(
				*params, r1Rc_Fw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
			// Driver wrapper for rRc_Bw
			EbwtRangeSourceDriver * dr1Rc_Bw = new EbwtRangeSourceDriver(
				*params, r1Rc_Bw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr1Rc_FwHalf = new EbwtRangeSourceDriver(
				*params, r1Rc_FwHalf, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, true, pool_, NULL);
			dr1RcVec->push_back(dr1Rc_Fw);
			dr1RcVec->push_back(dr1Rc_Bw);
			dr1RcVec->push_back(dr1Rc_FwHalf);
			if(!two_) {
				// Driver wrapper for rRc_Bw
				EbwtRangeSourceDriver * dr1Rc_BwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r1Rc_BwHalf, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, verbose_, quiet_, true, pool_, NULL);
				dr1RcVec->push_back(dr1Rc_BwHalf);
			}
		}

		TRangeSrcDrPtrVec *dr2FwVec;
		if(v1_) {
			dr2FwVec = new TRangeSrcDrPtrVec();
		} else {
			dr2FwVec = dr1FwVec;
		}

		if(do2Fw) {
			EbwtRangeSource *r2Fw_Bw = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, true,  verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r2Fw_Fw = new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r2Fw_BwHalf = new EbwtRangeSource(
				 ebwtBw_, true, 0xffffffff, false, verbose_, quiet_, 2, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r2Fw_FwHalf = two_ ? NULL : new EbwtRangeSource(
				&ebwtFw_, true, 0xffffffff, false, verbose_, quiet_, 3, seeded, maqPenalty_, qualOrder_);

			// Driver wrapper for rFw_Bw
			EbwtRangeSourceDriver * dr2Fw_Bw = new EbwtRangeSourceDriver(
				*params, r2Fw_Bw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, false, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr2Fw_Fw = new EbwtRangeSourceDriver(
				*params, r2Fw_Fw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, false, pool_, NULL);
			// Driver wrapper for rFw_Fw
			EbwtRangeSourceDriver * dr2Fw_BwHalf = new EbwtRangeSourceDriver(
				*params, r2Fw_BwHalf, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, false, pool_, NULL);
			dr2FwVec->push_back(dr2Fw_Bw);
			dr2FwVec->push_back(dr2Fw_Fw);
			dr2FwVec->push_back(dr2Fw_BwHalf);
			if(!two_) {
				EbwtRangeSourceDriver * dr2Fw_FwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r2Fw_FwHalf, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, verbose_, quiet_, false, pool_, NULL);
				dr2FwVec->push_back(dr2Fw_FwHalf);
			}
		}

		TRangeSrcDrPtrVec *dr2RcVec;
		if(v1_) {
			dr2RcVec = new TRangeSrcDrPtrVec();
		} else {
			dr2RcVec = dr1FwVec;
		}

		if(do2Rc) {
			EbwtRangeSource *r2Rc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, true,  verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r2Rc_Bw = new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, verbose_, quiet_, 0, seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r2Rc_FwHalf = new EbwtRangeSource(
				&ebwtFw_, false, 0xffffffff, false, verbose_, quiet_, 2,  seeded, maqPenalty_, qualOrder_);
			EbwtRangeSource *r2Rc_BwHalf = two_ ? NULL : new EbwtRangeSource(
				 ebwtBw_, false, 0xffffffff, false, verbose_, quiet_, 3,  seeded, maqPenalty_, qualOrder_);

			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr2Rc_Fw = new EbwtRangeSourceDriver(
				*params, r2Rc_Fw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, false, pool_, NULL);
			// Driver wrapper for rRc_Bw
			EbwtRangeSourceDriver * dr2Rc_Bw = new EbwtRangeSourceDriver(
				*params, r2Rc_Bw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				false,      // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_HI_HALF_EDGE, // right half is unrevisitable
				PIN_TO_HI_HALF_EDGE, // trumped by 0-mm
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, false, pool_, NULL);
			// Driver wrapper for rRc_Fw
			EbwtRangeSourceDriver * dr2Rc_FwHalf = new EbwtRangeSourceDriver(
				*params, r2Rc_FwHalf, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				0,          // seedLen (0 = whole read is seed)
				true,       // nudgeLeft (true for Fw index, false for Bw)
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				two_ ? PIN_TO_LEN : PIN_TO_HI_HALF_EDGE,
				PIN_TO_LEN,
				os_, verbose_, quiet_, false, pool_, NULL);
			dr2RcVec->push_back(dr2Rc_Fw);
			dr2RcVec->push_back(dr2Rc_Bw);
			dr2RcVec->push_back(dr2Rc_FwHalf);
			if(!two_) {
				EbwtRangeSourceDriver * dr2Rc_BwHalf = two_ ? NULL : new EbwtRangeSourceDriver(
					*params, r2Rc_BwHalf, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					0,          // seedLen (0 = whole read is seed)
					false,       // nudgeLeft (true for Fw index, false for Bw)
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_BEGINNING,
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_LEN,
					os_, verbose_, quiet_, false, pool_, NULL);
				dr2RcVec->push_back(dr2Rc_BwHalf);
			}
		}

		RefAligner<String<Dna5> >* refAligner;
		if(two_) {
			refAligner = new TwoMMRefAligner<String<Dna5> >(color_, verbose_, quiet_);
		} else {
			refAligner = new ThreeMMRefAligner<String<Dna5> >(color_, verbose_, quiet_);
		}

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		if(v1_) {
			PairedBWAlignerV1<EbwtRangeSource> *al = new PairedBWAlignerV1<EbwtRangeSource>(
				params,
				new TCostAwareRangeSrcDr(strandFix_, dr1FwVec, verbose_, quiet_, false),
				new TCostAwareRangeSrcDr(strandFix_, dr1RcVec, verbose_, quiet_, false),
				new TCostAwareRangeSrcDr(strandFix_, dr2FwVec, verbose_, quiet_, false),
				new TCostAwareRangeSrcDr(strandFix_, dr2RcVec, verbose_, quiet_, false),
				refAligner, rchase,
				sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
				peInner_, peOuter_, dontReconcile_, symCeil_, mixedThresh_,
				mixedAttemptLim_, refs_, rangeMode_, verbose_,
				quiet_, INT_MAX, pool_, NULL);
			delete dr1FwVec;
			delete dr1RcVec;
			delete dr2FwVec;
			delete dr2RcVec;
			return al;
		} else {
			PairedBWAlignerV2<EbwtRangeSource>* al = new PairedBWAlignerV2<EbwtRangeSource>(
				params, paramsSe1, paramsSe2,
				new TCostAwareRangeSrcDr(strandFix_, dr1FwVec, verbose_, quiet_, true),
				refAligner, rchase,
				sink_, sinkPtFactory_,
				sinkPt, sinkPtSe1, sinkPtSe2,
				mate1fw_, mate2fw_,
				peInner_, peOuter_,
				mixedAttemptLim_, refs_, rangeMode_,
				verbose_, quiet_, INT_MAX, pool_, NULL);
			delete dr1FwVec;
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

#endif /* ALIGNER_23MM_H_ */
