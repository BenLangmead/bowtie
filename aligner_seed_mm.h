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
#include "aligner_metrics.h"

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class UnpairedSeedAlignerFactory : public AlignerFactory {

	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;

public:
	UnpairedSeedAlignerFactory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool doFw,
			bool doRc,
			uint32_t seedMms,
			uint32_t seedLen,
			int qualCutoff,
			int maxBts,
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
			uint32_t seed,
			AlignerMetrics *metrics) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			doFw_(doFw), doRc_(doRc),
			seedMms_(seedMms),
			seedLen_(seedLen),
			qualCutoff_(qualCutoff),
			maxBts_(maxBts),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			cacheFw_(cacheFw),
			cacheBw_(cacheBw),
			cacheLimit_(cacheLimit),
			pool_(pool),
			refs_(refs),
			os_(os),
			strandFix_(strandFix),
			maqPenalty_(maqPenalty),
			qualOrder_(qualOrder),
			rangeMode_(rangeMode),
			verbose_(verbose),
			quiet_(quiet),
			metrics_(metrics)
	{
		assert(ebwtFw.isInMemory());
	}

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		HitSinkPerThread* sinkPt = sinkPtFactory_.create();
		EbwtSearchParams<String<Dna> >* params =
			new EbwtSearchParams<String<Dna> >(*sinkPt, os_);
		int *btCnt = new int[1];
		*btCnt = maxBts_;

		TRangeSrcDrPtrVec *drVec = new TRangeSrcDrPtrVec();
		if(seedMms_ == 0) {
			const int halfAndHalf = 0;
			bool mate1 = true;
			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, true,  qualCutoff_, true, verbose_, quiet_,
				 halfAndHalf, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, false, qualCutoff_, true, verbose_, quiet_,
				halfAndHalf, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSourceDriver * driverFw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, true, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriver * driverRc = new EbwtRangeSourceDriver(
				*params, rRc_Fw, false, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft (not applicable)
				PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				PIN_TO_SEED_EDGE, // "
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			if(doFw_) drVec->push_back(driverFw);
			if(doRc_) drVec->push_back(driverRc);

		} else if(seedMms_ == 1) {
			const int halfAndHalf = 0;
			bool fw = true;
			bool mate1 = true;

			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				 halfAndHalf, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				 halfAndHalf, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_,
				 halfAndHalf, true,  maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, fw, false, maqPenalty_, qualOrder_,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed, fw, false, maqPenalty_, qualOrder_,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen, fw, true, maqPenalty_,
				qualOrder_, sink_, sinkPt, seedLen_,
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);

			fw = false;

			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, true,  verbose_, quiet_,
				halfAndHalf, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
				&ebwtFw_, fw, qualCutoff_, true,  verbose_, quiet_,
				halfAndHalf, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
				 ebwtBw_, fw, qualCutoff_, false, verbose_, quiet_,
				halfAndHalf, true,  maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
				*params, rRc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
			EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, verbose_, quiet_, true);

			if(doFw_) {
				drVec->push_back(drFw_Bw);
				drVec->push_back(drFw_Seed);
			}
			if(doRc_) {
				drVec->push_back(drRc_Fw);
				drVec->push_back(drRc_Seed);
			}
		} else if(seedMms_ == 2) {

			bool fw = true;
			bool mate1 = true;

			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_, 0, true,  maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_, 2,  false, maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
			EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
				*params, rFw_BwHalf, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);

			fw = false;

			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
				&ebwtFw_, fw, qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
				 ebwtBw_, fw, qualCutoff_, false, verbose_, quiet_, 0, true,  maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
				&ebwtFw_, fw, qualCutoff_, false, verbose_, quiet_, 2,  false, maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
				*params, rRc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE, // no mismatches in hi half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // up to 2 in lo half
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE, // no mismatches in lo half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // up to 2 in hi half
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, verbose_, quiet_, true);
			EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
				*params, rRc_FwHalf, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);

			if(doFw_) {
				drVec->push_back(drFw_Bw);
				drVec->push_back(drFw_Seed);
				drVec->push_back(drFw_BwHalf);
			}
			if(doRc_) {
				drVec->push_back(drRc_Fw);
				drVec->push_back(drRc_Seed);
				drVec->push_back(drRc_FwHalf);
			}
		} else if(seedMms_ > 2) {

			bool fw = true;
			bool mate1 = true;

			EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				 0, false, maqPenalty_, qualOrder_, metrics_);

			// Partial and full aligners for alignments with 0
			// mismatches in the lo-half and up to 3 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rFw_BwSeed03 = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				 0, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rFw_FwSeedGen03 = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_,
				 0, true,  maqPenalty_, qualOrder_, metrics_);

			// Partial and full aligners for alignments with 1
			// mismatch in the lo-half and up to 2 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rFw_BwSeed12 = new EbwtRangeSourceFactory(
				 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				 0, false, maqPenalty_, qualOrder_, metrics_);
			// Note: the following is half-and-half (unlike the 03 version)
			EbwtRangeSource *rFw_FwSeedGen12 = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_,
				 3,  true,  maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSource *rFw_BwHalf12 = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_,
				 2,  false, maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
				*params, rFw_Bw, fw, false, maqPenalty_, qualOrder_,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE, // 0 mismatches in hi-half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // up to 3 mismatches in lo-half
				os_, verbose_, quiet_, mate1, pool_, btCnt);

			EbwtRangeSourceDriverFactory * drFw_BwSeed03 = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed03, fw, false, maqPenalty_, qualOrder_,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drFw_FwSeedGen03 = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen03, fw, true, maqPenalty_, qualOrder_,
				sink_, sinkPt, seedLen_,
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drFw_Seed03 = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed03, drFw_FwSeedGen03, fw, seedLen_, verbose_, quiet_, mate1);

			EbwtRangeSourceDriverFactory * drFw_BwSeed12 = new EbwtRangeSourceDriverFactory(
				*params, rFw_BwSeed12, fw, false, maqPenalty_, qualOrder_,
				sink_, sinkPt, seedLen_,
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drFw_FwSeedGen12 = new EbwtRangeSourceDriver(
				*params, rFw_FwSeedGen12, fw, true, maqPenalty_, qualOrder_,
				sink_, sinkPt, seedLen_,
				false,      // nudgeLeft
				PIN_TO_BEGINNING,
				PIN_TO_HI_HALF_EDGE, // 1-mismatch in lo-half
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,    // 1 or 2 mismatches in hi-half
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drFw_Seed12 = new EbwtSeededRangeSourceDriver(
				drFw_BwSeed12, drFw_FwSeedGen12, fw, seedLen_, verbose_, quiet_, mate1);

			EbwtRangeSourceDriver * drFw_BwHalf12 = new EbwtRangeSourceDriver(
				*params, rFw_BwHalf12, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);

			fw = false;

			EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				0, false, maqPenalty_, qualOrder_, metrics_);

			// Partial and full aligners for alignments with 0
			// mismatches in the lo-half and up to 3 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rRc_FwSeed03 = new EbwtRangeSourceFactory(
				&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				0, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rRc_BwSeedGen03 = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_,
				0, true,  maqPenalty_, qualOrder_, metrics_);

			// Partial and full aligners for alignments with 1
			// mismatch in the lo-half and up to 2 mismatches in the
			// hi-half
			EbwtRangeSourceFactory *rRc_FwSeed12 = new EbwtRangeSourceFactory(
				&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_,
				0, false, maqPenalty_, qualOrder_, metrics_);
			EbwtRangeSource *rRc_BwSeedGen12 = new EbwtRangeSource(
				 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_,
				3,  true,  maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSource *rRc_FwHalf12 = new EbwtRangeSource(
				&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_,
				2,  false, maqPenalty_, qualOrder_, metrics_);

			EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
				*params, rRc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);

			EbwtRangeSourceDriverFactory * drRc_FwSeed03 = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed03, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drRc_BwSeedGen03 = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen03, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drRc_Seed03 = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed03, drRc_BwSeedGen03, fw, seedLen_, verbose_, quiet_, mate1);

			EbwtRangeSourceDriverFactory * drRc_FwSeed12 = new EbwtRangeSourceDriverFactory(
				*params, rRc_FwSeed12, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtRangeSourceDriver * drRc_BwSeedGen12 = new EbwtRangeSourceDriver(
				*params, rRc_BwSeedGen12, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				false,      // nudgeLeft
				PIN_TO_BEGINNING,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);
			EbwtSeededRangeSourceDriver * drRc_Seed12 = new EbwtSeededRangeSourceDriver(
				drRc_FwSeed12, drRc_BwSeedGen12, fw, seedLen_, verbose_, quiet_, mate1);

			EbwtRangeSourceDriver * drRc_FwHalf12 = new EbwtRangeSourceDriver(
				*params, rRc_FwHalf12, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
				seedLen_,   // seedLen
				true,       // nudgeLeft
				PIN_TO_BEGINNING,    // nothing's unrevisitable
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_HI_HALF_EDGE,
				PIN_TO_SEED_EDGE,
				os_, verbose_, quiet_, mate1, pool_, btCnt);

			if(doFw_) {
				drVec->push_back(drFw_Bw);
				drVec->push_back(drFw_Seed03);
				drVec->push_back(drFw_Seed12);
				drVec->push_back(drFw_BwHalf12);
			}
			if(doRc_) {
				drVec->push_back(drRc_Fw);
				drVec->push_back(drRc_Seed03);
				drVec->push_back(drRc_Seed12);
				drVec->push_back(drRc_FwHalf12);
			}
		} else {
			cerr << "Unsupported --stateful mode: " << seedMms_ << endl;
		}
		TCostAwareRangeSrcDr* dr = new TCostAwareRangeSrcDr(strandFix_, drVec, verbose_, quiet_, false);
		delete drVec;

		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_, metrics_);

		return new UnpairedAlignerV2<EbwtRangeSource>(
			params, dr, rchase,
			sink_, sinkPtFactory_, sinkPt, os_, refs_,
			rangeMode_, verbose_, quiet_, maxBts_, pool_, btCnt,
			metrics_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	bool doFw_;
	bool doRc_;
	const uint32_t seedMms_;
	const uint32_t seedLen_;
	const int qualCutoff_;
	const int maxBts_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	RangeCache *cacheFw_;
	RangeCache *cacheBw_;
	const uint32_t cacheLimit_;
	ChunkPool *pool_;
	BitPairReference* refs_;
	vector<String<Dna5> >& os_;
	bool strandFix_;
	bool maqPenalty_;
	bool qualOrder_;
	bool rangeMode_;
	bool verbose_;
	bool quiet_;
	AlignerMetrics *metrics_;
};

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class PairedSeedAlignerFactory : public AlignerFactory {
	typedef RangeSourceDriver<EbwtRangeSource> TRangeSrcDr;
	typedef std::vector<TRangeSrcDr*> TRangeSrcDrPtrVec;
	typedef CostAwareRangeSourceDriver<EbwtRangeSource> TCostAwareRangeSrcDr;
public:
	PairedSeedAlignerFactory(
			Ebwt<String<Dna> >& ebwtFw,
			Ebwt<String<Dna> >* ebwtBw,
			bool color,
			bool v1,
			bool doFw,
			bool doRc,
			uint32_t seedMms,
			uint32_t seedLen,
			int qualCutoff,
			int maxBts,
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
			bool qualOrder,
			bool strandFix,
			bool rangeMode,
			bool verbose,
			bool quiet,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			color_(color),
			v1_(v1),
			doFw_(doFw),
			doRc_(doRc),
			seedMms_(seedMms),
			seedLen_(seedLen),
			qualCutoff_(qualCutoff),
			maxBts_(maxBts),
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
			quiet_(quiet)
	{
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
		RefAligner<String<Dna5> >* refAligner = NULL;
		int *btCnt = new int[1];
		*btCnt = maxBts_;
		if(seedMms_ == 0) {
			refAligner = new Seed0RefAligner<String<Dna5> >(color_, verbose_, quiet_, seedLen_, qualCutoff_, maqPenalty_);
		} else if(seedMms_ == 1) {
			refAligner = new Seed1RefAligner<String<Dna5> >(color_, verbose_, quiet_, seedLen_, qualCutoff_, maqPenalty_);
		} else if(seedMms_ == 2) {
			refAligner = new Seed2RefAligner<String<Dna5> >(color_, verbose_, quiet_, seedLen_, qualCutoff_, maqPenalty_);
		} else {
			refAligner = new Seed3RefAligner<String<Dna5> >(color_, verbose_, quiet_, seedLen_, qualCutoff_, maqPenalty_);
		}
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
		TRangeSrcDrPtrVec *dr1RcVec;
		TRangeSrcDrPtrVec *dr2FwVec;
		TRangeSrcDrPtrVec *dr2RcVec;
		if(v1_) {
			dr1RcVec = new TRangeSrcDrPtrVec();
			dr2FwVec = new TRangeSrcDrPtrVec();
			dr2RcVec = new TRangeSrcDrPtrVec();
		} else {
			dr1RcVec = dr1FwVec;
			dr2FwVec = dr1FwVec;
			dr2RcVec = dr1FwVec;
		}
		if(seedMms_ == 0) {
			const int halfAndHalf = 0;
			if(do1Fw) {
				bool mate1 = true;
				bool fw = true;
				EbwtRangeSource *r1Fw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true, verbose_, quiet_,
					 halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver *dr1Fw_Bw = new EbwtRangeSourceDriver(
					*params, r1Fw_Bw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr1FwVec->push_back(dr1Fw_Bw);
			}
			if(do2Fw) {
				bool mate1 = false;
				bool fw = true;
				EbwtRangeSource *r2Fw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true, verbose_, quiet_, halfAndHalf,
					 false, maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver *dr2Fw_Bw = new EbwtRangeSourceDriver(
					*params, r2Fw_Bw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr2FwVec->push_back(dr2Fw_Bw);
			}
			if(do1Rc) {
				bool mate1 = true;
				bool fw = false;
				EbwtRangeSource *r1Rc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true, verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver *dr1Rc_Fw = new EbwtRangeSourceDriver(
					*params, r1Rc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr1RcVec->push_back(dr1Rc_Fw);
			}
			if(do2Rc) {
				bool mate1 = false;
				bool fw = false;
				EbwtRangeSource *r2Rc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true, verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver *dr2Rc_Fw = new EbwtRangeSourceDriver(
					*params, r2Rc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft (not applicable)
					PIN_TO_SEED_EDGE, // whole alignment is unrevisitable
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					PIN_TO_SEED_EDGE, // "
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				dr2RcVec->push_back(dr2Rc_Fw);
			}
		} else if(seedMms_ == 1) {
			const int halfAndHalf = 0;
			if(do1Fw) {
				bool mate1 = true;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_,
					halfAndHalf, true,  maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				dr1FwVec->push_back(drFw_Bw);
				dr1FwVec->push_back(drFw_Seed);
			}
			if(do2Fw) {
				bool mate1 = false;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_,
					halfAndHalf, true,  maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				dr2FwVec->push_back(drFw_Bw);
				dr2FwVec->push_back(drFw_Seed);
			}
			if(do1Rc) {
				bool mate1 = true;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_,
					halfAndHalf, true,  maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				dr1RcVec->push_back(drRc_Fw);
				dr1RcVec->push_back(drRc_Seed);
			}
			if(do2Rc) {
				bool mate1 = false;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_,
					halfAndHalf, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_,
					halfAndHalf, true,  maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, NULL); // no backtrack limit for -n 1/2
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				dr2RcVec->push_back(drRc_Fw);
				dr2RcVec->push_back(drRc_Seed);
			}
		} else if(seedMms_ > 1) {
			bool two = seedMms_ == 2;
			if(do1Fw) {
				bool mate1 = true;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_, 0, true,  maqPenalty_, qualOrder_);
				EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_, 2,  false, maqPenalty_, qualOrder_);
				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				EbwtRangeSourceDriverFactory * drFw_BwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rFw_BwSeed12 = new EbwtRangeSourceFactory(
						 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
					drFw_BwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rFw_BwSeed12, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drFw_FwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rFw_FwSeedGen12 = new EbwtRangeSource(
						&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_, 3,  true,  maqPenalty_, qualOrder_);
					drFw_FwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rFw_FwSeedGen12, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drFw_Seed12 = NULL;
				if(!two) {
					drFw_Seed12 = new EbwtSeededRangeSourceDriver(
						drFw_BwSeed12, drFw_FwSeedGen12, fw, seedLen_, verbose_, quiet_, mate1);
				}
				EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
					*params, rFw_BwHalf, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				dr1FwVec->push_back(drFw_Bw);
				dr1FwVec->push_back(drFw_Seed);
				if(drFw_Seed12 != NULL) {
					dr1FwVec->push_back(drFw_Seed12);
				}
				dr1FwVec->push_back(drFw_BwHalf);
			}
			if(do2Fw) {
				bool mate1 = false;
				bool fw = true;
				EbwtRangeSource *rFw_Bw = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rFw_BwSeed = new EbwtRangeSourceFactory(
					 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rFw_FwSeedGen = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_, 0, true,  maqPenalty_, qualOrder_);
				EbwtRangeSource *rFw_BwHalf = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_, 2,  false, maqPenalty_, qualOrder_);

				EbwtRangeSourceDriver * drFw_Bw = new EbwtRangeSourceDriver(
					*params, rFw_Bw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drFw_BwSeed = new EbwtRangeSourceDriverFactory(
					*params, rFw_BwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drFw_FwSeedGen = new EbwtRangeSourceDriver(
					*params, rFw_FwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drFw_Seed = new EbwtSeededRangeSourceDriver(
					drFw_BwSeed, drFw_FwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				EbwtRangeSourceDriverFactory * drFw_BwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rFw_BwSeed12 = new EbwtRangeSourceFactory(
						 ebwtBw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
					drFw_BwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rFw_BwSeed12, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drFw_FwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rFw_FwSeedGen12 = new EbwtRangeSource(
						&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_, 3,  true,  maqPenalty_, qualOrder_);
					drFw_FwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rFw_FwSeedGen12, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drFw_Seed12 = NULL;
				if(!two) {
					drFw_Seed12 = new EbwtSeededRangeSourceDriver(
						drFw_BwSeed12, drFw_FwSeedGen12, fw, seedLen_, verbose_, quiet_, mate1);
				}
				EbwtRangeSourceDriver * drFw_BwHalf = new EbwtRangeSourceDriver(
					*params, rFw_BwHalf, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				dr2FwVec->push_back(drFw_Bw);
				dr2FwVec->push_back(drFw_Seed);
				if(drFw_Seed12 != NULL) {
					dr2FwVec->push_back(drFw_Seed12);
				}
				dr2FwVec->push_back(drFw_BwHalf);
			}
			if(do1Rc) {
				bool mate1 = true;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_, 0, true,  maqPenalty_, qualOrder_);
				EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_, 2,  false, maqPenalty_, qualOrder_);

				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				EbwtRangeSourceDriverFactory * drRc_FwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rRc_FwSeed12 = new EbwtRangeSourceFactory(
						&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
					drRc_FwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rRc_FwSeed12, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drRc_BwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rRc_BwSeedGen12 = new EbwtRangeSource(
						 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_, 3,  true,  maqPenalty_, qualOrder_);
					drRc_BwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rRc_BwSeedGen12, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drRc_Seed12 = NULL;
				if(!two) {
					drRc_Seed12 = new EbwtSeededRangeSourceDriver(
						drRc_FwSeed12, drRc_BwSeedGen12, fw, seedLen_, verbose_, quiet_, mate1);
				}
				EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
					*params, rRc_FwHalf, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				dr1RcVec->push_back(drRc_Fw);
				dr1RcVec->push_back(drRc_Seed);
				if(drRc_Seed12 != NULL) {
					dr1RcVec->push_back(drRc_Seed12);
				}
				dr1RcVec->push_back(drRc_FwHalf);
			}
			if(do2Rc) {
				bool mate1 = false;
				bool fw = false;
				EbwtRangeSource *rRc_Fw = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSourceFactory *rRc_FwSeed = new EbwtRangeSourceFactory(
					&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
				EbwtRangeSource *rRc_BwSeedGen = new EbwtRangeSource(
					 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_, 0, true,  maqPenalty_, qualOrder_);
				EbwtRangeSource *rRc_FwHalf = new EbwtRangeSource(
					&ebwtFw_, fw,  qualCutoff_, false, verbose_, quiet_, 2,  false, maqPenalty_, qualOrder_);

				EbwtRangeSourceDriver * drRc_Fw = new EbwtRangeSourceDriver(
					*params, rRc_Fw, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriverFactory * drRc_FwSeed = new EbwtRangeSourceDriverFactory(
					*params, rRc_FwSeed, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtRangeSourceDriver * drRc_BwSeedGen = new EbwtRangeSourceDriver(
					*params, rRc_BwSeedGen, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					false,      // nudgeLeft
					PIN_TO_HI_HALF_EDGE,
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				EbwtSeededRangeSourceDriver * drRc_Seed = new EbwtSeededRangeSourceDriver(
					drRc_FwSeed, drRc_BwSeedGen, fw, seedLen_, verbose_, quiet_, mate1);
				EbwtRangeSourceDriverFactory * drRc_FwSeed12 = NULL;
				if(!two) {
					EbwtRangeSourceFactory *rRc_FwSeed12 = new EbwtRangeSourceFactory(
						&ebwtFw_, fw,  qualCutoff_, true,  verbose_, quiet_, 0, false, maqPenalty_, qualOrder_);
					drRc_FwSeed12 = new EbwtRangeSourceDriverFactory(
						*params, rRc_FwSeed12, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						true,       // nudgeLeft
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtRangeSourceDriver * drRc_BwSeedGen12 = NULL;
				if(!two) {
					EbwtRangeSource *rRc_BwSeedGen12 = new EbwtRangeSource(
						 ebwtBw_, fw,  qualCutoff_, false, verbose_, quiet_, 3,  true,  maqPenalty_, qualOrder_);
					drRc_BwSeedGen12 = new EbwtRangeSourceDriver(
						*params, rRc_BwSeedGen12, fw, true, maqPenalty_, qualOrder_, sink_, sinkPt,
						seedLen_,   // seedLen
						false,      // nudgeLeft
						PIN_TO_BEGINNING,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_HI_HALF_EDGE,
						PIN_TO_SEED_EDGE,
						os_, verbose_, quiet_, mate1, pool_, btCnt);
				}
				EbwtSeededRangeSourceDriver * drRc_Seed12 = NULL;
				if(!two) {
					drRc_Seed12 = new EbwtSeededRangeSourceDriver(
						drRc_FwSeed12, drRc_BwSeedGen12, fw, seedLen_, verbose_, quiet_, mate1);
				}
				EbwtRangeSourceDriver * drRc_FwHalf = new EbwtRangeSourceDriver(
					*params, rRc_FwHalf, fw, false, maqPenalty_, qualOrder_, sink_, sinkPt,
					seedLen_,   // seedLen
					true,       // nudgeLeft
					PIN_TO_BEGINNING,    // nothing's unrevisitable
					PIN_TO_HI_HALF_EDGE,
					two ? PIN_TO_SEED_EDGE : PIN_TO_HI_HALF_EDGE,
					PIN_TO_SEED_EDGE,
					os_, verbose_, quiet_, mate1, pool_, btCnt);
				dr2RcVec->push_back(drRc_Fw);
				dr2RcVec->push_back(drRc_Seed);
				if(drRc_Seed12 != NULL) {
					dr2RcVec->push_back(drRc_Seed12);
				}
				dr2RcVec->push_back(drRc_FwHalf);
			}
		} else {
			cerr << "Unsupported --stateful mode: " << seedMms_ << endl;
		}
		// Set up a RangeChaser
		RangeChaser<String<Dna> > *rchase =
			new RangeChaser<String<Dna> >(cacheLimit_, cacheFw_, cacheBw_);

		if(v1_) {
			PairedBWAlignerV1<EbwtRangeSource>* al = new PairedBWAlignerV1<EbwtRangeSource>(
				params,
				new TCostAwareRangeSrcDr(strandFix_, dr1FwVec, verbose_, quiet_, false),
				new TCostAwareRangeSrcDr(strandFix_, dr1RcVec, verbose_, quiet_, false),
				new TCostAwareRangeSrcDr(strandFix_, dr2FwVec, verbose_, quiet_, false),
				new TCostAwareRangeSrcDr(strandFix_, dr2RcVec, verbose_, quiet_, false),
				refAligner, rchase, sink_, sinkPtFactory_, sinkPt,
				mate1fw_, mate2fw_, peInner_, peOuter_, dontReconcile_,
				symCeil_, mixedThresh_, mixedAttemptLim_, refs_,
				rangeMode_, verbose_, quiet_, maxBts_, pool_,
				btCnt);
			delete dr1FwVec;
			delete dr1RcVec;
			delete dr2FwVec;
			delete dr2RcVec;
			return al;
		} else {
			// We dumped all the drivers into dr1FwVec
			PairedBWAlignerV2<EbwtRangeSource>* al = new PairedBWAlignerV2<EbwtRangeSource>(
				params, paramsSe1, paramsSe2,
				new TCostAwareRangeSrcDr(strandFix_, dr1FwVec, verbose_, quiet_, true),
				refAligner, rchase, sink_, sinkPtFactory_, sinkPt,
				sinkPtSe1, sinkPtSe2, mate1fw_, mate2fw_, peInner_, peOuter_,
				mixedAttemptLim_, refs_, rangeMode_, verbose_,
				quiet_, maxBts_, pool_, btCnt);
			delete dr1FwVec;
			return al;
		}
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	bool color_;
	const bool v1_; // whether to use V1 PairedAligner
	const bool doFw_;
	const bool doRc_;
	const uint32_t seedMms_;
	const uint32_t seedLen_;
	const int qualCutoff_;
	const int maxBts_;
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
};

#endif /* ALIGNER_SEED_MM_H_ */
