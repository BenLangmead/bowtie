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
			vector<String<Dna5> >& os,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwtFw_(ebwtFw),
			ebwtBw_(ebwtBw),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
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

		return new UnpairedAlignerV1<GreedyDFSRangeSource, GreedyDFSContinuationManager>(
			params, driverFw, driverRc,
			sink_, sinkPtFactory_, sinkPt, os_, rangeMode_, verbose_,
			seed_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
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
			uint32_t symCeil,
			uint32_t mixedThresh,
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
			symCeil_(symCeil),
			mixedThresh_(mixedThresh),
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
		HitSinkPerThread* sinkPt = sinkPtFactory_.createMult(2, 0xffffffff);
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

		return new PairedBWAlignerV1<GreedyDFSRangeSource, GreedyDFSContinuationManager>(
			params,
			driver1Fw, driver1Rc, driver2Fw, driver2Rc, refAligner,
			sink_, sinkPtFactory_, sinkPt, mate1fw_, mate2fw_,
			peInner_, peOuter_, symCeil_, mixedThresh_, os_,
			rangeMode_, verbose_, seed_);
	}

private:
	Ebwt<String<Dna> >& ebwtFw_;
	Ebwt<String<Dna> >* ebwtBw_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	bool mate1fw_;
	bool mate2fw_;
	uint32_t peInner_;
	uint32_t peOuter_;
	uint32_t symCeil_;
	uint32_t mixedThresh_;
	vector<String<Dna5> >& os_;
	bool rangeMode_;
	bool verbose_;
	uint32_t seed_;
};

#endif /* ALIGNER_0MM_H_ */
