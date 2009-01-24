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
		return new UnpairedAlignerV1<GreedyDFSRangeSource, GreedyDFSContinuationManager>(
			ebwt_, sink_, sinkPtFactory_, os_, rangeMode_, verbose_, seed_);
	}

	/**
	 * Create a new vector of new UnpairedExactAlignerV1s.
	 */
	virtual std::vector<Aligner*>* create(uint32_t n) const {
		std::vector<Aligner*>* v = new std::vector<Aligner*>;
		for(uint32_t i = 0; i < n; i++) {
			v->push_back(new UnpairedAlignerV1<GreedyDFSRangeSource, GreedyDFSContinuationManager>(
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

/**
 * Concrete factory class for constructing unpaired exact aligners.
 */
class PairedExactAlignerV1Factory : public AlignerFactory {
public:
	PairedExactAlignerV1Factory(
			Ebwt<String<Dna> >& ebwt,
			HitSink& sink,
			const HitSinkPerThreadFactory& sinkPtFactory,
			bool mate1fw,
			bool mate2fw,
			uint32_t peInner,
			uint32_t peOuter,
			vector<String<Dna5> >& os,
			bool rangeMode,
			bool verbose,
			uint32_t seed) :
			ebwt_(ebwt),
			sink_(sink),
			sinkPtFactory_(sinkPtFactory),
			mate1fw_(mate1fw),
			mate2fw_(mate2fw),
			peInner_(peInner),
			peOuter_(peOuter),
			os_(os),
			rangeMode_(rangeMode),
			verbose_(verbose),
			seed_(seed)
	{ }

	/**
	 * Create a new UnpairedExactAlignerV1s.
	 */
	virtual Aligner* create() const {
		return new PairedAlignerV1<GreedyDFSRangeSource, GreedyDFSContinuationManager>(
			ebwt_, sink_, sinkPtFactory_, mate1fw_, mate2fw_,
			peInner_, peOuter_, os_, rangeMode_, verbose_, seed_);
	}

	/**
	 * Create a new vector of new UnpairedExactAlignerV1s.
	 */
	virtual std::vector<Aligner*>* create(uint32_t n) const {
		std::vector<Aligner*>* v = new std::vector<Aligner*>;
		for(uint32_t i = 0; i < n; i++) {
			v->push_back(new PairedAlignerV1<GreedyDFSRangeSource, GreedyDFSContinuationManager>(
				ebwt_, sink_, sinkPtFactory_, mate1fw_, mate2fw_,
				peInner_, peOuter_, os_, rangeMode_, verbose_, seed_));
			assert(v->back() != NULL);
		}
		return v;
	}

private:
	Ebwt<String<Dna> >& ebwt_;
	HitSink& sink_;
	const HitSinkPerThreadFactory& sinkPtFactory_;
	bool mate1fw_;
	bool mate2fw_;
	uint32_t peInner_;
	uint32_t peOuter_;
	vector<String<Dna5> >& os_;
	bool rangeMode_;
	bool verbose_;
	uint32_t seed_;
};

#endif /* ALIGNER_0MM_H_ */
