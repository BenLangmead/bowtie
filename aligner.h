/**
 * aligner.h
 *
 * A generic class providing a stateful way to find alignments.
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <iostream>
#include <stdint.h>
#include "seqan/sequence.h"
#include "assert_helpers.h"
#include "ebwt.h"
#include "pat.h"
#include "range.h"
#include "range_source.h"

/**
 * State machine for carrying out an alignment, which usually consists
 * of a series of phases that conduct different alignments using
 * different backtracking constraints.
 *
 * Each Aligner should have a dedicated PatternSourcePerThread.
 */
class Aligner {
public:
	Aligner(const Ebwt<String<Dna> >* ebwtFw,
	        const Ebwt<String<Dna> >* ebwtRc,
	        bool rangeMode,
	        uint32_t seed) :
		ebwtFw_(ebwtFw), ebwtRc_(ebwtRc), patsrc_(NULL),
		bufa_(NULL), bufb_(NULL), rangeMode_(rangeMode), seed_(seed)
	{ }

	virtual ~Aligner() { }
	/// Advance the range search by one memory op
	virtual bool advance() = 0;
	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() = 0;

	/// Prepare Aligner for the next read
	virtual void setQuery(PatternSourcePerThread *patsrc) {
		assert(patsrc != NULL);
		patsrc_ = patsrc;
		bufa_ = &patsrc->bufa();
		assert(bufa_ != NULL);
		bufb_ = &patsrc->bufb();
		alen_ = bufa_->length();
		blen_ = (bufb_ != NULL) ? bufb_->length() : 0;
		rand_.init(seed_ + genRandSeed(bufa_->patFw, bufa_->qualFw, bufa_->name));
	}

	/**
	 * Report hit(s) from a given range.
	 */
	bool reportSingleEndHitFromRange(Range& ra,
	                                 EbwtSearchParams<String<Dna> >& params,
	                                 bool fw)
	{
		assert_gt(ra.bot, ra.top);
		assert(bufa_ != NULL);
		assert(!seqan::empty(bufa_->patFw));
		assert(!seqan::empty(bufa_->qualFw));
		assert(!seqan::empty(bufa_->name));
//		if(stackDepth == 0 && !_reportExacts) {
//			// We are not reporting exact hits (usually because we've
//			// already reported them as part of a previous invocation
//			// of the backtracker)
//			return false;
//		}
		if(rangeMode_) {
			return ebwtFw_->report(
					fw ?  bufa_->patFw  :  bufa_->patRc,
					fw ? &bufa_->qualFw : &bufa_->qualRc,
					&bufa_->name,
                    ra.mms, ra.refcs, ra.numMms, 0,
                    ra.top, ra.bot, alen_,
                    ra.stratum, params);
		}
		uint32_t spread = ra.bot - ra.top;
		// Pick a random spot in the range to begin report
		uint32_t r = ra.top + (rand_.nextU32() % spread);
		for(uint32_t i = 0; i < spread; i++) {
			uint32_t ri = r + i;
			if(ri >= ra.bot) ri -= spread;
			// reportChaseOne takes the _mms[] list in terms of
			// their indices into the query string; not in terms
			// of their offset from the 3' or 5' end.
			if(ebwtFw_->reportChaseOne(fw ?  bufa_->patFw  :  bufa_->patRc,
			                           fw ? &bufa_->qualFw : &bufa_->qualRc,
			                           &bufa_->name,
			                           ra.mms, ra.refcs, ra.numMms, ri,
			                           ra.top, ra.bot, alen_,
			                           ra.stratum, params))
			{
				// Return value of true means that we can stop
				return true;
			}
			// Return value of false means that we should continue
			// searching.  This could happen if we the call to
			// reportChaseOne() reported a hit, but the user asked for
			// multiple hits and we haven't reached the ceiling yet.
			// This might also happen if the call to reportChaseOne()
			// didn't report a hit because the alignment was spurious
			// (i.e. overlapped some padding).
		}
		// All range elements were examined and we should keep going
		return false;
	}

protected:

	// Index
	const Ebwt<String<Dna> >* ebwtFw_;
	const Ebwt<String<Dna> >* ebwtRc_;
	// Current read pair
	PatternSourcePerThread* patsrc_;
	ReadBuf* bufa_;
	uint32_t alen_;
	ReadBuf* bufb_;
	uint32_t blen_;
	// RandomSource for choosing alignments to report from ranges
	bool rangeMode_;
	uint32_t seed_;
	RandomSource rand_;
};

/**
 * Abstract parent factory class for constructing aligners of all kinds.
 */
class AlignerFactory {
public:
	virtual ~AlignerFactory() { }
	virtual Aligner* create() const = 0;
	virtual std::vector<Aligner*>* create(uint32_t) const = 0;

	/// Free memory associated with the aligner
	virtual void destroy(Aligner* al) const {
		assert(al != NULL);
		// Free the Aligner
		delete al;
	}

	/// Free memory associated with an aligner list
	virtual void destroy(std::vector<Aligner*>* als) const {
		assert(als != NULL);
		// Free all of the Aligners
		for(size_t i = 0; i < als->size(); i++) {
			if((*als)[i] != NULL) {
				delete (*als)[i];
				(*als)[i] = NULL;
			}
		}
		// Free the vector
		delete als;
	}
};

/**
 * Coordinates multiple aligners.
 */
class MultiAligner {
public:
	MultiAligner(
			uint32_t n,
			uint32_t qUpto,
			const AlignerFactory& alignFact,
			const PatternSourcePerThreadFactory& patsrcFact) :
			n_(n), qUpto_(qUpto),
			alignFact_(alignFact), patsrcFact_(patsrcFact),
			aligners_(NULL), patsrcs_(NULL)
	{
		aligners_ = alignFact_.create(n_);
		assert(aligners_ != NULL);
		patsrcs_ = patsrcFact_.create(n_);
		assert(patsrcs_ != NULL);
	}

	/// Free memory associated with the aligners and their pattern sources.
	virtual ~MultiAligner() {
		alignFact_.destroy(aligners_);
		patsrcFact_.destroy(patsrcs_);
	}

	/**
	 * Advance an array of aligners in parallel, using prefetches to
	 * try to hide all the latency.
	 */
	void run() {
		bool done = false;
		while(!done) {
			done = true;
			for(uint32_t i = 0; i < n_; i++) {
				if(!(*aligners_)[i]->done()) {
					// Advance an aligner already in progress
					done = false;
					(*aligners_)[i]->advance();
				} else if(qUpto_ > 0) {
					// Get a new read and initialize an aligner with it
					(*patsrcs_)[i]->nextReadPair();
					if(!(*patsrcs_)[i]->empty()) {
						qUpto_--;
						(*aligners_)[i]->setQuery((*patsrcs_)[i]);
						assert(!(*aligners_)[i]->done());
						done = false;
					} else {
						// No more reads; if done == true, it remains
						// true
					}
				} else {
					// Past read limit; if done == true, it remains true
				}
			}
		}
	}

protected:
	uint32_t n_;     /// Number of aligners
	uint32_t qUpto_; /// Number of reads to align before stopping
	const AlignerFactory&                  alignFact_;
	const PatternSourcePerThreadFactory&   patsrcFact_;
	std::vector<Aligner *>*                aligners_;
	std::vector<PatternSourcePerThread *>* patsrcs_;
};

#endif /* ALIGNER_H_ */
