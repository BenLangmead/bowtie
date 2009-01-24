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
#include "range_chaser.h"

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
 * Coordinates multiple aligners of the same type (i.e. either all
 * single-end or all paired-end).
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

/**
 * Coordinates multiple single-end and paired-end aligners, routing
 * reads to one or the other type as appropriate.
 */
class MixedMultiAligner {
public:
	MixedMultiAligner(
			uint32_t n,
			uint32_t qUpto,
			const AlignerFactory& alignSEFact,
			const AlignerFactory& alignPEFact,
			const PatternSourcePerThreadFactory& patsrcFact) :
			n_(n), qUpto_(qUpto),
			alignSEFact_(alignSEFact),
			alignPEFact_(alignPEFact),
			patsrcFact_(patsrcFact),
			alignersSE_(NULL),
			alignersPE_(NULL),
			seOrPe_(NULL),
			patsrcs_(NULL)
	{
		// Instantiate all single-end aligners
		alignersSE_ = alignSEFact_.create(n_);
		assert(alignersSE_ != NULL);
		// Instantiate all paired-end aligners
		alignersPE_ = alignPEFact_.create(n_);
		assert(alignersPE_ != NULL);
		// Allocate array of boolean flags indicating whether each of
		// the slots is currently using the single-end or paired-end
		// aligner
		seOrPe_ = new bool[n_];
		// Instantiate all read sources
		patsrcs_ = patsrcFact_.create(n_);
		assert(patsrcs_ != NULL);
	}

	/// Free memory associated with the aligners and their pattern sources.
	virtual ~MixedMultiAligner() {
		alignSEFact_.destroy(alignersSE_);
		alignPEFact_.destroy(alignersPE_);
		patsrcFact_.destroy(patsrcs_);
		delete[] seOrPe_;
	}

	/**
	 * Advance an array of aligners in parallel, using prefetches to
	 * try to hide all the latency.
	 */
	void run(bool verbose = false) {
		bool done = false;
		while(!done) {
			done = true;
			for(uint32_t i = 0; i < n_; i++) {
				Aligner *al = seOrPe_[i] ? (*alignersSE_)[i] :
				                           (*alignersPE_)[i];
				if(!al->done()) {
					// Advance an aligner already in progress; this is
					// the common case
					done = false;
					al->advance();
				} else if(qUpto_ > 0) {
					// Feed a new read to a vacant aligner
					PatternSourcePerThread *ps = (*patsrcs_)[i];
					// Get a new read
					ps->nextReadPair();
					if(!ps->empty()) {
						qUpto_--;
						if(ps->paired()) {
							// Read currently in buffer is paired-end
							if(verbose) cout << "Paired input: " << ps->bufa().patFw << "," << ps->bufb().patFw << endl;
							(*alignersPE_)[i]->setQuery(ps);
							seOrPe_[i] = false; // false -> paired
						} else {
							// Read currently in buffer is single-end
							if(verbose) cout << "Unpaired input: " << ps->bufa().patFw << endl;
							(*alignersSE_)[i]->setQuery(ps);
							seOrPe_[i] = true; // true = unpaired
						}
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
	const AlignerFactory&                  alignSEFact_;
	const AlignerFactory&                  alignPEFact_;
	const PatternSourcePerThreadFactory&   patsrcFact_;
	std::vector<Aligner *>*                alignersSE_;
	std::vector<Aligner *>*                alignersPE_;
	bool *                                 seOrPe_;
	std::vector<PatternSourcePerThread *>* patsrcs_;
};

/**
 * An aligner for finding exact matches of unpaired reads.  Always
 * tries the forward-oriented version of the read before the reverse-
 * oriented read.
 */
template<typename TRangeSource, typename TContMan>
class UnpairedAlignerV1 : public Aligner {
	typedef RangeSourceDriver<TRangeSource, TContMan> TDriver;
public:
	UnpairedAlignerV1(
		const Ebwt<String<Dna> >& ebwt,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		vector<String<Dna5> >& os,
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(&ebwt, NULL, rangeMode, seed),
		doneFw_(true), done_(true), firstFw_(true), firstRc_(true),
		chaseFw_(false), chaseRc_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPtFactory.create()),
		params_(*sinkPt_, os, true, true, true, rangeMode),
		rchase_(ebwt, rand_),
		rFw_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		    false, NULL, NULL, verbose, seed, &os, false, false, false),
		cFw_(),
		rRc_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		    false, NULL, NULL, verbose, seed, &os, false, false, false),
		cRc_(),
		driverFw_(ebwt, params_, rFw_, cFw_, true,  sink, sinkPt_, os, verbose, seed),
		driverRc_(ebwt, params_, rRc_, cRc_, false, sink, sinkPt_, os, verbose, seed)
	{ }

	virtual ~UnpairedAlignerV1() {
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		Aligner::setQuery(patsrc); // set fields & random seed
		driverFw_.setQuery(patsrc);
		driverRc_.setQuery(patsrc);
		doneFw_  = false;
		done_    = false;
		firstFw_ = true;
		firstRc_ = true;
		chaseFw_ = false;
		chaseRc_ = false;
	}

	/**
	 * Helper for reporting an alignment.
	 */
	inline bool report(const Range& ra,
	                   uint32_t first,
	                   uint32_t second,
	                   uint32_t tlen,
	                   bool fw)
	{
		return params_.reportHit(
				fw ?  bufa_->patFw  :  bufa_->patRc,
				fw ? &bufa_->qualFw : &bufa_->qualRc,
				&bufa_->name,
				ra.mms,                   // mismatch positions
				ra.refcs,                 // reference characters for mms
				ra.numMms,                // # mismatches
				make_pair(first, second), // position
				make_pair(ra.top, ra.bot),// arrows
				tlen,                     // textlen
				alen_,                    // qlen
				ra.stratum,               // alignment stratum
				ra.bot - ra.top - 1);     // # other hits
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual bool advance() {
		assert(!done_);
		assert(!chaseFw_ || !chaseRc_);
		if(chaseFw_) {
			assert(!rangeMode_);
			assert(driverFw_.foundRange());
			const Range& ra = driverFw_.range();
			if(!rchase_.foundOff() && !rchase_.done()) {
				rchase_.advance();
				return false;
			}
			if(rchase_.foundOff()) {
				done_ = report(ra, rchase_.off().first, rchase_.off().second, rchase_.tlen(), true);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				chaseFw_ = false;
				done_ = driverFw_.done();
			}
		} else if(chaseRc_) {
			assert(!rangeMode_);
			assert(driverRc_.foundRange());
			const Range& ra = driverRc_.range();
			if(!rchase_.foundOff() && !rchase_.done()) {
				rchase_.advance();
				return false;
			}
			if(rchase_.foundOff()) {
				done_ = report(ra, rchase_.off().first, rchase_.off().second, rchase_.tlen(), false);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				chaseRc_ = false;
				done_ = driverRc_.done();
			}
		}
		// Still advancing a
		if(!done_ && !chaseFw_ && !chaseRc_) {
			if(!doneFw_) {
				driverFw_.advance();
				if(driverFw_.foundRange()) {
					const Range& ra = driverFw_.range();
					if(rangeMode_) {
						done_ = report(ra, ra.top, ra.bot, 0, true);
					} else {
						rchase_.setTopBot(ra.top, ra.bot, alen_);
						if(rchase_.foundOff()) {
							done_ = report(ra, rchase_.off().first, rchase_.off().second, rchase_.tlen(), true);
							rchase_.reset();
						}
						if(!rchase_.done()) {
							// Keep chasing this range
							chaseFw_ = true;
						}
					}
				}
				if(!done_ && !chaseFw_) {
					doneFw_ = driverFw_.done();
				}
			} else {
				driverRc_.advance();
				// Advance the RangeSource for the reverse-complement read
				if(driverRc_.foundRange()) {
					const Range& ra = driverRc_.range();
					if(rangeMode_) {
						done_ = report(ra, ra.top, ra.bot, 0, false);
					} else {
						rchase_.setTopBot(ra.top, ra.bot, alen_);
						if(rchase_.foundOff()) {
							done_ = report(ra, rchase_.off().first, rchase_.off().second, rchase_.tlen(), false);
							rchase_.reset();
						}
						if(!rchase_.done()) {
							// Keep chasing this range
							chaseRc_ = true;
						}
					}
				}
				if(!done_ && !chaseRc_) {
					done_ = driverRc_.done();
				}
			}
		}
		if(done_) {
			sinkPt_->finishRead(*patsrc_, true);
		}
		return done_;
	}

	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() {
		return done_;
	}

protected:
	// Progress state
	bool doneFw_;
	bool done_;
	bool firstFw_;
	bool firstRc_;
	bool chaseFw_;
	bool chaseRc_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> > params_;

	// State for getting alignments from ranges statefully
	RandomScanningRangeChaser<String<Dna> > rchase_;

	// Range-finding state
	TRangeSource rFw_;
	TContMan     cFw_;
	TRangeSource rRc_;
	TContMan     cRc_;
	TDriver      driverFw_;
	TDriver      driverRc_;
};

/**
 * An aligner for finding exact matches of paired reads.
 */
template<typename TRangeSource, typename TContMan>
class PairedAlignerV1 : public Aligner {

	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef std::vector<U32Pair> U32PairVec;
	typedef RangeSourceDriver<TRangeSource, TContMan> TDriver;

public:
	PairedAlignerV1(
		const Ebwt<String<Dna> >& ebwt,
		HitSink& sink,
		const HitSinkPerThreadFactory& sinkPtFactory,
		bool fw1, bool fw2,
		uint32_t minInsert,
		uint32_t maxInsert,
		vector<String<Dna5> >& os,
		bool rangeMode,
		bool verbose,
		uint32_t seed) :
		Aligner(&ebwt, NULL, rangeMode, seed),
		doneFw_(true), done_(true), firstFw_(true), firstRc_(true),
		chase1Fw_(false), chase1Rc_(false),
		chase2Fw_(false), chase2Rc_(false),
		delayedChase1Fw_(false), delayedChase1Rc_(false),
		delayedChase2Fw_(false), delayedChase2Rc_(false),
		sinkPtFactory_(sinkPtFactory),
		sinkPt_(sinkPtFactory.createMult(2)), // using Mult(2) allows us to print two records per paired alignment
		params_(*sinkPt_, os, true, true, true, rangeMode),
		minInsert_(minInsert),
		maxInsert_(maxInsert),
		fw1_(fw1), fw2_(fw2),
		rchase_(ebwt, rand_),
		r1Fw_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c1Fw_(),
		r1Rc_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c1Rc_(),
		driver1Fw_(ebwt, params_, r1Fw_, c1Fw_, fw1,  sink, sinkPt_, os, verbose, seed),
		driver1Rc_(ebwt, params_, r1Rc_, c1Rc_, !fw1, sink, sinkPt_, os, verbose, seed),
		offs1FwSz_(0), offs1RcSz_(0),
		r2Fw_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c2Fw_(),
		r2Rc_(&ebwt, params_, 0xffffffff, BacktrackLimits(), 0, true,
		      false, NULL, NULL, verbose, seed, &os, false, false, false),
		c2Rc_(),
		driver2Fw_(ebwt, params_, r2Fw_, c2Fw_, fw2,  sink, sinkPt_, os, verbose, seed),
		driver2Rc_(ebwt, params_, r2Rc_, c2Rc_, !fw2, sink, sinkPt_, os, verbose, seed),
		offs2FwSz_(0), offs2RcSz_(0)
	{ }

	virtual ~PairedAlignerV1() {
		sinkPtFactory_.destroy(sinkPt_); sinkPt_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc) {
		assert(!patsrc->bufa().empty());
		Aligner::setQuery(patsrc); // set fields & random seed
		assert(!patsrc->bufb().empty());
		// Give all of the drivers pointers to the relevant read info
		driver1Fw_.setQuery(patsrc, true  /* mate1 */);
		driver1Rc_.setQuery(patsrc, true  /* mate1 */);
		driver2Fw_.setQuery(patsrc, false /* mate2 */);
		driver2Rc_.setQuery(patsrc, false /* mate2 */);
		// Neither orientation is done
		doneFw_   = false;
		done_     = false;
		// No ranges are being chased yet
		chase1Fw_ = false;
		chase1Rc_ = false;
		chase2Fw_ = false;
		chase2Rc_ = false;
		delayedChase1Fw_ = false;
		delayedChase1Rc_ = false;
		delayedChase2Fw_ = false;
		delayedChase2Rc_ = false;
		// Clear all intermediate ranges
		offs1Fw_.clear(); offs1Rc_.clear();
		offs2Fw_.clear(); offs2Rc_.clear();
		offs1FwSz_ = offs1RcSz_ = offs2FwSz_ = offs2RcSz_ = 0;
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 *
	 * A call to this function does one of many things:
	 * 1. Advance a RangeSourceDriver and check if it found a new range
	 * 2. Advance a RowChaseDriver and check if it found a reference
	 *    offset for a an alignment in a range
	 */
	virtual bool advance() {
		assert(!done_);
		bool verbose = false;
		if(!doneFw_) {
			// Mate 1 aligns upstream of mate 2.  Expected orientations
			// are stored in fw1_, fw2_.
			advanceOrientation(chase1Fw_,        chase2Fw_,
	                           delayedChase1Fw_, delayedChase2Fw_,
			                   driver1Fw_,       driver2Fw_,
			                   offs1Fw_,         offs2Fw_,
			                   offs1FwSz_,       offs2FwSz_,
			                   doneFw_,
			                   fw1_, fw2_,
			                   true, verbose);
		} else {
			// In reverse-complement space, we expect mate 2 to align
			// upstream of mate 1.  Also, mates 1 and 2 have reversed
			// orientations from what we stored in fw1_ and fw2_.
			advanceOrientation(chase2Rc_,        chase1Rc_,
			                   delayedChase2Rc_, delayedChase1Rc_,
			                   driver2Rc_,       driver1Rc_,
			                   offs2Rc_,         offs1Rc_,
			                   offs2RcSz_,       offs1RcSz_,
			                   done_,
			                   !fw2_, !fw1_,
			                   false, verbose);
		}
		if(done_) {
			sinkPt_->finishRead(*patsrc_, true);
		}
		return done_;
	}

	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() {
		return done_;
	}

protected:

	/**
	 * Helper for reporting a pair of alignments.  As of now, we report
	 * a paired alignment by reporting two consecutive alignments, one
	 * for each mate.
	 */
	bool report(const Range& ra, // range for mate1
	            const Range& rb, // range for mate2
	            uint32_t first,  // ref idx
	            uint32_t upstreamOff,// offset for mate1
	            uint32_t dnstreamOff,// offset for mate2
	            uint32_t tlen,   // length of ref
	            bool fwa,        // whether mate1 is in fw orientation
	            bool fwb,        // whether mate2 is in fw orientation
	            bool pairFw)    // whether the pair is being mapped to fw strand
	{
		assert_lt(upstreamOff, dnstreamOff);
		uint32_t spreada = ra.bot - ra.top;
		uint32_t spreadb = rb.bot - rb.top;
		uint32_t oms = min(spreada, spreadb) - 1;
		bool ret;
		if(pairFw) {
			params_.setFw(fwa);
			ret = params_.reportHit(
					fwa ?  bufa_->patFw  :  bufa_->patRc,
					fwa ? &bufa_->qualFw : &bufa_->qualRc,
					&bufa_->name,
					ra.mms,                       // mismatch positions
					ra.refcs,                     // reference characters for mms
					ra.numMms,                    // # mismatches
					make_pair(first, upstreamOff),// position
					make_pair(ra.top, ra.bot),    // arrows
					tlen,                         // textlen
					alen_,                        // qlen
					ra.stratum,                   // alignment stratum
					oms);                         // # other hits
			assert(!ret);
			params_.setFw(fwb);
			ret = params_.reportHit(
					fwb ?  bufb_->patFw  :  bufb_->patRc,
					fwb ? &bufb_->qualFw : &bufb_->qualRc,
					&bufb_->name,
					rb.mms,                       // mismatch positions
					rb.refcs,                     // reference characters for mms
					rb.numMms,                    // # mismatches
					make_pair(first, dnstreamOff),// position
					make_pair(rb.top, rb.bot),    // arrows
					tlen,                         // textlen
					blen_,                        // qlen
					rb.stratum,                   // alignment stratum
					oms);                         // # other hits
		} else {
			params_.setFw(fwb);
			ret = params_.reportHit(
					fwb ?  bufb_->patFw  :  bufb_->patRc,
					fwb ? &bufb_->qualFw : &bufb_->qualRc,
					&bufb_->name,
					rb.mms,                       // mismatch positions
					rb.refcs,                     // reference characters for mms
					rb.numMms,                    // # mismatches
					make_pair(first, upstreamOff),// position
					make_pair(rb.top, rb.bot),    // arrows
					tlen,                         // textlen
					blen_,                        // qlen
					rb.stratum,                   // alignment stratum
					oms);                         // # other hits
			assert(!ret);
			params_.setFw(fwa);
			ret = params_.reportHit(
					fwa ?  bufa_->patFw  :  bufa_->patRc,
					fwa ? &bufa_->qualFw : &bufa_->qualRc,
					&bufa_->name,
					ra.mms,                       // mismatch positions
					ra.refcs,                     // reference characters for mms
					ra.numMms,                    // # mismatches
					make_pair(first, dnstreamOff),// position
					make_pair(ra.top, ra.bot),    // arrows
					tlen,                         // textlen
					alen_,                        // qlen
					ra.stratum,                   // alignment stratum
					oms);                         // # other hits
		}
		return ret;
	}

	/**
	 * Given a new reference position where one mate aligns, add it
	 * to the list of positions for that mate and reconcile it against
	 * all of the positions for the other mate, reporting paired
	 * alignments wherever the mating constraints are met.
	 *
	 * The parameters ending in 1 pertain to the upstream mate.
	 */
	bool reconcileAndAdd(const U32Pair& h,
	                     bool newFrom1,
	                     U32PairVec& offs1,
	                     const U32PairVec& offs2,
		                 TDriver& dr1,
		                 TDriver& dr2,
		                 bool fwa,
		                 bool fwb,
		                 bool pairFw,
		                 bool verbose = false)
	{
		// For each known hit for the other mate, check if this new
		// alignment can be mated with it.  If so, report the mates.
		size_t offsDstSz = newFrom1 ? offs2.size() : offs1.size();
		if(verbose) cout << "reconcileAndAdd called" << endl;
		if(offsDstSz > 0) {
			// Start in a random spot in the offs2 array and scan
			// linearly
			uint32_t rand = rand_.nextU32() % offsDstSz;
			for(size_t i = 0; i < offsDstSz; i++) {
				rand++;
				if(rand == offsDstSz) rand = 0;
				const U32Pair& h2 = newFrom1 ? offs2[rand] : offs1[rand];
				if(h.first == h2.first) {
					if(verbose) cout << "Found pair on some reference" << endl;
					// Incoming hit hits same reference as buffered
					uint32_t left  = newFrom1 ? h.second : h2.second;
					uint32_t right = newFrom1 ? h2.second : h.second;
					if(right > left) {
						if(verbose) cout << "...and in right orientation" << endl;
						uint32_t gap = right - left;
						if(gap >= minInsert_ && gap <= maxInsert_) {
							if(verbose) cout << "...and with the right sized gap; reporting" << endl;
							// Gap between the two alignments satisfies
							// the paired-end policy, so we can report
							// them
							if(report(
							       dr1.range(), dr2.range(),
								   h.first,
								   left, right,
								   this->ebwtFw_->_plen[h.first],
								   fwa, fwb,
								   pairFw)) return true;
						}
					}
				}
			}
		}
		offs1.push_back(h);
		return false;
	}

	/**
	 * Advance paired-end alignment where both reads are in their
	 * forward orientation.
	 */
	void advanceOrientation(bool& chase1,
	                        bool& chase2,
	                        bool& delayedChase1,
	                        bool& delayedChase2,
	                        TDriver& dr1,
	                        TDriver& dr2,
	                        U32PairVec& offs1,
	                        U32PairVec& offs2,
	                        uint32_t& offs1sz,
	                        uint32_t& offs2sz,
	                        bool& done,
	                        bool fw1,
	                        bool fw2,
	                        bool pairFw,
	                        bool verbose = false)
	{
		assert(!done_);
		assert(!done);
		assert(!chase1 || !chase2);
		if(chase1) {
			assert(!rangeMode_);
			assert(!delayedChase1);
			assert(dr1.foundRange());
			if(!rchase_.foundOff() && !rchase_.done()) {
				// Keep trying to resolve the reference loci for
				// alignments in this range
				rchase_.advance();
				return;
			} else if(rchase_.foundOff()) {
				// Resolve this against the reference loci
				// determined for the other mate
				done_ = reconcileAndAdd(rchase_.off(), true,
				                        offs1, offs2, dr1, dr2,
				                        fw1, fw2, pairFw, verbose);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				chase1 = false;
				if(verbose) cout << "Done with case for first mate" << endl;
				if(delayedChase2) {
					// Start chasing the delayed range
					if(verbose) cout << "Resuming delayed chase for second mate" << endl;
					assert(dr2.foundRange());
					rchase_.setTopBot(dr2.range().top, dr2.range().bot, dr2.qlen());
					chase2 = true;
					delayedChase2 = false;
				}
			}
		} else if(chase2) {
			assert(!rangeMode_);
			assert(!delayedChase2);
			assert(dr2.foundRange());
			if(!rchase_.foundOff() && !rchase_.done()) {
				// Keep trying to resolve the reference loci for
				// alignments in this range
				rchase_.advance();
				return;
			} else if(rchase_.foundOff()) {
				// Resolve this against the reference loci
				// determined for the other mate
				done_ = reconcileAndAdd(rchase_.off(), false,
				                        offs1, offs2, dr1, dr2,
				                        fw1, fw2, pairFw, verbose);
				rchase_.reset();
			} else {
				assert(rchase_.done());
				// Forget this range; keep looking for ranges
				chase2 = false;
				if(verbose) cout << "Done with case for second mate" << endl;
				if(delayedChase1) {
					// Start chasing the delayed range
					if(verbose) cout << "Resuming delayed chase for first mate" << endl;
					assert(dr1.foundRange());
					rchase_.setTopBot(dr1.range().top, dr1.range().bot, dr1.qlen());
					chase1 = true;
					delayedChase1 = false;
				}
			}
		}
		if(!done_ && !done && !chase1 && !chase2) {
			// Search for more ranges for whichever mate currently has
			// fewer ranges
			if((offs1sz < offs2sz || dr2.done()) && !dr1.done()) {
				// If there are no more ranges for the other mate and
				// there are no candidate alignments either, then we're
				// not going to find a paired alignment in this
				// orientation.
				assert(!delayedChase1);
				if(dr2.done() && offs2sz == 0) {
					// Give up on this orientation
					if(verbose) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 1" << endl;
					done = true;
					return;
				}
				dr1.advance();
				if(dr1.foundRange()) {
					// Add the size of this range to the total for this mate
					offs1sz += (dr1.range().bot - dr1.range().top);
					if(offs2sz == 0) {
						// Delay chasing this range; we delay to avoid
						// needlessly chasing rows in this range when
						// the other mate doesn't end up aligning
						// anywhere
						if(verbose) cout << "Delaying a chase for first mate" << endl;
						delayedChase1 = true;
					} else {
						// Start chasing this range
						if(verbose) cout << "Chasing a range for first mate" << endl;
						rchase_.setTopBot(dr1.range().top, dr1.range().bot, dr1.qlen());
						chase1 = true;
					}
				}
			} else if(!dr2.done()) {
				// If there are no more ranges for the other mate and
				// there are no candidate alignments either, then we're
				// not going to find a paired alignment in this
				// orientation.
				assert(!delayedChase2);
				if(dr1.done() && offs1sz == 0) {
					// Give up on this orientation
					if(verbose) cout << "Giving up on paired orientation " << (pairFw? "fw" : "rc") << " in mate 2" << endl;
					done = true;
					return;
				}
				dr2.advance();
				if(dr2.foundRange()) {
					// Add the size of this range to the total for this mate
					offs2sz += (dr2.range().bot - dr2.range().top);
					if(offs1sz == 0) {
						// Delay chasing this range; we delay to avoid
						// needlessly chasing rows in this range when
						// the other mate doesn't end up aligning
						// anywhere
						if(verbose) cout << "Delaying a chase for second mate" << endl;
						delayedChase2 = true;
					} else {
						// Start chasing this range
						if(verbose) cout << "Chasing a range for second mate" << endl;
						rchase_.setTopBot(dr2.range().top, dr2.range().bot, dr2.qlen());
						chase2 = true;
					}
				}
			} else {
				// Finished processing ranges for both mates
				assert(dr1.done() && dr2.done());
				done = true;
			}
		}
	}

	// Progress state
	bool doneFw_;   // finished with forward orientation of both mates?
	bool done_;
	bool firstFw_;
	bool firstRc_;

	bool chase1Fw_;
	bool chase1Rc_;
	bool chase2Fw_;
	bool chase2Rc_;

	bool delayedChase1Fw_;
	bool delayedChase1Rc_;
	bool delayedChase2Fw_;
	bool delayedChase2Rc_;

	// Temporary HitSink; to be deleted
	const HitSinkPerThreadFactory& sinkPtFactory_;
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> > params_;

	// Paired-end boundaries
	const uint32_t minInsert_;
	const uint32_t maxInsert_;

	// Orientation of upstream/downstream mates when aligning to
	// forward strand
	const bool fw1_;
	const bool fw2_;

	// State for getting alignments from ranges statefully
	RandomScanningRangeChaser<String<Dna> > rchase_;

	// Range-finding state for first mate
	TRangeSource r1Fw_;
	TContMan     c1Fw_;
	TRangeSource r1Rc_;
	TContMan     c1Rc_;
	TDriver      driver1Fw_;
	TDriver      driver1Rc_;
	U32PairVec   offs1Fw_;
	uint32_t     offs1FwSz_; // total size of all ranges found in this category
	U32PairVec   offs1Rc_;
	uint32_t     offs1RcSz_; // total size of all ranges found in this category

	// Range-finding state for second mate
	TRangeSource r2Fw_;
	TContMan     c2Fw_;
	TRangeSource r2Rc_;
	TContMan     c2Rc_;
	TDriver      driver2Fw_;
	TDriver      driver2Rc_;
	U32PairVec   offs2Fw_;
	uint32_t     offs2FwSz_; // total size of all ranges found in this category
	U32PairVec   offs2Rc_;
	uint32_t     offs2RcSz_; // total size of all ranges found in this category
};

#endif /* ALIGNER_H_ */
