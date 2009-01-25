/*
 * range_source.h
 *
 *  Created on: Jan 21, 2009
 *      Author: langmead
 */

#ifndef RANGE_SOURCE_H_
#define RANGE_SOURCE_H_

#include <stdint.h>
#include <vector>
#include "seqan/sequence.h"
#include "ebwt.h"
#include "range.h"

/**
 * Encapsulates an algorithm that navigates the Bowtie index to produce
 * candidate ranges of alignments in the Burrows-Wheeler matrix.  A
 * higher authority is responsible for reporting hits out of those
 * ranges, and stopping when the consumer is satisfied.
 */
template<typename TContinuationManager>
class RangeSource {
public:
	RangeSource()  { }
	virtual ~RangeSource() { }

	/// Set query to find ranges for
	virtual void setQuery(seqan::String<Dna5>* qry,
	                      seqan::String<char>* qual,
	                      seqan::String<char>* name) = 0;
	/// Set up the range search.
	virtual void initConts(TContinuationManager& conts, uint32_t ham) = 0;
	/// Advance the range search by one memory op
	virtual void advance(TContinuationManager& conts) = 0;
	/// Returns true iff the last call to advance yielded a range
	virtual bool foundRange() = 0;
	/// Return the last valid range found
	virtual Range& range() = 0;
	/// All searching w/r/t the current query is finished
	virtual bool done() = 0;
};

/**
 * Encapsulates an algorithm that navigates the Bowtie index to produce
 * candidate ranges of alignments in the Burrows-Wheeler matrix.  A
 * higher authority is responsible for reporting hits out of those
 * ranges, and stopping when the consumer is satisfied.
 */
template<typename TContinuationManager>
class ListRangeSource : public RangeSource<TContinuationManager> {

	typedef std::vector<RangeSource<TContinuationManager>*> TRangeSourcePtrVec;

public:

	ListRangeSource(const TRangeSourcePtrVec& rss) :
		cur_(0), done_(false), ham_(0), rss_(rss)
	{
		assert_gt(rss_.size(), 0);
	}

	virtual ~ListRangeSource() { }

	/// Set query to find ranges for
	virtual void setQuery(seqan::String<Dna5>* qry,
	                      seqan::String<char>* qual,
	                      seqan::String<char>* name)
	{
		cur_ = 0;
		done_ = false;
		rss_[0]->setQuery(qry, qual, name);
	}

	/// Set up the range search.
	virtual void initConts(TContinuationManager& conts, uint32_t ham) {
		assert_eq(0, cur_);
		ham_ = ham;
		rss_[0]->initConts(conts, ham);
	}

	/// Advance the range search by one memory op
	virtual void advance(TContinuationManager& conts) {
		assert(!done_);
		if(rss_[cur_]->done()) {
			assert(conts.empty());
			if(cur_ < rss_.size()-1) {
				rss_[++cur_]->initConts(conts, ham_);
			} else {
				done_ = true;
				return;
			}
		} else {
			rss_[cur_]->advance(conts);
		}
	}

	/// Returns true iff the last call to advance yielded a range
	virtual bool foundRange() { return rss_[cur_]->foundRange(); }
	/// Return the last valid range found
	virtual Range& range() { return rss_[cur_]->range(); }
	/// All searching w/r/t the current query is finished
	virtual bool done() { return done_; }

protected:

	uint32_t cur_;
	bool done_;
	uint32_t ham_;
	TRangeSourcePtrVec rss_;
};


template<typename TCont>
class ContinuationManager {
public:
	virtual ~ContinuationManager() { }
	virtual TCont& front() = 0;
	virtual TCont& back() = 0;
	virtual void   prep(const EbwtParams& ep, const uint8_t* ebwt) = 0;
	virtual void   pop() = 0;
	virtual void   expand1() = 0;
	virtual size_t size() const = 0;
	virtual bool   empty() const = 0;
	virtual void   clear() = 0;
};

/**
 * An aligner for finding exact matches of unpaired reads.  Always
 * tries the forward-oriented version of the read before the reverse-
 * oriented read.
 */
template<typename TRangeSource, typename TContMan>
class RangeSourceDriver {

	typedef Ebwt<String<Dna> > EbwtT;

public:
	RangeSourceDriver(
		const EbwtT& ebwtFw,
		const EbwtT* ebwtBw,
		EbwtSearchParams<String<Dna> >& params,
		TRangeSource& rs,
		TContMan& cm,
		bool fw,
		HitSink& sink,
		HitSinkPerThread* sinkPt,
		vector<String<Dna5> >& os,
		bool verbose,
		uint32_t seed) :
		done_(true), first_(true), len_(0),
		pat_(NULL), qual_(NULL), name_(NULL),
		sinkPt_(sinkPt), params_(params),
		fw_(fw), ebwtFw_(ebwtFw), ebwtBw_(ebwtBw), rs_(rs), cm_(cm)
	{
		assert(cm_.empty());
	}

	virtual ~RangeSourceDriver() { }

	/**
	 * Prepare this aligner for the next read.
	 */
	void setQuery(PatternSourcePerThread* patsrc, bool mate1 = true) {
		if(mate1) {
			if(fw_) {
				pat_  = &patsrc->bufa().patFw;
				qual_ = &patsrc->bufa().qualFw;
			} else {
				pat_  = &patsrc->bufa().patRc;
				qual_ = &patsrc->bufa().qualRc;
			}
			name_ = &patsrc->bufa().name;
		} else {
			if(fw_) {
				pat_  = &patsrc->bufb().patFw;
				qual_ = &patsrc->bufb().qualFw;
			} else {
				pat_  = &patsrc->bufb().patRc;
				qual_ = &patsrc->bufb().qualRc;
			}
			name_ = &patsrc->bufb().name;
		}
		assert(pat_ != NULL);
		assert(qual_ != NULL);
		assert(name_ != NULL);
		len_ = seqan::length(*pat_);
		assert_gt(len_, 0);
		done_ = false;
		first_ = true;
		cm_.clear();
		rs_.setQuery(pat_, qual_, name_);
		initRangeSource(rs_);
		rs_.initConts(cm_, 0); // set up initial continuation
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	void advance() {
		assert(!done_);
		assert(pat_ != NULL);
		params_.setFw(fw_);
		if(first_) {
			// Set up the RangeSource for the forward read
			first_ = false;
		} else {
			// Advance the RangeSource for the forward-oriented read
			rs_.advance(cm_);
		}
		if(!done_) {
			// Finished
			done_ = cm_.empty();
		}
		if(!done_) {
			// Hopefully, this will prefetch enough so that when we
			// resume, stuff is already in cache
			cm_.prep(ebwtFw_.eh(), ebwtFw_.ebwt());
		}
		return;
	}

	/**
	 * Return true iff we just found a range.
	 */
	bool foundRange() {
		return rs_.foundRange();
	}

	/**
	 * Return the range found.
	 */
	Range& range() {
		return rs_.range();
	}

	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	bool done() {
		return done_;
	}

	/// Return length of current query
	uint32_t qlen() {
		return len_;
	}

	/**
	 * Return a const ptr to whichever Ebwt (forward or backward) the
	 * current range came from.
	 */
	virtual const EbwtT *curEbwt() {
		return &ebwtFw_;
	}

protected:

	virtual void initRangeSource(TRangeSource& rs) = 0;

	// Progress state
	bool done_;
	bool first_;
	uint32_t len_;
	String<Dna5>* pat_;
	String<char>* qual_;
	String<char>* name_;

	// Temporary HitSink; to be deleted
	HitSinkPerThread* sinkPt_;

	// State for alignment
	EbwtSearchParams<String<Dna> >& params_;
	bool                            fw_;
	const EbwtT&                    ebwtFw_;
	const EbwtT*                    ebwtBw_;
	TRangeSource&                   rs_;
	TContMan&                       cm_;
};

#endif /* RANGE_SOURCE_H_ */
