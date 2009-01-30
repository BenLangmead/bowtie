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
	typedef Ebwt<String<Dna> > TEbwt;
public:
	RangeSource() : curEbwt_(NULL) { }
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
	/// Return ptr to index this RangeSource is currently getting ranges from
	const TEbwt *curEbwt() const { return curEbwt_; }
protected:
	/// ptr to index this RangeSource is currently getting ranges from
	const TEbwt *curEbwt_;
};

/**
 * Encapsulates an algorithm that navigates the Bowtie index to produce
 * candidate ranges of alignments in the Burrows-Wheeler matrix.  A
 * higher authority is responsible for reporting hits out of those
 * ranges, and stopping when the consumer is satisfied.
 */
template<typename TContinuationManager>
class ListRangeSource : public RangeSource<TContinuationManager> {

	typedef Ebwt<String<Dna> > TEbwt;
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
		cur_ = 0; // go back to first RangeSource in list
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
		assert_lt(cur_, rss_.size());
		if(rss_[cur_]->done()) {
			// Move on to next RangeSource
			assert(conts.empty());
			if(cur_ < rss_.size()-1) {
				conts.clear();
				rss_[++cur_]->initConts(conts, ham_);
			} else {
				// No RangeSources in list; done
				done_ = true;
				return;
			}
		} else {
			// Advance current RangeSource
			rss_[cur_]->advance(conts);
		}
	}

	/// Returns true iff the last call to advance yielded a range
	virtual bool foundRange() { assert(!done_); return rss_[cur_]->foundRange(); }
	/// Return the last valid range found
	virtual Range& range() { return rss_[cur_]->range(); }
	/// All searching w/r/t the current query is finished
	virtual bool done() { return done_; }
	/// Return curEbwt for the currently-active RangeSource
	const TEbwt* curEbwt() const {
		return rss_[cur_]->curEbwt();
	}

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
 * Abstract parent of RangeSourceDrivers
 */
template<typename TRangeSource, typename TContMan>
class RangeSourceDriver {
	typedef Ebwt<String<Dna> > TEbwt;
public:
	virtual ~RangeSourceDriver() { }
	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc, bool mate1 = true) = 0;
	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advance() = 0;
	/**
	 * Return true iff we just found a range.
	 */
	virtual bool foundRange() const = 0;
	/**
	 * Return the range found.
	 */
	virtual Range& range() = 0;
	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() const = 0;
	/// Return ptr to index this RangeSource is currently getting ranges from
	virtual const TEbwt* curEbwt() const = 0;

	virtual uint32_t qlen() const = 0;
protected:

	virtual void initRangeSource(TRangeSource& rs) = 0;
};

/**
 * An aligner for finding exact matches of unpaired reads.  Always
 * tries the forward-oriented version of the read before the reverse-
 * oriented read.
 */
template<typename TRangeSource, typename TContMan>
class SingleRangeSourceDriver : public RangeSourceDriver<TRangeSource, TContMan> {

	typedef Ebwt<String<Dna> > TEbwt;

public:
	SingleRangeSourceDriver(
		EbwtSearchParams<String<Dna> >& params,
		TRangeSource* rs,
		bool fw,
		HitSink& sink,
		HitSinkPerThread* sinkPt,
		vector<String<Dna5> >& os,
		bool verbose,
		uint32_t seed) :
		done_(true), first_(true), len_(0),
		pat_(NULL), qual_(NULL), name_(NULL), sinkPt_(sinkPt),
		params_(params),
		fw_(fw), rs_(rs), cm_()
	{
		assert(rs_ != NULL);
	}

	virtual ~SingleRangeSourceDriver() {
		delete rs_; rs_ = NULL;
	}

	/**
	 * Prepare this aligner for the next read.
	 */
	virtual void setQuery(PatternSourcePerThread* patsrc, bool mate1 = true) {
		if(mate1) {
			if(fw_) {
				pat_  = (curEbwt()->fw() ? &patsrc->bufa().patFw  : &patsrc->bufa().patFwRev);
				qual_ = (curEbwt()->fw() ? &patsrc->bufa().qualFw : &patsrc->bufa().qualFwRev);
			} else {
				pat_  = (curEbwt()->fw() ? &patsrc->bufa().patRc  : &patsrc->bufa().patRcRev);
				qual_ = (curEbwt()->fw() ? &patsrc->bufa().qualRc : &patsrc->bufa().qualRcRev);
			}
			name_ = &patsrc->bufa().name;
		} else {
			if(fw_) {
				pat_  = (curEbwt()->fw() ? &patsrc->bufb().patFw  : &patsrc->bufb().patFwRev);
				qual_ = (curEbwt()->fw() ? &patsrc->bufb().qualFw : &patsrc->bufb().qualFwRev);
			} else {
				pat_  = (curEbwt()->fw() ? &patsrc->bufb().patRc  : &patsrc->bufb().patRcRev);
				qual_ = (curEbwt()->fw() ? &patsrc->bufb().qualRc : &patsrc->bufb().qualRcRev);
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
		rs_->setQuery(pat_, qual_, name_);
		initRangeSource(*rs_);
		rs_->initConts(cm_, 0); // set up initial continuation
	}

	/**
	 * Advance the aligner by one memory op.  Return true iff we're
	 * done with this read.
	 */
	virtual void advance() {
		assert(!done_);
		assert(pat_ != NULL);
		params_.setFw(fw_);
		if(first_) {
			// Set up the RangeSource for the forward read
			first_ = false;
		} else {
			// Advance the RangeSource for the forward-oriented read
			rs_->advance(cm_);
		}
		if(!done_) {
			// Finished
			done_ = cm_.empty();
		}
		if(!done_) {
			// Hopefully, this will prefetch enough so that when we
			// resume, stuff is already in cache
			const TEbwt* ebwt = curEbwt();
			assert(ebwt != NULL);
			cm_.prep(ebwt->_eh, ebwt->_ebwt);
		}
		return;
	}

	/**
	 * Return true iff we just found a range.
	 */
	virtual bool foundRange() const {
		return rs_->foundRange();
	}

	/**
	 * Return the range found.
	 */
	virtual Range& range() {
		return rs_->range();
	}

	/**
	 * Returns true if all searching w/r/t the current query is
	 * finished or if there is no current query.
	 */
	virtual bool done() const {
		return done_;
	}

	/// Return length of current query
	virtual uint32_t qlen() const {
		return len_;
	}

	/// Return ptr to index this RangeSource is currently getting ranges from
	virtual const TEbwt* curEbwt() const {
		return rs_->curEbwt();
	}

	virtual void initRangeSource(TRangeSource& rs) = 0;

protected:

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
	TRangeSource*                   rs_; // delete this in destructor
	TContMan                        cm_;
};

/**
 * Encapsulates an algorithm that navigates the Bowtie index to produce
 * candidate ranges of alignments in the Burrows-Wheeler matrix.  A
 * higher authority is responsible for reporting hits out of those
 * ranges, and stopping when the consumer is satisfied.
 */
template<typename TRangeSource, typename TContMan>
class ListRangeSourceDriver : public RangeSourceDriver<TRangeSource, TContMan> {

	typedef Ebwt<String<Dna> > TEbwt;
	typedef std::vector<SingleRangeSourceDriver<TRangeSource, TContMan>*> TRangeSrcDrPtrVec;

public:

	ListRangeSourceDriver(const TRangeSrcDrPtrVec& rss) :
		RangeSourceDriver<TRangeSource, TContMan>(),
		cur_(0), done_(false), ham_(0), rss_(rss) /* copy */,
		patsrc_(NULL), mate1_(true)
	{
		assert_gt(rss_.size(), 0);
	}

	virtual ~ListRangeSourceDriver() {
		for(size_t i = 0; i < rss_.size(); i++) {
			delete rss_[i];
		}
	}

	/// Set query to find ranges for
	virtual void setQuery(PatternSourcePerThread* patsrc, bool mate1 = true) {
		cur_ = 0; // go back to first RangeSource in list
		done_ = false;
		rss_[0]->setQuery(patsrc, mate1);
		patsrc_ = patsrc; // so that we can call setQuery on the other elements later
		mate1_ = mate1;   // so that we can call setQuery on the other elements later
	}

	/// Advance the range search by one memory op
	virtual void advance() {
		assert(!done_);
		assert_lt(cur_, rss_.size());
		if(rss_[cur_]->done()) {
			// Move on to next RangeSourceDriver
			if(cur_ < rss_.size()-1) {
				rss_[++cur_]->setQuery(patsrc_, mate1_);
			} else {
				// No RangeSources in list; done
				done_ = true;
				return;
			}
		} else {
			// Advance current RangeSource
			rss_[cur_]->advance();
		}
	}

	/// Returns true iff the last call to advance yielded a range
	virtual bool foundRange() const { return rss_[cur_]->foundRange(); }

	/// Return the last valid range found
	virtual Range& range() { return rss_[cur_]->range(); }

	/// All searching w/r/t the current query is finished
	virtual bool done() const { return done_; }

	/// Return curEbwt for the currently-active RangeSource
	const TEbwt* curEbwt() const {
		return rss_[cur_]->curEbwt();
	}

	/// Return length of current query
	virtual uint32_t qlen() const {
		return rss_[cur_]->qlen();
	}

	/// Pass call to initRangeSource to current driver
	virtual void initRangeSource(TRangeSource& rs) {
		rss_[cur_]->initRangeSource(rs);
	}

protected:

	uint32_t cur_;
	bool done_;
	uint32_t ham_;
	TRangeSrcDrPtrVec rss_;
	PatternSourcePerThread* patsrc_;
	bool mate1_;
};



#endif /* RANGE_SOURCE_H_ */
