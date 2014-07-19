/*
 * range_chaser.h
 */

#ifndef RANGE_CHASER_H_
#define RANGE_CHASER_H_

#include <stdint.h>
#include "ebwt.h"
#include "random_source.h"
#include "row_chaser.h"
#include "range_cache.h"
#include "aligner_metrics.h"

/**
 * A class that statefully processes a range by picking one row
 * randomly and then linearly scanning forward through the range,
 * reporting reference offsets as we go.
 */
template<typename TStr>
class RangeChaser {

	typedef Ebwt<TStr> TEbwt;
	typedef std::pair<TIndexOffU,TIndexOffU> UPair;
	typedef std::vector<UPair> UPairVec;
	typedef RowChaser<TStr> TRowChaser;

public:
	RangeChaser(uint32_t cacheThresh,
	            RangeCache* cacheFw, RangeCache* cacheBw,
	            AlignerMetrics *metrics = NULL) :
		done(false),
		ebwt_(NULL),
		qlen_(0),
		cacheThresh_(cacheThresh),
		top_(OFF_MASK),
		bot_(OFF_MASK),
		irow_(OFF_MASK),
		row_(OFF_MASK),
		off_(make_pair(OFF_MASK, 0)),
		tlen_(0),
		chaser_(metrics),
		cached_(false),
		cacheFw_(cacheFw), cacheBw_(cacheBw),
		metrics_(metrics)
	{ }

	~RangeChaser() { }

	/**
	 * Convert a range to a vector of reference loci, where a locus is
	 * a u32 pair of <ref-idx, ref-offset>.
	 */
	static void toOffs(const TEbwt& ebwt,
			TIndexOffU qlen,
	                   RandomSource& rand,
	                   TIndexOffU top,
	                   TIndexOffU bot,
	                   UPairVec& dest)
	{
		RangeChaser rc(ebwt, rand);
		rc.setTopBot(top, bot, qlen);
		rc.prep();
		while(!rc.done) {
			rc.advance();
			if(rc.foundOff()) {
				dest.push_back(rc.off());
				rc.reset();
				assert(!rc.foundOff());
			}
			rc.prep();
		}
	}

	/**
	 * Set the row to chase
	 */
	void setRow(TIndexOffU row) {
		// Must be within bounds of range
		assert_lt(row, bot_);
		assert_geq(row, top_);
		row_ = row;
		while(true) {
			// First thing to try is the cache
			if(cached_) {
				assert(cacheEnt_.valid());
				TIndexOffU cached = cacheEnt_.get(row_ - top_);
				assert(cacheEnt_.valid());
				if(cached != RANGE_NOT_SET) {
					// Assert that it matches what we would have got...
					ASSERT_ONLY(uint32_t sanity = TRowChaser::toFlatRefOff(ebwt_, 1, row_));
					assert_eq(sanity, cached);
					// We have a cached result.  Cached result is in the
					// form of an offset into the joined reference string,
					// so now we have to convert it to a tidx/toff pair.
					ebwt_->joinedToTextOff(qlen_, cached, off_.first, off_.second, tlen_);
					// Note: tidx may be 0xffffffff, if alignment overlaps a
					// reference boundary
					if(off_.first != OFF_MASK) {
						// Bingo, we found a valid result using the cache
						assert(foundOff());
						return;
					}
				} else {
					// Wasn't in the cache; use the RowChaser
				}
			}
			// Second thing to try is the chaser
			chaser_.setRow(row_, qlen_, ebwt_);
			assert(chaser_.prepped_ || chaser_.done);
			// It might be done immediately...
			if(chaser_.done) {
				// We're done immediately
				off_ = chaser_.off();
				if(off_.first != OFF_MASK) {
					// This is a valid result
					if(cached_) {
						// Install the result in the cache
						assert(cacheEnt_.valid());
						cacheEnt_.install(row_ - top_, chaser_.flatOff());
						//if(ebwt_->fw()) assert(cacheFw_->repOk());
						//else            assert(cacheBw_->repOk());
					}
					tlen_ = chaser_.tlen();
					assert(foundOff());
					return; // found result
				}
			} else {
				// Pursue this row
				break;
			}
			// That row didn't have a valid result, move to the next
			row_++;
			if(row_ == bot_) {
				// Wrap back to top_
				row_ = top_;
			}
			if(row_ == irow_) {
				// Exhausted all possible rows
				done = true;
				assert_eq(OFF_MASK, off_.first);
				return;
			}
		}
		assert(chaser_.prepped_);
	}

	/**
	 * Set the next range for us to "chase" (i.e. convert row-by-row
	 * to reference loci).
	 */
	void setTopBot(TIndexOffU top,
					TIndexOffU bot,
					TIndexOffU qlen,
	               RandomSource& rand,
	               const TEbwt* ebwt)
	{
		assert_neq(OFF_MASK, top);
		assert_neq(OFF_MASK, bot);
		assert_gt(bot, top);
		assert_gt(qlen, 0);
		assert(ebwt != NULL);
		ebwt_ = ebwt;
		qlen_ = qlen;
		top_ = top;
		bot_ = bot;
		TIndexOffU spread = bot - top;
		irow_ = top + (rand.nextU32() % spread); // initial row
		done = false;
		cached_ = false;
		reset();
		if(cacheFw_ != NULL || cacheBw_ != NULL) {
			if(spread > cacheThresh_) {
				bool ret = false;
				if(ebwt->fw() && cacheFw_ != NULL) {
					ret = cacheFw_->lookup(top, bot, cacheEnt_);
					if(ret) assert(cacheEnt_.ebwt()->fw());
				} else if(!ebwt->fw() && cacheBw_ != NULL) {
					ret = cacheBw_->lookup(top, bot, cacheEnt_);
					if(ret) assert(!cacheEnt_.ebwt()->fw());
				} else {
					cacheEnt_.reset();
				}
				assert_eq(cacheEnt_.valid(), ret);
				cached_ = ret;
			} else {
				cacheEnt_.reset();
			}
		}
		setRow(irow_);
		assert(chaser_.prepped_ || foundOff() || done);
	}

	/**
	 * Advance the step-left process by one step.  Check if we're done.
	 */
	void advance() {
		assert(!done);
		assert(chaser_.prepped_ || chaser_.done);
		reset();
		if(chaser_.done) {
			// chaser finished with this row
			row_++;
			if(row_ == bot_) {
				row_ = top_;
			}
			if(row_ == irow_) {
				// Exhausted all possible rows
				done = true;
				assert_eq(OFF_MASK, off_.first);
				return;
			}
			setRow(row_);
			assert(chaser_.prepped_ || foundOff() || done);
		} else {
			chaser_.advance();
			assert(chaser_.prepped_ || chaser_.done);
			if(chaser_.done) {
				// We're done immediately
				off_ = chaser_.off();
				if(off_.first != OFF_MASK) {
					if(cached_) {
						// Install the result in the cache
						assert(cacheEnt_.valid());
						cacheEnt_.install(row_ - top_, chaser_.flatOff());
						//if(ebwt_->fw()) assert(cacheFw_->repOk());
						//else            assert(cacheBw_->repOk());
					}
					// Found a reference position
					tlen_ = chaser_.tlen();
					assert(foundOff());
				}
			}
		}
	}

	/**
	 * Prepare for the next call to advance() by prefetching relevant
	 * data.  In this case, 'chaser_' is doing this for us.
	 */
	void prep() {
		// nothing
	}

	/**
	 * Return true iff off_ contains a valid reference location for
	 * this range.
	 */
	bool foundOff() const {
		return off_.first != OFF_MASK;
	}

	/**
	 * Reset the chaser so that 'off_' does not hold a valid offset and
	 * foundOff() returns false.
	 */
	void reset() {
		off_.first = OFF_MASK;
	}

	/**
	 * Get the calculated offset.
	 */
	UPair off() const {
		return off_;
	}

	/**
	 * Get the length of the hit reference.
	 */
	TIndexOffU tlen() const {
		return tlen_;
	}

	bool done;             /// true = chase is done & answer is in off_

protected:

	const TEbwt* ebwt_;    /// index to resolve row in
	TIndexOffU qlen_;        /// length of read; needed to convert to ref. coordinates
	uint32_t cacheThresh_; /// ranges wider than thresh use cacheing
	TIndexOffU top_;         /// range top
	TIndexOffU bot_;         /// range bottom
	TIndexOffU irow_;        /// initial randomly-chosen row within range
	TIndexOffU row_;         /// current row within range
	UPair off_;          /// calculated offset (OFF_MASK if not done)
	TIndexOffU tlen_;        /// length of text hit
	TRowChaser chaser_;    /// stateful row chaser
	RangeCacheEntry cacheEnt_; /// current cache entry
	bool cached_;          /// cacheEnt is active for current range?
	RangeCache* cacheFw_; /// cache for the forward index
	RangeCache* cacheBw_; /// cache for the backward index
	AlignerMetrics *metrics_;
};

#endif /* RANGE_CHASER_H_ */
