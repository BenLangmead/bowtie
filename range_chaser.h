/*
 * range_chaser.h
 */

#ifndef RANGE_CHASER_H_
#define RANGE_CHASER_H_

#include <stdint.h>
#include "ebwt.h"
#include "random_source.h"
#include "row_chaser.h"

/**
 * Abstract parent for classes that report alignments from ranges in a
 * stateful manner.
 */
template<typename TStr>
class RangeChaser {
	typedef std::pair<uint32_t,uint32_t> U32Pair;
public:
	RangeChaser() : prepped_(false) { }
	virtual ~RangeChaser() { }
	virtual void setTopBot(uint32_t top, uint32_t bot, uint32_t qlen) = 0;
	virtual bool done() const = 0;
	virtual void advance() = 0;
	virtual void prep() = 0;
	virtual bool foundOff() const = 0;
	virtual U32Pair off() const = 0;
	bool prepped_; /// true = prefetch is issued and it's OK to call advance()
};

/**
 * A class that statefully processes a range by picking one row
 * randomly and then linearly scanning forward through the range,
 * reporting reference offsets as we go.
 */
template<typename TStr>
class RandomScanningRangeChaser : public RangeChaser<TStr> {

	typedef Ebwt<TStr> EbwtT;
	typedef std::pair<uint32_t,uint32_t> U32Pair;
	typedef std::vector<U32Pair> U32PairVec;

public:
	RandomScanningRangeChaser(const EbwtT& ebwt, RandomSource& rand) :
		RangeChaser<TStr>(),
		ebwt_(ebwt),
		qlen_(0),
		rand_(rand),
		top_(0xffffffff),
		bot_(0xffffffff),
		irow_(0xffffffff),
		row_(0xffffffff),
		done_(false),
		off_(make_pair(0xffffffff, 0)),
		tlen_(0),
		chaser_(ebwt)
	{ }

	virtual ~RandomScanningRangeChaser() { }

	/**
	 * Convert a range to a vector of reference loci, where a locus is
	 * a u32 pair of <ref-idx, ref-offset>.
	 */
	static void toOffs(const EbwtT& ebwt,
	                   uint32_t qlen,
	                   RandomSource& rand,
	                   uint32_t top,
	                   uint32_t bot,
	                   U32PairVec& dest)
	{
		RandomScanningRangeChaser rc(ebwt, rand);
		rc.setTopBot(top, bot, qlen);
		rc.prep();
		while(!rc.done()) {
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
	void setRow(uint32_t row) {
		// Must be within bounds of range
		assert_lt(row, bot_);
		assert_geq(row, top_);
		row_ = row;
		while(true) {
			// Set up the chaser
			chaser_.setRow(row_, qlen_);
			// It might be done immediately...
			if(chaser_.done()) {
				// We're done immediately
				off_ = chaser_.off();
				if(off_.first != 0xffffffff) {
					// This is a valid result
					tlen_ = chaser_.tlen();
					assert(foundOff());
					return; // found result
				}
				// That row didn't have a valid result, move to the next
				row_++;
				if(row_ == bot_) {
					// Wrap back to top_
					row_ = top_;
				}
				if(row_ == irow_) {
					// Exhausted all possible rows
					done_ = true;
					assert_eq(0xffffffff, off_.first);
					return;
				}
			} else {
				// Pursue this row
				break;
			}
		}
	}

	/**
	 * Set the next range for us to "chase" (i.e. convert row-by-row
	 * to reference loci).
	 */
	virtual void setTopBot(uint32_t top, uint32_t bot, uint32_t qlen) {
		assert_neq(0xffffffff, top);
		assert_neq(0xffffffff, bot);
		assert_gt(bot, top);
		assert_gt(qlen, 0);
		qlen_ = qlen;
		top_ = top;
		bot_ = bot;
		uint32_t spread = bot - top;
		irow_ = top + (rand_.nextU32() % spread); // initial row
		done_ = false;
		reset();
		setRow(irow_);
		assert(chaser_.prepped_);
	}

	/**
	 * Return true iff off_ now holds the reference location
	 * corresponding to the row last set with setRow().
	 */
	virtual bool done() const {
		return done_;
	}

	/**
	 * Advance the step-left process by one step.  Check if we're done.
	 */
	virtual void advance() {
		assert(!done_);
		assert(chaser_.prepped_);
		reset();
		if(chaser_.done()) {
			// chaser finished with this row
			row_++;
			if(row_ == bot_) {
				row_ = top_;
			}
			if(row_ == irow_) {
				// Exhausted all possible rows
				done_ = true;
				assert_eq(0xffffffff, off_.first);
				return;
			}
			setRow(row_);
			assert(chaser_.prepped_);
		} else {
			chaser_.advance();
			assert(chaser_.prepped_);
			if(chaser_.done()) {
				// We're done immediately
				off_ = chaser_.off();
				if(off_.first != 0xffffffff) {
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
	virtual void prep() {
		// nothing
	}

	/**
	 * Return true iff off_ contains a valid reference location for
	 * this range.
	 */
	virtual bool foundOff() const {
		return off_.first != 0xffffffff;
	}

	/**
	 * Reset the chaser so that 'off_' does not hold a valid offset and
	 * foundOff() returns false.
	 */
	virtual void reset() {
		off_.first = 0xffffffff;
	}

	/**
	 * Get the calculated offset.
	 */
	virtual U32Pair off() const {
		return off_;
	}

	/**
	 * Get the length of the hit reference.
	 */
	uint32_t tlen() const {
		return tlen_;
	}

protected:

	const EbwtT& ebwt_;  /// index to resolve row in
	uint32_t qlen_;      /// length of read; needed to convert to ref. coordinates
	RandomSource& rand_; /// pseudo-random number generator
	uint32_t top_;       /// range top
	uint32_t bot_;       /// range bottom
	uint32_t irow_;      /// initial randomly-chosen row within range
	uint32_t row_;       /// current row within range
	bool done_;          /// true = chase is done & answer is in off_
	U32Pair off_;        /// calculated offset (0xffffffff if not done)
	uint32_t tlen_;      /// length of text hit
	RowChaser<TStr> chaser_; /// stateful row chaser
};

#endif /* RANGE_CHASER_H_ */
