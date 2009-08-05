#include "hit.h"
#include "hit_set.h"

using namespace std;
using namespace seqan;

/// Sort by text-id then by text-offset
bool operator< (const Hit& a, const Hit& b) {
    return a.h < b.h;
}

/**
 * Report a batch of hits to a chaining file.
 */
void ChainingHitSink::reportHits(vector<Hit>& hs) {
	size_t hssz = hs.size();
	assert_gt(hssz, 0);
	assert_eq(0, hs[0].mate);
	// Convert vector<Hit> into HitSet
	{
		// Critical section for output stream 0
		HitSet s;
		Hit::toHitSet(hs, s, amap_);
		lock(0);
		s.serialize(out(0));
		unlock(0);
	}
	{
		// Global critical section
		mainlock();
		commitHits(hs); // Commit to recalibration table
		first_ = false;
		numReported_ += hssz;
		mainunlock();
	}
}

/**
 * Report a maxed-out read.  Typically we do nothing, but we might
 * want to print a placeholder when output is chained.
 */
void ChainingHitSink::reportMaxed(const vector<Hit>& hs, PatternSourcePerThread& p) {
	assert(!hs.empty());
	if(strata_) {
		// Get the stratum
		int8_t loStrat = hs.front().stratum;
#ifndef NDEBUG
		// Strata should be uniform across hits
		vector<Hit>::const_iterator it;
		for(it = hs.begin(); it != hs.end(); it++) {
			assert_eq(loStrat, it->stratum)
		}
#endif
		if(loStrat > 0) {
			// Critical section for output stream 0
			HitSet s;
			p.bufa().toHitSet(s); // grab read details from ReadBuf
			lock(0);
			s.serialize(out(0));
			unlock(0);
		} else {
			// We eliminated all possible strata, so this read is done.
			// Don't serialize anything.
			assert_eq(0, loStrat);
		}
	} else {
		// Maxed out in unstratified mode - this read is done
	}
}

/**
 * Report an unaligned read.  Typically we do nothing, but we might
 * want to print a placeholder when output is chained.
 */
void ChainingHitSink::reportUnaligned(PatternSourcePerThread& p) {
	// Read is unaligned; just report a huge starting stratum
	HitSet s;
	p.bufa().toHitSet(s);
	lock(0);
	s.serialize(out(0));
	unlock(0);
}
