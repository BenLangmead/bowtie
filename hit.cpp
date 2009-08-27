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
		numAligned_++;
		mainunlock();
	}
}

/**
 * Report a maxed-out read.  Typically we do nothing, but we might
 * want to print a placeholder when output is chained.
 */
void ChainingHitSink::reportMaxed(const vector<Hit>& hs, PatternSourcePerThread& p) {
	HitSink::reportMaxed(hs, p);
	assert(!hs.empty());
	int8_t loStrat = (strata_ ? hs.front().stratum : 0);
	HitSet s;
	p.bufa().toHitSet(s);
	s.maxedStratum = loStrat;
	lock(0);
	s.serialize(out(0));
	unlock(0);
}

/**
 * Report an unaligned read.  Typically we do nothing, but we might
 * want to print a placeholder when output is chained.
 */
void ChainingHitSink::reportUnaligned(PatternSourcePerThread& p) {
	HitSink::reportUnaligned(p);
	// Read is unaligned; just report a huge starting stratum
	HitSet s;
	p.bufa().toHitSet(s);
	lock(0);
	s.serialize(out(0));
	unlock(0);
}
