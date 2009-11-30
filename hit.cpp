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

/**
 * Report a maxed-out read.
 */
void VerboseHitSink::reportMaxed(const vector<Hit>& hs, PatternSourcePerThread& p) {
	HitSink::reportMaxed(hs, p);
	if(sampleMax_) {
		rand_.init(p.bufa().seed);
		assert_gt(hs.size(), 0);
		size_t num = 1;
		for(size_t i = 1; i < hs.size(); i++) {
			assert_geq(hs[i].stratum, hs[i-1].stratum);
			if(hs[i].stratum == hs[i-1].stratum) num++;
			else break;
		}
		assert_leq(num, hs.size());
		uint32_t ch = rand_.nextU32() % num;
		//hs[ch].oms = hs.size()-1;
		reportHit(hs[ch]);
	}
}
