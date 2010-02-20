#include "hit.h"
#include "hit_set.h"
#include "search_globals.h"

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
void ChainingHitSink::reportMaxed(vector<Hit>& hs, PatternSourcePerThread& p) {
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
void VerboseHitSink::reportMaxed(vector<Hit>& hs, PatternSourcePerThread& p) {
	HitSink::reportMaxed(hs, p);
	if(sampleMax_) {
		RandomSource rand;
		rand.init(p.bufa().seed);
		assert_gt(hs.size(), 0);
		bool paired = hs.front().mate > 0;
		size_t num = 1;
		if(paired) {
			num = 0;
			int bestStratum = 999;
			for(size_t i = 0; i < hs.size()-1; i += 2) {
				int strat = min(hs[i].stratum, hs[i+1].stratum);
				if(strat < bestStratum) {
					bestStratum = strat;
					num = 1;
				} else if(strat == bestStratum) {
					num++;
				}
			}
			assert_leq(num, hs.size());
			uint32_t r = rand.nextU32() % num;
			num = 0;
			for(size_t i = 0; i < hs.size()-1; i += 2) {
				int strat = min(hs[i].stratum, hs[i+1].stratum);
				if(strat == bestStratum) {
					if(num == r) {
						hs[i].oms = hs[i+1].oms = hs.size()/2;
						reportHits(hs, i, i+2);
						break;
					}
					num++;
				}
			}
			assert_eq(num, r);
		} else {
			for(size_t i = 1; i < hs.size(); i++) {
				assert_geq(hs[i].stratum, hs[i-1].stratum);
				if(hs[i].stratum == hs[i-1].stratum) num++;
				else break;
			}
			assert_leq(num, hs.size());
			uint32_t r = rand.nextU32() % num;
			Hit& h = hs[r];
			h.oms = hs.size();
			reportHit(h, false);
		}
	}
}

/**
 * Append a verbose, readable hit to the given output stream.
 */
void VerboseHitSink::append(ostream& ss,
                   const Hit& h,
                   const vector<string>* refnames,
                   ReferenceMap *rmap,
                   AnnotationMap *amap,
                   bool fullRef,
                   int partition,
                   int offBase,
                   bool colorSeq,
                   bool colorQual,
                   bool cost,
                   const Bitset& suppress)
{
	bool spill = false;
	int spillAmt = 0;
	uint32_t pdiv = 0xffffffff;
	uint32_t pmod = 0xffffffff;
	do {
		bool dospill = false;
		if(spill) {
			// The read spilled over a partition boundary and so
			// needs to be printed more than once
			spill = false;
			dospill = true;
			spillAmt++;
		}
		assert(!spill);
		size_t field = 0;
		bool firstfield = true;
		if(partition != 0) {
			int pospart = abs(partition);
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// Output a partitioning key
				// First component of the key is the reference index
				if(refnames != NULL && rmap != NULL) {
					printUptoWs(ss, rmap->getName(h.h.first), !fullRef);
				} else if(refnames != NULL && h.h.first < refnames->size()) {
					printUptoWs(ss, (*refnames)[h.h.first], !fullRef);
				} else {
					ss << h.h.first;
				}
			}
			ostringstream ss2, ss3;
			// Next component of the key is the partition id
			if(!dospill) {
				pdiv = (h.h.second + offBase) / pospart;
				pmod = (h.h.second + offBase) % pospart;
			}
			assert_neq(0xffffffff, pdiv);
			assert_neq(0xffffffff, pmod);
			if(dospill) assert_gt(spillAmt, 0);
			ss2 << (pdiv + (dospill ? spillAmt : 0));
			if(partition > 0 &&
			   (pmod + h.length()) >= ((uint32_t)pospart * (spillAmt + 1))) {
				// Spills into the next partition so we need to
				// output another alignment for that partition
				spill = true;
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// Print partition id with leading 0s so that Hadoop
				// can do lexicographical sort (modern Hadoop versions
				// seen to support numeric)
				string s2 = ss2.str();
				size_t partDigits = 1;
				if(pospart >= 10) partDigits++;
				if(pospart >= 100) partDigits++;
				if(pospart >= 1000) partDigits++;
				if(pospart >= 10000) partDigits++;
				if(pospart >= 100000) partDigits++;
				for(size_t i = s2.length(); i < (10-partDigits); i++) {
					ss << "0";
				}
				ss << s2.c_str();
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// Print offset with leading 0s
				ss3 << (h.h.second + offBase);
				string s3 = ss3.str();
				for(size_t i = s3.length(); i < 9; i++) {
					ss << "0";
				}
				ss << s3;
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.fw? "+":"-");
			}
			// end if(partition != 0)
		} else {
			assert(!dospill);
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.patName;
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.fw? '+' : '-');
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				// .first is text id, .second is offset
				if(refnames != NULL && rmap != NULL) {
					printUptoWs(ss, rmap->getName(h.h.first), !fullRef);
				} else if(refnames != NULL && h.h.first < refnames->size()) {
					printUptoWs(ss, (*refnames)[h.h.first], !fullRef);
				} else {
					ss << h.h.first;
				}
			}
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (h.h.second + offBase);
			}
			// end else clause of if(partition != 0)
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			const String<Dna5>* pat = &h.patSeq;
			if(h.color && colorSeq) pat = &h.colSeq;
			ss << *pat;
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			const String<char>* qual = &h.quals;
			if(h.color && colorQual) qual = &h.colQuals;
			ss << *qual;
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			ss << h.oms;
		}
		if(!suppress.test(field++)) {
			if(firstfield) firstfield = false;
			else ss << '\t';
			// Look for SNP annotations falling within the alignment
			map<int, char> snpAnnots;
			const size_t len = length(h.patSeq);
			if(amap != NULL) {
				AnnotationMap::Iter ai = amap->lower_bound(h.h);
				for(; ai != amap->end(); ai++) {
					assert_geq(ai->first.first, h.h.first);
					if(ai->first.first != h.h.first) {
						// Different chromosome
						break;
					}
					if(ai->first.second >= h.h.second + len) {
						// Doesn't fall into alignment
						break;
					}
					if(ai->second.first != 'S') {
						// Not a SNP annotation
						continue;
					}
					size_t off = ai->first.second - h.h.second;
					if(!h.fw) off = len - off - 1;
					snpAnnots[off] = ai->second.second;
				}
			}
			// Output mismatch column
			bool firstmm = true;
			for (unsigned int i = 0; i < len; ++ i) {
				if(h.mms.test(i)) {
					// There's a mismatch at this position
					if (!firstmm) ss << ",";
					ss << i; // position
					assert_gt(h.refcs.size(), i);
					char refChar = toupper(h.refcs[i]);
					char qryChar = (h.fw ? h.patSeq[i] : h.patSeq[length(h.patSeq)-i-1]);
					assert_neq(refChar, qryChar);
					ss << ":" << refChar << ">" << qryChar;
					firstmm = false;
				} else if(snpAnnots.find(i) != snpAnnots.end()) {
					if (!firstmm) ss << ",";
					ss << i; // position
					char qryChar = (h.fw ? h.patSeq[i] : h.patSeq[length(h.patSeq)-i-1]);
					ss << "S:" << snpAnnots[i] << ">" << qryChar;
					firstmm = false;
				}
			}
			if(partition != 0 && firstmm) ss << '-';
		}
		if(partition != 0) {
			// Fields addded as of Crossbow 0.1.4
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)h.mate;
			}
			// Print label, or whole read name if label isn't found
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				int labelOff = -1;
				// If LB: field is present, print its value
				for(int i = 0; i < (int)seqan::length(h.patName)-3; i++) {
					if(h.patName[i]   == 'L' &&
					   h.patName[i+1] == 'B' &&
					   h.patName[i+2] == ':' &&
					   ((i == 0) || h.patName[i-1] == ';'))
					{
						labelOff = i+3;
						for(int j = labelOff; j < (int)seqan::length(h.patName); j++) {
							if(h.patName[j] != ';') {
								ss << h.patName[j];
							} else {
								break;
							}
						}
					}
				}
				// Otherwise, print the whole read name
				if(labelOff == -1) ss << h.patName;
			}
		}
		if(cost) {
			// Stratum
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)h.stratum;
			}
			// Cost
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << (int)h.cost;
			}
		}
		if(showSeed) {
			// Seed
			if(!suppress.test(field++)) {
				if(firstfield) firstfield = false;
				else ss << '\t';
				ss << h.seed;
			}
		}
		ss << endl;
	} while(spill);
}
