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
 * Report a maxed-out read.
 */
void VerboseHitSink::reportMaxed(
	vector<Hit>& hs,
	size_t threadId,
	PatternSourcePerThread& p)
{
	HitSink::reportMaxed(hs, threadId, p);
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
						hs[i].oms = hs[i+1].oms = (uint32_t)(hs.size()/2);
						reportHits(NULL, &hs, i, i+2, threadId, 0, 0, true, p.rdid());
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
			h.oms = (uint32_t)hs.size();
			reportHits(&h, NULL, 0, 1, threadId, 0, 0, true, p.rdid());
		}
	}
}

/**
 * Append a verbose, readable hit to the given output stream.
 */
void VerboseHitSink::append(
	BTString& o,
	const Hit& h,
	const vector<string>* refnames,
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
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				// Output a partitioning key
				// First component of the key is the reference index
				if(refnames != NULL && h.h.first < refnames->size()) {
					printUptoWs(o, (*refnames)[h.h.first], !fullRef);
				} else {
					o << h.h.first;
				}
			}
			// Next component of the key is the partition id
			if(!dospill) {
				pdiv = (h.h.second + offBase) / pospart;
				pmod = (h.h.second + offBase) % pospart;
			}
			assert_neq(0xffffffff, pdiv);
			assert_neq(0xffffffff, pmod);
			if(dospill) assert_gt(spillAmt, 0);
			if(partition > 0 &&
			   (pmod + h.length()) >= ((uint32_t)pospart * (spillAmt + 1))) {
				// Spills into the next partition so we need to
				// output another alignment for that partition
				spill = true;
			}
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) {
					firstfield = false;
				} else {
					o << '\t';
				}
				// Print partition id with leading 0s so that Hadoop
				// can do lexicographical sort (modern Hadoop versions
				// seen to support numeric)
				int padding = 10;
				uint32_t part = (pdiv + (dospill ? spillAmt : 0));
				uint32_t parttmp = part;
				while(parttmp > 0) {
					padding--;
					parttmp /= 10;
				}
				assert_geq(padding, 0);
				for(int i = 0; i < padding; i++) {
					o << '0';
				}
				o << part;
			}
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) {
					firstfield = false;
				} else {
					o << '\t';
				}
				// Print offset with leading 0s
				int padding = 9;
				uint32_t off = h.h.second + offBase;
				uint32_t offtmp = off;
				while(offtmp > 0) {
					padding--;
					offtmp /= 10;
				}
				assert_geq(padding, 0);
				for(int i = 0; i < padding; i++) {
					o << '0';
				}
				o << off;
			}
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				o << (h.fw? "+":"-");
			}
			// end if(partition != 0)
		} else {
			assert(!dospill);
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				for(size_t i = 0; i < seqan::length(h.patName); i++) {
					o << (char)(h.patName[i]);
				}
			}
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				o << (h.fw? '+' : '-');
			}
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				// .first is text id, .second is offset
				if(refnames != NULL && h.h.first < refnames->size()) {
					printUptoWs(o, (*refnames)[h.h.first], !fullRef);
				} else {
					o << h.h.first;
				}
			}
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				o << (h.h.second + offBase);
			}
			// end else clause of if(partition != 0)
		}
		if(!suppress.test((uint32_t)field++)) {
			if(firstfield) firstfield = false;
			else o << '\t';
			const String<Dna5>* pat = &h.patSeq;
			if(h.color && colorSeq) pat = &h.colSeq;
			for(size_t i = 0; i < seqan::length(*pat); i++) {
				o << (char)((*pat)[i]);
			}
		}
		if(!suppress.test((uint32_t)field++)) {
			if(firstfield) firstfield = false;
			else o << '\t';
			const String<char>* qual = &h.quals;
			if(h.color && colorQual) qual = &h.colQuals;
			for(size_t i = 0; i < seqan::length(*qual); i++) {
				o << (char)((*qual)[i]);
			}
		}
		if(!suppress.test((uint32_t)field++)) {
			if(firstfield) firstfield = false;
			else o << '\t';
			o << h.oms;
		}
		if(!suppress.test((uint32_t)field++)) {
			if(firstfield) firstfield = false;
			else o << '\t';
			const size_t len = length(h.patSeq);
			// Output mismatch column
			bool firstmm = true;
			for (unsigned int i = 0; i < len; ++ i) {
				if(h.mms.test(i)) {
					// There's a mismatch at this position
					if (!firstmm) {
						o << ",";
					}
					o << i; // position
					assert_gt(h.refcs.size(), i);
					char refChar = toupper(h.refcs[i]);
					char qryChar = (h.fw ? h.patSeq[i] : h.patSeq[length(h.patSeq)-i-1]);
					assert_neq(refChar, qryChar);
					o << ":" << refChar << ">" << qryChar;
					firstmm = false;
				}
			}
			if(partition != 0 && firstmm) o << '-';
		}
		if(partition != 0) {
			// Fields addded as of Crossbow 0.1.4
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				o << (int)h.mate;
			}
			// Print label, or whole read name if label isn't found
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
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
								o << h.patName[j];
							} else {
								break;
							}
						}
					}
				}
				// Otherwise, print the whole read name
				if(labelOff == -1) {
					for(size_t i = 0; i < seqan::length(h.patName); i++) {
						o << (char)(h.patName[i]);
					}
				}
			}
		}
		if(cost) {
			// Stratum
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				o << (int)h.stratum;
			}
			// Cost
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				o << (int)h.cost;
			}
		}
		if(showSeed) {
			// Seed
			if(!suppress.test((uint32_t)field++)) {
				if(firstfield) firstfield = false;
				else o << '\t';
				o << h.seed;
			}
		}
		o << '\n';
	} while(spill);
}
