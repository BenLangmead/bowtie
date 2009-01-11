#include "ebwt_search_util.h"
#include "seqan/file.h"

using namespace std;
using namespace seqan;

/**
 * Print a hit along with information about the backtracking
 * regions constraining the hit.
 */
void printHit(const vector<String<Dna5> >& os,
			  const Hit& h,
			  const String<Dna5>& qry,
			  size_t qlen,
			  uint32_t unrevOff,
			  uint32_t oneRevOff,
			  uint32_t twoRevOff,
			  uint32_t threeRevOff,
			  bool ebwtFw)
{
	// Print pattern sequence
	cout << "  Pat:  " << qry << endl;
	// Print text sequence
	cout << "  Tseg: ";
	if(ebwtFw) {
		for(size_t i = 0; i < qlen; i++) {
			cout << os[h.h.first][h.h.second + i];
		}
	} else {
		for(int i = (int)qlen-1; i >= 0; i--) {
			cout << os[h.h.first][h.h.second + i];
		}
	}
	cout << endl;
	cout << "  Bt:   ";
	for(int i = (int)qlen-1; i >= 0; i--) {
		if     (i < (int)unrevOff)    cout << "0";
		else if(i < (int)oneRevOff)   cout << "1";
		else if(i < (int)twoRevOff)   cout << "2";
		else if(i < (int)threeRevOff) cout << "3";
		else cout << "X";
	}
	cout << endl;
}

/**
 * Naively search for the same hits that should be found by the
 * backtracking search.  This is used only if --orig and -s have
 * been specified for bowtie-debug.
 */
void naiveOracle(const vector<String<Dna5> >& os,
				 const String<Dna5>& qry,
				 uint32_t qlen,
				 const String<char>& qual,
				 const String<char>& name,
				 uint32_t patid,
				 vector<Hit>& hits,
				 vector<int>& strata,
				 uint32_t qualThresh,
				 uint32_t unrevOff,
				 uint32_t oneRevOff,
				 uint32_t twoRevOff,
				 uint32_t threeRevOff,
				 bool fw,
				 bool ebwtFw,
				 uint32_t iham = 0,
				 String<QueryMutation>* muts = NULL,
				 bool maqPenalty = true,
				 bool halfAndHalf = false,
				 bool reportExacts = true,
				 bool invert = false)
{
	bool fivePrimeOnLeft = (ebwtFw == fw);
	uint32_t plen = qlen;
	uint8_t *pstr = (uint8_t *)begin(qry, Standard());
	// For each text...
	for(size_t i = 0; i < os.size(); i++) {
		// For each text position...
		if(length(os[i]) < plen) continue;
		uint32_t olen = length(os[i]);
		uint8_t *ostr = (uint8_t *)begin(os[i], Standard());
		// For each possible alignment of pattern against text
		for(size_t j = 0; j <= olen - plen; j++) {
			size_t mms = 0; // mismatches observed over the whole thing
			size_t rev1mm  = 0; // mismatches observed in the 1-revisitable region
			size_t rev2mm  = 0; // mismatches observed in the 2-revisitable region
			size_t rev3mm  = 0; // mismatches observed in the 3-revisitable region
			uint32_t ham = iham; // weighted hamming distance so far
			FixedBitset<max_read_bp> diffs; // mismatch bitvector
			vector<char> refcs; // reference characters for mms
			refcs.resize(qlen, 0);
			// For each alignment column, from right to left
			bool success = true;
			int ok, okInc;
			if(ebwtFw) {
				ok = j+(int)plen-1;
				okInc = -1;
			} else {
				ok = olen-(j+((int)plen-1))-1;
				okInc = 1;
			}
			bool rejectN = false;
			for(int k = (int)plen-1; k >= 0; k--) {
				size_t kr = plen-1-k;
				if((int)ostr[ok] == 4) {
					rejectN = true;
					break;
				}
				if(pstr[k] != ostr[ok]) {
					mms++;
					ham += mmPenalty(maqPenalty, phredCharToPhredQual(qual[k]));
					if(ham > qualThresh) {
						// Alignment is invalid because it exceeds
						// our target weighted hamming distance
						// threshold
						success = false;
						break;
					}
					size_t koff = kr;
					if(invert) koff = (size_t)k;
					// What region does the mm fall into?
					if(koff < unrevOff) {
						// Alignment is invalid because it contains
						// a mismatch in the unrevisitable region
						success = false;
						break;
					} else if(koff < oneRevOff) {
						rev1mm++;
						if(rev1mm > 1) {
							// Alignment is invalid because it
							// contains more than 1 mismatch in the
							// 1-revisitable region
							success = false;
							break;
						}
					} else if(koff < twoRevOff) {
						rev2mm++;
						if(rev2mm > 2) {
							// Alignment is invalid because it
							// contains more than 2 mismatches in the
							// 2-revisitable region
							success = false;
							break;
						}
					} else if(koff < threeRevOff) {
						rev3mm++;
						if(rev3mm > 3) {
							// Alignment is invalid because it
							// contains more than 3 mismatches in the
							// 3-revisitable region
							success = false;
							break;
						}
					}
					if(halfAndHalf) {
						if(twoRevOff == threeRevOff) {
							// 1 and 1
							assert_eq(0, rev3mm);
							if(rev1mm > 1 || rev2mm > 1) {
								// Half-and-half alignment is invalid
								// because it contains more than 1 mismatch
								// in either one or the other half
								success = false;
								break;
							}
						} else {
							// 1 and 1,2
							assert_eq(unrevOff, oneRevOff);
							assert_eq(0, rev1mm);
							if(rev2mm > 2 || rev3mm > 2) {
								success = false;
								break;
							}
						}
					}
					// Update 'diffs' and 'refcs' to reflect this
					// mismatch
					if(fivePrimeOnLeft) {
						diffs.set(k);
						refcs[k] = (char)ostr[ok];
					} else {
						// The 3' end is on on the left end of the
						// pattern, but the diffs vector should
						// encode mismatches w/r/t the 5' end, so
						// we flip
						diffs.set(plen-k-1);
						refcs[plen-k-1] = (char)ostr[ok];
					}
				}
				ok += okInc;
			}
			if(rejectN) {
				// Rejected because the reference half of the
				// alignment contained one or more Ns
				continue;
			}
			if(halfAndHalf) {
				if(twoRevOff == threeRevOff) {
					if(rev1mm != 1 || rev2mm != 1) {
						success = false;
					}
				} else {
					if(rev2mm == 0 || rev3mm == 0 ||
					   rev2mm + rev3mm < 2 ||
					   rev2mm + rev3mm > 3)
					{
						success = false;
					}
				}
			}
			if(!reportExacts && mms == 0) {
				// Reject this exact alignment because we've been
				// instructed to ignore them (usually because a
				// previous invocation already reported them)
				success = false;
			}
			if(success) {
				// It's a hit
				uint32_t off = j;
				int stratum = (int)(rev1mm + rev2mm + rev3mm);
				if(!ebwtFw) {
					off = olen - off;
					off -= plen;
				}
				// Add in mismatches from _muts
				if(muts != NULL) {
					for(size_t i = 0; i < length(*muts); i++) {
						// Entries in _mms[] are in terms of offset into
						// _qry - not in terms of offset from 3' or 5' end
						if(fivePrimeOnLeft) {
							diffs.set((*muts)[i].pos);
							refcs[(*muts)[i].pos] = (*muts)[i].newBase;
						} else {
							diffs.set(plen - (*muts)[i].pos - 1);
							refcs[plen - (*muts)[i].pos - 1] = (*muts)[i].newBase;
						}
						stratum++;
					}
				}
				Hit h(make_pair(i, off),
					  patid,  // read id
					  name,   // read name
					  qry,    // read sequence
					  qual,   // read qualities
					  fw,     // forward/reverse-comp
					  diffs,  // mismatch bitvector
					  refcs);
				hits.push_back(h);
				strata.push_back(stratum);
			} // For each pattern character
		} // For each alignment over current text
	} // For each text
}
