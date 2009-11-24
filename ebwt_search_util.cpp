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
