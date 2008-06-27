#include <stdint.h>
#include <iostream>

using namespace std;

/**
 * Tricky, bit-bashing population count function.
 * 
 * I find it difficult to grok this population-count function, so I'm
 * testing it here.
 */
inline static int pop(uint64_t x) {
   x = x - ((x >> 1) & 0x5555555555555555llu);
   x = (x & 0x3333333333333333llu) + ((x >> 2) & 0x3333333333333333llu);
   x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Fllu;
   x = x + (x >> 8);
   x = x + (x >> 16);
   x = x + (x >> 32);
   return x & 0x3F;
} 

/// Silly population count
inline static int popNaive(uint64_t x) {
	int count = 0;
	for(int i = 0; i < 64; i++) {
		if(((x >> i) & 1) == 1) count++;
	}
	return count;
}

/**
 * For each value in the range delineated by the two command-line
 * arguments, calculate the population count using the bit-bashing
 * and naive approaches and make sure they're equal.  Scream if they're
 * not.
 */
int main(int argc, char **argv) {
	uint64_t iLo = atol(argv[1]);
	uint64_t iHi = atol(argv[2]);
	//cout << "i = (" << iLo << "," << iHi << ")" << endl;
	for(int i = iLo; i <= iHi; i++) {
		int p1 = pop(i);
		//cout << "pop(" << i << ") = " << pop(i) << endl;
		if(p1 != popNaive(i)) {
			cerr << "pop(" << i << ") (" << p1 << ") != popNaive(" << i << ") (" << popNaive(i) << ")" << endl;
			return 1;
		}
	}
	return 0;
}
