/*
 * random_test.cpp
 *
 *  Created on: Dec 11, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include <string.h>
#include "random_source.h"

using namespace std;

int main(void) {
	RandomSource rand;
	rand.init(0);
	uint32_t ts[32];
	memset(ts, 0, 32*sizeof(uint32_t));
	uint32_t r = rand.nextU32();
	cout << "Without reseeding:" << endl;
	for(int i = 0; i < 10000; i++) {
		uint32_t nr = rand.nextU32();
		for(int j = 0; j < 32; j++) {
			if(((r >> j) & 1) != ((nr >> j) & 1)) {
				ts[j]++;
			}
		}
	}
	for(int j = 0; j < 32; j++) {
		cout << ts[j] << endl;
	}
	memset(ts, 0, 32*sizeof(uint32_t));
	rand.init(0);
	r = rand.nextU32();
	cout << "With reseeding:" << endl;
	for(int i = 0; i < 10000; i++) {
		rand.init(i+1);
		uint32_t nr = rand.nextU32();
		for(int j = 0; j < 32; j++) {
			if(((r >> j) & 1) != ((nr >> j) & 1)) {
				ts[j]++;
			}
		}
	}
	for(int j = 0; j < 32; j++) {
		cout << ts[j] << endl;
	}
}
