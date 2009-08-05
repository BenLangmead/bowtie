/*
 * hit_set.cpp
 *
 *  Created on: Jul 31, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include <vector>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "hit_set.h"

using namespace std;
using namespace seqan;

/**
 * Report up to 'khits' hits from this HitSet.
 */
void HitSet::reportUpTo(ostream& os, int khits) {
	khits = min(khits, (int)size());
	String<Dna5> seqrc;
	String<char> qualr;
	for(int i = 0; i < khits; i++) {
		const HitSetEnt& h = ents[i];
		if(!h.fw && seqan::empty(seqrc)) {
			// Lazily initialize seqrc and qualr
			seqrc = seq;
			reverseComplementInPlace(seqrc, false);
			assert_eq(seqan::length(seqrc), seqan::length(seq));
			qualr = qual;
			reverseInPlace(qualr);
			assert_eq(seqan::length(qualr), seqan::length(qual));
		}
		os << name << '\t'
		   << (h.fw ? '+' : '-') << '\t'
		   << h.h.first << '\t'
		   << h.h.second << '\t'
		   << (h.fw ? seq : seqrc) << '\t'
		   << (h.fw ? qual : qualr) << '\t'
		   << h.oms << '\t';
		for(size_t i = 0; i < h.edits.size(); i++) {
			const Edit& e = h.edits[i];
			os << e.pos;
			if(e.type == EDIT_TYPE_SNP) os << "S";
			os << ":" << "ACGT"[e.chr] << ">" << seq[e.pos];
			if(i < h.edits.size()-1) os << ",";
		}
		os << endl;
	}
}

ostream& operator << (ostream& os, const HitSetEnt& hs) {
	os << "\t" << hs.h.first << ":" << hs.h.second;
	return os;
}

ostream& operator << (ostream& os, const HitSet& hs) {
	os << hs.name << ":" << hs.seq << ":" << hs.qual << endl;
	vector<HitSetEnt>::const_iterator it;
	for(it = hs.ents.begin(); it != hs.ents.end(); it++) {
		os << (*it);
	}
	return os;
}
