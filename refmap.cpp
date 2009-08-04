/*
 * refmap.cpp
 *
 *  Created on: Aug 3, 2009
 *      Author: Ben Langmead
 */

#include "refmap.h"

using namespace std;

/**
 * Given a refid,offset pair in the index space, transform it into the
 * reference coordinate space according to the reference mappings
 * provided by the user.
 */
void ReferenceMap::map(U32Pair& h) const {
	if(h.first >= map_.size()) {
		cerr << "Could not find a reference-map entry for reference "
				  << h.first << " in map file \"" << fname_ << "\""
				  << endl;
		exit(1);
	}
	h.second += map_[h.first].second;
	h.first = map_[h.first].first;
}

/**
 * Parse a reference-map file.
 */
void ReferenceMap::parse() {
	ifstream in(fname_);
	if(!in.good() && in.is_open()) {
		cerr << "Could not open reference map file " << fname_ << endl;
		exit(1);
	}
	while(!in.eof()) {
		map_.resize(map_.size()+1);
		uint32_t id, off;
		in >> id >> off;
		map_.back().first = id;
		map_.back().second = off;
	}
	in.close();
}
