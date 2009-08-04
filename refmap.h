/*
 * refmap.h
 *
 *  Created on: Aug 3, 2009
 *      Author: Ben Langmead
 */

#ifndef REFMAP_H_
#define REFMAP_H_

#include <stdint.h>
#include <vector>
#include <iostream>
#include <fstream>

class ReferenceMap {
	typedef std::pair<uint32_t, uint32_t> U32Pair;

public:
	ReferenceMap(const char *fname) {
		fname_ = fname;
		parse();
	}

	/**
	 * Give a reference coordinate in the index, translate it into a
	 * new reference coordinate via the reference map supplied by the
	 * user.
	 */
	void map(U32Pair& h) const;

protected:

	/**
	 * Parse a reference-map file.
	 */
	void parse();

	const char *fname_;
	std::vector<U32Pair> map_;
};

#endif /* REFMAP_H_ */
