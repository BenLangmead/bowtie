/*
 * annot.h
 *
 *  Created on: Aug 3, 2009
 *      Author: Ben Langmead
 */

#ifndef ANNOT_H_
#define ANNOT_H_

#include <stdint.h>
#include <map>
#include <iostream>
#include <fstream>

/**
 * Encapsulates a sorted list of reference positions that are annotated
 * somehow (e.g. as a SNP).
 */
class AnnotationMap {
	typedef std::pair<uint32_t, uint32_t> U32Pair;
	typedef std::map<U32Pair, char> AnnotMap;

public:
	AnnotationMap(const char *fname) {
		fname_ = fname;
		parse();
	}

	/**
	 * Give a reference coordinate in the index, translate it into a
	 * new reference coordinate via the reference map supplied by the
	 * user.
	 */
	AnnotMap::const_iterator lower_bound(U32Pair& h) const {
		return map_.lower_bound(h);
	}

protected:

	/**
	 * Parse an annotation-map file.
	 */
	void parse();

	/// filename of file containing the annotation map
	const char *fname_;
	/// maps reference positions to character annotations
	std::map<U32Pair, char> map_;
};

#endif /* ANNOT_H_ */
