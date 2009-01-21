/**
 * range.h
 */

#ifndef RANGE_H_
#define RANGE_H_

#include <stdint.h>

/**
 * A range along with the alignment it represents.
 */
struct Range {
	uint32_t top;     // top of range
	uint32_t bot;     // bottom of range
	uint32_t stratum; // stratum
	uint32_t numMms;  // # mismatches
	uint32_t *mms;    // list of positions with mismatches
	char     *refcs;  // reference characters at mismatch positions
};

#endif /* RANGE_H_ */
