/**
 * range.h
 */

#ifndef RANGE_H_
#define RANGE_H_

#include <vector>
#include <stdint.h>

/**
 * A range along with the alignment it represents.
 */
struct Range {
	Range() : top(0), bot(0), stratum(0), numMms(0), ebwt(NULL) {
		mms.clear();
		refcs.clear();
	}
	uint32_t top;     // top of range
	uint32_t bot;     // bottom of range
	uint32_t stratum; // stratum
	uint32_t numMms;  // # mismatches
	bool fw;          // the forward orientation of read aligned?
	std::vector<uint32_t> mms;   // list of positions with mismatches
	std::vector<uint8_t>  refcs; // reference characters at mismatch positions
	const Ebwt<seqan::String<seqan::Dna> > *ebwt;
};

#endif /* RANGE_H_ */
