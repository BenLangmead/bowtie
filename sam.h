/*
 * sam.h
 *
 *  Created on: Sep 23, 2009
 *      Author: Ben Langmead
 */

#ifndef SAM_H_
#define SAM_H_

#include "hit.h"

class ReferenceMap;
class AnnotationMap;
class PatternSourcePerThread;

enum {
	SAM_FLAG_PAIRED = 1,
	SAM_FLAG_MAPPED_PAIRED = 2,
	SAM_FLAG_UNMAPPED = 4,
	SAM_FLAG_MATE_UNMAPPED = 8,
	SAM_FLAG_QUERY_STRAND = 16,
	SAM_FLAG_MATE_STRAND = 32,
	SAM_FLAG_FIRST_IN_PAIR = 64,
	SAM_FLAG_SECOND_IN_PAIR = 128,
	SAM_FLAG_NOT_PRIMARY = 256,
	SAM_FLAG_FAILS_CHECKS = 512,
	SAM_FLAG_DUPLICATE = 1024
};

/**
 * Sink that prints lines in SAM format:
 */
class SAMHitSink : public HitSink {
public:
	/**
	 * Construct a single-stream VerboseHitSink (default)
	 */
	SAMHitSink(OutFileBuf* out,
	           int offBase,
	           ReferenceMap *rmap,
	           AnnotationMap *amap,
	           bool fullRef,
	           DECL_HIT_DUMPS2) :
	HitSink(out, PASS_HIT_DUMPS2),
	offBase_(offBase), rmap_(rmap), amap_(amap), fullRef_(fullRef) { }

	/**
	 * Construct a multi-stream VerboseHitSink with one stream per
	 * reference string (see --refout)
	 */
	SAMHitSink(size_t numOuts,
	           int offBase,
	           ReferenceMap *rmap,
	           AnnotationMap *amap,
	           bool fullRef,
	           DECL_HIT_DUMPS2) :
	HitSink(numOuts, PASS_HIT_DUMPS2),
	offBase_(offBase), rmap_(rmap), amap_(amap), fullRef_(fullRef) { }

	/**
	 * Append a SAM alignment to the given output stream.
	 */
	static void append(ostream& ss,
	                   const Hit& h,
	                   const vector<string>* refnames,
	                   ReferenceMap *rmap,
	                   AnnotationMap *amap,
	                   bool fullRef,
	                   int offBase);

	/**
	 * Append a SAM alignment for an aligned read to the given output
	 * stream.
	 */
	static void appendAligned(ostream& ss,
	                          const Hit& h,
	                          const vector<string>* refnames,
	                          ReferenceMap *rmap,
	                          AnnotationMap *amap,
	                          bool fullRef,
	                          int offBase);

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(ostream& ss, const Hit& h) {
		SAMHitSink::append(ss, h, _refnames, rmap_, amap_, fullRef_, offBase_);
	}

	/**
	 * Write the SAM header lines.
	 */
	void appendHeaders(OutFileBuf& os,
	                   size_t numRefs,
	                   const vector<string>& refnames,
	                   bool nosq,
	                   ReferenceMap *rmap,
	                   const uint32_t* plen,
	                   bool fullRef,
	                   const char *cmdline);

protected:

	/**
	 *
	 */
	void reportUnOrMax(PatternSourcePerThread& p,
	                   const vector<Hit>* hs,
	                   bool un);

	/**
	 * Report a verbose, human-readable alignment to the appropriate
	 * output stream.
	 */
	virtual void reportHit(const Hit& h);

	/**
	 * See sam.cpp
	 */
	virtual void reportMaxed(const vector<Hit>& hs, PatternSourcePerThread& p) {
		reportUnOrMax(p, &hs, false);
	}

	/**
	 * See sam.cpp
	 */
	virtual void reportUnaligned(PatternSourcePerThread& p) {
		reportUnOrMax(p, NULL, true);
	}

private:
	int      _partition;   /// partition size, or 0 if partitioning is disabled
	int      offBase_;     /// Add this to reference offsets before outputting.
	                       /// (An easy way to make things 1-based instead of
	                       /// 0-based)
    ReferenceMap *rmap_;   /// mapping to reference coordinate system.
	AnnotationMap *amap_;  ///
	bool fullRef_;         /// print full reference name, not just up to whitespace
};

#endif /* SAM_H_ */
