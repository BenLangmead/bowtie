/*
 * sam.h
 *
 *  Created on: Sep 23, 2009
 *      Author: Ben Langmead
 */

#ifndef SAM_H_
#define SAM_H_

#include "refmap.h"
#include "annot.h"
#include "pat.h"
#include "random_source.h"

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
	           bool sampleMax,
	           int defaultMapq,
	           DECL_HIT_DUMPS2) :
	HitSink(out, PASS_HIT_DUMPS2),
	offBase_(offBase), sampleMax_(sampleMax), defaultMapq_(defaultMapq),
	rmap_(rmap), amap_(amap), fullRef_(fullRef) { }

	/**
	 * Construct a multi-stream VerboseHitSink with one stream per
	 * reference string (see --refout)
	 */
	SAMHitSink(size_t numOuts,
	           int offBase,
	           ReferenceMap *rmap,
	           AnnotationMap *amap,
	           bool fullRef,
	           bool sampleMax,
	           int defaultMapq,
	           DECL_HIT_DUMPS2) :
	HitSink(numOuts, PASS_HIT_DUMPS2),
	offBase_(offBase),  sampleMax_(sampleMax), defaultMapq_(defaultMapq),
	rmap_(rmap), amap_(amap), fullRef_(fullRef) { }

	/**
	 * Append a SAM alignment to the given output stream.
	 */
	static void append(ostream& ss,
	                   const Hit& h,
	                   int mapq,
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
	                          int mapq,
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
		SAMHitSink::append(ss, h, defaultMapq_, _refnames, rmap_, amap_, fullRef_, offBase_);
	}

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(ostream& ss, const Hit& h, int mapq) {
		SAMHitSink::append(ss, h, mapq, _refnames, rmap_, amap_, fullRef_, offBase_);
	}

	/**
	 * Write the SAM header lines.
	 */
	void appendHeaders(OutFileBuf& os,
	                   size_t numRefs,
	                   const vector<string>& refnames,
	                   bool color,
	                   bool nosq,
	                   ReferenceMap *rmap,
	                   const uint32_t* plen,
	                   bool fullRef,
	                   const char *cmdline,
	                   const char *rgline);

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
	virtual void reportHit(const Hit& h) {
		reportHit(h, defaultMapq_);
	}

	/**
	 * Report a verbose, human-readable alignment to the appropriate
	 * output stream with the given mapping quality.
	 */
	virtual void reportHit(const Hit& h, int mapq);

	/**
	 * See sam.cpp
	 */
	virtual void reportMaxed(const vector<Hit>& hs, PatternSourcePerThread& p) {
		if(sampleMax_) {
			HitSink::reportMaxed(hs, p);
			rand_.init(p.bufa().seed);
			assert_gt(hs.size(), 0);
			size_t num = 1;
			for(size_t i = 1; i < hs.size(); i++) {
				assert_geq(hs[i].stratum, hs[i-1].stratum);
				if(hs[i].stratum == hs[i-1].stratum) num++;
				else break;
			}
			assert_leq(num, hs.size());
			reportHit(hs[rand_.nextU32() % num], 0 /* MAPQ=0 */);
		} else {
			reportUnOrMax(p, &hs, false);
		}
	}

	/**
	 * See sam.cpp
	 */
	virtual void reportUnaligned(PatternSourcePerThread& p) {
		reportUnOrMax(p, NULL, true);
	}

private:
	int  offBase_;        /// Add this to reference offsets before outputting.
	                      /// (An easy way to make things 1-based instead of
	                      /// 0-based)
	bool sampleMax_;      /// When reporting a maxed-out read, randomly report
	                      /// 1 of the alignments with mapping quality = 0
	int  defaultMapq_;    /// Default mapping quality to report when one is
	                      /// not specified
	ReferenceMap *rmap_;  /// mapping to reference coordinate system.
	AnnotationMap *amap_; ///
	bool fullRef_;        /// print full reference name, not just up to whitespace
	RandomSource rand_;   /// for pseudo-randoms
};

#endif /* SAM_H_ */
