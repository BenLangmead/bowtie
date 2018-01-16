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
#include "btypes.h"

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
	SAMHitSink(
		const std::string& ofn,
		size_t output_buffer_size,
		int offBase,
		ReferenceMap *rmap,
		AnnotationMap *amap,
		bool fullRef,
		bool noQnameTrunc,
		const std::string& dumpAl,
		const std::string& dumpUnal,
		const std::string& dumpMax,
		bool onePairFile,
		bool sampleMax,
		std::vector<std::string>* refnames,
		size_t nthreads,
		int perThreadBufSize,
		int nOutput) :
		HitSink(
			ofn,
			output_buffer_size,
			dumpAl,
			dumpUnal,
			dumpMax,
			onePairFile,
			sampleMax,
			refnames,
			nthreads,
			perThreadBufSize,
			nOutput),
		offBase_(offBase),
		rmap_(rmap),
		amap_(amap),
		fullRef_(fullRef),
		noQnameTrunc_(noQnameTrunc) { }

	/**
	 * Construct a multi-stream VerboseHitSink with one stream per
	 * reference string (see --refout)
	 */
	SAMHitSink(
		size_t numOuts,
		size_t output_buffer_size,
		int offBase,
		ReferenceMap *rmap,
		AnnotationMap *amap,
		bool fullRef,
		const std::string& dumpAl,
		const std::string& dumpUnal,
		const std::string& dumpMax,
		bool onePairFile,
		bool sampleMax,
		std::vector<std::string>* refnames,
		size_t nthreads,
		int perThreadBufSize,
		int nOutput) :
		HitSink(
			numOuts,
			output_buffer_size,
			dumpAl,
			dumpUnal,
			dumpMax,
			onePairFile,
			sampleMax,
			refnames,
			nthreads,
			perThreadBufSize,
			nOutput),
		offBase_(offBase),
		rmap_(rmap),
		amap_(amap),
		fullRef_(fullRef) { }

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	virtual void append(BTString& o, const Hit& h, int mapq, int xms);

	/**
	 * Write the SAM header lines.
	 */
	void appendHeaders(
		FILE *ofh,
		size_t numRefs,
		const vector<string>& refnames,
		bool color,
		bool nosq,
		ReferenceMap *rmap,
		const TIndexOffU* plen,
		bool fullRef,
		bool noQnameTrunc,
		const char *cmdline,
		const char *rgline);

protected:

	/**
	 * Both
	 */
	void reportUnOrMax(
		PatternSourcePerThread& p,
		vector<Hit>* hs, // could be NULL
		size_t threadId,
		bool un);

	/**
	 * See sam.cpp
	 */
	virtual void reportMaxed(
		vector<Hit>& hs,
		size_t threadId,
		PatternSourcePerThread& p);

	/**
	 * See sam.cpp
	 */
	virtual void reportUnaligned(
		size_t threadId,
		PatternSourcePerThread& p)
	{
		reportUnOrMax(p, NULL, threadId, true);
	}

private:
	int  offBase_;        /// Add this to reference offsets before outputting.
	                      /// (An easy way to make things 1-based instead of
	                      /// 0-based)
	ReferenceMap *rmap_;  /// mapping to reference coordinate system.
	AnnotationMap *amap_; ///
	bool fullRef_;        /// print full reference name, not just up to whitespace
	bool noQnameTrunc_;   /// true -> don't truncate QNAME at first whitespace
};

#endif /* SAM_H_ */
