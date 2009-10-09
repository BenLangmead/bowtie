/*
 * sam.cpp
 *
 *  Created on: Sep 23, 2009
 *      Author: Ben Langmead
 */

#include <vector>
#include <string>
#include <iostream>
#include "pat.h"
#include "hit.h"
#include "sam.h"

using namespace std;

/**
 * Write the SAM header lines.
 */
void SAMHitSink::appendHeaders(OutFileBuf& os,
                               size_t numRefs,
                               const vector<string>& refnames,
                               bool nosq,
                               ReferenceMap *rmap,
                               const uint32_t* plen,
                               bool fullRef,
                               const char *cmdline)
{
	ostringstream ss;
	ss << "@HD\tVN:1.0\tSO:unsorted" << endl;
	if(!nosq) {
		for(size_t i = 0; i < numRefs; i++) {
			// RNAME
			ss << "@SQ\tSN:";
			if(!refnames.empty() && rmap != NULL) {
				printUptoWs(ss, rmap->getName(i), !fullRef);
			} else if(i < refnames.size()) {
				printUptoWs(ss, refnames[i], !fullRef);
			} else {
				ss << i;
			}
			ss << "\tLN:" << plen[i] << endl;
		}
	}
	ss << "@PG\tID=Bowtie\tVN=" << BOWTIE_VERSION << "\tCL=\"" << cmdline << "\"" << endl;
	os.writeString(ss.str());
}

/**
 * Append a SAM output record for an unaligned read.
 */
void SAMHitSink::appendAligned(ostream& ss,
                               const Hit& h,
                               const vector<string>* refnames,
                               ReferenceMap *rmap,
                               AnnotationMap *amap,
                               bool fullRef,
                               int offBase)
{
	// QNAME
	ss << h.patName << "\t";
	// FLAG
	int flags = 0;
	if(h.mate == 1) {
		flags |= SAM_FLAG_PAIRED | SAM_FLAG_FIRST_IN_PAIR | SAM_FLAG_MAPPED_PAIRED;
	} else if(h.mate == 2) {
		flags |= SAM_FLAG_PAIRED | SAM_FLAG_SECOND_IN_PAIR | SAM_FLAG_MAPPED_PAIRED;
	}
	if(!h.fw) flags |= SAM_FLAG_QUERY_STRAND;
	if(h.mate > 0 && !h.mfw) flags |= SAM_FLAG_MATE_STRAND;
	ss << flags << "\t";
	// RNAME
	if(refnames != NULL && rmap != NULL) {
		printUptoWs(ss, rmap->getName(h.h.first), !fullRef);
	} else if(refnames != NULL && h.h.first < refnames->size()) {
		printUptoWs(ss, (*refnames)[h.h.first], !fullRef);
	} else {
		ss << h.h.first;
	}
	// POS
	ss << '\t' << (h.h.second + 1);
	// MAPQ
	ss << "\t255";
	// CIGAR
	ss << '\t' << h.length() << 'M';
	// MRNM
	if(h.mate > 0) {
		ss << "\t=";
	} else {
		ss << "\t*";
	}
	// MPOS
	if(h.mate > 0) {
		ss << '\t' << (h.mh.second + 1);
	} else {
		ss << "\t0";
	}
	// ISIZE
	ss << '\t';
	if(h.mate > 0) {
		assert_eq(h.h.first, h.mh.first);
		int64_t inslen = 0;
		if(h.h.second > h.mh.second) {
			inslen = (int64_t)h.h.second - (int64_t)h.mh.second + (int64_t)h.length();
			inslen = -inslen;
		} else {
			inslen = (int64_t)h.mh.second - (int64_t)h.h.second + (int64_t)h.mlen;
		}
		ss << inslen;
	} else {
		ss << '0';
	}
	// SEQ
	ss << '\t' << h.patSeq;
	// QUAL
	ss << '\t' << h.quals;
	//
	// Optional fields
	//
	// Always output stratum
	ss << "\tXA:i:" << (int)h.stratum;
	// Look for SNP annotations falling within the alignment
	// Output MD field
	const size_t len = length(h.patSeq);
	int nm = 0;
	int run = 0;
	ss << "\tMD:Z:";
	if(h.fw) {
		for (int i = 0; i < (int)len; ++ i) {
			if(h.mms.test(i)) {
				nm++;
				// There's a mismatch at this position
				assert_gt((int)h.refcs.size(), i);
				char refChar = toupper(h.refcs[i]);
				ASSERT_ONLY(char qryChar = (h.fw ? h.patSeq[i] : h.patSeq[length(h.patSeq)-i-1]));
				assert_neq(refChar, qryChar);
				ss << run << refChar;
				run = 0;
			} else {
				run++;
			}
		}
	} else {
		for (int i = len-1; i >= 0; -- i) {
			if(h.mms.test(i)) {
				nm++;
				// There's a mismatch at this position
				assert_gt((int)h.refcs.size(), i);
				char refChar = toupper(h.refcs[i]);
				ASSERT_ONLY(char qryChar = (h.fw ? h.patSeq[i] : h.patSeq[length(h.patSeq)-i-1]));
				assert_neq(refChar, qryChar);
				ss << run << refChar;
				run = 0;
			} else {
				run++;
			}
		}
	}
	ss << run;
	// Add optional edit distance field
	ss << "\tNM:i:" << nm;
	ss << endl;
}

/**
 * Report a verbose, human-readable alignment to the appropriate
 * output stream.
 */
void SAMHitSink::reportHit(const Hit& h) {
	HitSink::reportHit(h);
	ostringstream ss;
	append(ss, h);
	// Make sure to grab lock before writing to output stream
	lock(h.h.first);
	out(h.h.first).writeString(ss.str());
	unlock(h.h.first);
}

/**
 *
 */
void SAMHitSink::reportUnOrMax(PatternSourcePerThread& p,
                               const vector<Hit>* hs,
                               bool un)
{
	if(un) HitSink::reportUnaligned(p);
	else   HitSink::reportMaxed(*hs, p);
	ostringstream ss;
	bool paired = !p.bufb().empty();
	ss << p.bufa().name << "\t"
	   << (SAM_FLAG_UNMAPPED | (paired ? (SAM_FLAG_PAIRED | SAM_FLAG_FIRST_IN_PAIR | SAM_FLAG_MATE_UNMAPPED) : 0)) << "\t*"
	   << "\t0\t0\t*\t*\t0\t0\t"
	   << p.bufa().patFw << "\t" << p.bufa().qual << "\tXM:i:"
	   << (un ? '0' : '1') << endl;
	if(paired) {
		ss << p.bufb().name << "\t"
		   << (SAM_FLAG_UNMAPPED | (paired ? (SAM_FLAG_PAIRED | SAM_FLAG_SECOND_IN_PAIR | SAM_FLAG_MATE_UNMAPPED) : 0)) << "\t*"
		   << "\t0\t0\t*\t*\t0\t0\t"
		   << p.bufb().patFw << "\t" << p.bufb().qual << "\tXM:i:"
		   << (un ? '0' : '1') << endl;
	}
	lock(0);
	out(0).writeString(ss.str());
	unlock(0);
}

/**
 * Append a SAM alignment to the given output stream.
 */
void SAMHitSink::append(ostream& ss,
                        const Hit& h,
                        const vector<string>* refnames,
                        ReferenceMap *rmap,
                        AnnotationMap *amap,
                        bool fullRef,
                        int offBase)
{
	appendAligned(ss, h, refnames, rmap, amap, fullRef, offBase);
}
