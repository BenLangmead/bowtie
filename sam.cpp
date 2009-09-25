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
		flags |= SAM_FLAG_PAIRED | SAM_FLAG_SECOND_IN_PAIR;
	} else if(h.mate == 2) {
		flags |= SAM_FLAG_PAIRED | SAM_FLAG_FIRST_IN_PAIR;
	}
	if(!h.fw) flags |= SAM_FLAG_QUERY_STRAND;
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
	ss << "\t=";
	// MPOS
	ss << '\t' << (h.h.second + 1);
	// ISIZE
	ss << "\t0";
	// SEQ
	ss << '\t' << h.patSeq;
	// QUAL
	ss << '\t' << h.quals;
	//
	// Optional fields
	//
	// Always output stratum
	ss << "\tXS:i:" << (int)h.stratum;
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
void SAMHitSink::reportUnOrMax(PatternSourcePerThread& p, bool un) {
	HitSink::reportUnaligned(p);
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
