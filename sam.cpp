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
                               bool color,
                               bool nosq,
                               ReferenceMap *rmap,
                               const TIndexOffU* plen,
                               bool fullRef,
                               bool noQnameTrunc,
                               const char *cmdline,
                               const char *rgline)
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
			ss << "\tLN:" << (plen[i] + (color ? 1 : 0)) << endl;
		}
	}
	if(rgline != NULL) {
		ss << "@RG\t" << rgline << endl;
	}
	ss << "@PG\tID:Bowtie\tVN:" << BOWTIE_VERSION << "\tCL:\"" << cmdline << "\"" << endl;
	os.writeString(ss.str());
}

/**
 * Append a SAM output record for an unaligned read.
 */
void SAMHitSink::appendAligned(ostream& ss,
                               const Hit& h,
                               int mapq,
                               int xms, // value for XM:I field
                               const vector<string>* refnames,
                               ReferenceMap *rmap,
                               AnnotationMap *amap,
                               bool fullRef,
                               bool noQnameTrunc,
                               int offBase)
{
	// QNAME
	if(h.mate > 0) {
		// truncate final 2 chars
		for(int i = 0; i < (int)seqan::length(h.patName)-2; i++) {
			if(!noQnameTrunc && isspace((int)h.patName[i])) break;
			ss << h.patName[i];
		}
	} else {
		for(int i = 0; i < (int)seqan::length(h.patName); i++) {
			if(!noQnameTrunc && isspace((int)h.patName[i])) break;
			ss << h.patName[i];
		}
	}
	ss << '\t';
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
	ss << "\t" << mapq;
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
	// Always output cost
	//ss << "\tXC:i:" << (int)h.cost;
	// Look for SNP annotations falling within the alignment
	// Output MD field
	size_t len = length(h.patSeq);
	int nm = 0;
	int run = 0;
	ss << "\tMD:Z:";
	const FixedBitset<1024> *mms = &h.mms;
	ASSERT_ONLY(const String<Dna5>* pat = &h.patSeq);
	const vector<char>* refcs = &h.refcs;
	if(h.color && false) {
		// Disabled: print MD:Z string w/r/t to colors, not letters
		mms = &h.cmms;
		ASSERT_ONLY(pat = &h.colSeq);
		assert_eq(length(h.colSeq), len+1);
		len = length(h.colSeq);
		refcs = &h.crefcs;
	}
	if(h.fw) {
		for (int i = 0; i < (int)len; ++ i) {
			if(mms->test(i)) {
				nm++;
				// There's a mismatch at this position
				assert_gt((int)refcs->size(), i);
				char refChar = toupper((*refcs)[i]);
				ASSERT_ONLY(char qryChar = (h.fw ? (*pat)[i] : (*pat)[len-i-1]));
				assert_neq(refChar, qryChar);
				ss << run << refChar;
				run = 0;
			} else {
				run++;
			}
		}
	} else {
		for (int i = (int)len-1; i >= 0; -- i) {
			if(mms->test(i)) {
				nm++;
				// There's a mismatch at this position
				assert_gt((int)refcs->size(), i);
				char refChar = toupper((*refcs)[i]);
				ASSERT_ONLY(char qryChar = (h.fw ? (*pat)[i] : (*pat)[len-i-1]));
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
	if(h.color) ss << "\tCM:i:" << h.cmms.count();
	// Add optional fields reporting the primer base and the downstream color,
	// which, if they were present, were clipped when the read was read in
	if(h.color && gReportColorPrimer) {
		if(h.primer != '?') {
			ss << "\tZP:Z:" << h.primer;
			assert(isprint(h.primer));
		}
		if(h.trimc != '?') {
			ss << "\tZp:Z:" << h.trimc;
			assert(isprint(h.trimc));
		}
	}
	if(xms > 0)  ss << "\tXM:i:" << xms;
	ss << endl;
}

/**
 * Report a verbose, human-readable alignment to the appropriate
 * output stream.
 */
void SAMHitSink::reportSamHit(const Hit& h, int mapq, int xms) {
	if(xms == 0) {
		// Otherwise, this is actually a sampled read and belongs in
		// the same category as maxed reads
		HitSink::reportHit(h);
	}
	ostringstream ss;
	append(ss, h, mapq, xms);
	// Make sure to grab lock before writing to output stream
	lock(h.h.first);
	out(h.h.first).writeString(ss.str());
	unlock(h.h.first);
}

/**
 * Report a batch of hits from a vector, perhaps subsetting it.
 */
void SAMHitSink::reportSamHits(
	vector<Hit>& hs,
    size_t start,
    size_t end,
    int mapq,
    int xms)
{
	assert_geq(end, start);
	if(end-start == 0) return;
	assert_gt(hs[start].mate, 0);
	char buf[4096];
	lock(0);
	for(size_t i = start; i < end; i++) {
		ostringstream ss(ssmode_);
		ss.rdbuf()->pubsetbuf(buf, 4096);
		append(ss, hs[i], mapq, xms);
		out(0).writeChars(buf, ss.tellp());
	}
	unlock(0);
	mainlock();
	commitHits(hs);
	first_ = false;
	numAligned_++;
	numReportedPaired_ += (end-start);
	mainunlock();
}

/**
 * Report either an unaligned read or a read that exceeded the -m
 * ceiling.  We output placeholders for most of the fields in this
 * case.
 */
void SAMHitSink::reportUnOrMax(PatternSourcePerThread& p,
                               vector<Hit>* hs,
                               bool un) // lower bound on number of other hits
{
	if(un) HitSink::reportUnaligned(p);
	else   HitSink::reportMaxed(*hs, p);
	ostringstream ss;
	bool paired = !p.bufb().empty();
	assert(paired || p.bufa().mate == 0);
	assert(!paired || p.bufa().mate > 0);
	assert(un || hs->size() > 0);
	assert(!un || hs == NULL || hs->size() == 0);
	size_t hssz = 0;
	if(hs != NULL) hssz = hs->size();
	if(paired) {
		// truncate final 2 chars
		for(int i = 0; i < (int)seqan::length(p.bufa().name)-2; i++) {
			if(!noQnameTrunc_ && isspace((int)p.bufa().name[i])) break;
			ss << p.bufa().name[i];
		}
	} else {
		for(int i = 0; i < (int)seqan::length(p.bufa().name); i++) {
			if(!noQnameTrunc_ && isspace((int)p.bufa().name[i])) break;
			ss << p.bufa().name[i];
		}
	}
	ss << "\t"
	   << (SAM_FLAG_UNMAPPED | (paired ? (SAM_FLAG_PAIRED | SAM_FLAG_FIRST_IN_PAIR | SAM_FLAG_MATE_UNMAPPED) : 0)) << "\t*"
	   << "\t0\t0\t*\t*\t0\t0\t"
	   << p.bufa().patFw << "\t" << p.bufa().qual << "\tXM:i:"
	   << (paired ? (hssz+1)/2 : hssz);
	// Add optional fields reporting the primer base and the downstream color,
	// which, if they were present, were clipped when the read was read in
	if(p.bufa().color && gReportColorPrimer) {
		if(p.bufa().primer != '?') {
			ss << "\tZP:Z:" << p.bufa().primer;
			assert(isprint(p.bufa().primer));
		}
		if(p.bufa().trimc != '?') {
			ss << "\tZp:Z:" << p.bufa().trimc;
			assert(isprint(p.bufa().trimc));
		}
	}
	ss << endl;
	if(paired) {
		// truncate final 2 chars
		for(int i = 0; i < (int)seqan::length(p.bufb().name)-2; i++) {
			ss << p.bufb().name[i];
		}
		ss << "\t"
		   << (SAM_FLAG_UNMAPPED | (paired ? (SAM_FLAG_PAIRED | SAM_FLAG_SECOND_IN_PAIR | SAM_FLAG_MATE_UNMAPPED) : 0)) << "\t*"
		   << "\t0\t0\t*\t*\t0\t0\t"
		   << p.bufb().patFw << "\t" << p.bufb().qual << "\tXM:i:"
		   << (hssz+1)/2;
		// Add optional fields reporting the primer base and the downstream color,
		// which, if they were present, were clipped when the read was read in
		if(p.bufb().color && gReportColorPrimer) {
			if(p.bufb().primer != '?') {
				ss << "\tZP:Z:" << p.bufb().primer;
				assert(isprint(p.bufb().primer));
			}
			if(p.bufb().trimc != '?') {
				ss << "\tZp:Z:" << p.bufb().trimc;
				assert(isprint(p.bufb().trimc));
			}
		}
		ss << endl;
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
                        int mapq,
                        int xms,
                        const vector<string>* refnames,
                        ReferenceMap *rmap,
                        AnnotationMap *amap,
                        bool fullRef,
                        bool noQnameTrunc,
                        int offBase)
{
	appendAligned(ss, h, mapq, xms, refnames, rmap, amap, fullRef, noQnameTrunc, offBase);
}

/**
 * Report maxed-out read; if sampleMax_ is set, then report 1 alignment
 * at random.
 */
void SAMHitSink::reportMaxed(vector<Hit>& hs, PatternSourcePerThread& p) {
	if(sampleMax_) {
		HitSink::reportMaxed(hs, p);
		RandomSource rand;
		rand.init(p.bufa().seed);
		assert_gt(hs.size(), 0);
		bool paired = hs.front().mate > 0;
		size_t num = 1;
		if(paired) {
			num = 0;
			int bestStratum = 999;
			for(size_t i = 0; i < hs.size()-1; i += 2) {
				int strat = min(hs[i].stratum, hs[i+1].stratum);
				if(strat < bestStratum) {
					bestStratum = strat;
					num = 1;
				} else if(strat == bestStratum) {
					num++;
				}
			}
			assert_leq(num, hs.size());
			uint32_t r = rand.nextU32() % num;
			num = 0;
			for(size_t i = 0; i < hs.size()-1; i += 2) {
				int strat = min(hs[i].stratum, hs[i+1].stratum);
				if(strat == bestStratum) {
					if(num == r) {
						reportSamHits(hs, i, i+2, 0, (int)(hs.size()/2)+1);
						break;
					}
					num++;
				}
			}
			assert_eq(num, r);
		} else {
			for(size_t i = 1; i < hs.size(); i++) {
				assert_geq(hs[i].stratum, hs[i-1].stratum);
				if(hs[i].stratum == hs[i-1].stratum) num++;
				else break;
			}
			assert_leq(num, hs.size());
			uint32_t r = rand.nextU32() % num;
			reportSamHit(hs[r], /*MAPQ*/0, /*XM:I*/(int)hs.size()+1);
		}
	} else {
		reportUnOrMax(p, &hs, false);
	}
}
