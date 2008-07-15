#include "ref_read.h"

/**
 * Reads past the next sequence from the given FASTA FILE* and returns
 * the size of the passed sequence.  Does not do anything with the
 * sequence characters themselves; this is purely for counting lengths.
 * Does not mutate refparams.
 */
size_t fastaRefReadSize(istream& in,
                        RefReadInParams& refparams, 
                        bool first = false)
{
	int c;
	assert_neq(refparams.baseCutoff, 0);
	assert_neq(refparams.numSeqCutoff, 0);
	size_t seqCharsRead = 0;
	// Pick off the first carat
	if(first) {
		c = in.get(); if(in.eof()) return seqCharsRead;
		assert(c == '>' || c == '#');
	}
	// Skip to the end of the id line; if the next line is either
	// another id line or a comment line, keep skipping
	// TODO: grab name here, stick in *name
	do {
		if((c = skipLine(in)) == -1) return seqCharsRead;
	} while (c == '>' || c == '#');
	// in now points just past the first character of a sequence
	// line, and c holds the first character
	while(c != '>' && c != '#') {
		// Note: can't have a comment in the middle of a sequence,
		// though a comment can end a sequence
		if(isalpha(c)) {
			seqCharsRead++;
			if((int64_t)seqCharsRead >= refparams.baseCutoff) {
				return seqCharsRead;
			}
		}
		c = in.get();
		if(in.eof()) break;
	}
	return seqCharsRead;
}

/**
 * Calculate a vector containing the sizes of all of the patterns in
 * all of the given input files, in order.  Returns the total size of
 * all references combined.  Rewinds each istream before returning.
 * Does not mutate refparams.
 */
size_t fastaRefReadSizes(vector<istream*>& in,
                         vector<uint32_t>& szs,
                         RefReadInParams& refparams)
{
	size_t tot = 0;
	RefReadInParams rpcp = refparams;
	assert_gt(in.size(), 0);
	// For each input istream
	for(size_t i = 0; i < in.size(); i++) {
		bool first = true;
		assert(in[i]->good());
		streampos pos = in[i]->tellg();
		assert_geq(rpcp.baseCutoff, -1);
		assert_neq(rpcp.baseCutoff, 0);
		assert_geq(rpcp.numSeqCutoff, -1);
		assert_neq(rpcp.numSeqCutoff, 0);
		// For each pattern in this istream
		while(in[i]->good() && rpcp.baseCutoff != 0 && rpcp.numSeqCutoff != 0) {
			uint32_t sz = fastaRefReadSize(*in[i], refparams, first);
			assert_leq(sz, rpcp.baseCutoff);
			if(rpcp.baseCutoff != -1)   rpcp.baseCutoff -= sz;
			if(rpcp.numSeqCutoff != -1) rpcp.numSeqCutoff--;
			assert_geq(rpcp.baseCutoff, -1);
			assert_geq(rpcp.numSeqCutoff, -1);
			tot += sz;
			first = false;
			if(sz == 0) continue;
			szs.push_back(sz);
		}
		in[i]->clear();
		in[i]->seekg(pos);
		assert(!in[i]->bad());
		assert(!in[i]->fail());
		in[i]->clear();
		assert(in[i]->good());
		assert(!in[i]->eof());
		#ifndef NDEBUG
		int c = in[i]->get();
		assert_eq('>', c);
		assert(in[i]->good());
		assert(!in[i]->eof());
		in[i]->seekg(pos);
		in[i]->clear();
		assert(in[i]->good());
		assert(!in[i]->eof());
		#endif
	}
	assert_gt(tot, 0);
	return tot;
}
