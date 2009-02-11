#include "ref_read.h"

/**
 * Reads past the next sequence from the given FASTA file and returns
 * the size of the passed sequence.  Does not do anything with the
 * sequence characters themselves; this is purely for counting lengths.
 */
RefRecord fastaRefReadSize(FileBuf& in,
                           const RefReadInParams& refparams,
                           bool first,
                           BitpairOutFileBuf* bpout)
{
	int c;
	static int lastc = '>';
	assert_neq(refparams.baseCutoff, 0);
	assert_neq(refparams.numSeqCutoff, 0);

	// RefRecord params
	size_t seqCharsRead = 0;
	size_t seqOff = 0;
	bool seqFirst = true;

	// Pick off the first carat and any preceding whitespace
	if(first) {
		assert(!in.eof());
		lastc = '>';
		c = skipWhitespace(in);
		if(in.eof()) {
			// Got eof right away; emit warning
			cerr << "Warning: Empty input file" << endl;
			lastc = -1;
			return RefRecord(seqOff, seqCharsRead, seqFirst);
		}
		assert(c == '>' || c == '#');
	}
	// Skip to the end of the id line; if the next line is either
	// another id line or a comment line, keep skipping
	if(lastc == '>') {
		// Skip to the end of the name line
		do {
			if((c = skipLine(in)) == -1) {
				cerr << "Warning: Encountered empty reference sequence" << endl;
				lastc = -1;
				return RefRecord(seqOff, seqCharsRead, seqFirst);
			}
			if(c == '>') {
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
		} while (c == '>' || c == '#');
	} else {
		// Skip until we arrive at a '>', in which case let lastc = '>'
		// and return, or until we arrive at a DNA char
		seqFirst = false;
		seqOff = 1; // The gap has already been consumed, so count it
		c = in.get(); // Get next char
		if(c == -1) {
			// Don't emit a warning, since this might legitimately be
			// a gap on the end of the final sequence in the file
			lastc = -1;
			return RefRecord(seqOff, seqCharsRead, seqFirst);
		}
	}

	// Now skip to the first DNA character, counting gap characters
	// as we go
	while(true) {
		if(refparams.nsToAs && dna4Cat[c] == 2) {
			// Turn this 'N' (or other ambiguous char) into an 'A'
			c = 'A';
		}
		int cat = dna4Cat[c];
		ASSERT_ONLY(int cc = toupper(c));
		if(cat == 1) {
			// This is a DNA character
			assert(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T');
			break; // to read-in loop
		} else if(cat == 2) {
			assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
			seqOff++; // skip over gap character and increment
		} else if(c == '>') {
			if(seqOff > 0 && lastc == '>') {
				cerr << "Warning: Encountered reference sequence with only gaps" << endl;
			} else if(lastc == '>') {
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			lastc = '>';
			return RefRecord(seqOff, seqCharsRead, seqFirst);
		}
		c = in.get();
		if(c == -1) {
			// End-of-file
			if(seqOff > 0 && lastc == '>') {
				cerr << "Warning: Encountered reference sequence with only gaps" << endl;
			} else if(lastc == '>') {
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			lastc = -1;
			return RefRecord(seqOff, seqCharsRead, seqFirst);
		}
	}
	assert_eq(1, dna4Cat[c]);

	// in now points just past the first character of a sequence
	// line, and c holds the first character
	while(c != -1 && c != '>') {
		if(refparams.nsToAs && dna4Cat[c] == 2) {
			// Turn this 'N' (or other ambiguous char) into an 'A'
			c = 'A';
		}
		uint8_t cat = dna4Cat[c];
		ASSERT_ONLY(int cc = toupper(c));
		if(cat == 1) {
			// It's a DNA character
			assert(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T');
			// Consume it
			seqCharsRead++;
			// Output it
			if(bpout) bpout->write(charToDna5[c]);
			if(refparams.baseCutoff != -1 &&
			   (int64_t)seqCharsRead >= refparams.baseCutoff)
			{
				lastc = -1;
				return RefRecord(seqOff, seqCharsRead, seqFirst);
			}
		} else if(cat == 2) {
			// It's an N or a gap
			lastc = c;
			assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
			return RefRecord(seqOff, seqCharsRead, seqFirst);
		} else {
			// Not DNA and not a gap, ignore it
#ifndef NDEBUG
			if(!isspace(c)) {
				cerr << "Unexpected character in sequence: ";
				if(isprint(c)) {
					cerr << ((char)c) << endl;
				} else {
					cerr << "(" << c << ")" << endl;
				}
			}
#endif
		}
		c = in.get();
	}
	lastc = c;
	return RefRecord(seqOff, seqCharsRead, seqFirst);
}

/**
 * Calculate a vector containing the sizes of all of the patterns in
 * all of the given input files, in order.  Returns the total size of
 * all references combined.  Rewinds each istream before returning.
 */
size_t fastaRefReadSizes(vector<FileBuf*>& in,
                         vector<RefRecord>& recs,
                         const RefReadInParams& refparams,
                         BitpairOutFileBuf* bpout)
{
	uint32_t tot = 0;
	RefReadInParams rpcp = refparams;
	assert_gt(in.size(), 0);
	// For each input istream
	for(size_t i = 0; i < in.size(); i++) {
		bool first = true;
		assert(!in[i]->eof());
		assert_geq(rpcp.baseCutoff, -1);
		assert_neq(rpcp.baseCutoff, 0);
		assert_geq(rpcp.numSeqCutoff, -1);
		assert_neq(rpcp.numSeqCutoff, 0);
		// For each pattern in this istream
		while(!in[i]->eof() && rpcp.baseCutoff != 0 && rpcp.numSeqCutoff != 0) {
			RefRecord rec = fastaRefReadSize(*in[i], refparams, first, bpout);
			if(rpcp.baseCutoff > 0) assert_leq((int64_t)rec.len, rpcp.baseCutoff);
			if(rpcp.baseCutoff != -1)   rpcp.baseCutoff -= rec.len;
			if(rpcp.numSeqCutoff != -1) rpcp.numSeqCutoff--;
			assert_geq(rpcp.baseCutoff, -1);
			assert_geq(rpcp.numSeqCutoff, -1);
			if((tot + rec.len) < tot) {
				cerr << "Error: Reference sequence has more than 2^32-1 characters!  Please divide the" << endl
				     << "reference into batches or chunks of about 3.6 billion characters or less each" << endl
				     << "and index each independently." << endl;
				exit(1);
			}
			tot += rec.len;
			first = false;
			if(rec.len == 0 && rec.off == 0 && !rec.first) continue;
			recs.push_back(rec);
		}
		in[i]->reset();
		assert(!in[i]->eof());
		#ifndef NDEBUG
		int c = in[i]->get();
		assert_eq('>', c);
		in[i]->reset();
		assert(!in[i]->eof());
		#endif
	}
	assert_gt(tot, 0);
	return tot;
}
