#include "ref_read.h"

/**
 * Reads past the next ambiguous or unambiguous stretch of sequence
 * from the given FASTA file and returns its length.  Does not do
 * anything with the sequence characters themselves; this is purely for
 * measuring lengths.
 */
RefRecord fastaRefReadSize(FileBuf& in,
                           const RefReadInParams& rparms,
                           bool first,
                           BitpairOutFileBuf* bpout)
{
	int c;
	static int lastc = '>'; // last character seen
	assert_neq(rparms.baseCutoff, 0); // should have stopped
	assert_neq(rparms.numSeqCutoff, 0); // should have stopped

	// RefRecord params
	size_t len = 0; // 'len' counts toward total length
	// 'off' counts number of ambiguous characters before first
	// unambiguous character
	size_t off = 0;

	// Pick off the first carat and any preceding whitespace
	if(first) {
		assert(!in.eof());
		lastc = '>';
		c = in.getPastWhitespace();
		if(in.eof()) {
			// Got eof right away; emit warning
			cerr << "Warning: Empty input file" << endl;
			lastc = -1;
			return RefRecord(0, 0, true);
		}
		// TODO: having # could yield spurious warnings (below)
		assert(c == '>' || c == '#');
	}

	first = true;
	// Skip to the end of the id line; if the next line is either
	// another id line or a comment line, keep skipping
	if(lastc == '>') {
		// Skip to the end of the name line
		do {
			if((c = in.getPastNewline()) == -1) {
				// No more input
				cerr << "Warning: Encountered empty reference sequence" << endl;
				lastc = -1;
				return RefRecord(0, 0, true);
			}
			if(c == '>') {
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			// continue until a non-name, non-comment line
		} while (c == '>' || c == '#');
	} else {
		first = false; // not the first in a sequence
		off = 1; // The gap has already been consumed, so count it
		if((c = in.get()) == -1) {
			// Don't emit a warning, since this might legitimately be
			// a gap on the end of the final sequence in the file
			lastc = -1;
			return RefRecord(off, len, first);
		}
	}

	// Now skip to the first DNA character, counting gap characters
	// as we go
	int lc = -1; // last-DNA char variable for color conversion
	while(true) {
		int cat = dna4Cat[c];
		if(rparms.nsToAs && cat == 2) {
			// Turn this 'N' (or other ambiguous char) into an 'A'
			c = 'A';
		}
		if(cat == 1) {
			// This is a DNA character
			if(rparms.color) {
				if(lc != -1) {
					// Got two consecutive unambiguous DNAs
					break; // to read-in loop
				}
				// Keep going; we need two consecutive unambiguous DNAs
				lc = charToDna5[(int)c];
				// The 'if(off > 0)' takes care of the case where
				// the reference is entirely unambiguous and we don't
				// want to incorrectly increment off.
				if(off > 0) off++;
			} else {
				break; // to read-in loop
			}
		} else if(cat == 2) {
			lc = -1;
			off++; // skip over gap character and increment
		} else if(c == '>') {
			if(off > 0 && lastc == '>') {
				cerr << "Warning: Encountered reference sequence with only gaps" << endl;
			} else if(lastc == '>') {
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			lastc = '>';
			return RefRecord(0, 0, false);
		}
		c = in.get();
		if(c == -1) {
			// End-of-file
			if(off > 0 && lastc == '>') {
				cerr << "Warning: Encountered reference sequence with only gaps" << endl;
			} else if(lastc == '>') {
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			lastc = -1;
			return RefRecord(0, 0, false);
		}
	}
	assert(!rparms.color || (lc != -1));
	assert_eq(1, dna4Cat[c]); // C must be unambiguous base
	if(off > 0 && rparms.color && first) {
		// Handle the case where the first record has ambiguous
		// characters but we're in color space; one of those counts is
		// spurious
		off--;
	}

	// in now points just past the first character of a sequence
	// line, and c holds the first character
	while(c != -1 && c != '>') {
		if(rparms.nsToAs && dna4Cat[c] == 2) {
			// Turn this 'N' (or other ambiguous char) into an 'A'
			c = 'A';
		}
		uint8_t cat = dna4Cat[c];
		int cc = toupper(c);
		if(rparms.bisulfite && cc == 'C') c = cc = 'T';
		if(cat == 1) {
			// It's a DNA character
			assert(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T');
			// Consume it
			len++;
			// Output it
			if(bpout != NULL) {
				if(rparms.color) {
					// output color
					bpout->write(dinuc2color[charToDna5[(int)c]][lc]);
				} else if(!rparms.color) {
					// output nucleotide
					bpout->write(charToDna5[c]);
				}
			}
			if(rparms.baseCutoff != -1 &&
			   (int64_t)len >= rparms.baseCutoff)
			{
				lastc = -1;
				return RefRecord(off, len, first);
			}
			lc = charToDna5[(int)c];
		} else if(cat == 2) {
			// It's an N or a gap
			lastc = c;
			assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
			return RefRecord(off, len, first);
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
	return RefRecord(off, len, first);
}

/**
 * Calculate a vector containing the sizes of all of the patterns in
 * all of the given input files, in order.  Returns the total size of
 * all references combined.  Rewinds each istream before returning.
 */
std::pair<size_t, size_t>
fastaRefReadSizes(vector<FileBuf*>& in,
                  vector<RefRecord>& recs,
                  const RefReadInParams& rparms,
                  BitpairOutFileBuf* bpout,
                  int& numSeqs)
{
	uint32_t unambigTot = 0;
	uint32_t bothTot = 0;
	RefReadInParams rpcp = rparms;
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
			RefRecord rec = fastaRefReadSize(*in[i], rparms, first, bpout);
			if(rpcp.baseCutoff > 0) assert_leq((int64_t)rec.len, rpcp.baseCutoff);
			if(rpcp.baseCutoff != -1)   rpcp.baseCutoff -= rec.len;
			if(rpcp.numSeqCutoff != -1) rpcp.numSeqCutoff--;
			assert_geq(rpcp.baseCutoff, -1);
			assert_geq(rpcp.numSeqCutoff, -1);
			if((unambigTot + rec.len) < unambigTot) {
				cerr << "Error: Reference sequence has more than 2^32-1 characters!  Please divide the" << endl
				     << "reference into batches or chunks of about 3.6 billion characters or less each" << endl
				     << "and index each independently." << endl;
				throw 1;
			}
			// Add the length of this record.
			if(rec.first) numSeqs++;
			unambigTot += rec.len;
			bothTot += rec.len;
			bothTot += rec.off;
			first = false;
			if(rec.len == 0 && rec.off == 0 && !rec.first) continue;
			recs.push_back(rec);
		}
		// Reset the input stream
		in[i]->reset();
		assert(!in[i]->eof());
#ifndef NDEBUG
		// Check that it's really reset
		int c = in[i]->get();
		assert_eq('>', c);
		in[i]->reset();
		assert(!in[i]->eof());
#endif
	}
	assert_gt(bothTot, 0);
	assert_gt(unambigTot, 0);
	return make_pair(
		unambigTot, // total number of unambiguous DNA characters read
		bothTot); // total number of DNA characters read, incl. ambiguous ones
}
