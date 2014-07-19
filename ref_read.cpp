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

	// RefRecord params
	TIndexOffU len = 0; // 'len' counts toward total length
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
		assert(c == '>');
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
		} while (c == '>');
	} else {
		first = false; // not the first in a sequence
		off = 1; // The gap has already been consumed, so count it
		if((c = in.get()) == -1) {
			// Don't emit a warning, since this might legitimately be
			// a gap on the end of the final sequence in the file
			lastc = -1;
			return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
		}
	}

	// Now skip to the first DNA character, counting gap characters
	// as we go
	int lc = -1; // last-DNA char variable for color conversion
	while(true) {
		int cat = dna4Cat[c];
		if(rparms.nsToAs && cat == 2) c = 'A';
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
			if(lc != -1 && off == 0) off++;
			lc = -1;
			off++; // skip over gap character and increment
		} else if(c == '>') {
			if(off > 0 && lastc == '>') {
				cerr << "Warning: Encountered reference sequence with only gaps" << endl;
			} else if(lastc == '>') {
				cerr << "Warning: Encountered empty reference sequence" << endl;
			}
			lastc = '>';
			return RefRecord((TIndexOffU)off, 0, first);
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
			return RefRecord((TIndexOffU)off, 0, first);
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
		if(rparms.nsToAs && dna4Cat[c] == 2) c = 'A';
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
			lc = charToDna5[(int)c];
		} else if(cat == 2) {
			// It's an N or a gap
			lastc = c;
			assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
			return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
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
	return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
}

static void
printRecords(ostream& os, const vector<RefRecord>& l) {
	for(size_t i = 0; i < l.size(); i++) {
		os << l[i].first << ", " << l[i].off << ", " << l[i].len << endl;
	}
}

/**
 * Reverse the 'src' list of RefRecords into the 'dst' list.  Don't
 * modify 'src'.
 */
void reverseRefRecords(const vector<RefRecord>& src,
					   vector<RefRecord>& dst,
					   bool recursive,
					   bool verbose)
{
	dst.clear();
	{
		vector<RefRecord> cur;
		for(int64_t i = (int64_t)src.size()-1; i >= 0; i--) {
			bool first = (i == (int)src.size()-1 || src[i+1].first);
			if(src[i].len) {
				cur.push_back(RefRecord(0, src[i].len, first));
				first = false;
			}
			if(src[i].off) cur.push_back(RefRecord(src[i].off, 0, first));
		}
		for(int64_t i = 0; i < (int64_t)cur.size(); i++) {
			assert(cur[i].off == 0 || cur[i].len == 0);
			if(i < (int64_t)cur.size()-1 && cur[i].off != 0 && !cur[i+1].first) {
				dst.push_back(RefRecord(cur[i].off, cur[i+1].len, cur[i].first));
				i++;
			} else {
				dst.push_back(cur[i]);
			}
		}
	}
	if(verbose) {
		cout << "Source: " << endl;
		printRecords(cout, src);
		cout << "Dest: " << endl;
		printRecords(cout, dst);
	}
#ifndef NDEBUG
	if(!recursive) {
		vector<RefRecord> tmp;
		reverseRefRecords(dst, tmp, true);
		assert_eq(tmp.size(), src.size());
		for(size_t i = 0; i < src.size(); i++) {
			assert_eq(src[i].len, tmp[i].len);
			assert_eq(src[i].off, tmp[i].off);
			assert_eq(src[i].first, tmp[i].first);
		}
	}
#endif
}

/**
 * Calculate a vector containing the sizes of all of the patterns in
 * all of the given input files, in order.  Returns the total size of
 * all references combined.  Rewinds each istream before returning.
 */
std::pair<size_t, size_t>
fastaRefReadSizes(vector<FileBuf*>& in,
                  vector<RefRecord>& recs,
                  vector<uint32_t>& plens,
                  const RefReadInParams& rparms,
                  BitpairOutFileBuf* bpout,
                  TIndexOff& numSeqs)
{
	TIndexOffU unambigTot = 0;
	size_t bothTot = 0;
	assert_gt(in.size(), 0);
	uint32_t both = 0, unambig = 0;
	// For each input istream
	for(size_t i = 0; i < in.size(); i++) {
		bool first = true;
		assert(!in[i]->eof());
		// For each pattern in this istream
		while(!in[i]->eof()) {
			RefRecord rec = fastaRefReadSize(*in[i], rparms, first, bpout);
			// Update plens
			if(rec.first) {
				if(unambig > 0) {
					plens.push_back(both);
				}
				both = 0;
				unambig = 0;
			}
#ifndef ACCOUNT_FOR_ALL_GAP_REFS
			if(rec.len == 0) rec.first = false;
#endif
			if((unambigTot + rec.len) < unambigTot) {
#ifdef BOWTIE_64BIT_INDEX
				cerr << "Error: Reference sequence has more than 2^32-1 characters!  Please divide the" << endl
				     << "reference into smaller chunks and index each independently." << endl;
#else
				cerr << "Error: Reference sequence has more than 2^32-1 characters!  Please try to" << endl
				     << "build a large index instead using the appropiate options." << endl;
#endif
				throw 1;
			}
			// Add the length of this record.
			if(rec.first) numSeqs++;
			unambigTot += rec.len; unambig += rec.len;
			bothTot += rec.len;    both += rec.len;
			bothTot += rec.off;    both += rec.off;
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
	assert_geq(bothTot, 0);
	assert_geq(unambigTot, 0);
	if(unambig > 0) {
		plens.push_back(both);
	}
	return make_pair(
		unambigTot, // total number of unambiguous DNA characters read
		bothTot); // total number of DNA characters read, incl. ambiguous ones
}
