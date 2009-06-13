#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "endian_swap.h"

/**
 * Concrete reference representation that bulk-loads the reference from
 * the bit-pair-compacted binary file and stores it in memory also in
 * bit-pair-compacted format.  The user may request reference
 * characters either on a per-character bases or by "stretch" using
 * getBase(...) and getStretch(...) respectively.
 *
 * Most of the complexity in this class is due to the fact that we want
 * to represent references with ambiguous (non-A/C/G/T) characters but
 * we don't want to use more than two bits per base.  This means we
 * need a way to encode the ambiguous stretches of the reference in a
 * way that is external to the bitpair sequence.  To accomplish this,
 * we use the RefRecords vector, which is stored in the .3.ebwt index
 * file.  The bitpairs themselves are stored in the .4.ebwt index file.
 *
 * Once it has been loaded, a BitPairReference is read-only, and is
 * safe for many threads to access at once.
 */
class BitPairReference {

public:
	/**
	 * Load from .3.ebwt/.4.ebwt Bowtie index files.
	 */
	BitPairReference(const string& in,
	                 bool sanity = false,
	                 std::vector<string>* infiles = NULL,
	                 std::vector<String<Dna5> >* origs = NULL,
	                 bool infilesSeq = false,
	                 bool useMm = false,
	                 bool verbose = false) :
	buf_(NULL),
	sanityBuf_(NULL),
	loaded_(true),
	sanity_(sanity),
	useMm_(useMm),
	verbose_(verbose)
	{
		string s3 = in + ".3.ebwt";
		string s4 = in + ".4.ebwt";
		int f3, f4;
		if((f3 = open(s3.c_str(), O_RDONLY)) < 0) {
			cerr << "Could not open reference-string index file " << s3 << " for reading." << endl;
			cerr << "This is most likely because your index was built with an older version" << endl
			     << "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
			     << "index (or download one from the Bowtie website) and try again." << endl;
			loaded_ = false;
			return;
		}
		if((f4 = open(s4.c_str(), O_RDONLY)) < 0) {
			cerr << "Could not open reference-string index file " << s4 << " for reading." << endl;
			loaded_ = false;
			return;
		}

		char *mmFile = NULL;
#ifdef BOWTIE_MM
		if(useMm_) {
			if(verbose_) {
				cout << "  Memory-mapping reference index file " << s4 << endl;
			}
			struct stat sbuf;
			if (stat(s4.c_str(), &sbuf) == -1) {
				perror("stat");
				cerr << "Error: Could not stat index file " << s4.c_str() << " prior to memory-mapping" << endl;
				exit(1);
			}
			mmFile = (char*)mmap((void *)0, sbuf.st_size,
			                     PROT_READ, MAP_SHARED, f4, 0);
			if(mmFile == (void *)(-1) || mmFile == NULL) {
				perror("mmap");
				cerr << "Error: Could not memory-map the index file " << s4.c_str() << endl;
				exit(1);
			}
		}
#endif

		// Read endianness sentinel, set 'swap'
		uint32_t one;
		bool swap = false;
		one = readU32(f3, swap);
		if(one != 1) {
			if(useMm_) {
				cerr << "Error: Can't use memory-mapped files when the index is the opposite endianness" << endl;
				exit(1);
			}
			assert_eq(0x1000000, one);
			swap = true; // have to endian swap U32s
		}

		// Read # records
		uint32_t sz;
		sz = readU32(f3, swap);
		if(sz == 0) {
			cerr << "Error: number of reference records is 0 in " << s3 << endl;
			exit(1);
		}

		// Read records
		nrefs_ = 0;

		// Cumulative count of all unambiguous characters on a per-
		// stretch 8-bit alignment (i.e. count of bytes we need to
		// allocate in buf_)
		uint32_t cumsz = 0;
		uint32_t cumlen = 0;
		// For each unambiguous stretch...
		for(uint32_t i = 0; i < sz; i++) {
			recs_.push_back(RefRecord(f3, swap));
			if(recs_.back().first) {
				// Remember that this is the first record for this
				// reference sequence (and the last record for the one
				// before)
				refRecOffs_.push_back(recs_.size()-1);
				refOffs_.push_back(cumsz);
				if(nrefs_ > 0) {
					refLens_.push_back(cumlen);
				}
				cumlen = 0;
				nrefs_++;
			} else if(i == 0) {
				cerr << "First record in reference index file was not marked as 'first'" << endl;
				exit(1);
			}
			cumsz += recs_.back().len;
			cumlen += recs_.back().off;
			cumlen += recs_.back().len;
		}
		if(verbose_) {
			cout << "Read " << nrefs_ << " reference strings from " << sz << " records" << endl;
		}
		// Store a cap entry for the end of the last reference seq
		refRecOffs_.push_back(recs_.size());
		refOffs_.push_back(cumsz);
		refLens_.push_back(cumlen);
		bufSz_ = cumsz;
		assert_eq(nrefs_, refLens_.size());
		assert_eq(sz, recs_.size());
		close(f3); // done with .3.ebwt file
		// Round cumsz up to nearest byte boundary
		if((cumsz & 3) != 0) {
			cumsz += (4 - (cumsz & 3));
		}
		bufAllocSz_ = cumsz >> 2;
		assert_eq(0, cumsz & 3); // should be rounded up to nearest 4
		if(useMm_) {
			buf_ = (uint8_t*)mmFile;
			if(sanity_) {
				FILE *ftmp = fopen(s4.c_str(), "rb");
				sanityBuf_ = new uint8_t[cumsz >> 2];
				size_t ret = fread(sanityBuf_, 1, cumsz >> 2, ftmp);
				if(ret != (cumsz >> 2)) {
					cerr << "Only read " << ret << " bytes (out of " << (cumsz >> 2) << ") from reference index file " << s4 << endl;
					exit(1);
				}
				fclose(ftmp);
				for(size_t i = 0; i < (cumsz >> 2); i++) {
					assert_eq(sanityBuf_[i], buf_[i]);
				}
			}
		} else {
			// Allocate a buffer to hold the reference string
			try {
				buf_ = new uint8_t[cumsz >> 2];
				if(buf_ == NULL) throw std::bad_alloc();
			} catch(std::bad_alloc& e) {
				cerr << "Error: Ran out of memory allocating space for the bitpacked reference.  Please" << endl
				     << "re-run on a computer with more memory." << endl;
				exit(1);
			}
			// Open the bitpair-encoded reference file
			FILE *f4 = fopen(s4.c_str(), "rb");
			if(f4 == NULL) {
				cerr << "Could not open reference-string index file " << s4 << " for reading." << endl;
				cerr << "This is most likely because your index was built with an older version" << endl
					 << "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
					 << "index (or download one from the Bowtie website) and try again." << endl;
				loaded_ = false;
				return;
			}
			// Read the whole thing in
			size_t ret = fread(buf_, 1, cumsz >> 2, f4);
			// Didn't read all of it?
			if(ret != (cumsz >> 2)) {
				cerr << "Only read " << ret << " bytes (out of " << (cumsz >> 2) << ") from reference index file " << s4 << endl;
				exit(1);
			}
			// Make sure there's no more
			char c;
			ret = fread(&c, 1, 1, f4);
			assert_eq(0, ret); // should have failed
			fclose(f4);
		}

		// Populate byteToU32_
		bool big = currentlyBigEndian();
		for(int i = 0; i < 256; i++) {
			uint32_t word = 0;
			if(big) {
				word |= ((i >> 0) & 3) << 24;
				word |= ((i >> 2) & 3) << 16;
				word |= ((i >> 4) & 3) << 8;
				word |= ((i >> 6) & 3) << 0;
			} else {
				word |= ((i >> 0) & 3) << 0;
				word |= ((i >> 2) & 3) << 8;
				word |= ((i >> 4) & 3) << 16;
				word |= ((i >> 6) & 3) << 24;
			}
			byteToU32_[i] = word;
		}

#ifndef NDEBUG
		if(sanity_) {
			// Compare the sequence we just read from the compact index
			// file to the true reference sequence.
			std::vector<seqan::String<seqan::Dna5> > *os; // for holding references
			std::vector<seqan::String<seqan::Dna5> > osv; // for holding references
			if(infiles != NULL) {
				if(infilesSeq) {
					for(size_t i = 0; i < infiles->size(); i++) {
						// Remove initial backslash; that's almost
						// certainly being used to protect the first
						// character of the sequence from getopts (e.g.,
						// when the first char is -)
						if((*infiles)[i].at(0) == '\\') {
							(*infiles)[i].erase(0, 1);
						}
						osv.push_back(String<Dna5>((*infiles)[i]));
					}
				} else {
					readSequenceFiles<seqan::String<seqan::Dna5>, seqan::Fasta>(*infiles, osv);
				}
				os = &osv;
			} else {
				assert(origs != NULL);
				os = origs;
			}
			for(size_t i = 0; i < os->size(); i++) {
				size_t olen = seqan::length((*os)[i]);
				size_t olenU32 = (olen + 12) / 4;
				uint32_t *buf = new uint32_t[olenU32];
				uint8_t *bufadj = (uint8_t*)buf;
				bufadj += getStretch(buf, i, 0, olen);
				for(size_t j = 0; j < olen; j++) {
					assert_eq((*os)[i][j], bufadj[j]);
					assert_eq((int)(*os)[i][j], getBase(i, j));
				}
				delete[] buf;
			}
		}
#endif
	}

	~BitPairReference() {
		if(buf_ != NULL && !useMm_) delete[] buf_;
		if(buf_ != NULL) delete[] sanityBuf_;
	}

	/**
	 * Return a single base of the reference.  Calling this repeatedly
	 * is not an efficient way to retrieve bases from the reference;
	 * use loadStretch() instead.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getBase(uint32_t tidx, uint32_t toff) const {
		uint32_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint32_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint32_t bufOff = refOffs_[tidx];
		uint32_t off = 0;
		// For all records pertaining to the target reference sequence...
		for(uint32_t i = reci; i < recf; i++) {
			assert_geq(toff, off);
			off += recs_[i].off;
			if(toff < off) {
				return 4;
			}
			assert_geq(toff, off);
			uint32_t recOff = off + recs_[i].len;
			if(toff < recOff) {
				toff -= off;
				bufOff += toff;
				assert_lt(bufOff, bufSz_);
				const uint32_t bufElt = (bufOff) >> 2;
				const uint32_t shift = (bufOff & 3) << 1;
				return ((buf_[bufElt] >> shift) & 3);
			}
			bufOff += recs_[i].len;
			off = recOff;
			assert_geq(toff, off);
		} // end for loop over records
		return 4;
	}

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getStretchNaive(uint32_t *destU32,
	                    uint32_t tidx,
	                    uint32_t toff,
	                    uint32_t count) const
	{
		uint8_t *dest = (uint8_t*)destU32;
		uint32_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint32_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint32_t cur = 0;
		uint32_t bufOff = refOffs_[tidx];
		uint32_t off = 0;
		// For all records pertaining to the target reference sequence...
		for(uint32_t i = reci; i < recf; i++) {
			assert_geq(toff, off);
			off += recs_[i].off;
			for(; toff < off && count > 0; toff++) {
				dest[cur++] = 4;
				count--;
			}
			if(count == 0) break;
			assert_geq(toff, off);
			if(toff < off + recs_[i].len) {
				bufOff += (toff - off); // move bufOff pointer forward
			} else {
				bufOff += recs_[i].len;
			}
			off += recs_[i].len;
			for(; toff < off && count > 0; toff++) {
				assert_lt(bufOff, bufSz_);
				const uint32_t bufElt = (bufOff) >> 2;
				const uint32_t shift = (bufOff & 3) << 1;
				dest[cur++] = (buf_[bufElt] >> shift) & 3;
				bufOff++;
				count--;
			}
			if(count == 0) break;
			assert_geq(toff, off);
		} // end for loop over records
		// In any chars are left after scanning all the records,
		// they must be ambiguous
		while(count > 0) {
			count--;
			dest[cur++] = 4;
		}
		assert_eq(0, count);
		return 0;
	}

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getStretch(uint32_t *destU32,
	               uint32_t tidx,
	               uint32_t toff,
	               uint32_t count) const
	{
		ASSERT_ONLY(uint32_t origCount = count);
		ASSERT_ONLY(uint32_t origToff = toff);
		if(count == 0) return 0;
		uint8_t *dest = (uint8_t*)destU32;
#ifndef NDEBUG
		uint32_t *destU32_2 = new uint32_t[(origCount >> 2) + 2];
		int off2 = getStretchNaive(destU32_2, tidx, origToff, origCount);
		uint8_t *dest_2 = ((uint8_t*)destU32_2) + off2;
#endif
		destU32[0] = 0x04040404; // Add Ns, which we might end up using later
		uint32_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint32_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint32_t cur = 4; // keep a cushion of 4 bases at the beginning
		uint32_t bufOff = refOffs_[tidx];
		uint32_t off = 0;
		int offset = 4;
		bool firstStretch = true;
		// For all records pertaining to the target reference sequence...
		for(uint32_t i = reci; i < recf; i++) {
			ASSERT_ONLY(uint32_t origBufOff = bufOff);
			assert_geq(toff, off);
			off += recs_[i].off;
			assert_gt(count, 0);
			if(toff < off) {
				uint32_t cpycnt = min(off - toff, count);
				memset(&dest[cur], 4, cpycnt);
				count -= cpycnt;
				toff += cpycnt;
				cur += cpycnt;
				if(count == 0) break;
			}
			assert_geq(toff, off);
			if(toff < off + recs_[i].len) {
				bufOff += (toff - off); // move bufOff pointer forward
			} else {
				bufOff += recs_[i].len;
			}
			off += recs_[i].len;
			if(toff < off) {
				if(firstStretch) {
					if(toff + 8 < off && count > 8) {
						// We already added some Ns, so we have to do
						// a fixup at the beginning of the buffer so
						// that we can start clobbering at cur >> 2
						if(cur & 3) {
							offset -= (cur & 3);
						}
						uint32_t curU32 = cur >> 2;
						// Do the initial few bases
						if(bufOff & 3) {
							const uint32_t bufElt = (bufOff) >> 2;
							const int low2 = bufOff & 3;
							destU32[curU32] = byteToU32_[buf_[bufElt]];
							for(int j = 0; j < low2; j++) {
								((char *)(&destU32[curU32]))[j] = 4;
							}
							curU32++;
							offset += low2;
							const int chars = 4 - low2;
							count -= chars;
							bufOff += chars;
							toff += chars;
						}
						assert_eq(0, bufOff & 3);
						uint32_t bufOffU32 = bufOff >> 2;
						uint32_t countLim = count >> 2;
						uint32_t offLim = (off - (toff + 4)) >> 2;
						uint32_t lim = min(countLim, offLim);
						// Do the fast thing for as far as possible
						for(uint32_t j = 0; j < lim; j++) {
							destU32[curU32] = byteToU32_[buf_[bufOffU32++]];
							assert_eq(dest[(curU32 << 2) + 0], dest_2[(curU32 << 2) - offset + 0]);
							assert_eq(dest[(curU32 << 2) + 1], dest_2[(curU32 << 2) - offset + 1]);
							assert_eq(dest[(curU32 << 2) + 2], dest_2[(curU32 << 2) - offset + 2]);
							assert_eq(dest[(curU32 << 2) + 3], dest_2[(curU32 << 2) - offset + 3]);
							curU32++;
						}
						toff += (lim << 2);
						assert_leq(toff, off);
						assert_leq((lim << 2), count);
						count -= (lim << 2);
						bufOff = bufOffU32 << 2;
						cur = curU32 << 2;
					}
					// Do the slow thing for the rest
					for(; toff < off && count > 0; toff++) {
						assert_lt(bufOff, bufSz_);
						const uint32_t bufElt = (bufOff) >> 2;
						const uint32_t shift = (bufOff & 3) << 1;
						dest[cur++] = (buf_[bufElt] >> shift) & 3;
						bufOff++;
						count--;
					}
					firstStretch = false;
				} else {
					// Do the slow thing
					for(; toff < off && count > 0; toff++) {
						assert_lt(bufOff, bufSz_);
						const uint32_t bufElt = (bufOff) >> 2;
						const uint32_t shift = (bufOff & 3) << 1;
						dest[cur++] = (buf_[bufElt] >> shift) & 3;
						bufOff++;
						count--;
					}
				}
			}
			if(count == 0) break;
			assert_eq(recs_[i].len, bufOff - origBufOff);
			assert_geq(toff, off);
		} // end for loop over records
		// In any chars are left after scanning all the records,
		// they must be ambiguous
		while(count > 0) {
			count--;
			dest[cur++] = 4;
		}
		assert_eq(0, count);
#ifndef NDEBUG
		delete[] destU32_2;
#endif
		return offset;
	}

	/// Return the number of reference sequences.
	uint32_t numRefs() const {
		return nrefs_;
	}

	/// Return the number of reference sequences.
	uint32_t approxLen(uint32_t elt) const {
		assert_lt(elt, nrefs_);
		return refLens_[elt];
	}

	/// Return true iff buf_ and all the vectors are populated.
	bool loaded() const {
		return loaded_;
	}

protected:

	uint32_t byteToU32_[256];

	std::vector<RefRecord> recs_;       /// records describing unambiguous stretches
	std::vector<uint32_t>  refLens_;    /// approx lens of ref seqs (excludes trailing ambig chars)
	std::vector<uint32_t>  refOffs_;    /// buf_ begin offsets per ref seq
	std::vector<uint32_t>  refRecOffs_; /// record begin/end offsets per ref seq
	uint8_t *buf_;      /// the whole reference as a big bitpacked byte array
	uint8_t *sanityBuf_;/// for sanity-checking buf_
	uint32_t bufSz_;    /// size of buf_
	uint32_t bufAllocSz_;
	uint32_t nrefs_;    /// the number of reference sequences
	bool     loaded_;   /// whether it's loaded
	bool     sanity_;   /// do sanity checking
	bool     useMm_;    /// load the reference as a memory-mapped file
	bool     verbose_;
};

#endif
