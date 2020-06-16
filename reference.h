#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <stdexcept>

#include "btypes.h"
#include "ds.h"
#include "endian_swap.h"
#include "mm.h"
#include "ref_read.h"
#include "sequence_io.h"
#include "shmem.h"
#include "sstring.h"
#include "timer.h"
#include "word_io.h"

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
	                 bool sanity,
	                 EList<string>* infiles,
	                 EList<BTRefString >* origs,
	                 bool infilesSeq,
	                 bool loadSequence, // as opposed to just records
	                 bool useMm,
	                 bool useShmem,
	                 bool mmSweep,
	                 bool verbose,
	                 bool startVerbose) :
	buf_(NULL),
	sanityBuf_(NULL),
	loaded_(true),
	sanity_(sanity),
	useMm_(useMm),
	useShmem_(useShmem),
	verbose_(verbose)
	{
		string s3 = in + ".3." + gEbwt_ext;
		string s4 = in + ".4." + gEbwt_ext;

		FILE *f3, *f4;
		if((f3 = fopen(s3.c_str(), "rb")) == NULL) {
			cerr << "Could not open reference-string index file " << s3 << " for reading." << endl;
			cerr << "This is most likely because your index was built with an older version" << endl
			     << "(<= 0.9.8.1) of bowtie-build.  Please re-run bowtie-build to generate a new" << endl
			     << "index (or download one from the Bowtie website) and try again." << endl;
			loaded_ = false;
			return;
		}
		if((f4 = fopen(s4.c_str(), "rb")) ==NULL) {
			cerr << "Could not open reference-string index file " << s4 << " for reading." << endl;
			loaded_ = false;
			return;
		}
#ifdef BOWTIE_MM
		char *mmFile = NULL;
		if(useMm_) {
			if(verbose_ || startVerbose) {
				cerr << "  Memory-mapping reference index file " << s4 << ": ";
				logTime(cerr);
			}
			struct stat sbuf;
			if (stat(s4.c_str(), &sbuf) == -1) {
				perror("stat");
				cerr << "Error: Could not stat index file " << s4.c_str() << " prior to memory-mapping" << endl;
				throw 1;
			}
			mmFile = (char*)mmap((void *)0, sbuf.st_size,
			                     PROT_READ, MAP_SHARED, fileno(f4), 0);
			if(mmFile == (void *)(-1) || mmFile == NULL) {
				perror("mmap");
				cerr << "Error: Could not memory-map the index file " << s4.c_str() << endl;
				throw 1;
			}
			if(mmSweep) {
				TIndexOff sum = 0;
				for(off_t i = 0; i < sbuf.st_size; i += 1024) {
					sum += (TIndexOff) mmFile[i];
				}
				if(startVerbose) {
					cerr << "  Swept the memory-mapped ref index file; checksum: " << sum << ": ";
					logTime(cerr);
				}
			}
		}
#endif

		// Read endianness sentinel, set 'swap'
		uint32_t one;
		bool swap = false;
		one = readU<int32_t>(f3, swap);
		if(one != 1) {
			if(useMm_) {
				cerr << "Error: Can't use memory-mapped files when the index is the opposite endianness" << endl;
				throw 1;
			}
			assert_eq(0x1000000, one);
			swap = true; // have to endian swap U32s
		}

		// Read # records
		TIndexOffU sz;
		sz = readU<TIndexOffU>(f3, swap);
		if(sz == 0) {
			cerr << "Error: number of reference records is 0 in " << s3 << endl;
			throw 1;
		}

		// Read records
		nrefs_ = 0;
		nNoGapRefs_ = 0;

		// Cumulative count of all unambiguous characters on a per-
		// stretch 8-bit alignment (i.e. count of bytes we need to
		// allocate in buf_)
		TIndexOffU cumsz = 0;
		TIndexOffU cumlen = 0;
		TIndexOffU unambiglen = 0;
		TIndexOffU maxlen = 0;
		// For each unambiguous stretch...
		for(TIndexOffU i = 0; i < sz; i++) {
			recs_.push_back(RefRecord(f3, swap));
		}
		for(TIndexOffU i = 0; i < sz; i++) {
			if(recs_[i].first) {
				if(nrefs_ > 0) {
					refLens_.push_back(cumlen);
				}
				// Stupid hack to get around the fact that, for
				// colorspace Ebwts, not only do we omit reference
				// sequences that are *all* gaps from
				// nPat/plen/refnames, but we also omit reference
				// sequences that never have a stretch of more than 1
				// unambiguous character.
				if(unambiglen > 0 && maxlen > 1) {
					refApproxLens_.push_back(cumlen);
				}
				// More hackery to detect references that won't be
				// in the Ebwt even though they have non-zero length
				bool willBeInEbwt = true;
				if(recs_[i].len > 0 && willBeInEbwt) {
					// Remember that this is the first record for this
					// reference sequence (and the last record for the one
					// before)
					refRecOffs_.push_back(i);
					refOffs_.push_back(cumsz);
					expandIdx_.push_back(nrefs_);
					shrinkIdx_.push_back(nNoGapRefs_);
					nNoGapRefs_++;
					isGaps_.push_back(false);
				} else {
					shrinkIdx_.push_back(nNoGapRefs_);
					isGaps_.push_back(true);
				}
				cumlen = 0;
				unambiglen = 0;
				maxlen = 0;
				nrefs_++;
				assert_eq(nNoGapRefs_, expandIdx_.size());
				assert_eq(nrefs_, shrinkIdx_.size());
			} else if(i == 0) {
				//cerr << "First record in reference index file was not marked as 'first'" << endl;
				//throw 1;
			}
			cumUnambig_.push_back(cumsz);
			cumRefOff_.push_back(cumlen);
			cumsz += recs_[i].len;
#ifdef ACCOUNT_FOR_ALL_GAP_REFS
			cumlen += recs_[i].off;
			cumlen += recs_[i].len;
#else
			if(recs_[i].len > 0) {
				cumlen += recs_[i].off;
				cumlen += recs_[i].len;
			}
#endif
			unambiglen += recs_[i].len;
			if(recs_[i].len > maxlen) maxlen = recs_[i].len;
		}
		if(verbose_ || startVerbose) {
			cerr << "Read " << nrefs_ << " reference strings (" << nNoGapRefs_ << " non-empty) from " << sz << " records: ";
			logTime(cerr);
		}
		// Store a cap entry for the end of the last reference seq
		refRecOffs_.push_back((TIndexOffU)recs_.size());
		refOffs_.push_back(cumsz);
		if(unambiglen > 0 && maxlen > 1) {
			refApproxLens_.push_back(cumlen);
		}
		refLens_.push_back(cumlen);
		cumUnambig_.push_back(cumsz);
		cumRefOff_.push_back(cumlen);
		bufSz_ = cumsz;
		assert_eq(nNoGapRefs_, refApproxLens_.size());
		assert_eq(sz, recs_.size());
		if (f3 != NULL) fclose(f3); // done with .3.ebwt file
		// Round cumsz up to nearest byte boundary
		if((cumsz & 3) != 0) {
			cumsz += (4 - (cumsz & 3));
		}
		bufAllocSz_ = cumsz >> 2;
		assert_eq(0, cumsz & 3); // should be rounded up to nearest 4
		if(!loadSequence) return;
		if(useMm_) {
#ifdef BOWTIE_MM
			buf_ = (uint8_t*)mmFile;
			if(sanity_) {
				FILE *ftmp = fopen(s4.c_str(), "rb");
				sanityBuf_ = new uint8_t[cumsz >> 2];
				size_t ret = fread(sanityBuf_, 1, cumsz >> 2, ftmp);
				if(ret != (cumsz >> 2)) {
					cerr << "Only read " << ret << " bytes (out of " << (cumsz >> 2) << ") from reference index file " << s4 << endl;
					throw 1;
				}
				fclose(ftmp);
				for(size_t i = 0; i < (cumsz >> 2); i++) {
					assert_eq(sanityBuf_[i], buf_[i]);
				}
			}
#else
			cerr << "Shouldn't be at " << __FILE__ << ":" << __LINE__ << " without BOWTIE_MM defined" << endl;
			throw 1;
#endif
		} else {
			bool shmemLeader = true;
			if(!useShmem_) {
				// Allocate a buffer to hold the reference string
				try {
					buf_ = new uint8_t[cumsz >> 2];
					if(buf_ == NULL) throw std::bad_alloc();
				} catch(std::bad_alloc& e) {
					cerr << "Error: Ran out of memory allocating space for the bitpacked reference.  Please" << endl
						 << "re-run on a computer with more memory." << endl;
					throw 1;
				}
			} else {
				shmemLeader = ALLOC_SHARED_U8(
					(s4 + "[ref]"), (cumsz >> 2), &buf_,
					"ref", (verbose_ || startVerbose));
			}
			if(shmemLeader) {
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
					throw 1;
				}
				// Make sure there's no more
				char c;
				ret = fread(&c, 1, 1, f4);
				assert_eq(0, ret); // should have failed
				fclose(f4);
				if(useShmem_) NOTIFY_SHARED(buf_, (cumsz >> 2));
			} else {
				if(useShmem_) WAIT_SHARED(buf_, (cumsz >> 2));
			}
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
			EList<BTRefString > *os; // for holding references
			EList<BTRefString > osv; // for holding references
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
						osv.push_back(BTRefString((*infiles)[i].c_str(), true));
					}
				} else {
					readSequenceFiles(*infiles, osv);
				}
				os = &osv;
			} else {
				assert(origs != NULL);
				os = origs;
			}

			// Go through the loaded reference files base-by-base and
			// sanity check against what we get by calling getBase and
			// getStretch
			size_t refi = 0;
			int longestStretch = 0;
			int curStretch = 0;
			for(size_t i = 0; i < os->size(); i++) {
				size_t olen = (*os)[i].length();
				for(size_t j = 0; j < olen; j++) {
					if((int)(*os)[i][j] < 4) {
						curStretch++;
						if(curStretch > longestStretch) longestStretch = curStretch;
					} else {
						curStretch = 0;
					}
				}
				if(longestStretch == 0) {
					continue;
				}
				longestStretch = 0;
				size_t olenU32 = (olen + 12) / 4;
				uint32_t *buf = new uint32_t[olenU32];
				uint8_t *bufadj = (uint8_t*)buf;
				bufadj += getStretch(buf, refi, 0, olen);
				for(size_t j = 0; j < olen; j++) {
					assert_eq((int)(*os)[i][j], (int)bufadj[j]);
					assert_eq((int)(*os)[i][j], (int)getBase(refi, j));
				}
				refi++;
				delete[] buf;
			}
		}
#endif
	}

	~BitPairReference() {
		if(buf_ != NULL && !useMm_ && !useShmem_) delete[] buf_;
		if(sanityBuf_ != NULL) delete[] sanityBuf_;
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
	int getBase(size_t tidx, size_t toff) const {
		uint64_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint64_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint64_t bufOff = refOffs_[tidx];
		uint64_t off = 0;
		// For all records pertaining to the target reference sequence...
		for(uint64_t i = reci; i < recf; i++) {
			assert_geq(toff, off);
			off += recs_[i].off;
			if(toff < off) {
				return 4;
			}
			assert_geq(toff, off);
			uint64_t recOff = off + recs_[i].len;
			if(toff < recOff) {
				toff -= off;
				bufOff += toff;
				assert_lt(bufOff, bufSz_);
				const uint64_t bufElt = (bufOff) >> 2;
				const uint64_t shift = (bufOff & 3) << 1;
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
	                    size_t tidx,
	                    size_t toff,
	                    size_t count) const
	{
		uint8_t *dest = (uint8_t*)destU32;
		uint64_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint64_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint64_t cur = 0;
		uint64_t bufOff = refOffs_[tidx];
		uint64_t off = 0;
		// For all records pertaining to the target reference sequence...
		for(uint64_t i = reci; i < recf; i++) {
			assert_geq(toff, off);
			off += recs_[i].off;
			for(; toff < off && count > 0; toff++) {
				dest[cur++] = 4;
				count--;
			}
			if(count == 0) break;
			assert_geq(toff, off);
			if(toff < off + recs_[i].len) {
				bufOff += (TIndexOffU)(toff - off); // move bufOff pointer forward
			} else {
				bufOff += recs_[i].len;
			}
			off += recs_[i].len;
			for(; toff < off && count > 0; toff++) {
				assert_lt(bufOff, bufSz_);
				const uint64_t bufElt = (bufOff) >> 2;
				const uint64_t shift = (bufOff & 3) << 1;
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
	               size_t tidx,
	               size_t toff,
	               size_t count) const
	{
		ASSERT_ONLY(size_t origCount = count);
		ASSERT_ONLY(size_t origToff = toff);
		if(count == 0) return 0;
		uint8_t *dest = (uint8_t*)destU32;
#ifndef NDEBUG
		uint32_t *destU32_2 = new uint32_t[(origCount >> 2) + 2];
		int off2 = getStretchNaive(destU32_2, tidx, origToff, origCount);
		uint8_t *dest_2 = ((uint8_t*)destU32_2) + off2;
#endif
		destU32[0] = 0x04040404; // Add Ns, which we might end up using later
		uint64_t reci = refRecOffs_[tidx];   // first record for target reference sequence
		uint64_t recf = refRecOffs_[tidx+1]; // last record (exclusive) for target seq
		assert_gt(recf, reci);
		uint64_t cur = 4; // keep a cushion of 4 bases at the beginning
		uint64_t bufOff = refOffs_[tidx];
		uint64_t off = 0;
		int64_t offset = 4;
		bool firstStretch = true;
		ASSERT_ONLY(bool binarySearched = false);
		uint64_t left  = reci;
		uint64_t right = recf;
		uint64_t mid   = 0;
		// For all records pertaining to the target reference sequence...
		for(uint64_t i = reci; i < recf; i++) {
#ifndef NDEBUG
			uint64_t origBufOff = bufOff;
#endif
			assert_geq(toff, off);
			if (firstStretch && recf > reci + 16){
				// binary search finds smallest i s.t. toff >= cumRefOff_[i]
				while (left < right-1) {
					mid = left + ((right - left) >> 1);
					if (cumRefOff_[mid] <= toff)
						left = mid;
					else
						right = mid;
				}
				off = cumRefOff_[left];
				bufOff = cumUnambig_[left];
#ifndef NDEBUG
				origBufOff = bufOff;
#endif
				i = left;
				assert(cumRefOff_[i+1] == 0 || cumRefOff_[i+1] > toff);
				ASSERT_ONLY(binarySearched = true);
			}
			off += recs_[i].off; // skip Ns at beginning of stretch
			assert_gt(count, 0);
			if(toff < off) {
				size_t cpycnt = min((size_t)(off - toff), count);
				memset(&dest[cur], 4, cpycnt);
				count -= cpycnt;
				toff += cpycnt;
				cur += cpycnt;
				if(count == 0) break;
			}
			assert_geq(toff, off);
			if(toff < off + recs_[i].len) {
				bufOff += toff - off; // move bufOff pointer forward
			} else {
				bufOff += recs_[i].len;
			}
			off += recs_[i].len;
			assert(off == cumRefOff_[i+1] || cumRefOff_[i+1] == 0);
			assert(!binarySearched || toff < off);
			if(toff < off) {
				if(firstStretch) {
					if(toff + 8 < off && count > 8) {
						// We already added some Ns, so we have to do
						// a fixup at the beginning of the buffer so
						// that we can start clobbering at cur >> 2
						if(cur & 3) {
							offset -= (cur & 3);
						}
						uint64_t curU32 = cur >> 2;
						// Do the initial few bases
						if(bufOff & 3) {
							const uint64_t bufElt = (bufOff) >> 2;
							const int64_t low2 = bufOff & 3;
							// Lots of cache misses on the following line
							destU32[curU32] = byteToU32_[buf_[bufElt]];
							for(int j = 0; j < low2; j++) {
								((char *)(&destU32[curU32]))[j] = 4;
							}
							curU32++;
							offset += low2;
							const int64_t chars = 4 - low2;
							count -= chars;
							bufOff += chars;
							toff += chars;
						}
						assert_eq(0, bufOff & 3);
						uint64_t bufOffU32 = bufOff >> 2;
						uint64_t countLim = count >> 2;
						uint64_t offLim = ((off - (toff + 4)) >> 2);
						uint64_t lim = min(countLim, offLim);
						// Do the fast thing for as far as possible
						for(uint64_t j = 0; j < lim; j++) {
							// Lots of cache misses on the following line
							destU32[curU32] = byteToU32_[buf_[bufOffU32++]];
#ifndef NDEBUG
							if(dest_2 != NULL) {
								assert_eq(dest[(curU32 << 2) + 0], dest_2[(curU32 << 2) - offset + 0]);
								assert_eq(dest[(curU32 << 2) + 1], dest_2[(curU32 << 2) - offset + 1]);
								assert_eq(dest[(curU32 << 2) + 2], dest_2[(curU32 << 2) - offset + 2]);
								assert_eq(dest[(curU32 << 2) + 3], dest_2[(curU32 << 2) - offset + 3]);
							}
#endif
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
						const uint64_t bufElt = (bufOff) >> 2;
						const uint64_t shift = (bufOff & 3) << 1;
						dest[cur++] = (buf_[bufElt] >> shift) & 3;
						bufOff++;
						count--;
					}
					firstStretch = false;
				} else {
					// Do the slow thing
					for(; toff < off && count > 0; toff++) {
						assert_lt(bufOff, bufSz_);
						const uint64_t bufElt = (bufOff) >> 2;
						const uint64_t shift = (bufOff & 3) << 1;
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
		return (int)offset;
	}

	/// Return the number of reference sequences.
	TIndexOffU numRefs() const {
		return nrefs_;
	}

	/// Return the number of reference sequences that don't consist
	/// entirely of stretches of ambiguous characters.
	uint32_t numNonGapRefs() const {
		return nNoGapRefs_;
	}

	/**
	 *
	 */
	uint32_t shrinkIdx(uint32_t idx) const {
		assert_lt(idx, shrinkIdx_.size());
		return shrinkIdx_[idx];
	}

	/**
	 *
	 */
	uint32_t expandIdx(uint32_t idx) const {
		assert_lt(idx, expandIdx_.size());
		return expandIdx_[idx];
	}

	/// Return the lengths of reference sequences.
	TIndexOffU approxLen(TIndexOffU elt) const {
		assert_lt(elt, nrefs_);
		return refApproxLens_[elt];
	}

	/// Return the lengths of reference sequences.
	uint32_t len(uint32_t elt) const {
		assert_lt(elt, nrefs_);
		return refLens_[elt];
	}

	/// Return true iff ref 'elt' is all gaps
	bool isAllGaps(uint32_t elt) const {
		assert_lt(elt, nrefs_);
		assert_eq(isGaps_.size(), nrefs_);
		return isGaps_[elt];
	}

	/// Return true iff buf_ and all the vectors are populated.
	bool loaded() const {
		return loaded_;
	}

	/**
	 * Return constant reference to the RefRecord list.
	 */
	const EList<RefRecord>& refRecords() const { return recs_; }

protected:

	uint32_t byteToU32_[256];

	EList<RefRecord> recs_;       /// records describing unambiguous stretches
	EList<uint32_t>  refApproxLens_; /// approx lens of ref seqs (excludes trailing ambig chars)
	EList<TIndexOffU>  refLens_;    /// approx lens of ref seqs (excludes trailing ambig chars)
	EList<TIndexOffU>  refOffs_;    /// buf_ begin offsets per ref seq
	EList<TIndexOffU>  cumUnambig_;    /// # unambig ref chars up to each record
	EList<TIndexOffU>  cumRefOff_;    /// # ref chars up to each record
	EList<TIndexOffU>  refRecOffs_; /// record begin/end offsets per ref seq
	EList<uint32_t>  expandIdx_; /// map from small idxs (e.g. w/r/t plen) to large ones (w/r/t refnames)
	EList<uint32_t>  shrinkIdx_; /// map from large idxs to small
	EList<bool>      isGaps_;    /// ref i is all gaps?
	uint8_t *buf_;      /// the whole reference as a big bitpacked byte array
	uint8_t *sanityBuf_;/// for sanity-checking buf_
	TIndexOffU bufSz_;    /// size of buf_
	TIndexOffU bufAllocSz_;
	TIndexOffU nrefs_;      /// the number of reference sequences
	uint32_t nNoGapRefs_; /// the number of reference sequences that aren't totally ambiguous
	bool     loaded_;   /// whether it's loaded
	bool     sanity_;   /// do sanity checking
	bool     useMm_;    /// load the reference as a memory-mapped file
	bool     useShmem_; /// load the reference into shared memory
	bool     verbose_;
};

#endif
