/*
 * ref_aligner.h
 */

#ifndef REF_ALIGNER_H_
#define REF_ALIGNER_H_

#include <stdint.h>
#include <iostream>
#include <vector>
#include "seqan/sequence.h"
#include "range.h"
#include "reference.h"

// Let the reference-aligner buffer size be 16K by default.  If more
// room is required, a new buffer must be allocated from the heap.
const static int REF_ALIGNER_BUFSZ = 16 * 1024;

/**
 * Abstract parent class for classes that look for alignments by
 * matching against the reference sequence directly.  This is useful
 * both for sanity-checking results from the Bowtie index and for
 * finding mates when the reference location of the opposite mate is
 * known.
 */
template<typename TStr>
class RefAligner {

	typedef seqan::String<seqan::Dna5> TDna5Str;
	typedef seqan::String<char> TCharStr;
	typedef std::vector<uint32_t> TU32Vec;
	typedef std::vector<Range> TRangeVec;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	RefAligner(uint32_t seedLen) :
		seedLen_(seedLen), refbuf_(buf_), refbufSz_(REF_ALIGNER_BUFSZ),
		freeRefbuf_(false) { }

	/**
	 * Free the reference-space alignment buffer if this object
	 * allocated it.
	 */
	virtual ~RefAligner() {
		if(freeRefbuf_) {
			delete[] refbuf_;
		}
	}

	/**
	 * Find one alignment of qry:quals in the range begin-end in
	 * reference string ref.  Store the alignment details in range.
	 */
	virtual void find(uint32_t numToFind,
	                  const uint32_t tidx,
	                  const BitPairReference *refs,
	                  const TDna5Str& qry,
	                  const TCharStr& quals,
	                  uint32_t begin,
	                  uint32_t end,
	                  TRangeVec& ranges,
	                  TU32Vec& results,
	                  TSetPairs* pairs = NULL,
	                  uint32_t aoff = 0xffffffff,
	                  bool low = false)
	{
		assert_gt(numToFind, 0);
		uint32_t spread = end - begin;
		uint32_t spreadPlus = spread + 12;
		// Make sure the buffer is large enough to accommodate the spread
		if(spreadPlus > this->refbufSz_) {
			this->newBuf(spreadPlus);
		}
		// Read in the relevant stretch of the reference string
		int offset = refs->getStretch(this->refbuf_, tidx, begin, spread);
		uint8_t *buf = ((uint8_t*)this->refbuf_) + offset;
		// Look for alignments
		return seed64Find(numToFind, tidx, buf, qry, quals, begin,
		                  end, ranges, results, pairs, aoff, low);
	}

	virtual void seed64Find(uint32_t numToFind,
	                uint32_t tidx,
	                uint8_t* ref,
                    const TDna5Str& qry,
                    const TCharStr& quals,
                    uint32_t begin,
                    uint32_t end,
	                TRangeVec& ranges,
	                TU32Vec& results,
	                TSetPairs* pairs = NULL,
	                uint32_t aoff = 0xffffffff,
	                bool low = false) const = 0;

	/**
	 * Set a new reference-sequence buffer.
	 */
	void setBuf(uint32_t *newbuf, uint32_t newsz) {
		if(freeRefbuf_) {
			delete[] refbuf_;
			freeRefbuf_ = false;
		}
		refbuf_ = newbuf;
		refbufSz_ = newsz;
	}

	/**
	 * Set a new reference-sequence buffer.
	 */
	void newBuf(uint32_t newsz) {
		if(freeRefbuf_) {
			delete[] refbuf_;
		}
		try {
			refbuf_ = new uint32_t[(newsz + 3) / 4];
			if(refbuf_ == NULL) throw std::bad_alloc();
		} catch(std::bad_alloc& e) {
			cerr << "Error: Could not allocate reference-space alignment buffer of " << newsz << "B" << endl;
			exit(1);
		}
		refbufSz_ = newsz;
		freeRefbuf_ = true;
	}

protected:
	uint32_t  seedLen_;   /// length of seed region for read
	uint32_t *refbuf_;    /// pointer to current reference buffer
	uint32_t  refbufSz_;  /// size of current reference buffer
	uint32_t  buf_[REF_ALIGNER_BUFSZ / 4]; /// built-in reference buffer (may be superseded)
	bool      freeRefbuf_; /// whether refbuf_ points to something we should delete
};

/**
 * Concrete RefAligner for finding nearby exact hits given an anchor
 * hit.
 */
template<typename TStr>
class ExactRefAligner : public RefAligner<TStr> {

	typedef seqan::String<seqan::Dna5> TDna5Str;
	typedef seqan::String<char> TCharStr;
	typedef std::vector<uint32_t> TU32Vec;
	typedef std::vector<Range> TRangeVec;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	ExactRefAligner(uint32_t seedLen) : RefAligner<TStr>(seedLen) { }
	virtual ~ExactRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
	               uint32_t tidx,
	               uint8_t* ref,
                   const TDna5Str& qry,
                   const TCharStr& quals,
                   uint32_t begin,
                   uint32_t end,
	               TRangeVec& ranges,
	               TU32Vec& results,
	               TSetPairs* pairs,
	               uint32_t aoff,
	               bool low) const
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
					// Mismatch
					match = false;
					break;
				}
				// Match; continue
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
					// Mismatch
					match = false;
					break;
				}
				// Match; continue
#endif
			}
			if(match) {
				ranges.resize(ranges.size()+1);
				Range& range = ranges.back();
				range.stratum = 0;
				range.numMms = 0;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				results.push_back(ri);
			}
		}
		return; // no match
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void seed64Find(uint32_t numToFind,
                    uint32_t tidx,
	                uint8_t *ref,
                    const TDna5Str& qry,
                    const TCharStr& quals,
                    uint32_t begin,
                    uint32_t end,
 	                TRangeVec& ranges,
	                TU32Vec& results,
	                TSetPairs* pairs,
	                uint32_t aoff, // offset of anchor mate
	                bool low) const
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t rangesInitSz = ranges.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get all naive hits
		TRangeVec r2; TU32Vec re2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          re2, pairs, aoff, low);
#endif
		const uint32_t seedBitPairs = min<int>(qlen, 32);
		const uint32_t seedOverhang = qlen <= 32 ? 0 : qlen - 32;
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t seed = 0llu;
		uint64_t buffw = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the seed area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(seedBitPairs < 32) {
			clearMask >>= ((32-seedBitPairs) << 1);
			useMask = true;
		}
		const int lhsShift = ((seedBitPairs - 1) << 1);
		// Build the contents of the 'seed' dword and the initial
		// contents of the 'buffw' dword.  If there are fewer than 32
		// seedBitPairs, the content will be packed into the least
		// significant bits of the word.
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		for(uint32_t i = 0; i < seedBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			assert_leq(c, 4);
			if(c & 4) {
				assert_eq(r2.size(), ranges.size() - rangesInitSz);
				return; // can't match if query has Ns
			}
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, seedBitPairs - i);
			}
			assert_lt(r, 4);
			assert_lt(c, 4);
			seed  = ((seed  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns outside of the seed
		// region
		for(uint32_t i = seedBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				assert_eq(r2.size(), ranges.size() - rangesInitSz);
				return; // can't match if query has Ns
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the seed along until
		// it's 'seedOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'seedOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiSeed = rirHi + seedBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++; rirHi++; rirHiSeed++;
				r = (int)ref[rirHiSeed];
				if(r & 4) {
					r = 0;
					skipLeftToRights = seedBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				if(buffw != seed) {
					continue;
				}
			} else {
				hi = true;
				// Moving right-to-left
				riLo--; rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				if(bufbw != seed) {
					continue;
				}
			}
			// Seed hit!
			bool foundHit = true;
			uint32_t ri = hi ? riLo : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			if(seedOverhang > 0) {
				// Does the non-seed part of the alignment (the
				// "overhang") ruin it?
				bool skipCandidate = false;
				for(uint32_t j = 0; j < seedOverhang; j++) {
					assert_lt(ri + seedBitPairs + j, end);
					int rc = (int)ref[rir + seedBitPairs + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = seedOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = seedBitPairs + j;
						}
						skipCandidate = true;
						break; // Skip this candidate
					}
					if((int)qry[32 + j] != rc) {
						// Yes, overhang ruins it
						foundHit = false;
						break;
					}
				}
				if(skipCandidate) continue;
			}
			if(foundHit) {
				if(pairs != NULL) {
					TU64Pair p;
					if(ri < aoff) {
						// By convention, the upstream ("low") mate's
						// coordinates go in the 'first' field
						p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
						p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
					} else {
						p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
						p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					}
					if(pairs->find(p) != pairs->end()) {
						// We already found this hit!  Continue.
						ASSERT_ONLY(duplicates++);
						ASSERT_ONLY(r2i++);
						continue;
					} else {
						// Record this hit
						pairs->insert(p);
					}
				}
				assert_lt(r2i, r2.size());
				assert_eq(re2[r2i], ri);
				ranges.resize(ranges.size()+1);
				Range& range = ranges.back();
				range.stratum = 0;
				range.numMms = 0;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				ASSERT_ONLY(r2i++);
				results.push_back(ri);
				if(--numToFind == 0) return;
			} else {
				// Keep scanning
			}
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, ranges.size() - rangesInitSz);
		return; // no hit
	}
};

/**
 * Defined in ref_aligner.cpp.  Maps an octet representing the XOR of
 * two two-bit-per-base-encoded DNA sequences to the number of bases
 * that mismatch between the two.
 */
extern unsigned char u8toMms[];

/**
 * Concrete RefAligner for finding nearby 1-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class OneMMRefAligner : public RefAligner<TStr> {

	typedef seqan::String<seqan::Dna5> TDna5Str;
	typedef seqan::String<char> TCharStr;
	typedef std::vector<Range> TRangeVec;
	typedef std::vector<uint32_t> TU32Vec;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	OneMMRefAligner(uint32_t seedLen) : RefAligner<TStr>(seedLen) { }
	virtual ~OneMMRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
	               uint32_t tidx,
	               uint8_t* ref,
                   const TDna5Str& qry,
                   const TCharStr& quals,
                   uint32_t begin,
                   uint32_t end,
	               TRangeVec& ranges,
	               TU32Vec& results,
	               TSetPairs* pairs,
	               uint32_t aoff,
	               bool low) const
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc = -1;
			uint32_t mmOff = 0xffffffff;
			int mms = 0;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					if(++mms > 1) {
						// Too many; reject this alignment
						match = false;
						break;
					} else {
						// First one; remember offset and ref char
						refc = "ACGTN"[r];
						mmOff = j;
					}
				}
	 		}
			if(match) {
				assert_leq(mms, 1);
				ranges.resize(ranges.size()+1);
				Range& range = ranges.back();
				range.stratum = mms;
				range.numMms = mms;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				if(mms == 1) {
					assert_lt(mmOff, qlen);
					range.mms.push_back(mmOff);
					range.refcs.push_back(refc);
				}
				results.push_back(ri);
			}
		}
		return;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void seed64Find(uint32_t numToFind,
	                uint32_t tidx,
	                uint8_t* ref,
                    const TDna5Str& qry,
                    const TCharStr& quals,
                    uint32_t begin,
                    uint32_t end,
	                TRangeVec& ranges,
	                TU32Vec& results,
	                TSetPairs* pairs = NULL,
	                uint32_t aoff = 0xffffffff,
	                bool low = false) const
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t rangesInitSz = ranges.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		TRangeVec r2; TU32Vec re2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
		          re2, pairs, aoff, low);
#endif
		const uint32_t seedBitPairs = min<int>(qlen, 32);
		const int lhsShift = ((seedBitPairs - 1) << 1);
		const uint32_t seedCushion  = 32 - seedBitPairs;
		const uint32_t seedOverhang = (qlen <= 32 ? 0 : (qlen - 32));
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t seed = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the seed area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(seedBitPairs < 32) {
			clearMask >>= ((32-seedBitPairs) << 1);
			useMask = true;
		}
		int nsInSeed = 0;
		int nPos = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		// Construct the 'seed' 64-bit buffer so that it holds all of
		// the first 'seedBitPairs' bit pairs of the query.
		for(uint32_t i = 0; i < seedBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, seedBitPairs - i);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInSeed > 1) {
					// More than one 'N' in the seed region; can't
					// possibly have a 1-mismatch hit anywhere
					assert_eq(r2.size(), ranges.size() - rangesInitSz);
					return;   // can't match if query has Ns
				}
				nPos = (int)i;
				// Make it look like an 'A' in the seed
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			seed  = ((seed  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		// Check whether read is disqualified by Ns outside of the seed
		// region
		for(uint32_t i = seedBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				if(++nsInSeed > 1) {
					assert_eq(r2.size(), ranges.size() - rangesInitSz);
					return; // can't match if query has Ns
				}
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the seed along until
		// it's 'seedOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'seedOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiSeed = rirHi + seedBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiSeed++;
				r = (int)ref[rirHiSeed];
				if(r & 4) {
					r = 0;
					skipLeftToRights = seedBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ seed) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ seed) | diffMask;
			}
			if((diff & 0xffffffff00000000llu) &&
			   (diff & 0x00000000ffffffffllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the seed
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 1) continue;
			diffs += u8toMms[(int)diff8[1]] +
			         u8toMms[(int)diff8[2]] +
			         u8toMms[(int)diff8[3]] +
			         u8toMms[(int)diff8[4]] +
			         u8toMms[(int)diff8[5]] +
			         u8toMms[(int)diff8[6]];
			uint32_t mmpos = 0xffffffff;
			int refc = -1;
			if(diffs > 1) {
				// Too many differences
				continue;
			} else if(diffs == 1 && nPos != -1) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos = nPos;
				refc = "ACGT"[(int)ref[rir + nPos]];
			} else if(diffs == 1) {
				// Figure out which position mismatched
				mmpos = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos--; }
				assert_neq(0, diff);
				assert_geq(mmpos, 0);
				assert_lt(mmpos, 32);
				mmpos -= seedCushion;
				refc = "ACGT"[(int)ref[rir + mmpos]];
			}
			// Now extend the seed into a longer alignment
			bool foundHit = true;
			if(seedOverhang > 0) {
				assert_leq(ri + seedBitPairs + seedOverhang, end);
				bool skipCandidate = false;
				for(uint32_t j = 0; j < seedOverhang; j++) {
					int rc = (int)ref[rir + 32 + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = seedOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = seedBitPairs + j;
						}
						skipCandidate = true; // Skip this candidate
						break;
					}
					if((int)qry[32 + j] != rc) {
						if(++diffs > 1) {
							foundHit = false;
							break;
						} else {
							mmpos = 32 + j;
							refc = "ACGT"[(int)ref[rir + 32 + j]];
						}
					}
				}
				if(skipCandidate) continue;
			}
			if(!foundHit) continue;
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream ("low") mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			assert_leq(diffs, 1);
			assert_lt(r2i, r2.size());
			assert_eq(re2[r2i], ri);
			ranges.resize(ranges.size()+1);
			Range& range = ranges.back();
			assert_eq(diffs, r2[r2i].numMms);
			range.stratum = diffs;
			range.numMms = diffs;
			assert_eq(0, range.mms.size());
			assert_eq(0, range.refcs.size());
			if(diffs == 1) {
				assert_neq(mmpos, 0xffffffff);
				assert_eq(mmpos, r2[r2i].mms[0]);
				assert_neq(-1, refc);
				assert_eq(refc, r2[r2i].refcs[0]);
				range.mms.push_back(mmpos);
				range.refcs.push_back(refc);
			}
			ASSERT_ONLY(r2i++);
			results.push_back(ri);
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, ranges.size() - rangesInitSz);
		return; // no hit
	}
};

/**
 * Concrete RefAligner for finding nearby 2-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class TwoMMRefAligner : public RefAligner<TStr> {

	typedef seqan::String<seqan::Dna5> TDna5Str;
	typedef seqan::String<char> TCharStr;
	typedef std::vector<Range> TRangeVec;
	typedef std::vector<uint32_t> TU32Vec;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	TwoMMRefAligner(uint32_t seedLen) : RefAligner<TStr>(seedLen) { }
	virtual ~TwoMMRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
				   uint32_t tidx,
				   uint8_t* ref,
				   const TDna5Str& qry,
				   const TCharStr& quals,
				   uint32_t begin,
				   uint32_t end,
				   TRangeVec& ranges,
				   TU32Vec& results,
				   TSetPairs* pairs,
				   uint32_t aoff,
				   bool low) const
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc1 = -1;
			uint32_t mmOff1 = 0xffffffff;
			int refc2 = -1;
			uint32_t mmOff2 = 0xffffffff;
			int mms = 0;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					if(++mms > 2) {
						// Too many; reject this alignment
						match = false;
						break;
					} else if(mms == 2) {
						// Second one; remember offset and ref char
						refc2 = "ACGTN"[r];
						mmOff2 = j;
					} else {
						assert_eq(1, mms);
						// First one; remember offset and ref char
						refc1 = "ACGTN"[r];
						mmOff1 = j;
					}
				}
			}
			if(match) {
				assert_leq(mms, 2);
				ranges.resize(ranges.size()+1);
				Range& range = ranges.back();
				range.stratum = mms;
				range.numMms = mms;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				if(mms > 0) {
					assert_lt(mmOff1, qlen);
					assert(refc1 == 'A' || refc1 == 'C' || refc1 == 'G' || refc1 == 'T');
					range.mms.push_back(mmOff1);
					range.refcs.push_back(refc1);
					if(mms > 1) {
						assert_eq(2, mms);
						assert_lt(mmOff2, qlen);
						assert(refc2 == 'A' || refc2 == 'C' || refc2 == 'G' || refc2 == 'T');
						range.mms.push_back(mmOff2);
						range.refcs.push_back(refc2);
					}
				}
				results.push_back(ri);
			}
		}
		return;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void seed64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					TRangeVec& ranges,
					TU32Vec& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool low = false) const
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t rangesInitSz = ranges.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		TRangeVec r2; TU32Vec re2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
				  re2, pairs, aoff, low);
#endif
		const uint32_t seedBitPairs = min<int>(qlen, 32);
		const int lhsShift = ((seedBitPairs - 1) << 1);
		const uint32_t seedCushion  = 32 - seedBitPairs;
		const uint32_t seedOverhang = (qlen <= 32 ? 0 : (qlen - 32));
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t seed = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the seed area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(seedBitPairs < 32) {
			clearMask >>= ((32-seedBitPairs) << 1);
			useMask = true;
		}
		int nsInSeed = 0;
		uint32_t nPoss = 0;
		int nPos1 = -1;
		int nPos2 = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		// Construct the 'seed' 64-bit buffer so that it holds all of
		// the first 'seedBitPairs' bit pairs of the query.
		for(uint32_t i = 0; i < seedBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, seedBitPairs - i);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInSeed > 2) {
					// More than two 'N's in the seed region; can't
					// possibly have a 2-mismatch hit anywhere
					return;   // can't match if query has Ns
				} else if(nsInSeed == 2) {
					nPos2 = (int)i;
					nPoss++;
				} else {
					assert_eq(1, nsInSeed);
					nPos1 = (int)i;
					nPoss++;
				}
				// Make it look like an 'A' in the seed
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			seed  = ((seed  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		assert_leq(nPoss, 2);
		// Check whether read is disqualified by Ns outside of the seed
		// region
		for(uint32_t i = seedBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				if(++nsInSeed > 2) {
					return; // can't match if query has Ns
				}
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the seed along until
		// it's 'seedOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'seedOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiSeed = rirHi + seedBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		uint32_t i;
		for(i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiSeed++;
				r = (int)ref[rirHiSeed];
				if(r & 4) {
					r = 0;
					skipLeftToRights = seedBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ seed) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ seed) | diffMask;
			}
			if((diff & 0xfffff00000000000llu) &&
			   (diff & 0x00000ffffff00000llu) &&
			   (diff & 0x00000000000fffffllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the seed
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 2) continue;
			diffs += u8toMms[(int)diff8[1]] +
					 u8toMms[(int)diff8[2]] +
					 u8toMms[(int)diff8[3]] +
					 u8toMms[(int)diff8[4]] +
					 u8toMms[(int)diff8[5]] +
					 u8toMms[(int)diff8[6]];
			uint32_t mmpos1 = 0xffffffff;
			int refc1 = -1;
			uint32_t mmpos2 = 0xffffffff;
			int refc2 = -1;
			if(diffs > 2) {
				// Too many differences
				continue;
			} else if(nPoss > 1 && diffs == nPoss) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos1 = nPos1;
				refc1 = "ACGT"[(int)ref[rir + nPos1]];
				if(nPoss == 2) {
					mmpos2 = nPos2;
					refc2 = "ACGT"[(int)ref[rir + nPos2]];
				}
			} else if(diffs > 0) {
				// Figure out which position mismatched
				uint64_t diff2 = diff;
				mmpos1 = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos1 -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos1 -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos1 -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos1 -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos1--; }
				assert_neq(0, diff);
				assert_geq(mmpos1, 0);
				assert_lt(mmpos1, 32);
				mmpos1 -= seedCushion;
				refc1 = "ACGT"[(int)ref[rir + mmpos1]];
				if(diffs > 1) {
					// Figure out the second mismatched position
					ASSERT_ONLY(uint64_t origDiff2 = diff2);
					diff2 &= ~(0xc000000000000000llu >> (uint64_t)((mmpos1+seedCushion) << 1));
					assert_neq(diff2, origDiff2);
					mmpos2 = 31;
					if((diff2 & 0xffffffffllu) == 0) { diff2 >>= 32llu; mmpos2 -= 16; }
					assert_neq(0, diff2);
					if((diff2 & 0xffffllu) == 0)     { diff2 >>= 16llu; mmpos2 -=  8; }
					assert_neq(0, diff2);
					if((diff2 & 0xffllu) == 0)       { diff2 >>= 8llu;  mmpos2 -=  4; }
					assert_neq(0, diff2);
					if((diff2 & 0xfllu) == 0)        { diff2 >>= 4llu;  mmpos2 -=  2; }
					assert_neq(0, diff2);
					if((diff2 & 0x3llu) == 0)        { mmpos2--; }
					assert_neq(0, diff2);
					assert_geq(mmpos2, 0);
					assert_lt(mmpos2, 32);
					mmpos2 -= seedCushion;
					assert_neq(mmpos1, mmpos2);
					refc2 = "ACGT"[(int)ref[rir + mmpos2]];
					if(mmpos2 < mmpos1) {
						uint32_t mmtmp = mmpos1;
						mmpos1 = mmpos2;
						mmpos2 = mmtmp;
						int refctmp = refc1;
						refc1 = refc2;
						refc2 = refctmp;
					}
					assert_lt(mmpos1, mmpos2);
				}
			}
			// Now extend the seed into a longer alignment
			bool foundHit = true;
			if(seedOverhang > 0) {
				assert_leq(ri + seedBitPairs + seedOverhang, end);
				bool skipCandidate = false;
				for(uint32_t j = 0; j < seedOverhang; j++) {
					int rc = (int)ref[rir + 32 + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = seedOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = seedBitPairs + j;
						}
						skipCandidate = true; // Skip this candidate
						break;
					}
					if((int)qry[32 + j] != rc) {
						if(++diffs > 2) {
							foundHit = false;
							break;
						} else if(diffs == 2) {
							mmpos2 = 32 + j;
							refc2 = "ACGT"[(int)ref[rir + 32 + j]];
						} else {
							assert_eq(1, diffs);
							mmpos1 = 32 + j;
							refc1 = "ACGT"[(int)ref[rir + 32 + j]];
						}
					}
				}
				if(skipCandidate) continue;
			}
			if(!foundHit) continue;
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream ("low") mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			assert_leq(diffs, 2);
			assert_lt(r2i, r2.size());
			assert_eq(re2[r2i], ri);
			ranges.resize(ranges.size()+1);
			Range& range = ranges.back();
			assert_eq(diffs, r2[r2i].numMms);
			range.stratum = diffs;
			range.numMms = diffs;
			assert_eq(0, range.mms.size());
			assert_eq(0, range.refcs.size());
			if(diffs > 0) {
				assert_neq(mmpos1, 0xffffffff);
				assert_eq(mmpos1, r2[r2i].mms[0]);
				assert_neq(-1, refc1);
				assert_eq(refc1, r2[r2i].refcs[0]);
				range.mms.push_back(mmpos1);
				range.refcs.push_back(refc1);
				if(diffs > 1) {
					assert_neq(mmpos2, 0xffffffff);
					assert_eq(mmpos2, r2[r2i].mms[1]);
					assert_neq(-1, refc2);
					assert_eq(refc2, r2[r2i].refcs[1]);
					range.mms.push_back(mmpos2);
					range.refcs.push_back(refc2);
				}
			}
			ASSERT_ONLY(r2i++);
			results.push_back(ri);
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, ranges.size() - rangesInitSz);
		return; // no hit
	}
};

/**
 * Concrete RefAligner for finding nearby 2-mismatch hits given an
 * anchor hit.
 */
template<typename TStr>
class ThreeMMRefAligner : public RefAligner<TStr> {

	typedef seqan::String<seqan::Dna5> TDna5Str;
	typedef seqan::String<char> TCharStr;
	typedef std::vector<Range> TRangeVec;
	typedef std::vector<uint32_t> TU32Vec;
	typedef std::pair<uint64_t, uint64_t> TU64Pair;
	typedef std::set<TU64Pair> TSetPairs;

public:
	ThreeMMRefAligner(uint32_t seedLen) : RefAligner<TStr>(seedLen) { }
	virtual ~ThreeMMRefAligner() { }

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void naiveFind(uint32_t numToFind,
				   uint32_t tidx,
				   uint8_t* ref,
				   const TDna5Str& qry,
				   const TCharStr& quals,
				   uint32_t begin,
				   uint32_t end,
				   TRangeVec& ranges,
				   TU32Vec& results,
				   TSetPairs* pairs,
				   uint32_t aoff,
				   bool low) const
	{
		assert_gt(numToFind, 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		uint32_t qend = end - qlen;
		uint32_t lim = qend - begin;
		uint32_t halfway = begin + (lim >> 1);
		bool hi = false;
		for(uint32_t i = 1; i <= lim+1; i++) {
			uint32_t ri;  // leftmost position in candidate alignment
			uint32_t rir; // same, minus begin; for indexing into ref[]
			if(hi) {
				ri = halfway + (i >> 1); rir = ri - begin;
				assert_leq(ri, qend);
			} else {
				ri = halfway - (i >> 1); rir = ri - begin;
				assert_geq(ri, begin);
			}
			hi = !hi;
			// Do the naive comparison
			bool match = true;
			int refc1 = -1;
			uint32_t mmOff1 = 0xffffffff;
			int refc2 = -1;
			uint32_t mmOff2 = 0xffffffff;
			int refc3 = -1;
			uint32_t mmOff3 = 0xffffffff;
			int mms = 0;
			for(uint32_t j = 0; j < qlen; j++) {
#if 0
				// Count Ns in the reference as mismatches
				const int q = (int)qry[j];
				const int r = (int)ref[rir + j];
				assert_leq(q, 4);
				assert_leq(r, 4);
				if(q == 4 || r == 4 || q != r) {
#else
				// Disallow alignments that involve an N in the
				// reference
				const int r = (int)ref[rir + j];
				if(r & 4) {
					// N in reference; bail
					match = false;
					break;
				}
				const int q = (int)qry[j];
				assert_leq(q, 4);
				assert_lt(r, 4);
				if(q != r) {
#endif
					// Mismatch!
					if(++mms > 3) {
						// Too many; reject this alignment
						match = false;
						break;
					} else if(mms == 3) {
						// Second one; remember offset and ref char
						refc3 = "ACGTN"[r];
						mmOff3 = j;
					} else if(mms == 2) {
						// Second one; remember offset and ref char
						refc2 = "ACGTN"[r];
						mmOff2 = j;
					} else {
						assert_eq(1, mms);
						// First one; remember offset and ref char
						refc1 = "ACGTN"[r];
						mmOff1 = j;
					}
				}
			}
			if(match) {
				assert_leq(mms, 3);
				ranges.resize(ranges.size()+1);
				Range& range = ranges.back();
				range.stratum = mms;
				range.numMms = mms;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				if(mms > 0) {
					assert_lt(mmOff1, qlen);
					assert(refc1 == 'A' || refc1 == 'C' || refc1 == 'G' || refc1 == 'T');
					range.mms.push_back(mmOff1);
					range.refcs.push_back(refc1);
					if(mms > 1) {
						assert_lt(mmOff2, qlen);
						assert(refc2 == 'A' || refc2 == 'C' || refc2 == 'G' || refc2 == 'T');
						range.mms.push_back(mmOff2);
						range.refcs.push_back(refc2);
						if(mms > 2) {
							assert_eq(3, mms);
							assert_lt(mmOff3, qlen);
							assert(refc3 == 'A' || refc3 == 'C' || refc3 == 'G' || refc3 == 'T');
							range.mms.push_back(mmOff3);
							range.refcs.push_back(refc3);
						}
					}
				}
				results.push_back(ri);
			}
		}
		return;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	virtual void seed64Find(uint32_t numToFind,
					uint32_t tidx,
					uint8_t* ref,
					const TDna5Str& qry,
					const TCharStr& quals,
					uint32_t begin,
					uint32_t end,
					TRangeVec& ranges,
					TU32Vec& results,
					TSetPairs* pairs = NULL,
					uint32_t aoff = 0xffffffff,
					bool low = false) const
	{
		assert_gt(numToFind, 0);
		ASSERT_ONLY(const uint32_t rangesInitSz = ranges.size());
		ASSERT_ONLY(uint32_t duplicates = 0);
		ASSERT_ONLY(uint32_t r2i = 0);
		const uint32_t qlen = seqan::length(qry);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		// Get results from the naive matcher for sanity-checking
		TRangeVec r2; TU32Vec re2;
		naiveFind(numToFind, tidx, ref, qry, quals, begin, end, r2,
				  re2, pairs, aoff, low);
#endif
		const uint32_t seedBitPairs = min<int>(qlen, 32);
		const int lhsShift = ((seedBitPairs - 1) << 1);
		const uint32_t seedCushion  = 32 - seedBitPairs;
		const uint32_t seedOverhang = (qlen <= 32 ? 0 : (qlen - 32));
		const uint32_t lim = end - qlen - begin;
		const uint32_t halfway = begin + (lim >> 1);
		uint64_t seed = 0llu;
		uint64_t buffw = 0llu;  // rotating ref sequence buffer
		// OR the 'diff' buffer with this so that we can always count
		// 'N's as mismatches
		uint64_t diffMask = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the seed area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(seedBitPairs < 32) {
			clearMask >>= ((32-seedBitPairs) << 1);
			useMask = true;
		}
		int nsInSeed = 0;
		uint32_t nPoss = 0;
		int nPos1 = -1;
		int nPos2 = -1;
		int nPos3 = -1;
		uint32_t skipLeftToRights = 0;
		uint32_t skipRightToLefts = 0;
		// Construct the 'seed' 64-bit buffer so that it holds all of
		// the first 'seedBitPairs' bit pairs of the query.
		for(uint32_t i = 0; i < seedBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r & 4) {
				r = 0;
				// The right-to-left direction absorbs the candidate
				// alignment based at 'halfway'; so skipLeftToRights is
				// i, not i+1
				skipLeftToRights = max(skipLeftToRights, i);
				skipRightToLefts = max(skipRightToLefts, seedBitPairs - i);
			}
			assert_leq(c, 4);
			assert_lt(r, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInSeed > 3) {
					// More than two 'N's in the seed region; can't
					// possibly have a 2-mismatch hit anywhere
					assert_eq(r2.size(), ranges.size() - rangesInitSz);
					return;   // can't match if query has Ns
				} else if(nsInSeed == 3) {
					nPos3 = (int)i;
					nPoss++;
				} else if(nsInSeed == 2) {
					nPos2 = (int)i;
					nPoss++;
				} else {
					assert_eq(1, nsInSeed);
					nPos1 = (int)i;
					nPoss++;
				}
				// Make it look like an 'A' in the seed
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			seed  = ((seed  << 2llu) | c);
			buffw = ((buffw << 2llu) | r);
		}
		assert_leq(nPoss, 3);
		// Check whether read is disqualified by Ns outside of the seed
		// region
		for(uint32_t i = seedBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				if(++nsInSeed > 3) {
					assert_eq(r2.size(), ranges.size() - rangesInitSz);
					return; // can't match if query has Ns
				}
			}
		}
		uint64_t bufbw = buffw;
		// We're moving the right-hand edge of the seed along until
		// it's 'seedOverhang' chars from the end of the target region.
		// Note that we're not making a 3'/5' distinction here; if we
		// were, we might need to make the 'seedOverhang' adjustment on
		// the left end of the range rather than the right.
		bool hi = false;
		uint32_t riHi  = halfway;
		uint32_t rirHi = halfway - begin;
		uint32_t rirHiSeed = rirHi + seedBitPairs - 1;
		uint32_t riLo  = halfway + 1;
		uint32_t rirLo = halfway - begin + 1;
		for(uint32_t i = 1; i <= lim + 1; i++) {
			int r;       // new reference char
			uint64_t diff;
			assert_lt(skipLeftToRights, qlen);
			assert_leq(skipRightToLefts, qlen);
			if(hi) {
				hi = false;
				// Moving left-to-right
				riHi++;
				rirHi++;
				rirHiSeed++;
				r = (int)ref[rirHiSeed];
				if(r & 4) {
					r = 0;
					skipLeftToRights = seedBitPairs;
				}
				assert_lt(r, 4);
				// Bring in new base pair at the least significant
				// position
				buffw = ((buffw << 2llu) | r);
				if(useMask) buffw &= clearMask;
				if(skipLeftToRights > 0) {
					skipLeftToRights--;
					continue;
				}
				diff = (buffw ^ seed) | diffMask;
			} else {
				hi = true;
				// Moving right-to-left
				riLo--;
				rirLo--;
				r = (int)ref[rirLo];
				if(r & 4) {
					r = 0;
					skipRightToLefts = qlen;
				}
				assert_lt(r, 4);
				if(i >= 2) {
					bufbw >>= 2llu;
					// Bring in new base pair at the most significant
					// position
					bufbw |= ((uint64_t)r << lhsShift);
				}
				if(skipRightToLefts > 0) {
					skipRightToLefts--;
					continue;
				}
				diff = (bufbw ^ seed) | diffMask;
			}
			if((diff & 0xffff000000000000llu) &&
			   (diff & 0x0000ffff00000000llu) &&
			   (diff & 0x00000000ffff0000llu) &&
			   (diff & 0x000000000000ffffllu)) continue;
			uint32_t ri  = hi ? riLo  : riHi;
			uint32_t rir = hi ? rirLo : rirHi;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the seed
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 3) continue;
			diffs += u8toMms[(int)diff8[1]] +
					 u8toMms[(int)diff8[2]] +
					 u8toMms[(int)diff8[3]] +
					 u8toMms[(int)diff8[4]] +
					 u8toMms[(int)diff8[5]] +
					 u8toMms[(int)diff8[6]];
			uint32_t mmpos1 = 0xffffffff;
			int refc1 = -1;
			uint32_t mmpos2 = 0xffffffff;
			int refc2 = -1;
			uint32_t mmpos3 = 0xffffffff;
			int refc3 = -1;
			if(diffs > 3) {
				// Too many differences
				continue;
			} else if(nPoss > 1 && diffs == nPoss) {
				// There was one difference, but there was also one N,
				// so we already know where the difference is
				mmpos1 = nPos1;
				refc1 = "ACGT"[(int)ref[rir + nPos1]];
				if(nPoss > 1) {
					mmpos2 = nPos2;
					refc2 = "ACGT"[(int)ref[rir + nPos2]];
					if(nPoss > 2) {
						mmpos3 = nPos3;
						refc3 = "ACGT"[(int)ref[rir + nPos3]];
					}
				}
			} else if(diffs > 0) {
				// Figure out which position mismatched
				uint64_t diff2 = diff;
				mmpos1 = 31;
				if((diff & 0xffffffffllu) == 0) { diff >>= 32llu; mmpos1 -= 16; }
				assert_neq(0, diff);
				if((diff & 0xffffllu) == 0)     { diff >>= 16llu; mmpos1 -=  8; }
				assert_neq(0, diff);
				if((diff & 0xffllu) == 0)       { diff >>= 8llu;  mmpos1 -=  4; }
				assert_neq(0, diff);
				if((diff & 0xfllu) == 0)        { diff >>= 4llu;  mmpos1 -=  2; }
				assert_neq(0, diff);
				if((diff & 0x3llu) == 0)        { mmpos1--; }
				assert_neq(0, diff);
				assert_geq(mmpos1, 0);
				assert_lt(mmpos1, 32);
				mmpos1 -= seedCushion;
				refc1 = "ACGT"[(int)ref[rir + mmpos1]];
				if(diffs > 1) {
					// Figure out the second mismatched position
					diff2 &= ~(0xc000000000000000llu >> (uint64_t)((mmpos1 + seedCushion) << 1));
					uint64_t diff3 = diff2;
					mmpos2 = 31;
					if((diff2 & 0xffffffffllu) == 0) { diff2 >>= 32llu; mmpos2 -= 16; }
					assert_neq(0, diff2);
					if((diff2 & 0xffffllu) == 0)     { diff2 >>= 16llu; mmpos2 -=  8; }
					assert_neq(0, diff2);
					if((diff2 & 0xffllu) == 0)       { diff2 >>= 8llu;  mmpos2 -=  4; }
					assert_neq(0, diff2);
					if((diff2 & 0xfllu) == 0)        { diff2 >>= 4llu;  mmpos2 -=  2; }
					assert_neq(0, diff2);
					if((diff2 & 0x3llu) == 0)        { mmpos2--; }
					assert_neq(0, diff2);
					mmpos2 -= seedCushion;
					assert_geq(mmpos2, 0);
					assert_lt(mmpos2, 32);
					assert_neq(mmpos1, mmpos2);
					refc2 = "ACGT"[(int)ref[rir + mmpos2]];
					uint32_t mmpos2orig = mmpos2;
					if(mmpos2 < mmpos1) {
						uint32_t mmtmp = mmpos1;
						mmpos1 = mmpos2;
						mmpos2 = mmtmp;
						int refctmp = refc1;
						refc1 = refc2;
						refc2 = refctmp;
					}
					assert_lt(mmpos1, mmpos2);
					if(diffs > 2) {
						// Figure out the second mismatched position
						diff3 &= ~(0xc000000000000000llu >> (uint64_t)((mmpos2orig + seedCushion) << 1));
						mmpos3 = 31;
						if((diff3 & 0xffffffffllu) == 0) { diff3 >>= 32llu; mmpos3 -= 16; }
						assert_neq(0, diff3);
						if((diff3 & 0xffffllu) == 0)     { diff3 >>= 16llu; mmpos3 -=  8; }
						assert_neq(0, diff3);
						if((diff3 & 0xffllu) == 0)       { diff3 >>= 8llu;  mmpos3 -=  4; }
						assert_neq(0, diff3);
						if((diff3 & 0xfllu) == 0)        { diff3 >>= 4llu;  mmpos3 -=  2; }
						assert_neq(0, diff3);
						if((diff3 & 0x3llu) == 0)        { mmpos3--; }
						assert_neq(0, diff3);
						mmpos3 -= seedCushion;
						assert_geq(mmpos3, 0);
						assert_lt(mmpos3, 32);
						assert_neq(mmpos1, mmpos3);
						assert_neq(mmpos2, mmpos3);
						refc3 = "ACGT"[(int)ref[rir + mmpos3]];
						if(mmpos3 < mmpos1) {
							uint32_t mmtmp = mmpos1;
							mmpos1 = mmpos3;
							mmpos3 = mmpos2;
							mmpos2 = mmtmp;
							int refctmp = refc1;
							refc1 = refc3;
							refc3 = refc2;
							refc2 = refctmp;
						} else if(mmpos3 < mmpos2) {
							uint32_t mmtmp = mmpos2;
							mmpos2 = mmpos3;
							mmpos3 = mmtmp;
							int refctmp = refc2;
							refc2 = refc3;
							refc3 = refctmp;
						}
						assert_lt(mmpos1, mmpos2);
						assert_lt(mmpos2, mmpos3);
					}
				}
			}
			// Now extend the seed into a longer alignment
			bool foundHit = true;
			if(seedOverhang > 0) {
				assert_leq(ri + seedBitPairs + seedOverhang, end);
				bool skipCandidate = false;
				for(uint32_t j = 0; j < seedOverhang; j++) {
					int rc = (int)ref[rir + 32 + j];
					if(rc == 4) {
						// Oops, encountered an N in the reference in
						// the overhang portion of the candidate
						// alignment
						// (Note that we inverted hi earlier)
						if(hi) {
							// Right-to-left
							skipRightToLefts = seedOverhang - j - 1;
						} else {
							// Left-to-right
							skipLeftToRights = seedBitPairs + j;
						}
						skipCandidate = true; // Skip this candidate
						break;
					}
					if((int)qry[32 + j] != rc) {
						if(++diffs > 3) {
							foundHit = false;
							break;
						} else if(diffs == 3) {
							mmpos3 = 32 + j;
							refc3 = "ACGT"[(int)ref[rir + 32 + j]];
						} else if(diffs == 2) {
							mmpos2 = 32 + j;
							refc2 = "ACGT"[(int)ref[rir + 32 + j]];
						} else {
							assert_eq(1, diffs);
							mmpos1 = 32 + j;
							refc1 = "ACGT"[(int)ref[rir + 32 + j]];
						}
					}
				}
				if(skipCandidate) continue;
			}
			if(!foundHit) continue;
			if(pairs != NULL) {
				TU64Pair p;
				if(ri < aoff) {
					// By convention, the upstream ("low") mate's
					// coordinates go in the 'first' field
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.second = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				} else {
					p.second = ((uint64_t)tidx << 32) | (uint64_t)ri;
					p.first  = ((uint64_t)tidx << 32) | (uint64_t)aoff;
				}
				if(pairs->find(p) != pairs->end()) {
					// We already found this hit!  Continue.
					ASSERT_ONLY(duplicates++);
					ASSERT_ONLY(r2i++);
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			assert_leq(diffs, 3);
			assert_lt(r2i, r2.size());
			assert_eq(re2[r2i], ri);
			ranges.resize(ranges.size()+1);
			Range& range = ranges.back();
			assert_eq(diffs, r2[r2i].numMms);
			range.stratum = diffs;
			range.numMms = diffs;
			assert_eq(0, range.mms.size());
			assert_eq(0, range.refcs.size());
			if(diffs > 0) {
				assert_neq(mmpos1, 0xffffffff);
				assert_eq(mmpos1, r2[r2i].mms[0]);
				assert_neq(-1, refc1);
				assert_eq(refc1, r2[r2i].refcs[0]);
				range.mms.push_back(mmpos1);
				range.refcs.push_back(refc1);
				if(diffs > 1) {
					assert_neq(mmpos2, 0xffffffff);
					assert_eq(mmpos2, r2[r2i].mms[1]);
					assert_neq(-1, refc2);
					assert_eq(refc2, r2[r2i].refcs[1]);
					range.mms.push_back(mmpos2);
					range.refcs.push_back(refc2);
					if(diffs > 2) {
						assert_neq(mmpos3, 0xffffffff);
						assert_eq(mmpos3, r2[r2i].mms[2]);
						assert_neq(-1, refc3);
						assert_eq(refc3, r2[r2i].refcs[2]);
						range.mms.push_back(mmpos3);
						range.refcs.push_back(refc3);
					}
				}
			}
			ASSERT_ONLY(r2i++);
			results.push_back(ri);
			if(--numToFind == 0) return;
		}
		assert_leq(duplicates, r2.size());
		assert_geq(r2.size() - duplicates, ranges.size() - rangesInitSz);
		return; // no hit
	}
};

#endif /* REF_ALIGNER_H_ */
