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
	 * Find all new alignments and store in results vector.
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
	                  uint32_t aoff = 0,
	                  bool low = false) = 0;

	/**
	 * Set a new reference-sequence buffer.
	 */
	void setBuf(uint8_t *newbuf, uint32_t newsz) {
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
			refbuf_ = new uint8_t[newsz];
			if(refbuf_ == NULL) throw std::bad_alloc();
		} catch(std::bad_alloc& e) {
			cerr << "Error: Could not allocate reference-space alignment buffer of " << newsz << "B" << endl;
			exit(1);
		}
		refbufSz_ = newsz;
		freeRefbuf_ = true;
	}

protected:
	uint32_t seedLen_;   /// length of seed region for read
	uint8_t *refbuf_;    /// pointer to current reference buffer
	uint32_t refbufSz_;  /// size of current reference buffer
	uint8_t  buf_[REF_ALIGNER_BUFSZ]; /// built-in reference buffer (may be superseded)
	bool     freeRefbuf_; /// whether refbuf_ points to something we should delete
};

/**
 * Concrete RefAligner for finding exact hits.
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
	                  uint32_t apos = 0,
	                  bool low = false)
	{
		assert_gt(numToFind, 0);
		uint32_t spread = end - begin;
		// Make sure the buffer is large enough to accommodate the spread
		if(spread > this->refbufSz_) {
			this->newBuf(spread);
		}
		// Read in the relevant stretch of the reference string
		refs->getStretch(this->refbuf_, tidx, begin, spread);
		// Look for alignments
		return seed64Find(numToFind, tidx, this->refbuf_, qry, quals,
		                  begin, end, ranges, results, pairs, apos,
		                  low);
	}

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
	               uint32_t aoff, // offset of anchor
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
				if(r == 4) {
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
				if(pairs != NULL) {
					TU64Pair p;
					if(low) {
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
						continue;
					} else {
						// *Don't* record this hit - this is just a
						// sanity check
						//pairs->insert(p);
					}
				}
				ranges.resize(ranges.size());
				Range& range = ranges.back();
				range.stratum = 0;
				range.numMms = 0;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				results.push_back(ri);
				if(--numToFind == 0) return;
			}
		}
		return; // no match
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void seed64Find(uint32_t numToFind,
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
			if(c == 4) {
				assert_eq(r2.size(), ranges.size() - rangesInitSz);
				return; // can't match if query has Ns
			}
			int r = (int)ref[halfway - begin + i]; // next reference character
			if(r == 4) {
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
		for(uint32_t i = 1; i <= lim + 1; i++) {
			uint32_t ri; // leftmost position in proposed alignment
			uint32_t rir;// same, but minus 'begin'
			int r;       // new reference char
			assert_leq(skipLeftToRights, seedBitPairs);
			assert_leq(skipRightToLefts, seedBitPairs);
			if(hi) {
				hi = false;
				// Moving left-to-right
				ri = halfway + (i >> 1); rir = ri - begin;
				r = (int)ref[rir + seedBitPairs - 1];
				if(r == 4) {
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
				ri = halfway - (i >> 1); rir = ri - begin;
				r = (int)ref[rir];
				if(r == 4) {
					r = 0;
					skipRightToLefts = seedBitPairs;
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
			if(seedOverhang > 0) {
				// Does the non-seed part of the alignment (the
				// "overhang") ruin it?
				for(uint32_t j = 0; j < seedOverhang; j++) {
					assert_lt(ri + seedBitPairs + j, end);
					if((int)qry[32 + j] != (int)ref[rir + seedBitPairs + j]) {
						// Yes, overhang ruins it
						foundHit = false;
						break;
					}
				}
			}
			if(foundHit) {
				if(pairs != NULL) {
					TU64Pair p;
					if(low) {
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
						continue;
					} else {
						// Record this hit
						pairs->insert(p);
					}
				}
				ASSERT_ONLY(uint32_t added = ranges.size() - rangesInitSz);
				assert_lt(added, r2.size());
				assert_eq(re2[added], ri);
				ranges.resize(ranges.size()+1);
				Range& range = ranges.back();
				range.stratum = 0;
				range.numMms = 0;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				results.push_back(ri);
				if(--numToFind == 0) return;
			} else {
				// Keep scanning
			}
		}
		assert_eq(r2.size(), ranges.size() - rangesInitSz);
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
 * Concrete RefAligner for finding exact hits.
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
		// Make sure the buffer is large enough to accommodate the spread
		if(spread > this->refbufSz_) {
			this->newBuf(spread);
		}
		// Read in the relevant stretch of the reference string
		refs->getStretch(this->refbuf_, tidx, begin, spread);
		// Look for alignments
		return seed64Find(numToFind, tidx, this->refbuf_, qry, quals,
		                  begin, end, ranges, results, pairs, aoff, low);
	}

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
				if(r == 4) {
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
				if(pairs != NULL) {
					TU64Pair p;
					if(low) {
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
						continue;
					} else {
						// *Don't* record this hit - this is just a
						// sanity check
						// pairs->insert(p);
					}
				}
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
				if(--numToFind == 0) return;
			}
		}
		return;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	void seed64Find(uint32_t numToFind,
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
			if(r == 4) {
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
		for(uint32_t i = 1; i <= lim + 1; i++) {
			uint32_t ri; // leftmost position in proposed alignment
			uint32_t rir;// same, but minus 'begin'
			int r;       // new reference char
			uint64_t diff;
			assert_leq(skipLeftToRights, seedBitPairs);
			assert_leq(skipRightToLefts, seedBitPairs);
			if(hi) {
				hi = false;
				// Moving left-to-right
				ri = halfway + (i >> 1); rir = ri - begin;
				r = (int)ref[rir + seedBitPairs - 1];
				if(r == 4) {
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
				ri = halfway - (i >> 1); rir = ri - begin;
				r = (int)ref[rir];
				if(r == 4) {
					r = 0;
					skipRightToLefts = seedBitPairs;
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
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			uint32_t mmpos = 0xffffffff;
			int refc = -1;
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the seed
			uint32_t diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 1) {
				continue;
			}
			diffs += u8toMms[(int)diff8[1]] +
			         u8toMms[(int)diff8[2]] +
			         u8toMms[(int)diff8[3]] +
			         u8toMms[(int)diff8[4]] +
			         u8toMms[(int)diff8[5]] +
			         u8toMms[(int)diff8[6]];
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
				for(uint32_t j = 0; j < seedOverhang; j++) {
					if((int)qry[32 + j] != (int)ref[rir + 32 + j]) {
						if(++diffs > 1) {
							foundHit = false;
							break;
						} else {
							mmpos = 32 + j;
							refc = "ACGT"[(int)ref[rir + 32 + j]];
						}
					}
				}
			}
			if(!foundHit) continue;
			if(pairs != NULL) {
				TU64Pair p;
				if(low) {
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
					continue;
				} else {
					// Record this hit
					pairs->insert(p);
				}
			}
			assert(diffs <= 1);
			ASSERT_ONLY(uint32_t added = ranges.size() - rangesInitSz);
			assert_lt(added, r2.size());
			assert_eq(re2[added], ri);
			ranges.resize(ranges.size()+1);
			Range& range = ranges.back();
			assert_eq(diffs, r2[added].numMms);
			range.stratum = diffs;
			range.numMms = diffs;
			assert_eq(0, range.mms.size());
			assert_eq(0, range.refcs.size());
			if(diffs == 1) {
				assert_neq(mmpos, 0xffffffff);
				assert_eq(mmpos, r2[added].mms[0]);
				assert_neq(-1, refc);
				assert_eq(refc, r2[added].refcs[0]);
				range.mms.push_back(mmpos);
				range.refcs.push_back(refc);
			}
			results.push_back(ri);
			if(--numToFind == 0) return;
		}
		assert_eq(r2.size(), ranges.size() - rangesInitSz);
		return; // no hit
	}
};

#endif /* REF_ALIGNER_H_ */
