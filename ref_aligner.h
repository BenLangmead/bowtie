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
public:
	RefAligner(uint32_t seedLen) : seedLen_(seedLen) { }
	virtual ~RefAligner() { }
	virtual uint32_t findOne(const TStr& ref,
	                         const TDna5Str& qry,
	                         const TCharStr& quals,
	                         uint32_t begin,
	                         uint32_t end,
	                         Range& range) const = 0;
	virtual void findAll    (const TStr& ref,
	                         const TDna5Str& qry,
	                         const TCharStr& quals,
	                         uint32_t begin,
	                         uint32_t end,
	                         std::vector<Range>& results) const = 0;
protected:
	uint32_t seedLen_;
};

/**
 * Abstract factory parent class for RefAligners.
 */
template<typename TStr>
class RefAlignerFactory {
public:
	RefAlignerFactory(uint32_t seedLen) : seedLen_(seedLen) { }
	virtual ~RefAlignerFactory() { }
	virtual RefAligner<TStr>* create() const = 0;
protected:
	uint32_t seedLen_;
};

/**
 * Concrete RefAligner for finding exact hits.
 */
template<typename TStr>
class ExactRefAligner : public RefAligner<TStr> {
	typedef seqan::String<seqan::Dna5> TDna5Str;
	typedef seqan::String<char> TCharStr;
public:
	ExactRefAligner(uint32_t seedLen) : RefAligner<TStr>(seedLen) { }
	virtual ~ExactRefAligner() { }


	/**
	 *
	 */
	virtual uint32_t findOne(const TStr& ref,
	                         const TDna5Str& qry,
	                         const TCharStr& quals,
	                         uint32_t begin,
	                         uint32_t end,
	                         Range& range) const
	{
		return seed64FindOne(ref, qry, quals, begin, end, range);
	}

	/**
	 *
	 */
	virtual void findAll(const TStr& ref,
	                     const TDna5Str& qry,
	                     const TCharStr& quals,
	                     uint32_t begin,
	                     uint32_t end,
	                     std::vector<Range>& results) const
	{
		throw;
	}

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	uint32_t naiveFindOne( const TStr& ref,
                           const TDna5Str& qry,
                           const TCharStr& quals,
                           uint32_t begin,
                           uint32_t end,
                           Range& range) const
	{
		const uint32_t qlen = seqan::length(qry);
		const uint32_t rlen = seqan::length(ref);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_leq(end, rlen);
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		for(uint32_t i = begin; i <= end - qlen; i++) {
			bool match = true;
			for(uint32_t j = 0; j < qlen; j++) {
				assert_lt((int)ref[i+j], 4);
				if(qry[j] != ref[i+j]) {
					match = false;
					break;
				}
	 		}
			if(match) {
				range.stratum = 0;
				range.numMms = 0;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				return i;
			}
		}
		return 0xffffffff;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	uint32_t seed64FindOne(const TStr& ref,
                           const TDna5Str& qry,
                           const TCharStr& quals,
                           uint32_t begin,
                           uint32_t end,
                           Range& range) const
	{
		const uint32_t qlen = seqan::length(qry);
		ASSERT_ONLY(const uint32_t rlen = seqan::length(ref));
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_leq(end, rlen);
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		uint32_t naive;
		{
			Range r2;
			naive = naiveFindOne(ref, qry, quals, begin, end, r2);
		}
#endif
		uint32_t seedBitPairs = min<int>(qlen, 32);
		uint32_t seedCushion = qlen <= 32 ? 0 : qlen - 32;
		uint64_t seed = 0llu;
		uint64_t buf = 0llu;
		// Set up a mask that we'll apply to the two bufs every round
		// to discard bits that were rotated out of the seed area
		uint64_t clearMask = 0xffffffffffffffffllu;
		bool useMask = false;
		if(seedBitPairs < 32) {
			clearMask >>= ((32-seedBitPairs) << 1);
			useMask = true;
		}
		for(uint32_t i = 0; i < seedBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[begin + i]; // next reference character
			if(c == 4) {
				assert_eq(0xffffffff, naive);
				return 0xffffffff;   // can't match if query has Ns
			}
			assert_lt(c, 4);
			seed = ((seed << 2llu) | c);
			buf  = ((buf  << 2llu) | r);
		}
		buf >>= 2;
		for(uint32_t i = seedBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				assert_eq(0xffffffff, naive);
				return 0xffffffff; // can't match if query has Ns
			}
		}
		const uint32_t lim = end - seedCushion;
		// We're moving the right-hand edge of the seed along
		for(uint32_t i = begin + (seedBitPairs-1); i < lim; i++) {
			const int r = (int)ref[i];
			assert_lt(r, 4);
			buf = ((buf << 2llu) | r);
			if(useMask) buf &= clearMask;
			if(buf != seed) {
				assert_neq(i - seedBitPairs + 1, naive);
				continue;
			}
			// Seed hit!
			if(seedCushion == 0) {
				uint32_t ret = i - seedBitPairs + 1;
				assert_eq(ret, naive);
				range.stratum = 0;
				range.numMms = 0;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				return ret;
			}
			bool foundHit = true;
			for(uint32_t j = 0; j < seedCushion; j++) {
				assert_lt(i+1+j, end);
				assert_lt(i+1+j, rlen);
				if((int)qry[32+j] != (int)ref[i+1+j]) {
					foundHit = false;
					break;
				}
			}
			if(foundHit) {
				uint32_t ret = i - seedBitPairs + 1;
				assert_eq(ret, naive);
				range.stratum = 0;
				range.numMms = 0;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				return ret;
			}
		}
		assert_eq(0xffffffff, naive);
		return 0xffffffff; // no hit
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
public:
	OneMMRefAligner(uint32_t seedLen) : RefAligner<TStr>(seedLen) { }
	virtual ~OneMMRefAligner() { }


	/**
	 *
	 */
	virtual uint32_t findOne(const TStr& ref,
	                         const TDna5Str& qry,
	                         const TCharStr& quals,
	                         uint32_t begin,
	                         uint32_t end,
	                         Range& range) const
	{
		return seed64FindOne(ref, qry, quals, begin, end, range);
	}

	/**
	 *
	 */
	virtual void findAll(const TStr& ref,
	                     const TDna5Str& qry,
	                     const TCharStr& quals,
	                     uint32_t begin,
	                     uint32_t end,
	                     std::vector<Range>& results) const
	{
		throw;
	}

protected:
	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	uint32_t naiveFindOne( const TStr& ref,
                           const TDna5Str& qry,
                           const TCharStr& quals,
                           uint32_t begin,
                           uint32_t end,
                           Range& range) const
	{
		const uint32_t qlen = seqan::length(qry);
		const uint32_t rlen = seqan::length(ref);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_leq(end, rlen);
		assert_gt(end, begin);
		assert_gt(qlen, 0);
		for(uint32_t i = begin; i <= end - qlen; i++) {
			bool match = true;
			int mms = 0;
			int mmOff = -1;
			int refc = -1;
			for(uint32_t j = 0; j < qlen; j++) {
				assert_lt((int)ref[i+j], 4);
				if(qry[j] != ref[i+j]) {
					if(++mms > 1) {
						match = false;
						break;
					} else {
						refc = "ACGT"[(int)ref[i+j]];
						mmOff = j;
					}
				}
	 		}
			if(match) {
				assert_leq(mms, 1);
				range.stratum = mms;
				range.numMms = mms;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				if(mms == 1) {
					range.mms.push_back(mmOff);
					range.refcs.push_back(refc);
				}
				return i;
			}
		}
		return 0xffffffff;
	}

	/**
	 * Because we're doing end-to-end exact, we don't care which end of
	 * 'qry' is the 5' end.
	 */
	uint32_t seed64FindOne(const TStr& ref,
                           const TDna5Str& qry,
                           const TCharStr& quals,
                           uint32_t begin,
                           uint32_t end,
                           Range& range) const
	{
		const uint32_t qlen = seqan::length(qry);
		ASSERT_ONLY(const uint32_t rlen = seqan::length(ref));
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_leq(end, rlen);
		assert_gt(end, begin);
		assert_gt(qlen, 0);
#ifndef NDEBUG
		uint32_t naive;
		int naiveNumMms;
		int naiveMmsOff = -1;
		int naiveRefc = -1;
		{
			Range r2;
			naive = naiveFindOne(ref, qry, quals, begin, end, r2);
			naiveNumMms = r2.numMms;
			assert_leq(naiveNumMms, 1);
			if(naiveNumMms == 1) {
				naiveMmsOff = r2.mms[0];
				naiveRefc = r2.refcs[0];
			}
		}
#endif
		uint32_t seedBitPairs = min<int>(qlen, 32);
		uint32_t seedCushion = qlen <= 32 ? 0 : qlen - 32;
		uint64_t seed = 0llu;
		uint64_t buf = 0llu;
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
		// Construct the 'seed' 64-bit buffer so that it holds all of
		// the first 'seedBitPairs' bit pairs of the query.
		for(uint32_t i = 0; i < seedBitPairs; i++) {
			int c = (int)qry[i]; // next query character
			int r = (int)ref[begin + i]; // next reference character
			assert_leq(c, 4);
			// Special case: query has an 'N'
			if(c == 4) {
				if(++nsInSeed > 1) {
					// More than one 'N' in the seed region; can't
					// possibly have a 1-mismatch hit anywhere
					assert_eq(0xffffffff, naive);
					return 0xffffffff;   // can't match if query has Ns
				}
				nPos = (int)i;
				// Make it look like an 'A' in the seed
				c = 0;
				diffMask = (diffMask << 2llu) | 1llu;
			} else {
				diffMask <<= 2llu;
			}
			seed = ((seed << 2llu) | c);
			buf  = ((buf  << 2llu) | r);
		}
		buf >>= 2;
		for(uint32_t i = seedBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
				if(++nsInSeed > 1) {
					assert_eq(0xffffffff, naive);
					return 0xffffffff; // can't match if query has Ns
				}
			}
		}
		const uint32_t lim = end - seedCushion;
		// We're moving the right-hand edge of the seed along
		for(uint32_t i = begin + (seedBitPairs-1); i < lim; i++) {
			const int r = (int)ref[i];
			assert_lt(r, 4);
			buf = ((buf << 2llu) | r);
			if(useMask) buf &= clearMask;
			uint64_t diff = (buf ^ seed) | diffMask;
			// Could use pop count
			uint8_t *diff8 = reinterpret_cast<uint8_t*>(&diff);
			int mmpos = -1;
			int refc = -1;
			// As a first cut, see if there are too many mismatches in
			// the first and last parts of the seed
			int diffs = u8toMms[(int)diff8[0]] + u8toMms[(int)diff8[7]];
			if(diffs > 1) {
				assert_neq(i - seedBitPairs + 1, naive);
				continue;
			}
			diffs += u8toMms[(int)diff8[1]] +
			         u8toMms[(int)diff8[2]] +
			         u8toMms[(int)diff8[3]] +
			         u8toMms[(int)diff8[4]] +
			         u8toMms[(int)diff8[5]] +
			         u8toMms[(int)diff8[6]];
			if(diffs > 1) {
				assert_neq(i - seedBitPairs + 1, naive);
				continue;
			} else if(diffs == 1 && nPos != -1) {
				mmpos = nPos;
				refc = "ACGT"[(int)ref[i - seedBitPairs + nPos + 1]];
#ifndef NDEBUG
				if(naiveNumMms == 1) {
					assert_eq(mmpos, naiveMmsOff);
					assert_eq(refc, naiveRefc);
				}
#endif
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
				refc = "ACGT"[(int)ref[i - seedBitPairs + mmpos + 1]];
#ifndef NDEBUG
				if(naiveNumMms == 1) {
					assert_eq(mmpos, naiveMmsOff);
					assert_eq(refc, naiveRefc);
				}
#endif
			}
			// There's nothing else to match beyond the seed, so this
			// is a seed hit!
			if(seedCushion == 0) {
				uint32_t ret = i - seedBitPairs + 1;
				assert(diffs <= 1);
				assert_eq(ret, naive);
				assert_eq(diffs, naiveNumMms);
				range.stratum = diffs;
				range.numMms = diffs;
				assert_eq(0, range.mms.size());
				assert_eq(0, range.refcs.size());
				if(diffs == 1) {
					assert_geq(mmpos, 0);
					assert_eq(mmpos, naiveMmsOff);
					assert_neq(-1, refc);
					assert_eq(refc, naiveRefc);
					range.mms.push_back(mmpos);
					range.refcs.push_back(refc);
				}
				return ret;
			}
			// Now extend the seed into a longer alignment
			bool foundHit = true;
			for(uint32_t j = 0; j < seedCushion; j++) {
				assert_lt(i+1+j, end);
				assert_lt(i+1+j, rlen);
				if((int)qry[32+j] != (int)ref[i+1+j]) {
					if(++diffs > 1) {
						foundHit = false;
						break;
					} else {
						mmpos = 32+j;
						refc = "ACGT"[(int)ref[i+1+j]];
					}
				}
			}
			if(!foundHit) continue;
			uint32_t ret = i - seedBitPairs + 1;
			assert(diffs <= 1);
			assert_eq(ret, naive);
			assert_eq(diffs, naiveNumMms);
			range.stratum = diffs;
			range.numMms = diffs;
			assert_eq(0, range.mms.size());
			assert_eq(0, range.refcs.size());
			if(diffs == 1) {
				assert_geq(mmpos, 0);
				assert_eq(mmpos, naiveMmsOff);
				assert_neq(-1, refc);
				assert_eq(refc, naiveRefc);
				range.mms.push_back(mmpos);
				range.refcs.push_back(refc);
			}
			return ret;
		}
		assert_eq(0xffffffff, naive);
		return 0xffffffff; // no hit
	}
};

/**
 * Concrete factory for ExactRefAligners.
 */
template<typename TStr>
class ExactRefAlignerFactory : public RefAlignerFactory<TStr> {
public:
	ExactRefAlignerFactory(uint32_t seedLen) :
		ExactRefAligner<TStr>(seedLen) { }
	virtual RefAligner<TStr>* create() const {
		return new ExactRefAligner<TStr>(this->seedLen_);
	}
};

#endif /* REF_ALIGNER_H_ */
