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
				if(qry[j] != ref[i+j]) {
					match = false;
					break;
				}
	 		}
			if(match) return i;
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
		const uint32_t rlen = seqan::length(ref);
		assert_geq(end - begin, qlen); // caller should have checked this
		assert_leq(end, rlen);
		assert_gt(end, begin);
		assert_gt(qlen, 0);
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
#ifndef NDEBUG
				Range r2;
				assert_eq(0xffffffff, naiveFindOne(ref, qry, quals, begin, end, r2));
#endif
				return 0xffffffff;   // can't match if query has Ns
			}
			assert_lt(c, 4);
			seed = ((seed << 2llu) | c);
			buf  = ((buf  << 2llu) | r);
		}
		buf >>= 1;
		if(seed == buf)
		for(uint32_t i = seedBitPairs; i < qlen; i++) {
			if((int)qry[i] == 4) {
#ifndef NDEBUG
				Range r2;
				assert_eq(0xffffffff, naiveFindOne(ref, qry, quals, begin, end, r2));
#endif
				return 0xffffffff; // can't match if query has Ns
			}
		}
		const uint32_t lim = end - seedCushion - qlen;
		for(uint32_t i = begin + (seedBitPairs-1); i < lim; i++) {
			const int r = (int)ref[i];
			buf = ((buf << 2llu) | r);
			if(useMask) buf &= clearMask;
			if(buf != seed) continue;
			// Seed hit!
			if(seedCushion == 0) {
#ifndef NDEBUG
				Range r2;
				assert_eq(i, naiveFindOne(ref, qry, quals, begin, end, r2));
#endif
				return i;
			}
			bool foundHit = true;
			for(uint32_t j = 0; j < seedCushion; j++) {
				if((int)qry[32+j] != (int)ref[i+1+j]) {
					foundHit = false;
					break;
				}
			}
			if(foundHit) {
#ifndef NDEBUG
				Range r2;
				assert_eq(i, naiveFindOne(ref, qry, quals, begin, end, r2));
#endif
				return i;
			}
		}
#ifndef NDEBUG
		Range r2;
		assert_eq(0xffffffff, naiveFindOne(ref, qry, quals, begin, end, r2));
#endif
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
