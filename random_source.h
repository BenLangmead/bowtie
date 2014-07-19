#ifndef RANDOM_GEN_H_
#define RANDOM_GEN_H_

#include "assert_helpers.h"

/**
 * Simple pseudo-random linear congruential generator, a la Numerical
 * Recipes.
 */
class RandomSource {
public:
	static const uint32_t DEFUALT_A = 1664525;
	static const uint32_t DEFUALT_C = 1013904223;

	RandomSource() :
		a(DEFUALT_A), c(DEFUALT_C), inited_(false) { }
	RandomSource(uint32_t _a, uint32_t _c) :
		a(_a), c(_c), inited_(false) { }

	void init(uint32_t seed = 0) {
		last = seed;
		inited_ = true;
	}

	uint32_t nextU32() {
		assert(inited_);
		uint32_t ret;
		last = a * last + c;
		ret = last >> 16;
		last = a * last + c;
		ret ^= last;
		lastOff = 0;
		return ret;
	}

    uint64_t nextU64() {
		assert(inited_);
		uint64_t first = nextU32();
		first = first << 32;
		uint64_t second = nextU32();
		return first | second;
	}


	size_t nextSizeT() {
		if(sizeof(size_t) == 4) {
			return nextU32();
		} else {
			return nextU64();
		}
	}

	template <typename T>
	T nextU() {
		if(sizeof(T)>4)
			return nextU64();
		return nextU32();
	}

	uint32_t nextU2() {
		assert(inited_);
		if(lastOff > 30) {
			nextU32();
		}
		uint32_t ret = (last >> lastOff) & 3;
		lastOff += 2;
		return ret;
	}

	static uint32_t nextU32(uint32_t last,
	                        uint32_t a = DEFUALT_A,
	                        uint32_t c = DEFUALT_C)
	{
		return (a * last) + c;
	}

private:
	uint32_t a;
	uint32_t c;
	uint32_t last;
	uint32_t lastOff;
	bool inited_;
};

#endif /*RANDOM_GEN_H_*/
