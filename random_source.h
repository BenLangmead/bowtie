#ifndef RANDOM_GEN_H_
#define RANDOM_GEN_H_

/**
 * Simple pseudo-random linear congruential generator, a la Numerical
 * Recipes.
 */
class RandomSource {
public:
	static const uint32_t DEFUALT_A = 1664525;
	static const uint32_t DEFUALT_C = 1013904223;

	RandomSource(uint32_t seed = 0) :
		a(DEFUALT_A), c(DEFUALT_C), last(a * seed + c), lastOff(0) { }
	RandomSource(uint32_t _a, uint32_t _c, uint32_t seed = 0) :
		a(_a), c(_c), last(seed), lastOff(0) { }

	void init(uint32_t seed = 0) {
		last = seed;
	}

	uint32_t nextU32() {
		last = a * last + c;
		lastOff = 0;
		return last;
	}

	uint32_t nextU2() {
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
};

#endif /*RANDOM_GEN_H_*/
