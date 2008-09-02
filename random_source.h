#ifndef RANDOM_GEN_H_
#define RANDOM_GEN_H_

/**
 * Simple pseudo-random linear congruential generator, a la Numerical
 * Recipes.
 */
class RandomSource {
public:
	RandomSource(uint32_t seed = 0) : a(1664525), c(1013904223), last(seed) { }

	void init(uint32_t seed = 0) {
		last = seed;
	}

	uint32_t nextU32() {     // 0.97% in profile
		last = a * last + c; // 0.95% in profile
		return last;
	}

private:
	uint32_t a;
	uint32_t c;
	uint32_t last;
};

#endif /*RANDOM_GEN_H_*/
