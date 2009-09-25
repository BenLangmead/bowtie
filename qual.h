#ifndef QUAL_H_
#define QUAL_H_

#include <stdexcept>

extern unsigned char qualRounds[];
extern unsigned char solToPhred[];

/// Translate a Phred-encoded ASCII character into a Phred quality
static inline uint8_t phredCharToPhredQual(char c) {
	return ((uint8_t)c >= 33 ? ((uint8_t)c - 33) : 0);
}

/**
 * Convert a Solexa-scaled quality value into a Phred-scale quality
 * value.
 *
 * p = probability that base is miscalled
 * Qphred = -10 * log10 (p)
 * Qsolexa = -10 * log10 (p / (1 - p))
 * See: http://en.wikipedia.org/wiki/FASTQ_format
 *
 */
static inline uint8_t solexaToPhred(int sol) {
	assert_lt(sol, 256);
	if(sol < -10) return 0;
	return solToPhred[sol+10];
}

class SimplePhredPenalty {
public:
	static uint8_t mmPenalty (uint8_t qual) {
		return qual;
	}
	static uint8_t delPenalty(uint8_t qual) {
		return qual;
	}
	static uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right) {
		return std::max(qual_left, qual_right);
	}
};

class MaqPhredPenalty {
public:
	static uint8_t mmPenalty (uint8_t qual) {
		return qualRounds[qual];
	}
	static uint8_t delPenalty(uint8_t qual) {
		return qualRounds[qual];
	}
	static uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right) {
		return qualRounds[std::max(qual_left, qual_right)];
	}
};

static inline uint8_t mmPenalty(bool maq, uint8_t qual) {
	if(maq) {
		return MaqPhredPenalty::mmPenalty(qual);
	} else {
		return SimplePhredPenalty::mmPenalty(qual);
	}
}

static inline uint8_t delPenalty(bool maq, uint8_t qual) {
	if(maq) {
		return MaqPhredPenalty::delPenalty(qual);
	} else {
		return SimplePhredPenalty::delPenalty(qual);
	}
}

static inline uint8_t insPenalty(bool maq, uint8_t qual_left, uint8_t qual_right) {
	if(maq) {
		return MaqPhredPenalty::insPenalty(qual_left, qual_right);
	} else {
		return SimplePhredPenalty::insPenalty(qual_left, qual_right);
	}
}

/**
 * Take an ASCII-encoded quality value and convert it to a Phred33
 * ASCII char.
 */
inline static char charToPhred33(char c, bool solQuals, bool phred64Quals) {
	if(c == ' ') {
		cerr << "Saw a space but expected an ASCII-encoded quality value." << endl
		     << "Are quality values formatted as integers?  If so, try --integer-quals." << endl;
		throw 1;
	}
	if (solQuals) {
		// Convert solexa-scaled chars to phred
		// http://maq.sourceforge.net/fastq.shtml
		char cc = solexaToPhred((int)c - 64) + 33;
		if (cc < 33) {
			cerr << "Saw ASCII character "
			     << ((int)c)
			     << " but expected 64-based Solexa qual (converts to " << (int)cc << ")." << endl
			     << "Try not specifying --solexa-quals." << endl;
			throw 1;
		}
		c = cc;
	}
	else if(phred64Quals) {
		if (c < 64) {
			cerr << "Saw ASCII character "
			     << ((int)c)
			     << " but expected 64-based Phred qual." << endl
			     << "Try not specifying --solexa1.3-quals/--phred64-quals." << endl;
			throw 1;
		}
		// Convert to 33-based phred
		c -= (64-33);
	}
	else {
		// Keep the phred quality
		if (c < 33) {
			cerr << "Saw ASCII character "
			     << ((int)c)
			     << " but expected 33-based Phred qual." << endl;
			throw 1;
		}
	}
	return c;
}

/**
 * Take an integer quality value and convert it to a Phred33 ASCII
 * char.
 */
inline static char intToPhred33(int iQ, bool solQuals) {
	int pQ;
	if (solQuals) {
		// Convert from solexa quality to phred
		// quality and translate to ASCII
		// http://maq.sourceforge.net/qual.shtml
		pQ = solexaToPhred((int)iQ) + 33;
	} else {
		// Keep the phred quality and translate
		// to ASCII
		pQ = (iQ <= 93 ? iQ : 93) + 33;
	}
	if (pQ < 33) {
		cerr << "Saw negative Phred quality " << ((int)pQ-33) << "." << endl;
		throw 1;
	}
	assert_geq(pQ, 0);
	return (int)pQ;
}

/**
 * Fill the q[] array with the penalties that are determined by
 * subtracting the quality values of the alternate basecalls from
 * the quality of the primary basecall.
 */
inline static uint8_t penaltiesAt(size_t off, uint8_t *q,
                                  int alts,
                                  const String<char>& qual,
                                  const String<Dna5>* altQry,
                                  const String<char>* altQual)
{
	uint8_t primQ = qual[off]; // qual of primary call
	uint8_t bestPenalty = primQ - 33;
	q[0] = q[1] = q[2] = q[3] = bestPenalty;
	for(int i = 0; i < alts; i++) {
		uint8_t altQ = altQual[i][off]; // qual of alt call
		if(altQ == 33) break; // no alt call
		assert_leq(altQ, primQ);
		if(primQ - altQ < bestPenalty) {
			bestPenalty = primQ - altQ;
		}
		// Get the base
		int altC = (int)altQry[i][off];
		assert_lt(altC, 4);
		q[altC] = primQ - altQ;
	}
	return bestPenalty;
}

/**
 * Fill the q[] array with the penalties that are determined by
 * subtracting the quality values of the alternate basecalls from
 * the quality of the primary basecall.
 */
inline static uint8_t loPenaltyAt(size_t off, int alts,
                                  const String<char>& qual,
                                  const String<char>* altQual)
{
	uint8_t primQ = qual[off]; // qual of primary call
	uint8_t bestPenalty = primQ - 33;
	for(int i = 0; i < alts; i++) {
		uint8_t altQ = altQual[i][off]; // qual of alt call
		if(altQ == 33) break; // no more alt calls at this position
		assert_leq(altQ, primQ);
		if(primQ - altQ < bestPenalty) {
			bestPenalty = primQ - altQ;
		}
	}
	return bestPenalty;
}

#endif /*QUAL_H_*/
