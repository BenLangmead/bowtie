#ifndef QUAL_H_
#define QUAL_H_

extern unsigned char qualRounds[];

/// Translate a Phred-encoded ASCII character into a Phred quality
static inline uint8_t phredCharToPhredQual(char c) {
	return ((uint8_t)c >= 33 ? ((uint8_t)c - 33) : 0);
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

#endif /*QUAL_H_*/
