#ifndef QUAL_H_
#define QUAL_H_

extern unsigned char qualRounds[];

class Penalty {
public:
	virtual ~Penalty() { }
	virtual uint8_t mmPenalty (uint8_t qual) const = 0;
	virtual uint8_t delPenalty(uint8_t qual) const = 0;
	virtual uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right) const = 0;
};

class SimplePhredPenalty : public Penalty {
public:
	virtual ~SimplePhredPenalty() { }
	virtual uint8_t mmPenalty (uint8_t qual) const {
		return qual;
	}
	virtual uint8_t delPenalty(uint8_t qual) const {
		return qual;
	}
	virtual uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right) const {
		return std::max(qual_left, qual_right);
	}
};

class MaqPhredPenalty : public Penalty {
public:
	virtual ~MaqPhredPenalty() { }
	virtual uint8_t mmPenalty (uint8_t qual) const {
		return qualRounds[qual];
	}
	virtual uint8_t delPenalty(uint8_t qual) const {
		return qualRounds[qual];
	}
	virtual uint8_t insPenalty(uint8_t qual_left, uint8_t qual_right) const {
		return qualRounds[std::max(qual_left, qual_right)];
	}
};

#endif /*QUAL_H_*/
