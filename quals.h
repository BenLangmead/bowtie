#ifndef QUALS_H_
#define QUALS_H_

/// An array that transforms Phred qualities into their maq-like
/// equivalents by dividing by ten and rounding to the nearest 1,
/// but saturating at 3.
extern unsigned char qualRounds[];

#endif /*QUALS_H_*/
