#ifndef WORD_IO_H_
#define WORD_IO_H_

#include <stdint.h>
#include <iostream>

// Write
extern void     writeU32(std::ostream& out, uint32_t x, bool toBigEndian);
extern void     writeU32(std::ostream& out, uint32_t x);
extern void     writeI32(std::ostream& out, int32_t x, bool toBigEndian);
extern void     writeI32(std::ostream& out, int32_t x);

// Read
extern uint32_t readU32 (std::istream& in, bool toBigEndian);
extern int32_t  readI32 (std::istream& in, bool toBigEndian);

#endif /*WORD_IO_H_*/
