#ifndef ENDIAN_H
#define ENDIAN_H

bool currentlyBigEndian();
uint32_t endianSwapU32(uint32_t u);
uint64_t endianSwapU64(uint64_t u);
int32_t  endianSwapI32(int32_t i);
uint32_t endianizeU32(uint32_t u, bool toBig);
int32_t  endianizeI32(int32_t i, bool toBig);

#endif
