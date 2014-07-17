

#ifndef BOWTIE_INDEX_TYPES_H
#define	BOWTIE_INDEX_TYPES_H

#ifdef BOWTIE_64BIT_INDEX
#define OFF_MASK 0xffffffffffffffff
#define OFF_LEN_MASK 0xc000000000000000
#define LS_SIZE 0x100000000000000
#define OFF_SIZE 8
#define CACHE_WRAPPER_BIT 0x8000000000000000

typedef uint64_t TIndexOffU;
typedef int64_t TIndexOff;

#else
#define OFF_MASK 0xffffffff
#define OFF_LEN_MASK 0xc0000000
#define LS_SIZE 0x10000000
#define OFF_SIZE 4
#define CACHE_WRAPPER_BIT 0x80000000

typedef uint32_t TIndexOffU;
typedef int TIndexOff;

#endif /* BOWTIE_64BIT_INDEX */

extern const std::string gEbwt_ext;

#endif	/* BOWTIE_INDEX_TYPES_H */

