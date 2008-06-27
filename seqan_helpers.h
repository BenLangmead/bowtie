#ifndef SEQAN_HELPERS_H_
#define SEQAN_HELPERS_H_

#include <iostream>
#include <seqan/sequence.h>

#ifdef NDEBUG
#define appendNE(t, s) { \
	seqan::append(t, s); \
}

#define appendValueNE(t, s) { \
	seqan::appendValue(t, s); \
}
#else
#define appendNE(t, s) { \
	size_t cap = seqan::capacity(t); \
	seqan::append(t, s); \
	assert_eq(seqan::capacity(t), cap); \
}

#define appendValueNE(t, s) { \
	size_t cap = seqan::capacity(t); \
	seqan::appendValue(t, s); \
	assert_eq(seqan::capacity(t), cap); \
}
#endif

#endif /*SEQAN_HELPERS_H_*/
