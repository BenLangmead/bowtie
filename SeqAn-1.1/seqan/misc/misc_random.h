 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: misc_random.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MISC_RANDOM_H
#define SEQAN_HEADER_MISC_RANDOM_H

#include <cstdlib>
#include <ctime>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Mersenne Twister Random Number Generator
// Implementation by Michael Brundage;
// with some modifications
// http://www.qbrundage.com/michaelb/pubs/essays/random_number_generation.html
//////////////////////////////////////////////////////////////////////////////


#define SEQAN_MERSENNE_MT_LEN			624
#define SEQAN_MERSENNE_MT_IA			397
#define SEQAN_MERSENNE_MT_IB			(SEQAN_MERSENNE_MT_LEN - SEQAN_MERSENNE_MT_IA)
#define SEQAN_MERSENNE_UPPER_MASK      0x80000000
#define SEQAN_MERSENNE_LOWER_MASK      0x7FFFFFFF
#define SEQAN_MERSENNE_MATRIX_A        0x9908B0DF
#define SEQAN_MERSENNE_TWIST(b,i,j)    ((b)[i] & SEQAN_MERSENNE_UPPER_MASK) | ((b)[j] & SEQAN_MERSENNE_LOWER_MASK)
#define SEQAN_MERSENNE_MAGIC(s)        (((s)&1)*SEQAN_MERSENNE_MATRIX_A)

//////////////////////////////////////////////////////////////////////////////
//forward declaration

inline unsigned long mtRand(); 
inline void mtRandInit();

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _MersenneBuffer
{
    static unsigned long buffer[SEQAN_MERSENNE_MT_LEN];
	static int index;
	static bool is_initialized;

};
//____________________________________________________________________________

template <typename T>
unsigned long _MersenneBuffer<T>::buffer[SEQAN_MERSENNE_MT_LEN];

//____________________________________________________________________________

template <typename T>
int _MersenneBuffer<T>::index = 0;

//____________________________________________________________________________

template <typename T>
bool _MersenneBuffer<T>::is_initialized = false;

//////////////////////////////////////////////////////////////////////////////
// NOTE: mtRandInit() must have been called at least once before mtRand() is used.
// Can also be called several times since it is protected against multiple initalizations. 

inline void 
mtRandInit(bool 
#ifndef SEQAN_NOSRAN
		   _doSRand
#endif
		   )
{
	// test whether mtRandInit was already initialized
	// return immediately if this is the case
	if (_MersenneBuffer<>::is_initialized) return;
	_MersenneBuffer<>::is_initialized = true;

#ifndef SEQAN_NOSRAN
	if (_doSRand)
		::std::srand((unsigned) ::std::time(0));
#endif

	int i;
	for (i = 0; i < SEQAN_MERSENNE_MT_LEN; i++)
		_MersenneBuffer<>::buffer[i] = ::std::rand();

	mtRand(); //pop the first number, since it is not as "random" as we like it
}

inline void 
mtRandInit()
{
	mtRandInit(true);
}

//////////////////////////////////////////////////////////////////////////////
// NOTE: mtRandInit() must be called once before mtRand() is used.

inline unsigned long 
mtRand()
{
	unsigned long * b = _MersenneBuffer<>::buffer;
	int idx = _MersenneBuffer<>::index;
	unsigned long s;
	int i;
	
	if (idx == SEQAN_MERSENNE_MT_LEN*sizeof(unsigned long))
	{
		idx = 0;
		i = 0;
		for (; i < SEQAN_MERSENNE_MT_IB; i++) {
			s = SEQAN_MERSENNE_TWIST(b, i, i+1);
			b[i] = b[i + SEQAN_MERSENNE_MT_IA] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
		}
		for (; i < SEQAN_MERSENNE_MT_LEN-1; i++) {
			s = SEQAN_MERSENNE_TWIST(b, i, i+1);
			b[i] = b[i - SEQAN_MERSENNE_MT_IB] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
		}
        
		s = SEQAN_MERSENNE_TWIST(b, SEQAN_MERSENNE_MT_LEN-1, 0);
		b[SEQAN_MERSENNE_MT_LEN-1] = b[SEQAN_MERSENNE_MT_IA-1] ^ (s >> 1) ^ SEQAN_MERSENNE_MAGIC(s);
	}
	_MersenneBuffer<>::index = idx + sizeof(unsigned long);
	return *(unsigned long *)((unsigned char *)b + idx);
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Geometrical Distribution Heuristics
// from Numerical Recipes in C
//////////////////////////////////////////////////////////////////////////////


#define SEQAN_RG_IB1 1
#define SEQAN_RG_IB2 2
#define SEQAN_RG_IB5 16
#define SEQAN_RG_IB18 131072
#define SEQAN_RG_MASK ( SEQAN_RG_IB1 + SEQAN_RG_IB2 + SEQAN_RG_IB5 )

template <typename T>
inline T
geomRand()
{
	static unsigned long seed = ::std::rand();
	T value = 0;
	while ( true )
	{
		if( ( seed & SEQAN_RG_IB18 ) )
		{
			seed = ( ( seed ^ SEQAN_RG_MASK ) << 1 ) | SEQAN_RG_IB1;
			++value;
		}
		else 
		{
			seed <<= 1;
			break;
		}
	}
	return value;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
