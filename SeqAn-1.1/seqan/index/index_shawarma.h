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
  $Id: index_shawarma.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_SHAWARMA_H
#define SEQAN_HEADER_INDEX_SHAWARMA_H


namespace SEQAN_NAMESPACE_MAIN {

	namespace shawarma {

		extern void ds_ssort(unsigned char *t, int *sa, int n);
		extern int init_ds_ssort(int adist, int bs_ratio);

	}

//////////////////////////////////////////////////////////////////////////////
// SeqAn interface

	template <typename TSpec>
	struct Shawarma {};

	struct MSufSort{};			// MSufSort
	struct DivSufSort{};		// DivSufSort
	struct DeepShallow{};		// Deep-Shallow sort
	struct QSufSort{};			// QSufSort

	// WARNING:
	// 1. text value must have char size
	// 2. SA value must have int size
	// 3. Deep-Shallow sort expects overshoot bytes behind the text
	// 4. SA must be contiguous
	//

    template < typename TSA,
               typename TText >
    void createSuffixArray(
		TSA &SA,
		TText &s,
		Shawarma<DeepShallow> const)
	{
		typedef typename Value<TText>::Type	TValue;
		typedef typename Value<TSA>::Type	TSAValue;

		SEQAN_ASSERT(sizeof(TValue) == sizeof(unsigned char));
		SEQAN_ASSERT(sizeof(TSAValue) == sizeof(int));
		SEQAN_ASSERT(IsContiguous<TSA>::VALUE);

		int overshoot = shawarma::init_ds_ssort(500, 2000);

		SEQAN_ASSERT(overshoot > 0);
		reserve(s, length(s) + overshoot);
		shawarma::ds_ssort(
			(unsigned char*)toCString(s),		// text
			(int*)begin(SA, Standard()),		// SA
			length(s));							// n
	}


} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_H
