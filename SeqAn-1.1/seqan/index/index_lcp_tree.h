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
  $Id: index_lcp_tree.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_LCP_TREE_H
#define SEQAN_HEADER_INDEX_LCP_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	template <
		class LCPFwdIt,		// lcp table input iterator
		class FlatOutIt >	// flat tree output iterator
	inline FlatOutIt createLCPBinTree(
		LCPFwdIt _First, LCPFwdIt _Last,
		FlatOutIt _Dest)
	{
        typedef typename Value<LCPFwdIt>::Type  TValue;
        typedef typename Size<LCPFwdIt>::Type   TSize;

        TSize size = difference(_First, _Last);
        if (size <= 1) return _Dest;
		--size;

        // calculate the depth of the lcp tree
		unsigned treeLevels = 1;
		TSize _xSize = 1;
		for(; size > _xSize; _xSize *= 2, ++treeLevels);
	
		// get output iterators for every level in the flat tree
		FlatOutIt *level = new FlatOutIt[treeLevels];
		for(unsigned i = treeLevels - 1; _xSize; --i, _xSize /= 2) {
			level[i] = _Dest;
			goFurther(_Dest, (size + _xSize - 1) / _xSize);
		}

		// fields to keep track of minimum elements and state
		TValue *min = new TValue[treeLevels];
		bool *half = new bool[treeLevels];
		for(unsigned i = 0; i < treeLevels; ++i)
			half[i] = false;

		// it works like a binary counter of half[n]...half[1]half[0]
		for(TSize j = 0; j < size; ++j, ++_First) {
			*(level[0]) = min[0] = *_First;
			++(level[0]);
			for(unsigned i = 1; i < treeLevels; ++i) {
				if (half[i]) {
					if (min[i-1] < min[i]) min[i] = min[i-1];
					*(level[i]) = min[i];	// min[i] is the minimum of last 2 values in min[i-1]
					++(level[i]);
					half[i] = false;
				} else {
					min[i] = min[i-1];
					half[i] = true;
					break;
				}
			}
		}

		// complete half filled nodes
		bool carry = false;
		for(unsigned i = 1; i < treeLevels; ++i)
			if (half[i] || carry) {
				if (half[i]) {
					if (min[i-1] < min[i]) min[i] = min[i-1];
				} else
					min[i] = min[i-1];
				*(level[i]) = min[i];
				++(level[i]);
				carry = true;
			}

		// trailing zero
		*_Dest = 0;
		++_Dest;

		delete[] half;
		delete[] min;
		delete[] level;

		return _Dest;
    }



	template < typename TSize >
	inline TSize sizeofLCPE(TSize n)
	{
		if (n < 2) return n;	// 0 -> 0, 1 -> 1, 2 -> 2, 3 -> 4
		--n;
		TSize size = 2;
		for(TSize _xSize = 1; _xSize < n; _xSize *= 2)
			size += (n + _xSize - 1) / _xSize;
		return size;
	}

	template < typename TSize >
	inline TSize sizeofLCPH(TSize n)
	{
		return sizeofLCPE(n);
	}

	template <
		class LCPFwdIt,		// lcp table input iterator
		typename TSize >
	inline void sizeofLCPE(LCPFwdIt _First, LCPFwdIt _Last, TSize &_Size)
	{
		_Size = sizeofLCPE(difference(_First, _Last));
	}

	template <
		class LCPFwdIt,		// lcp table input iterator
		typename TSize >
	inline void sizeofLCPH(LCPFwdIt _First, LCPFwdIt _Last, TSize &_Size)
	{
		sizeofLCPE(_First, _Last, _Size);
		return;
	}


    template < typename TLCPE, typename TLCP >
    inline void createLCPBinTree(TLCPE &lcp_enhanced, TLCP &lcp) {
        createLCPBinTree(begin(lcp), end(lcp), begin(lcp_enhanced));
    }


	template < typename TSize >
    inline unsigned _treeLevels(TSize lcpSize)
	{
		unsigned treeLevels = 1;
		--lcpSize;
		TSize _xSize = 1;
		for(; lcpSize > _xSize; _xSize *= 2, ++treeLevels);
        return treeLevels;
    }

    template < typename TValue, typename TConfig, typename TLCP >
    inline void createLCPBinTree(String<TValue, External<TConfig> > &lcp_enhanced, TLCP &lcp) {
        unsigned writeHeads = _treeLevels(length(lcp)) + 1;   // plus 1 write back buffer
        if (lcp_enhanced.cache.size() < writeHeads)
            lcp_enhanced.resizeCache(writeHeads);
        createLCPBinTree(begin(lcp), end(lcp), begin(lcp_enhanced));
    }

}

#endif
