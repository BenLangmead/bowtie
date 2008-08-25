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
  $Id: index_sa_btree.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_SA_BTREE_H
#define SEQAN_HEADER_INDEX_SA_BTREE_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	template <
		class SAFwdIt,		// suffix array input iterator
		class FlatOutIt >	// flat tree output iterator
	inline FlatOutIt createSABTree(
		SAFwdIt _First, SAFwdIt _Last,
		FlatOutIt _Dest, unsigned BlockSize)
	{
        typedef typename Value<SAFwdIt>::Type	TSize;

        TSize size = difference(_First, _Last);
		if (!size) return _Dest;

        // calculate the depth of the sa b-tree
		TSize BlockElements = BlockSize - 1;

		unsigned treeLevels = 1;
		TSize _xSize;
		for(_xSize = 1; _xSize * BlockSize <= size; _xSize *= BlockSize, ++treeLevels);
	
		// get output iterators for every level in the flat tree
		FlatOutIt *level = new FlatOutIt[treeLevels];
		for(int i = treeLevels - 1; _xSize; --i, _xSize /= BlockSize) {
			level[i] = _Dest;
			goFurther(_Dest, ((size / _xSize + BlockSize - 1) / BlockSize) * BlockSize);
		}

		// counter for each b-tree level
		TSize *cnt = new TSize[treeLevels];
		for(int i = 0; i < treeLevels; ++i)
			cnt[i] = 0;

		// distribute to responsible levels
		for(TSize j = 0; j < size; ++j, ++_First)
			for(int i = 0; i < treeLevels; ++i) {
				*(level[i]) = *_First;
				++(level[i]);
				if (cnt[i] != BlockElements) {
					++cnt[i];
					break;
				} else
					cnt[i] = 0;
			}

		delete[] cnt;
		delete[] level;

		return _Dest;
    }


	template < typename TSize >
	inline TSize sizeofSAB(TSize n, unsigned BlockSize)
	{
		TSize size = 0;
		for(TSize _xSize = 1; _xSize <= n; _xSize *= BlockSize)
			size += (n / _xSize + BlockSize - 1) / BlockSize;
		return size * BlockSize;
	}

	template <
		class SAFwdIt,		// suffix array input iterator
		typename TSize >
	inline void sizeofSAB(SAFwdIt _First, SAFwdIt _Last, TSize &_Size, unsigned BlockSize)
	{
		_Size = sizeofSAB(difference(_First, _Last), BlockSize);
	}


    template < typename TSAB, typename TSA >
    inline void createSABTree(TSAB &sa_btree, TSA &sa, unsigned BlockSize) {
        createSABTree(begin(sa), end(sa), begin(sa_btree), BlockSize);
    }


	template < typename TSize >
    inline unsigned treeLevelsSAB(TSize saSize, unsigned BlockSize)
	{
		unsigned treeLevels = 1;
		for(TSize _xSize = 1; _xSize <= saSize; _xSize *= BlockSize, ++treeLevels);
        return treeLevels;
    }

    template < typename TValue, typename TConfig, typename TSA >
    inline void createSABTree(String<TValue, External<TConfig> > &sa_btree, TSA &sa, unsigned BlockSize) {
        int writeHeads = treeLevelsSAB(length(sa)) + 1;   // plus 1 write back buffer
        if (sa_btree.cache.size() < writeHeads)
            sa_btree.resizeCache(writeHeads);
        createSABTree(begin(sa), end(sa), begin(sa_btree), BlockSize);
    }

}

#endif
