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
  $Id: misc_map.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_MISC_MAP_H
#define SEQAN_HEADER_MISC_MAP_H

#include <algorithm>
#include "misc_base.h"


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TPair>
	inline typename TPair::T1 & keyOf(TPair &pair) {
		return getValueI1(pair);
	}
	template <typename TPair>
	inline typename TPair::T1 const & keyOf(TPair const &pair) {
		return getValueI1(pair);
	}
	template <typename TPair>
	inline typename TPair::T2 & objectOf(TPair &pair) {
		return getValueI2(pair);
	}
	template <typename TPair>
	inline typename TPair::T2 const & objectOf(TPair const &pair) {
		return getValueI2(pair);
	}


	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TString, typename TLess = ::std::less<typename Value<TString>::Type::T1> >
	struct SequenceMap;


	template <typename TString>
	struct Value< SequenceMap<TString> > {
		typedef typename Value<TString>::Type Type;
	};
	template <typename TString>
	struct Size< SequenceMap<TString> > {
		typedef typename Size<TString>::Type Type;
	};

	template <typename TString>
	inline typename Size< SequenceMap<TString> >::Type length(SequenceMap<TString> const &set) {
		return length(set.string);
	}



	template <typename TKey, typename TObject, typename TSpec, typename TLess>
	struct SequenceMap< String< Pair<TKey, TObject>, TSpec >, TLess > {
		typedef Pair<TKey, TObject>						TValue;
		typedef String< Pair<TKey, TObject>, TSpec >	TSequence;
		typedef typename TVector::size_type				TSize;

		TKey				maxKey;
		TLess				comp;
		Holder<TSequence>	string;

		SequenceMap():
			maxKey(0)	{}

		SequenceMap(TLess const &_comp):
			maxKey(0),
			comp(_comp) {}
	};


	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TKey>
	struct Map {
//		typedef ::std::map<TKey> Type;
	};
	template <typename TKey, typename TObject>
	struct Map< Pair<TKey, TObject> > {
//		typedef ::std::set< Pair<TKey, TObject>, SetLess< Pair<TKey, TObject> > > Type;
	};



	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TString, typename TLess>
	struct Iterator< SequenceMap<TString, TLess> > {
		typedef Iterator<TString>::Type	Type;
	}

	template <typename TString, typename TLess>
	struct Iterator< SequenceMap<TString, TLess> const > {
		typedef Iterator<TString const>::Type	Type;
	}


	//////////////////////////////////////////////////////////////////////////////
	//

	template <typename TString, typename TLess>
	inline void clear(SequenceMap<TString, TLess> &map) {
		clear(map.string);
		map.maxKey = 0;
	}

	template <typename TKey, typename TString, typename TLess>
	inline typename Iterator< SequenceMap<TString, TLess> >::Type find(TKey const &key, SequenceMap<TString, TLess> &map) {

		typedef typename Size< SequenceMap<TString, TLess> >::Type		TSize;
		typedef typename Iterator< SequenceMap<TString, TLess> >::Type	TIter;
		
		// accelerate binary search
		if (map.comp(map.maxKey, key))
			return end(map.string);

		TIter _First = begin(map.string);
		TSize _Count = end(map.string) - _First;

		for (; 0 < _Count; )
		{	// divide and conquer, find half that contains answer
			TSize _Count2 = _Count / 2;
			TIter _Mid = _First + _Count2;

			if (map.comp((*_Mid).i1, key))
				_First = ++_Mid, _Count -= _Count2 + 1;
			else
				_Count = _Count2;
		}
		return (_First);
	}

	template <typename TPair, typename TString, typename TLess>
	inline void insert(TPair &pair, SequenceMap<TString> &map)
	{
		typedef typename Iterator< SequenceMap<TString, TLess> >::Type	TIter;

		// accelerate binary search
		if (map.comp(map.maxKey, pair.i1)) {
// whatever
			map.maxKey = pair.i1;
			return append(map.string);
		}

		TIter iter = find(pair.i1, map);
// whatever
		insert(pair, map);
	}




	template <typename TKey, typename TString>
	inline void erase(TKey const &key, SequenceMap<TString> &map)
	{
		if (!max.comp(map.maxKey, key) && !max.comp(key, map.maxKey)) {
			// delete last element
			if (length(map.string))
				map.maxKey = *(map.string.end() - 1).i2;
			else
				map.maxKey = 0;
		}

		typedef typename Iterator< SequenceMap<TString, TLess> >::Type	TIter;
		TIter iter = find(pair.i1, map);
		delete iter;
	}

}

#endif

