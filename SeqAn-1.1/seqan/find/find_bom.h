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
  $Id: find_bom.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_BOM_H
#define SEQAN_HEADER_FIND_BOM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// BomAlgo
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BomAlgo:
..summary: Backward Oracle Matching algorithm. Exact string matching using a factor oracle.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, BomAlgo>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The types of the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.BomAlgo

struct _BomAlgo;
typedef Tag<_BomAlgo> BomAlgo;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, BomAlgo> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef typename Size<TNeedle>::Type TSize;
	Holder<TNeedle> data_needle;
	TSize needleLength;		
	TSize haystackLength;
	TSize step;
	Graph<Automaton<TAlphabet> > oracle;

//____________________________________________________________________________

	Pattern() {	
	}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
SEQAN_CHECKPOINT
		setHost(*this, ndl);
	}

	~Pattern() {
		SEQAN_CHECKPOINT
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, BomAlgo> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, BomAlgo> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BomAlgo> & me, TNeedle2 const& needle) 
{
	SEQAN_CHECKPOINT
	me.needleLength = length(needle);
	clear(me.oracle);
	createOracleOnReverse(me.oracle,needle);
	assignRoot(me.oracle,0);
	setValue(me.data_needle, needle);
}

template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, BomAlgo> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, BomAlgo> & me) 
{
SEQAN_CHECKPOINT
	me.step = 0;
}


//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BomAlgo>const>::Type & 
host(Pattern<TNeedle, BomAlgo> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BomAlgo>const>::Type & 
host(Pattern<TNeedle, BomAlgo> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool 
find(TFinder & finder, Pattern<TNeedle, BomAlgo> & me) 
{
	SEQAN_CHECKPOINT
	
	if (empty(finder)) {
		_patternInit(me);
		_finderSetNonEmpty(finder);
		me.haystackLength = length(container(finder));
	} else
		finder+=me.step;

	if (me.haystackLength < me.needleLength) return false;
	typedef typename Value<TNeedle>::Type TAlphabet;
	typedef Graph<Automaton<TAlphabet> > TOracle;
	typedef typename Size<TNeedle>::Type TSize;
	typedef typename VertexDescriptor<TOracle>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TOracle>::Type TEdgeDescriptor;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	while (position(finder) <= me.haystackLength - me.needleLength) {
		TVertexDescriptor current = getRoot(me.oracle);
		TSize j = me.needleLength;
		while ((j>0) &&
				(current != nilVal))
		{
			TAlphabet c = *(finder+(j-1));
			current = targetVertex(me.oracle, findEdge(me.oracle, current, c));
			--j;
		}
		if (current != nilVal) {
			me.step = j + 1;
			return true;
		}
		finder += j + 1;
	}
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
