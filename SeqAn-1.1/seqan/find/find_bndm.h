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
  $Id: find_bndm.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_BNDMALGO_H
#define SEQAN_HEADER_FIND_BNDMALGO_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// BndmAlgo
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BndmAlgo:
..summary: Backward Nondeterministic Dawg Matching algorithm. Exact string matching using bit parallelism.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, BndmAlgo>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The types of the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.BndmAlgo

struct _BndmAlgo;
typedef Tag<_BndmAlgo> BndmAlgo;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, BndmAlgo> {
//____________________________________________________________________________
public:
	typedef unsigned int TWord;

	Holder<TNeedle> data_needle;
	String<TWord> table;			// Look up table for each character in the alphabet (called B in "Navarro")
	String<TWord> activeFactors;	// The active factors in the pattern (called D in "Navarro")
	TWord needleLength;				// e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
	TWord haystackLength;			// Length of haystack
	TWord blockCount;				// #unsigned ints required to store needle
	TWord last;

//____________________________________________________________________________

	Pattern() {}

	template <typename TNeedle2>
	Pattern(TNeedle2 const & ndl)
	{
		setHost(*this, ndl);
	}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, BndmAlgo> & me, TNeedle2 const& needle) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TNeedle>::Type TValue;

	me.needleLength = length(needle);
	if (me.needleLength<1) me.blockCount=1;
	else me.blockCount=((me.needleLength-1) / BitsPerValue<TWord>::VALUE)+1;
			
	clear(me.table);
	fill(me.table, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

	for (TWord j = 0; j < me.needleLength; ++j) {
		// Determine character position in array table
		TWord pos = convert<TWord>(getValue(needle,j));
		me.table[me.blockCount*pos + j / BitsPerValue<TWord>::VALUE] |= (1<<(j%BitsPerValue<TWord>::VALUE));
	}

	setValue(me.data_needle, needle);
	/*
	// Debug code
	std::cout << "Alphabet size: " << ValueSize<TValue>::VALUE << ::std::endl;
	std::cout << "Needle length: " << me.needleLength << ::std::endl;
	std::cout << "Block count: " << me.blockCount << ::std::endl;

	for(unsigned int i=0;i<ValueSize<TValue>::VALUE;++i) {
		if ((i<97) || (i>122)) continue;
		std::cout << static_cast<char>(i) << ": ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.table[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
	}
	*/
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, BndmAlgo> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, BndmAlgo> & me) 
{
SEQAN_CHECKPOINT
	clear(me.activeFactors);
	fill(me.activeFactors, me.blockCount, 0, Exact());
	me.last = 0;
}

//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BndmAlgo>const>::Type & 
host(Pattern<TNeedle, BndmAlgo> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, BndmAlgo>const>::Type & 
host(Pattern<TNeedle, BndmAlgo> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool _findBndm_SmallNeedle(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	if (me.haystackLength < me.needleLength) return false;
	while (position(finder) <= me.haystackLength - me.needleLength) {
		me.last=me.needleLength;
		TWord j=me.needleLength;
		me.activeFactors[0]=~0;
		while (me.activeFactors[0]!=0) {
			TWord pos = convert<TWord>(*(finder+j-1));
			me.activeFactors[0] = (me.activeFactors[0] & me.table[me.blockCount*pos]);
			j--;
			if (me.activeFactors[0] & 1 != 0) {
				if (j>0) me.last=j;
				else return true;
			}
			me.activeFactors[0] = me.activeFactors[0] >> 1;
		}
		finder+=me.last;
	}
	return false;
}

template <typename TFinder, typename TNeedle>
inline bool _findBndm_LargeNeedle(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	TWord carryPattern = (1<< (BitsPerValue<TWord>::VALUE - 1));
	while (position(finder) <= me.haystackLength - me.needleLength) {
		me.last=me.needleLength;
		TWord j=me.needleLength;
		for(TWord block=0;block<me.blockCount;++block) me.activeFactors[block]=~0;
		bool zeroPrefSufMatch = false;
		while (!zeroPrefSufMatch) {
			TWord pos = convert<TWord>(*(finder+j-1));

			/*	
			// Debug code
			std::cout << "   ";
			for(int j=0;j<me.blockCount;++j) {
				for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
					std::cout << ((me.activeFactors[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
				}
			}
			std::cout << ::std::endl;
			*/

			for(TWord block=0;block<me.blockCount;++block) me.activeFactors[block] &= me.table[me.blockCount*pos+block];
			j--;
			if (me.activeFactors[0] & 1 != 0) {
				if (j>0) me.last=j;
				else return true;
			}
			bool carry=0;
			zeroPrefSufMatch=true;
			for(int block=me.blockCount-1;block>=0;--block) {
				bool newCarry=((me.activeFactors[block] & 1)!=0); 
				me.activeFactors[block]>>=1;
				if (carry) me.activeFactors[block]|=carryPattern;
				carry=newCarry;
				if (me.activeFactors[block]!=0) zeroPrefSufMatch=false;
			}
		}
		finder+=me.last;
	}
	return false;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, BndmAlgo> & me) {
	SEQAN_CHECKPOINT

	if (empty(finder)) {
		_patternInit(me);
		_finderSetNonEmpty(finder);
		me.haystackLength = length(container(finder));
	} else
		finder+=me.last;

	// Fast algorithm for needles < machine word?
	if (me.blockCount == 1) {
		return _findBndm_SmallNeedle(finder, me);
	} else {
		return _findBndm_LargeNeedle(finder, me);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
