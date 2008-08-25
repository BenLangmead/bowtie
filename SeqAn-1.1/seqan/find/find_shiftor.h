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
  $Id: find_shiftor.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_SHIFTOR_H
#define SEQAN_HEADER_FIND_SHIFTOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ShiftOr
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ShiftOr:
..summary: Exact string matching using bit parallelism. The Shift-Or algorithm is applicable to search small patterns in texts using a small alphabet.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, ShiftOr>
..param.TNeedle:The needle type.
...type:Class.String
..remarks.text:The types of the needle and the haystack have to match.
*/

///.Class.Pattern.param.TSpec.type:Spec.ShiftOr

struct _ShiftOr;
typedef Tag<_ShiftOr> ShiftOr;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, ShiftOr> {
//____________________________________________________________________________
public:
	typedef unsigned int TWord;
	Holder<TNeedle> data_needle;
	String<TWord> table;			// Look up table for each character in the alphabet (called B in "Navarro")
	String<TWord> prefSufMatch;		// Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
	TWord needleLength;				// e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
	TWord blockCount;				// #unsigned ints required to store needle	

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
void setHost (Pattern<TNeedle, ShiftOr> & me, TNeedle2 const & needle) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TNeedle>::Type TValue;
	
	me.needleLength = length(needle);
	if (me.needleLength<1) me.blockCount=1;
	else me.blockCount=((me.needleLength-1) / BitsPerValue<TWord>::VALUE)+1;
	
	clear(me.table);
	fill(me.table, me.blockCount * ValueSize<TValue>::VALUE, (TWord) ~0, Exact());

	for (TWord j = 0; j < me.needleLength; ++j) {
		// Determine character position in array table
		TWord pos = convert<TWord>(getValue(needle,j));
		me.table[me.blockCount*pos + j / BitsPerValue<TWord>::VALUE] ^= (1<<(j%BitsPerValue<TWord>::VALUE));
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
void setHost (Pattern<TNeedle, ShiftOr> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, ShiftOr> & me) 
{
SEQAN_CHECKPOINT
	typedef unsigned int TWord;

	clear(me.prefSufMatch);
	fill(me.prefSufMatch, me.blockCount, (TWord) ~0, Exact());
}

//____________________________________________________________________________


template <typename TNeedle>
inline typename Host<Pattern<TNeedle, ShiftOr>const>::Type & 
host(Pattern<TNeedle, ShiftOr> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, ShiftOr>const>::Type & 
host(Pattern<TNeedle, ShiftOr> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________

/*
template <typename TFinder, typename TNeedle>
bool _findShiftOr_SmallNeedle(TFinder & finder, Pattern<TNeedle, ShiftOr> & me) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	TWord compare= (~(1 << (me.needleLength-1)));
	while (!atEnd(finder)) {
		TWord pos = convert<TWord>(*finder);
		me.prefSufMatch[0] = (me.prefSufMatch[0] << 1) | me.table[me.blockCount*pos];
		if ((me.prefSufMatch[0] | compare) != (TWord) ~0) {
			finder-=(me.needleLength-1);
			return true; 
		}
		goNext(finder);
	}
	return false;
}
*/
template <typename TFinder, typename TNeedle>
bool _findShiftOr_SmallNeedle(TFinder & finder, Pattern<TNeedle, ShiftOr> & me) 
{
SEQAN_CHECKPOINT
	typedef typename Haystack<TFinder>::Type THaystack;
	THaystack & hstk = haystack(finder);

	typedef typename Iterator<THaystack, Standard>::Type THaystackIterator;
	THaystackIterator hayit = iter(hstk, position(finder));
	THaystackIterator hayit_end = end(hstk, Standard());

	typedef unsigned int TWord;
	TWord mask = (1 << (me.needleLength-1));
	TWord pref_suf_match = me.prefSufMatch[0];

	for (;hayit < hayit_end; ++hayit)
	{
		unsigned int pos = *hayit; //conversion alphabet to unsigned int
		pref_suf_match <<= 1;				//shift...
		pref_suf_match |= me.table[pos];	//...or
		if (pref_suf_match & mask) 
		{
			continue;
		}

		//found a hit!
		//set finder to start position
		setPosition(finder, (hayit - begin(hstk, Standard())) - me.needleLength + 1); 
		//save machine state
		me.prefSufMatch[0] = pref_suf_match; 
		return true;
	}
	return false;
}

template <typename TFinder, typename TNeedle>
bool _findShiftOr_LargeNeedle(TFinder & finder, Pattern<TNeedle, ShiftOr> & me) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;

	TWord compare = (~(1 << ((me.needleLength-1) % BitsPerValue<TWord>::VALUE)));
	while (!atEnd(finder)) {
		TWord pos = convert<TWord>(*finder);
		TWord carry = 0;
		for(TWord block=0;block<me.blockCount;++block) {
			bool newCarry = ((me.prefSufMatch[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0); 
			me.prefSufMatch[block]<<=1;
			me.prefSufMatch[block]|=carry;
			carry=newCarry;
		}
		for(TWord block=0;block<me.blockCount;++block) me.prefSufMatch[block] |= me.table[me.blockCount*pos+block];
		if ((me.prefSufMatch[me.blockCount-1] | compare) != (TWord) ~0) {
			finder-=(me.needleLength-1);
			return true;  
		}

		/*
		// Debug code
		std::cout << "   ";
		for(int j=0;j<me.blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
				std::cout << ((me.prefSufMatch[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		*/
		goNext(finder);
	}
	return false;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, ShiftOr> & me) {
	SEQAN_CHECKPOINT

	if (empty(finder)) {
		_patternInit(me);
		_finderSetNonEmpty(finder);
	} else
		finder += me.needleLength;

	// Fast algorithm for needles < machine word?
	if (me.blockCount == 1) {
		return _findShiftOr_SmallNeedle(finder, me);
	} else {
		return _findShiftOr_LargeNeedle(finder, me);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTOR_H
