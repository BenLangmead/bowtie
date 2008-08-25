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
  $Id: find_quasar.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_QUASAR_H
#define SEQAN_HEADER_FIND_QUASAR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Quasar
//////////////////////////////////////////////////////////////////////////////

//DISABLED .Class.Pattern.param.TSpec.type:Spec.Quasar

struct _Quasar;
typedef Tag<_Quasar> Quasar;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, Quasar> {
//____________________________________________________________________________
private:
	Pattern(Pattern const& other);
	Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
	Holder<TNeedle> data_needle;


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
struct Host< Pattern<TNeedle, Quasar> >
{
	typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, Quasar> const>
{
	typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, Quasar> & me, TNeedle2 const& needle) 
{
	SEQAN_CHECKPOINT
	setValue(me.data_needle, needle);
}

template <typename TNeedle, typename TNeedle2>
inline void 
setHost (Pattern<TNeedle, Quasar> & me, TNeedle2 & needle)
{
	setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, Quasar> & me) 
{
SEQAN_CHECKPOINT
}


//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Quasar>const>::Type & 
host(Pattern<TNeedle, Quasar> & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, Quasar>const>::Type & 
host(Pattern<TNeedle, Quasar> const & me)
{
SEQAN_CHECKPOINT
	return value(me.data_needle);
}

//____________________________________________________________________________


template <typename TFinder, typename TNeedle>
inline bool 
find(TFinder & finder, Pattern<TNeedle, Quasar> & me) 
{
	SEQAN_CHECKPOINT
	
	if (empty(finder)) {
		_patternInit(me);
		_finderSetNonEmpty(finder);
	} else {
		finder+=4;
	}

	return true;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
