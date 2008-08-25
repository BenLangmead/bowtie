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
  $Id: find_wumanber.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_DDDOC: do not generate documentation for this file

#ifndef SEQAN_HEADER_FIND_WUMANBER_H
#define SEQAN_HEADER_FIND_WUMANBER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// WuManber
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.WuManber:
..general:Class.Pattern
..cat:Searching
..summary:A fast online-algorithm for multi-pattern search
..signature:Pattern<TNeedle, WuManber<HashType> >
..param.TNeedle:The needle type.
...type:Class.String
..param.HashType: Defines the used Hashing-Function 
...default:@Class.DefaultHash@
*/

///.Class.Pattern.param.TSpec.type:Spec.WuManber

//////////////////////////////////////////////////////////////////////////////

/**
.Class.DefaultHash:
..cat:Miscellaneous
..summary:A simple Hash-Object
..summary:Hashing is based on bit-shift
..signature:DefaultHash
..remarks: The shift-width is defined by the used Alphabet
..remarks: The shift-width could also be set with the ctor 
..see:Spec.WuManber
*/

///.Spec.WuManber.param.HashType.type:Class.DefaultHash

class DefaultHash
{
private:
	unsigned int shift;
	Size<String<size_t> >::Type h_size;

	void * _CompareTypeId;

public:
	DefaultHash()
		: shift(0),
		h_size(0),
		_CompareTypeId(0)
	{
	}

	DefaultHash(DefaultHash const & _other)
		: shift(_shift_width(_other)),
		h_size(_hash_size(_other)),
		_CompareTypeId(_other._CompareTypeId)
	{
	}

	DefaultHash & 
	operator = (DefaultHash const & other_)
	{
SEQAN_CHECKPOINT
		shift = other_.shift;
		h_size = other_.h_size;
		_CompareTypeId = other_._CompareTypeId;
		return *this;
	}

/**
.DISABLED.Function.init:
..summary:Initializes the Hash-Object
..signature:init(Hash,TNeedle [, _shift_width])
..param.DefaultHash: Reference to the Hash-Object that should be Initialized
..param.TNeedle: Needle for which the Hash-Object should work
..param._shift_width: Specifiers a user defined shift with
*/

///.DISABLED.Function.init.param.Hash.type:Class.DefaultHash

	template <typename TNeedle>
	friend inline void
	init(DefaultHash & me, TNeedle & ndl)
	{
SEQAN_CHECKPOINT
		typedef typename Value<TNeedle>::Type TTargetAlphabet;
		me._CompareTypeId = _ClassIdentifier<TTargetAlphabet>::getID();

		if(_ClassIdentifier<TTargetAlphabet>::getID() == _ClassIdentifier<char>::getID())
			_setShift_Width(me,5);
		else
			_setShift_Width(me,BitsPerValue<TTargetAlphabet>::VALUE);

		// bestimmen der hash bzw. shift-tabellen größe
		typename Size<TNeedle>::Type ndl_size = length(ndl);

		int i = 1;
		size_t lmin;
		lmin = length(ndl[0]);
		size_t w_width;
		while(i < ndl_size)
		{
			if(lmin > length(ndl[i]))
				lmin = length(ndl[i]);
			++i;
		}
		w_width = static_cast<unsigned int>(ceil(std::log(static_cast<double>(2*lmin*ndl_size))/std::log(static_cast<double>(ValueSize<TTargetAlphabet>::VALUE))));

		me.h_size = 1 << me.shift * w_width;
	}

	template <typename TNeedle,typename THaystack>
	friend inline void
	init(DefaultHash & me,TNeedle & ndl,THaystack & finder,unsigned int _shift_width)
	{
SEQAN_CHECKPOINT
		typedef typename Value<TNeedle>::Type TTargetAlphabet;
		me._CompareTypeId = _ClassIdentifier<TTargetAlphabet>::getID();
		_setShift_Width(me,_shift_width);

		// bestimmen der hash bzw. shift-tabellen größe
		typename Size<TNeedle>::Type ndl_size = length(ndl);

		int i = 1;
		size_t lmin;
		lmin = length(ndl[0]);
		size_t w_width;
		while(i < ndl_size)
		{
			if(lmin > length(ndl[i]))
				lmin = length(ndl[i]);
			++i;
		}
		w_width = static_cast<unsigned int>(ceil(std::log(static_cast<double>(2*lmin*ndl_size))/std::log(static_cast<double>(ValueSize<TTargetAlphabet>::VALUE))));

		me.h_size = 1 << me.shift * w_width;
	}

/**
.Internal._shift_width:
..summary:Returns the Shift-Width that was defined by the Hash-Object used by @Spec.WuManber@ Pattern
..signature:shift_width(Hash)
..param.Hash: Reference to the Hash-Object
*/

///.Internal._shift_width.param.Hash.type:Class.DefaultHash
	friend inline unsigned int &
	_shift_width(DefaultHash & me)
	{
SEQAN_CHECKPOINT
		return me.shift;
	}

	friend inline unsigned int const &
	_shift_width(DefaultHash const & me)
	{
SEQAN_CHECKPOINT
		return me.shift;
	}

	friend inline void 
	_setShift_Width(DefaultHash & me, unsigned int _shift)
	{
SEQAN_CHECKPOINT
		me.shift = _shift;
	}

/**
.Internal._hash_size:
..summary:Returns the Size of the HASH-Tables used by @Spec.WuManber@ Pattern
..signature:hash_size(Hash)
..param.Hash: Reference to the Hash-Object
*/

///.Internal._hash_size.param.Hash.type:Class.DefaultHash
	friend inline Size<String<size_t> >::Type &
	_hash_size(DefaultHash & me)
	{
SEQAN_CHECKPOINT
		return me.h_size;
	}

	friend inline Size<String<size_t> >::Type const &
	_hash_size(DefaultHash const & me)
	{
SEQAN_CHECKPOINT
		return me.h_size;
	}

	friend inline Size<String<size_t> >::Type &
	_shift_size(DefaultHash & me)
	{
SEQAN_CHECKPOINT
		return me.h_size;
	}

	friend inline Size<String<size_t> >::Type const &
	_shift_size(DefaultHash const & me)
	{
SEQAN_CHECKPOINT
		return me.h_size;
	}

};





//////////////////////////////////////////////////////////////////////////////

template <typename THash = DefaultHash>
class WuManber {};

//////////////////////////////////////////////////////////////////////////////

template<typename TNeedle, typename THash>
class Pattern<TNeedle, WuManber<THash> >
{
private:
	unsigned int			found_in;
	String<size_t>			SHIFT;
	String<String<size_t> >	HASH;
	unsigned int			q;
	size_t					lmin;

	Holder<TNeedle>			data_needle;
	Holder<THash>			_hash;
	bool					_initialized;

public:

	Pattern():
		_initialized(false)
	{
SEQAN_CHECKPOINT
		// erzeuge neues Hashobjekt
		THash h;
		setHash(*this,h);
		create(_hash);

		q = 0;
		lmin = 0;
		found_in = -1;
	}

	Pattern(Pattern const & other_):
		found_in(other_.found_in),
		SHIFT(other_.SHIFT),
		HASH(other_.HASH),
		q(window_width(other_)),
		lmin(minLength(other_)),
		_initialized(other_._initialized)
	{
SEQAN_CHECKPOINT
		assignHash(*this,hash(other_));
	}

	~Pattern()
	{
		clear(_hash);
	}

	Pattern & 
	operator = (Pattern const & other_)
	{
SEQAN_CHECKPOINT
		found_in = needle(other_);
		SHIFT = other_.SHIFT;
		HASH = other_.HASH;
		q = other_.q;
		lmin = other_.lmin;
		assignHash(*this,hash(other_));
		_initialized = other_._initialized;
		return *this;
	}
//____________________________________________________________________________

	friend inline typename Host<Pattern>::Type & 
	host(Pattern & me)
	{
SEQAN_CHECKPOINT
		return value(me.data_needle);
	}

	friend inline typename Host<Pattern const>::Type & 
	host(Pattern const & me)
	{
SEQAN_CHECKPOINT
		return value(me.data_needle);
	}

//____________________________________________________________________________

	friend inline unsigned int &
	needle(Pattern & me)
	{
SEQAN_CHECKPOINT
		return me.found_in;
	}

	friend inline unsigned int const &
	needle(Pattern const & me)
	{
SEQAN_CHECKPOINT
		return me.found_in;
	}

/**
.Function.hash:
..cat:Searching
..summary:Returns the Hashobject, that is used by the Finderobject
..signature:hash(Pattern)
..param.Pattern: Reference to the Pattern-Object
*/

///.Function.hash.param.Pattern.type:Spec.WuManber


	friend inline THash &
	hash(Pattern & me)
	{
SEQAN_CHECKPOINT
		return value(me._hash);
	}

	friend inline THash const &
	hash(Pattern const & me)
	{
SEQAN_CHECKPOINT
		return value(me._hash);
	}


/**
.Function.setHash:
..cat:Searching
..summary:Sets the Hashobject to a new Value
..signature:setHash(Pattern,THash)
..param.Pattern: Reference to the Pattern-Object
..param.THash: The new Hashobject
*/

///.Function.setHash.param.Pattern.type:Spec.WuManber
///.Function.setHash.param.THash.type:Class.DefaultHash

	friend inline void 
	setHash(Pattern & me,THash & value)
	{
SEQAN_CHECKPOINT
		setValue(me._hash,value);
	}

/**
.Function.assignHash:
..cat:Searching
..summary: Sets the Hashobject to a new Value which is independent of the original Object
..signature:assignHash(Pattern,THash)
..param.Pattern: Reference to the Pattern-Object
..param.THash: The new Hashobject
*/

///.Function.assignHash.param.Pattern.type:Spec.WuManber
///.Function.assignHash.param.THash.type:Class.DefaultHash


	friend inline void
	assignHash(Pattern &  me,THash const & value)
	{
SEQAN_CHECKPOINT
		assignValue(me._hash,value);
	}
//____________________________________________________________________________


	friend inline void
	setNeedle(Pattern & me, unsigned int const needleIndex_)
	{
SEQAN_CHECKPOINT
		me.found_in = needleIndex_;
	}

/**
.Function.window_with:
..cat:Searching
..summary:Returns the window-with which is used to search the Haystack
..signature:window_with(Pattern)
..param.Pattern: Reference to the Pattern-Object
*/

///.Function.window_with.param.Pattern.type:Spec.WuManber


	friend inline unsigned int &
	window_width(Pattern & me)
	{
SEQAN_CHECKPOINT
		return me.q;
	}	

	friend inline unsigned int const &
	window_width(Pattern const & me)
	{
SEQAN_CHECKPOINT
		return me.q;
	}	

/**
.Function.setWindow_with:
..cat:Searching
..summary:Sets the window-with to a new value
..signature:window_with(Pattern,value_)
..param.Pattern: Reference to the Pattern-Object
..param.value_: The new window-with
*/

///.Function.setWindow_with.param.Pattern.type:Spec.WuManber

	friend inline void
	setWindow_width(Pattern & me,unsigned int value_)
	{
SEQAN_CHECKPOINT
		me.q = value_;
	}

/**
.Function.minLength:
..cat:Searching
..summary:Returns the length of the shortest Searchpattern
..signature:minLength(Pattern,)
..param.Pattern: Reference to the Pattern-Object
*/

///.Function.setWindow_with.param.Pattern.type:Spec.WuManber


	friend inline size_t &
	minLength(Pattern & me)
	{
SEQAN_CHECKPOINT
		return me.lmin;
	}

	friend inline size_t const &
	minLength(Pattern const & me)
	{
SEQAN_CHECKPOINT
		return me.lmin;
	}

/**
.Function.setMinLength:
..cat:Searching
..summary:Sets the length of the shortest Searchpattern
..signature:minLength(Pattern,value_)
..param.Pattern: Reference to the Pattern-Object
..param.value_: The new min. Length
*/

///.Function.setWindow_with.param.Pattern.type:Spec.WuManber

	friend inline void
	setMinLength(Pattern & me,size_t value_)
	{
SEQAN_CHECKPOINT
		me.lmin = value_;
	}
//____________________________________________________________________________

/**
.Function.setNeedle:
..signature:setNeedle(pattern, needle[, param_init])
..param.param_init: For @Spec.WuManber@ only. If true (default) the $window_with$ of the WuManber-Pattern is set by @Function.setNeedle@, otherwise you got to set $window_with$ manually.
..remarks:If $param_init$ is $false$ the $window_with$ must be defined before @Function.init@ is called, otherwise @Function.find@ won't work.
*/

	template <typename TNeedle2>
	friend inline void
	setHost(Pattern & me, TNeedle2 & ndl, bool param_init = true)
	{
SEQAN_CHECKPOINT

		setNeedle(me,0);
		clear(me.SHIFT);
		clear(me.HASH);

		typename Size<TNeedle>::Type ndl_size = length(ndl);
		typedef typename Value<TNeedle>::Type TNdl;
		typedef typename Value<TNdl>::Type TAlphabet;
		typedef typename Host<THaystack>::Type THaystackHost;
		THaystackHost & haystack_host = container(finder);

		// bestimmen von lmin
		int i = 1;
		setMinLength(me,length(ndl[0]));
		while(i < ndl_size)
		{
SEQAN_CHECKPOINT
			if(minLength(me) > length(ndl[i]))
				setMinLength(me,length(ndl[i]));
			++i;
		}

		// wennn param_init gesetzt wird, dann wird die fensterweite angepasst
		// ist dies nicht der fall muss diese schon vorher gesetzt worden sein
		if(param_init)
		{
SEQAN_CHECKPOINT
			setWindow_width(me,static_cast<unsigned int>(ceil(std::log(static_cast<double>(minLength(me)*ndl_size))/std::log(static_cast<double>(ValueSize<TAlphabet>::VALUE)))));
			if(window_width(me) > minLength(me))
				setWindow_width(me,minLength(me));
		}

		// SHIFT Tabelle aufbauen
		resize(me.SHIFT,_shift_size(hash(me)));
		arrayFill(begin(me.SHIFT),length(me.SHIFT),minLength(me) - window_width(me) + 1);

		// HASH-Tabelle aufbauen
		resize(me.HASH,_hash_size(hash(me)));
		i = 0;
		int j;
		int h_hash;
		int s_hash;
		int n_length;

		while(i < ndl_size)
		{
SEQAN_CHECKPOINT
			// für alle needles
			j = window_width(me) - 1;
			n_length = length(ndl[i]);

			while(j < n_length)
			{
SEQAN_CHECKPOINT
				s_hash = _compute_SHIFT_Hash(hash(me),infix(value(ndl,i),j - window_width(me) + 1,j),haystack_host,window_width(me));
				me.SHIFT[s_hash] = ::std::min(static_cast<size_t>(n_length) - (j + 1),static_cast<size_t>(me.SHIFT[s_hash]));
				++j;
			}

			h_hash = _compute_HASH_Hash(hash(me),infix(value(ndl,i),j - window_width(me),j),haystack_host,window_width(me));
			resize(me.HASH[h_hash],length(me.HASH[h_hash]) + 1);
			(me.HASH[h_hash])[length(me.HASH[h_hash]) - 1] =  i;
			++i;
		}

		data_needle = ndl;

#ifdef SEQAN_DEBUG
		std::cout << "initialized WuManber with following parameters" << ::std::endl;

		std::cout << "lmin " << minLength(me) << ::std::endl;
		std::cout << "q " << window_width(me) << ::std::endl;
		std::cout << "number of needles " << ndl_size << ::std::endl;
#endif
	}

//____________________________________________________________________________

/**

	template <typename TFinder>
	friend inline bool find(TFinder & finder, Pattern & me)
	{
SEQAN_CHECKPOINT
SEQAN_ASSERT2(length(ndl) > 0, "need to search for at least one needle in haystack")		
		// Parameter berechnen
		typedef typename Host<THaystack>::Type THaystackHost;
		THaystackHost & haystack_host = container(finder);

		typename Position<THaystack>::Type pos;
		typename Size<THaystackHost>::Type hstk_size = length(haystack_host);
		typename Size<TNeedle>::Type ndl_size = length(ndl);

		// Typ der einzelnen Suchwörter
		typedef typename Value<TNeedle>::Type TNdl;
		typedef typename Value<TNdl>::Type TAlphabet;

		TNdl ndl_support = value(ndl,0);

		int i = 1;

		// infix wird einmal berechnet und dann zwischengespeichert um nicht
		// mehrmals die infix()-Function aufzurufen
		size_t shift = 0;
		unsigned int h_hash = 0;
		unsigned int _needle = 0;
		unsigned int h_length;
		int infix_start;

SEQAN_ASSERT2(hstk_size > minLength(me),"nothing to find haystack is smaller then needle")		
		// Matchs berechnen
		if(empty(finder))
		{
SEQAN_CHECKPOINT
			goBegin(finder);
			init(hash(me), host(me), finder);
			finder += minLength(me);
			pos = position(finder);
		}
		else
		{
SEQAN_CHECKPOINT
			finder += length(host(me)[needle(me)]);
			pos =  position(finder);

			shift = me.SHIFT[_compute_SHIFT_Hash(hash(me),infix(haystack_host,pos - window_width(me),pos),ndl_support,window_width(me))];
			if(shift == 0)
			{
SEQAN_CHECKPOINT
				h_hash = _compute_HASH_Hash(hash(me),infix(haystack_host,pos - window_width(me),pos),ndl_support,window_width(me));
				//Überprüfen, ob es noch weitere Funde an dieser Position gibt
				i = 0;
				h_length = length(me.HASH[h_hash]);
				while(i < h_length)
				{
SEQAN_CHECKPOINT
					// verhindert doppelte funde
					_needle = (me.HASH[h_hash])[i];
					if(_needle > needle(me))
					{
SEQAN_CHECKPOINT
						infix_start = static_cast<int>(pos - length(ndl[_needle]));
						if(infix_start < 0)
						{
							++i;
							continue;
						}
						// TODO: Alternative Lösung
						if(ndl[_needle] == infix(haystack_host,infix_start,pos))
						{
SEQAN_CHECKPOINT
							setNeedle(me,_needle);
							// setzt den finder-Iterator auf den Anfang des gefundenen patterns
							finder -= length(ndl[_needle]);
							return true;
						}
					}
					++i;
				}
			}
			++finder;
			++pos;
		}

		while(pos <= hstk_size)
		{
SEQAN_CHECKPOINT
			shift = me.SHIFT[_compute_SHIFT_Hash(hash(me),infix(haystack_host,pos - window_width(me),pos),ndl_support,window_width(me))];
			if(shift == 0)
			{
SEQAN_CHECKPOINT
				//Überprüfen, ob es sich um einen Fund handelt
				h_hash = _compute_HASH_Hash(hash(me),infix(haystack_host,pos - window_width(me),pos),ndl_support,window_width(me));
				i = 0;
				h_length = length(me.HASH[h_hash]);
				while(i < h_length)
				{
SEQAN_CHECKPOINT
					_needle = (me.HASH[h_hash])[i];
					infix_start = static_cast<int>(pos - length(ndl[_needle]));
					if(infix_start < 0)
					{
						++i;
						continue;
					}
					if(ndl[_needle] == infix(haystack_host,infix_start,pos))
					{
SEQAN_CHECKPOINT
						setNeedle(me,_needle);
						// setzt den finder-Iterator auf den Anfang des gefundenen patterns
						finder -= length(ndl[_needle]);
						return true;
					}
					++i;
				}
				//Falls nichts an dieser Stelle gefunden wird, dann rücke um 1 weiter
				++finder;
				++pos;
			}
			else
			{
SEQAN_CHECKPOINT
				finder += shift;
				pos += shift;
			}

		}
		return false;
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Internal._compute_SHIFT_Hash:
..cat:Functions
..summary:Calculates a Hash-Value a specific part of the Sequence in the SHIFT-Table.
..signature:_compute_SHIFT_Hash(Hash & ,Sequence &,start, window_width)
..param.Hash: Reference to the Pattern-Object that should be Initialized
..param.Sequence: Needle for which the Pattern-Object should search
..param.start: Position of the first element of the segment, that should be hashed
..param.window_width: length of the infix
..see:Spec.WuManber
*/

///.Internal._compute_SHIFT_Hash.param.Hash.type:Class.DefaultHash
	template <typename TSequence, typename TSupport>
	friend inline unsigned int
	_compute_SHIFT_Hash(DefaultHash & me,TSequence const & sequence,TSupport & supp, unsigned int length)
	{
SEQAN_CHECKPOINT
		typedef typename CompareType< Value< Host < TSequence >::Type >::Type, Value< Container< TSupport >::Type >::Type >::Type TTargetAlphabet;

		unsigned int ret = 0;
		for(int i = 0;i < length;++i)
		{
SEQAN_CHECKPOINT
			ret <<= _shift_width(me);
			ret += static_cast<unsigned int>(static_cast<TTargetAlphabet>(value(sequence,i)));
		}
		return ret % _hash_size(me);
	}

/**
.Internal._compute_HASH_Hash:
..cat:Functions
..summary:Calculates a Hash-Value a specific part of the Sequence in the HASH-Table.
..signature:_compute_HASH_Hash(Hash & ,Sequence &,start, window_width)
..param.Hash: Reference to the Pattern-Object that should be Initialized
..param.Sequence: Needle for which the Pattern-Object should search
..param.start: Position of the first element of the segment, that should be hashed
..param.window_width: length of the infix
..see:Spec.WuManber
*/

///.Internal._compute_HASH_Hash.param.Hash.type:Class.DefaultHash

	template <typename TSequence, typename TSupport>
	friend inline unsigned int
	_compute_HASH_Hash(DefaultHash & me,TSequence const & sequence,TSupport & supp, unsigned int length)
	{
SEQAN_CHECKPOINT
		return _compute_SHIFT_Hash(me,sequence,supp,length);
	}
////////////////////////////////////////////////////////////////////////////////

};


//////////////////////////////////////////////////////////////////////////////
// Host
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename THash>
struct Host< Pattern<TNeedle, WuManber<THash> > >
{
	typedef TNeedle Type;
};

template <typename TNeedle, typename THash>
struct Host< Pattern<TNeedle, WuManber<THash> > const>
{
	typedef TNeedle const Type;
};


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
