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
  $Id: modifier_reverse.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_REVERSE_H
#define SEQAN_HEADER_MODIFIER_REVERSE_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ModReverse:
..summary:Mirrors the characters from begin to end.
..cat:Modifier
..general:Class.ModifiedIterator
..general:Class.ModifiedString
..signature:ModifiedIterator<THost, ModReverse>
..signature:ModifiedString<THost, ModReverse>
..param.THost:Original string/iterator.
...type:Concept.Iterator
*/

	struct ModReverse {};

//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// reverse iterator
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost>
	struct Cargo< ModifiedIterator<THost, ModReverse> > {
		typedef Cargo Type;		// to reduce namespace pollution
		bool _atEnd;
		Cargo(): _atEnd(false) {}
	};

	template <typename THost>
	class ModifiedIterator<THost, ModReverse> {
	public:
		Holder<THost, Simple>					data_host;
		typename Cargo<ModifiedIterator>::Type	data_cargo;

		ModifiedIterator(ModifiedIterator &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		ModifiedIterator(ModifiedIterator const &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		template <typename T>
		ModifiedIterator(T & _origin) {
			assign(*this, _origin);
		}

		template <typename T>
		ModifiedIterator(T const & _origin) {
			assign(*this, _origin);
		}
//____________________________________________________________________________

		template <typename T>
		inline ModifiedIterator const &
		operator = (T & _origin) {
			assign(*this, _origin);
			return *this;
		}

		template <typename T>
		inline ModifiedIterator const &
		operator = (T const & _origin) {
			assign(*this, _origin);
			return *this;
		}
	};

	//////////////////////////////////////////////////////////////////////////////
	// operator ++
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline void
	goNext(ModifiedIterator<THost, ModReverse> & me)
	{
	SEQAN_CHECKPOINT
		if (atBegin(host(me)))
			cargo(me)._atEnd = true;
		else
			goPrevious(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator --
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline void
	goPrevious(ModifiedIterator<THost, ModReverse> & me)
	{
	SEQAN_CHECKPOINT
		if (cargo(me)._atEnd)
			cargo(me)._atEnd = false;
		else
			goNext(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator +
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TDelta>
	inline ModifiedIterator<THost, ModReverse> &
	operator += (ModifiedIterator<THost, ModReverse> & me, TDelta delta_) 
	{
		typedef ModifiedIterator<THost, ModReverse> TIterator;
		typedef typename Position<TIterator>::Type TPosition;
		TPosition delta = delta_;

		if (delta == 0)
		{
			return me;
		}
		if (delta > 0)
		{
			if (position(host(me)) < delta) {
				cargo(me)._atEnd = true;
				--delta;
			}
			host(me) -= delta;
		}
		else
		{
			if (cargo(me)._atEnd) {
				cargo(me)._atEnd = false;
				++delta;
			}
			host(me) -= delta;
		} 
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator -
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TDelta>
	inline ModifiedIterator<THost, ModReverse> &
	operator -= (ModifiedIterator<THost, ModReverse> & me, TDelta delta) {
		if (delta > 0) {
			if (cargo(me)._atEnd) {
				cargo(me)._atEnd = false;
				--delta;
			}
			host(me) += delta;
		} else {
			if (position(host(me)) < -delta) {
				cargo(me)._atEnd = true;
				++delta;
			}
			host(me) -= -delta;
		}
		return me;
	}

	template <typename THost>
	inline typename Difference< ModifiedIterator<THost, ModReverse> >::Type
	operator - (ModifiedIterator<THost, ModReverse> const & a, ModifiedIterator<THost, ModReverse> const & b) {
		typename Difference< ModifiedIterator<THost, ModReverse> >::Type diff = host(b) - host(a);
		if (cargo(a)._atEnd) ++diff;
		if (cargo(b)._atEnd) --diff;
		return diff;
	}

	//////////////////////////////////////////////////////////////////////////////
	// position
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type 
	position(ModifiedIterator<THost, ModReverse> const & me)
	{
	SEQAN_CHECKPOINT
		if (cargo(me)._atEnd)
			return length(container(host(me)));
		else
			return length(container(host(me))) - 1 - position(host(me));
	}

	template <typename THost, typename TContainer>
	inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type 
	position(ModifiedIterator<THost, ModReverse> const & me, TContainer const &cont)
	{
	SEQAN_CHECKPOINT
		if (cargo(me)._atEnd)
			return length(cont);
		else
			return length(cont) - 1 - position(host(me), cont);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator ==
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost>
	inline bool
	operator == (ModifiedIterator<THost, ModReverse> const & a, ModifiedIterator<THost, ModReverse> const & b) {
		return cargo(a)._atEnd == cargo(b)._atEnd && host(a) == host(b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator <
	//////////////////////////////////////////////////////////////////////////////

	// redefinition candidate
	template <typename THost>
	inline bool
	operator < (ModifiedIterator<THost, ModReverse> const & a, ModifiedIterator<THost, ModReverse> const & b) {
		return (!cargo(a)._atEnd && cargo(b)._atEnd) ||
			   (!cargo(a)._atEnd && !cargo(b)._atEnd && host(a) > host(b));
	}

	//////////////////////////////////////////////////////////////////////////////
	// atBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TContainer>
	inline bool
	atBegin(ModifiedIterator<THost, ModReverse> const & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		return position(me, container) == 0;
	}

	template <typename THost>
	inline bool
	atBegin(ModifiedIterator<THost, ModReverse> const & me)
	{
	SEQAN_CHECKPOINT
		return position(me) == 0;
	}

	//////////////////////////////////////////////////////////////////////////////
	// atEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TContainer>
	inline bool
	atEnd(ModifiedIterator<THost, ModReverse> const & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		return cargo(me)._atEnd;
	}

	template <typename THost>
	inline bool
	atEnd(ModifiedIterator<THost, ModReverse> const & me)
	{
	SEQAN_CHECKPOINT
		return cargo(me)._atEnd;
	}


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// reverse string
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost>
	class ModifiedString<THost, ModReverse> {
	public:
		Holder<THost>							data_host;
		typename Cargo<ModifiedString>::Type	data_cargo;

		ModifiedString() {}

		ModifiedString(ModifiedString &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		ModifiedString(ModifiedString const &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		template <typename THostHost, typename THostSpec>
		ModifiedString(ModifiedString<THostHost, THostSpec> &_origin):
			data_host(_origin.data_host) {}

		ModifiedString(THost &_origin) {
			setHost(*this, _origin);
		}

		template <typename T>
		ModifiedString(T & _origin) {
			setValue(*this, _origin);
		}

		template <typename T>
		ModifiedString(T const & _origin) {
			setValue(*this, _origin);
		}

		template <typename T>
		inline ModifiedString const &
		operator = (T & _origin) {
			assign(*this, _origin);
			return *this;
		}

		template <typename TPos>
		inline typename Reference<ModifiedString>::Type 
		operator [] (TPos pos)
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<ModifiedString const>::Type 
		operator [] (TPos pos) const
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};


	template <typename THost>
	struct Iterator< ModifiedString<THost, ModReverse>, Standard > {
		typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, ModReverse> Type;
	};

	template <typename THost>
	struct Iterator< ModifiedString<THost, ModReverse> const, Standard > {
		typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, ModReverse> Type;
	};

	template <typename THost>
	struct DefaultIteratorSpec< ModifiedString<THost, ModReverse> >
	{
		typedef Rooted Type;
	};




	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TPos>
	inline typename Reference<ModifiedString<THost, ModReverse> >::Type 
	value(ModifiedString<THost, ModReverse> & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return value(host(me), (length(host(me)) - 1) - pos);
	}

	template <typename THost, typename TPos>
	inline typename Reference<ModifiedString<THost, ModReverse> const>::Type 
	value(ModifiedString<THost, ModReverse> const & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return value(host(me), (length(host(me)) - 1) - pos);
	}

	//////////////////////////////////////////////////////////////////////////////
	// begin
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TTag >
	inline typename Iterator< ModifiedString<THost, ModReverse> const >::Type 
	begin(ModifiedString<THost, ModReverse> const & me) {
		typename Iterator< ModifiedString<THost, ModReverse> const >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost >
	inline typename Iterator< ModifiedString<THost, ModReverse> >::Type 
	begin(ModifiedString<THost, ModReverse> & me) {
		typename Iterator< ModifiedString<THost, ModReverse> >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type temp_(end(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// end
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost >
	inline typename Iterator< ModifiedString<THost, ModReverse> const >::Type 
	end(ModifiedString<THost, ModReverse> const & me) {
		typename Iterator< ModifiedString<THost, ModReverse> const >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost >
	inline typename Iterator< ModifiedString<THost, ModReverse> >::Type 
	end(ModifiedString<THost, ModReverse> & me) {
		typename Iterator< ModifiedString<THost, ModReverse> >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}

	template < typename THost, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const) {
		typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
		_copyCargo(temp_, me);
		goNext(temp_);
		return temp_;
	}


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// reverseInPlace
	//////////////////////////////////////////////////////////////////////////////

	template < typename TSequence >
	inline void
	reverseInPlace(TSequence & sequence) 
	{
		typedef typename Iterator<TSequence, Standard>::Type	TIter;
		typedef typename Value<TSequence>::Type					TValue;

		TIter it1 = begin(sequence, Standard());
		TIter it2 = it1 + length(sequence) - 1;
		TIter itMid = it1 + length(sequence) / 2;

		for(; it1 != itMid; ++it1, --it2) {
			TValue tmp = *it1;
			*it1 = *it2;
			*it2 = tmp;
		}
	}

	template < typename TSequence, typename TSpec >
	inline void
	reverseInPlace(StringSet<TSequence, TSpec> & stringSet) 
	{
		unsigned seqCount = length(stringSet);
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
			reverseComplementInPlace(stringSet[seqNo]);
	}

}

#endif
