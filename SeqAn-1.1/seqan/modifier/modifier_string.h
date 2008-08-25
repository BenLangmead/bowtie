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
  $Id: modifier_string.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_STRING_H
#define SEQAN_HEADER_MODIFIER_STRING_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Class.ModifiedString:
..summary:Allows to modify arbitrary strings by specializing what differs from an origin.
..cat:Modifier
..signature:ModifiedString<THost[, TSpec]>
..param.THost:Original sequence type.
...type:Concept.Container
..param.TSpec:The modifier type.
...metafunction:Metafunction.Spec
..implements:Concept.Container
..remarks:$THost$ can also be a modified string, so you can create custom strings by combining predefined ones.
*/

	template < typename THost, typename TSpec = void >
	class ModifiedString {
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
	};

	template < typename THost, typename TSpec >
	struct Spec< ModifiedString<THost, TSpec> > {
		typedef TSpec Type;
	};

	template < typename THost, typename TSpec >
	struct Spec< ModifiedString<THost, TSpec> const > {
		typedef TSpec Type;
	};


	// use Value, GetValue, Reference, Size, ... from corresponding iterator
	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> >:
		Value< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> const >:
		Value< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};


	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> >:
		GetValue< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> const >:
		GetValue< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};


	template < typename THost, typename TSpec >
	struct Reference< ModifiedString<THost, TSpec> >:
		Reference< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

	template < typename THost, typename TSpec >
	struct Reference< ModifiedString<THost, TSpec> const >:
		Reference< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};

///.Metafunction.Size.param.T.type:Class.ModifiedString

	template < typename THost, typename TSpec >
	struct Size< ModifiedString<THost, TSpec> >:
		Size< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

	template < typename THost, typename TSpec >
	struct Position< ModifiedString<THost, TSpec> >:
		Position< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

	template < typename THost, typename TSpec >
	struct Difference< ModifiedString<THost, TSpec> >:
		Difference< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};


///.Metafunction.Iterator.param.T.type:Class.ModifiedString

	template <typename THost, typename TSpec>
	struct Iterator< ModifiedString<THost, TSpec>, Standard > {
		typedef ModifiedIterator<typename Iterator<THost, Standard>::Type, TSpec> Type;
	};

	template <typename THost, typename TSpec >
	struct Iterator< ModifiedString<THost, TSpec> const, Standard > {
		typedef ModifiedIterator<typename Iterator<THost const, Standard>::Type, TSpec> Type;
	};

	template <typename THost, typename TSpec>
	struct Iterator< ModifiedString<THost, TSpec>, Rooted > {
		typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, TSpec> Type;
	};

	template <typename THost, typename TSpec >
	struct Iterator< ModifiedString<THost, TSpec> const, Rooted > {
		typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, TSpec> Type;
	};


///.Metafunction.Host.param.T.type:Class.ModifiedString

	template < typename THost, typename TSpec >
	struct Host< ModifiedString<THost, TSpec> > {
		typedef THost Type;
	};

	template < typename THost, typename TSpec >
	struct Host< ModifiedString<THost, TSpec> const > {
		typedef THost const Type;
	};


///.Metafunction.IsSequence.param.T.type:Class.ModifiedString

	template < typename THost, typename TSpec >
	struct IsSequence< ModifiedString<THost, TSpec> > {
		typedef True Type;
		enum { VALUE = true };
	};

	template < typename THost, typename TSpec >
	struct AllowsFastRandomAccess< ModifiedString<THost, TSpec> >:
		AllowsFastRandomAccess<THost> {};

	//////////////////////////////////////////////////////////////////////////////
	// host interface
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline Holder<THost> &
	_dataHost(ModifiedString<THost, TSpec> & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_host;
	}
	
	template <typename THost, typename TSpec>
	inline Holder<THost> const &
	_dataHost(ModifiedString<THost, TSpec> const & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_host;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedString<THost, TSpec> >::Type >::Type
	cargo(ModifiedString<THost, TSpec> & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_cargo;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedString<THost, TSpec> const>::Type >::Type
	cargo(ModifiedString<THost, TSpec> const & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_cargo;
	}

	//////////////////////////////////////////////////////////////////////////////
	// assign
	//////////////////////////////////////////////////////////////////////////////

	template <typename TDest, typename TSource>
	inline void _copyCargoImpl(TDest &, TSource &, False const) {}
    
	template <typename TDest, typename TSource>
	inline void _copyCargoImpl(TDest & me, TSource & _origin, True const) {
		cargo(me) = cargo(_origin);
	}
    
	template <typename TDest, typename TSource>
	inline void _copyCargo(TDest & me, TSource & _origin) {
		_copyCargoImpl(me, _origin, typename _IsSameType<
				typename Cargo<TDest>::Type, 
				typename Cargo<TSource>::Type >::Type());
	}
    
    
	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, ModifiedString<THost2, TSpec> & _origin) {
		host(me) = host(_origin);
		_copyCargo(me, _origin);
		return me;
	}

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, ModifiedString<THost2, TSpec> const & _origin) {
		host(me) = host(_origin);
		_copyCargo(me, _origin);
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, T & _origin) {
		host(me) = _origin;
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, T const & _origin) {
		host(me) = _origin;
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////

    template <typename THost, typename TSpec, typename TPos>
	inline typename Reference<ModifiedString<THost, TSpec> >::Type 
	value(ModifiedString<THost, TSpec> & me, TPos pos)
	{
		return value(begin(me, Standard()) + pos);
	}

	template <typename THost, typename TSpec, typename TPos>
	inline typename Reference<ModifiedString<THost, TSpec> const >::Type 
	value(ModifiedString<THost, TSpec> const & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return value(begin(me, Standard()) + pos);
	}

	//////////////////////////////////////////////////////////////////////////////
	// setValue
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline ModifiedString<THost, TSpec> const &
	setValue(ModifiedString<THost, TSpec> & me, ModifiedString<THost, TSpec> const & _origin) {
		setHost(me, host(_origin));
		_copyCargo(me, _origin);
		return me;
	}

	template <typename THost, typename TSpec>
	inline ModifiedString<THost, TSpec> const &
	setValue(ModifiedString<THost, TSpec> & me, ModifiedString<THost, TSpec> & _origin) {
		setHost(me, host(_origin));
		_copyCargo(me, _origin);
		return me;
	}

	// pass _origin to parent modifier
	template <typename THost, typename THostSpec, typename TSpec, typename THost2>
	inline ModifiedString< ModifiedString<THost, THostSpec>, TSpec> const &
	setValue(
		ModifiedString< ModifiedString<THost, THostSpec>, TSpec> & me, 
		THost2 const & _origin) 
	{
		setValue(host(me), _origin);
		return me;
	}

	template <typename THost, typename THostSpec, typename TSpec, typename THost2>
	inline ModifiedString< ModifiedString<THost, THostSpec>, TSpec> const &
	setValue(
		ModifiedString< ModifiedString<THost, THostSpec>, TSpec> & me, 
		THost2 & _origin) 
	{
		setValue(host(me), _origin);
		return me;
	}

	// set _origin at the innermost modifier
	template <typename THost, typename TSpec>
	inline ModifiedString<THost, TSpec> const &
	setValue(ModifiedString<THost, TSpec> & me, THost const & _origin) {
		assignHost(me, _origin);
		return me;
	}

	template <typename THost, typename TSpec>
	inline ModifiedString<THost, TSpec> const &
	setValue(ModifiedString<THost, TSpec> & me, THost & _origin) {
		setHost(me, _origin);
		return me;
	}

	// allow _origin conversion at the innermost modifier
	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedString<THost, TSpec> const &
	setValue(ModifiedString<THost, TSpec> & me, THost2 & _origin) {
		assignHost(me, _origin);
		return me;
	}

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedString<THost, TSpec> const &
	setValue(ModifiedString<THost, TSpec> & me, THost2 const & _origin) {
		assignHost(me, _origin);
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// length
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec >
	inline typename Size< ModifiedString<THost, TSpec> >::Type 
	length(ModifiedString<THost, TSpec> const & me) {
		return length(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// begin
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
	begin(ModifiedString<THost, TSpec> const & me) {
		typename Iterator< ModifiedString<THost, TSpec> const >::Type temp_(begin(host(me)));
		_copyCargo(temp_, me);
		return temp_;
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
	begin(ModifiedString<THost, TSpec> & me) {
		typename Iterator< ModifiedString<THost, TSpec> >::Type temp_(begin(host(me)));
		_copyCargo(temp_, me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_) {
		typename Iterator< 
			ModifiedString<THost, TSpec> const, 
			Tag<TTagSpec> const 
		>::Type temp_(begin(host(me), tag_));
		_copyCargo(temp_, me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_) {
		typename Iterator< 
			ModifiedString<THost, TSpec>, 
			Tag<TTagSpec> const 
		>::Type temp_(begin(host(me), tag_));
		_copyCargo(temp_, me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// end
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
	end(ModifiedString<THost, TSpec> const & me) {
		typename Iterator< ModifiedString<THost, TSpec> const >::Type temp_(end(host(me)));
		_copyCargo(temp_, me);
		return temp_;
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
	end(ModifiedString<THost, TSpec> & me) {
		typename Iterator< ModifiedString<THost, TSpec> >::Type temp_(end(host(me)));
		_copyCargo(temp_, me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_) {
		typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type temp_(end(host(me), tag_));
		_copyCargo(temp_, me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_) {
		typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type temp_(end(host(me), tag_));
		_copyCargo(temp_, me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// stream operators
	//////////////////////////////////////////////////////////////////////////////

	template < typename TStream, typename THost, typename TSpec >
	inline TStream &
	operator << (TStream & target, ModifiedString<THost, TSpec> const & source)
	{
	SEQAN_CHECKPOINT
		write(target, source);
		return target;
	}

	//////////////////////////////////////////////////////////////////////////////

	template < typename TStream, typename THost, typename TSpec >
	inline TStream &
	operator >> (TStream & source, ModifiedString<THost, TSpec> & target)
	{
	SEQAN_CHECKPOINT
		read(source, target);
		return source;
	}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
inline void const *
id(ModifiedString<THost, TSpec> & me) 
{
SEQAN_CHECKPOINT
	return id(host(me));
}
template <typename THost, typename TSpec>
inline void const *
id(ModifiedString<THost, TSpec> const & me) 
{
SEQAN_CHECKPOINT
	return id(host(me));
}

}

#endif
