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
  $Id: modifier_view.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_VIEW_H
#define SEQAN_HEADER_MODIFIER_VIEW_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ModView:
..summary:Transforms the characters of the $THost$ string/iterator using a custom function.
..cat:Modifier
..general:Class.ModifiedIterator
..general:Class.ModifiedString
..signature:ModifiedIterator<THost, ModView<TFunctor> >
..signature:ModifiedString<THost, ModView<TFunctor> >
..param.THost:Original string/iterator.
...type:Concept.Iterator
..param.TFunctor:A unary function (see STL's $unary_function$).
...remarks:The argument type of $TFunctor$ must be $VALUE<THost>::Type$.
..remarks:The @Metafunction.Value@ type of this modifier is the result type of $TFunctor$.
*/

	template <typename TFunctor>
	struct ModView {};

	template <typename TFunctor>
	struct ModViewCargo {
		TFunctor	func;
	};


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// view iterator
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost, typename TFunctor>
	struct Cargo< ModifiedIterator<THost, ModView<TFunctor> > > {
		typedef ModViewCargo<TFunctor>	Type;
	};

	template <typename THost, typename TFunctor>
	class ModifiedIterator<THost, ModView<TFunctor> > {
	public:
		Holder<THost, Simple>					data_host;
		typename Cargo<ModifiedIterator>::Type	data_cargo;

		mutable typename Value<ModifiedIterator>::Type	tmp_value;

		ModifiedIterator() {}

		explicit ModifiedIterator(TFunctor &_func) {
			assignModViewFunctor(*this, _func);
		}

		explicit ModifiedIterator(TFunctor const &_func) {
			assignModViewFunctor(*this, _func);
		}

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

	template <typename THost, typename TFunctor>
	struct Value< ModifiedIterator<THost, ModView<TFunctor> > > {
		typedef typename TFunctor::result_type	Type;
	};

	template <typename THost, typename TFunctor>
	struct GetValue< ModifiedIterator<THost, ModView<TFunctor> > >:
		Value< ModifiedIterator<THost, ModView<TFunctor> > > {};

	template <typename THost, typename TFunctor>
	struct Reference< ModifiedIterator<THost, ModView<TFunctor> > > {
		typedef typename Value< ModifiedIterator<THost, ModView<TFunctor> > >::Type & Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TFunctor>
	inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > >::Type 
	value(ModifiedIterator<THost, ModView<TFunctor> > & me)
	{
	SEQAN_CHECKPOINT
		me.tmp_value = cargo(me).func(value(host(me)));
		return me.tmp_value;
	}

	template <typename THost, typename TFunctor>
	inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > const>::Type 
	value(ModifiedIterator<THost, ModView<TFunctor> > const & me)
	{
	SEQAN_CHECKPOINT
		me.tmp_value = cargo(me).func(value(host(me)));
		return me.tmp_value;
	}

	//////////////////////////////////////////////////////////////////////////////
	// assignModViewFunctor
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TFunctor>
	inline void
	assignModViewFunctor(ModifiedIterator<THost, ModView<TFunctor> > & me, TFunctor const & _func) 
	{
	SEQAN_CHECKPOINT
		cargo(me).func = _func;
	}


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// view string
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost, typename TFunctor>
	struct Cargo< ModifiedString<THost, ModView<TFunctor> > > {
		typedef ModViewCargo<TFunctor>	Type;
	};

	template <typename THost, typename TFunctor>
	class ModifiedString<THost, ModView<TFunctor> > {
	public:
		Holder<THost>							data_host;
		typename Cargo<ModifiedString>::Type	data_cargo;

		mutable typename Value<ModifiedString>::Type	tmp_value;

		ModifiedString() {}

		explicit ModifiedString(TFunctor &_func) {
			cargo(*this).func = _func;
		}

		explicit ModifiedString(TFunctor const &_func) {
			cargo(*this).func = _func;
		}

		explicit ModifiedString(ModifiedString const &_origin, TFunctor const &_func):
			data_host(_origin.data_host)
		{
			cargo(*this).func = _func;
		}

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


	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TFunctor, typename TPos>
	inline typename Reference<ModifiedString<THost, ModView<TFunctor> > >::Type 
	value(ModifiedString<THost, ModView<TFunctor> > & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		me.tmp_value = cargo(me).func(value(host(me), pos));
		return me.tmp_value;
	}

	template <typename THost, typename TFunctor, typename TPos>
	inline typename Reference<ModifiedString<THost, ModView<TFunctor> > const>::Type 
	value(ModifiedString<THost, ModView<TFunctor> > const & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		me.tmp_value = cargo(me).func(value(host(me), pos));
		return me.tmp_value;
	}

	//////////////////////////////////////////////////////////////////////////////
	// assignModViewFunctor
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TFunctor>
	inline void
	assignModViewFunctor(ModifiedString<THost, ModView<TFunctor> > & me, TFunctor const & _func)
	{
	SEQAN_CHECKPOINT
		cargo(me).func = _func;
	}


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// convertInPlace
	//////////////////////////////////////////////////////////////////////////////

	template < typename TSequence, typename TFunctor >
	inline void
	convertInPlace(TSequence & sequence, TFunctor const &F)
	{
		typedef typename Iterator<TSequence, Standard>::Type	TIter;

		TIter it = begin(sequence, Standard());
		TIter itEnd = end(sequence, Standard());
		for(; it != itEnd; ++it)
			*it = F(*it);
	}

}

#endif
