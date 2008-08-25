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
  $Id: pipe_iterator.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_ITERATOR_H
#define SEQAN_HEADER_PIPE_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	//////////////////////////////////////////////////////////////////////////////
	// input pipe iterator
	// interface between pipe modules and algorithms working with iterators

    template <typename TInput>
    struct IPipeIterator
    {
		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef IPipeIterator					iterator;
		typedef ::std::forward_iterator_tag		iterator_category;
		typedef typename Value<TInput>::Type	value_type;
		typedef const value_type			    const_value_type;
		typedef typename Size<TInput>::Type		difference_type;
		typedef value_type*						pointer;
		typedef value_type&						reference;
		typedef const value_type&				const_reference;
		
        TInput*			in;
		difference_type	rest;
        
	    IPipeIterator():
			in(NULL),
			rest(0) {}

	    IPipeIterator(TInput &_in):
			in(&_in),
			rest(length(_in))
		{
			beginRead(*in);
		}

		IPipeIterator(const iterator &I):
			in(I.in),
			rest(I.rest) {}

		const_value_type operator* () const {
			return **in;
		}
    
		//const_reference operator* () const {
		//	return **in;
		//}
    
		iterator& operator++ () {
			++(*in); 
			if (!--rest) endRead(*in);
			return *this;
		}

		iterator operator++ (int) {
			iterator before = *this;
			++*this;
			return before;
		}

		iterator operator+ (difference_type delta) {
			for(; delta != 0; --delta)
				++*this;
			return *this;
		};
		
        difference_type operator- (const iterator &I) const {
			return I.rest - rest;
		};
		
		bool _Equal(const iterator &I) const {
			return rest == I.rest;
		}
	};
    
    template <typename TInput>
	bool operator==(const IPipeIterator<TInput>& _Left, const IPipeIterator<TInput>& _Right)
	{
		return _Left._Equal(_Right);
	}

    template <typename TInput>
	bool operator!=(const IPipeIterator<TInput>& _Left, const IPipeIterator<TInput>& _Right)
	{
		return !(_Left == _Right);
	}



	//////////////////////////////////////////////////////////////////////////////
	// output pipe iterator
	// interface between algorithms working with iterators and pipe modules 

    template <typename TOutput>
    struct OPipeIterator
    {        
		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef OPipeIterator					iterator;
		typedef ::std::forward_iterator_tag		iterator_category;
		typedef typename Value<TOutput>::Type	value_type;
		typedef typename Size<TOutput>::Type	difference_type;
		typedef iterator*						pointer;
		typedef iterator&						reference;
		typedef value_type const &				const_reference;
		
        TOutput*		out;
		difference_type	rest;
        
	    OPipeIterator():
			out(NULL),
			rest(0) {}

	    OPipeIterator(TOutput &_out):
			out(&_out),
			rest(length(_out))
		{
			beginWrite(*out);
		}

		OPipeIterator(const iterator &I):
			out(I.out),
			rest(I.rest) {}

		reference operator* () {
			return *this;
		}

		iterator operator= (const_reference _v) {
			out->push(_v);
		}
    
		iterator& operator++ () {
			++(*out);
			if (!--rest) endWrite(*out);
			return *this;
		}

		iterator operator++ (int) {
			iterator before = *this;
			++*this;
			return before;
		}

        difference_type operator- (const iterator &I) const {
			return I.rest - rest;
		};
		
		bool _Equal(const iterator &I) const {
			return rest == I.rest;
		}
	};

    template <typename TOutput>
	bool operator==(const OPipeIterator<TOutput>& _Left, const OPipeIterator<TOutput>& _Right)
	{
		return _Left._Equal(_Right);
	}

    template <typename TOutput>
	bool operator!=(const OPipeIterator<TOutput>& _Left, const OPipeIterator<TOutput>& _Right)
	{
		return !(_Left == _Right);
	}


    template < typename TInput, typename TSpec, typename TIteratorSpec >
	struct Iterator< Pipe< TInput, TSpec >, TIteratorSpec> {
		typedef IPipeIterator< Pipe< TInput, TSpec > > Type;
	};

	template < typename TInput >
	struct Value< IPipeIterator< TInput > > {
		typedef typename Value<TInput>::Type Type;
	};

	template < typename TInput >
	struct Size< IPipeIterator< TInput > > {
		typedef typename Size<TInput>::Type Type;
	};

	template < typename TOutput >
	struct Value< OPipeIterator< TOutput > > {
		typedef typename Value<TOutput>::Type Type;
	};

	template < typename TOutput >
	struct Size< OPipeIterator< TOutput > > {
		typedef typename Size<TOutput>::Type Type;
	};


    template < typename TInput, typename TSpec >
	IPipeIterator< Pipe< TInput, TSpec > >
	begin(Pipe< TInput, TSpec > &pipe) {
		return IPipeIterator< Pipe< TInput, TSpec > >(pipe);
	}

    template < typename TInput, typename TSpec >
	IPipeIterator< Pipe< TInput, TSpec > >
	end(Pipe< TInput, TSpec > &pipe) {
		return IPipeIterator< Pipe< TInput, TSpec > >();
	}


	template < typename TInput >
    inline typename Difference<TInput>::Type
    difference(IPipeIterator<TInput> first, IPipeIterator<TInput> last) {
        return last - first;
    }

	template < typename TOutput >
    inline typename Difference<TOutput>::Type
    difference(OPipeIterator<TOutput> first, OPipeIterator<TOutput> last) {
        return last - first;
    }

//}

}

#endif
