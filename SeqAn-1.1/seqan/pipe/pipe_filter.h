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
  $Id: pipe_filter.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_FILTER_H
#define SEQAN_HEADER_PIPE_FILTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{
    
    template <typename InType, typename Result = typename InType::T1>
    struct filterI1 : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const { return x.i1; }
    };

    template <typename InType, typename Result = typename InType::T2>
    struct filterI2 : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const { return x.i2; }
    };

    template <typename InType, typename Result = typename InType::T3>
    struct filterI3 : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const { return x.i3; }
    };


    template < typename TFunctor >
    struct Filter;

	template < typename TInput, typename TFunctor >
    struct Value< Pipe< TInput, Filter<TFunctor> > > {
		typedef typename TFunctor::result_type Type;
	};


/**
.Spec.Filter:
..cat:Pipelining
..general:Class.Pipe
..summary:Applies a specific function to the input stream.
..signature:Pipe<TInput, Filter<TFunctor> >
..param.TInput:The type of the pipeline module this module reads from.
..param.TFunctor:A unary function (see STL's $unary_function$).
...remarks:The argument type of $TFunctor$ must be $VALUE<TInput>::Type$.
..remarks: The output type of this pipe is the result type of $TFunctor$.
*/

	//////////////////////////////////////////////////////////////////////////////
    // filter class
    template <typename TInput, typename TFunctor >
    struct Pipe< TInput, Filter<TFunctor> >
    {
		TInput      &in;
        TFunctor    F;
        
/**
.Memfunc.Filter#Pipe:
..class:Spec.Filter
..summary:Constructor
..signature:Pipe<TInput, Filter<TFunctor> > (in)
..signature:Pipe<TInput, Filter<TFunctor> > (in, func)
..param.in:Reference to an input pipe.
..param.func:A TFunctor object (copy constructor).
*/
        Pipe(TInput& _in):
            in(_in) {}
        
        Pipe(TInput& _in, const TFunctor& _F) :
            in(_in),
            F(_F) {}
        
        inline typename Value<Pipe>::Type const operator*() const {
            return F(*in);
        }

        Pipe& operator++() {
            ++in;
            return *this;
        }
                
    };
    
//}

}

#endif
