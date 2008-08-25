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
  $Id: pipe_counter.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_COUNTER_H
#define SEQAN_HEADER_PIPE_COUNTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct Counter;

	template < typename TInput >
    struct Value< Pipe< TInput, Counter > > {
		typedef Pair<
			typename Value<TInput>::Type,
			typename Size<TInput>::Type,
			Compressed
		> Type;
	};


/**
.Spec.Counter:
..cat:Pipelining
..general:Class.Pipe
..summary:Extends the input stream by a second field which enumerates the elements.
..signature:Pipe<TInput, Counter>
..param.TInput:The type of the pipeline module this module reads from.
..remarks:The output type is a @Class.Pair@ of input type and size type (i.e. $Pair<Value<TInput>::Type, Size<TInput>::Type>$).
..remarks:The first output field is the original input stream.
..remarks:The second output field begins with 0 and increases by 1 per element.
*/

    //////////////////////////////////////////////////////////////////////////////
    // counter class
    template < typename TInput >
    struct Pipe< TInput, Counter >
    {
		TInput                      &in;
        typename Value<Pipe>::Type	tmp;
        
        Pipe(TInput& _in):
            in(_in) {}
        
        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = *in;
            return tmp;
        }
        
        inline Pipe& operator++() {
            ++in;
            ++tmp.i2;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
	inline bool control(Pipe< TInput, Counter > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.tmp.i2 = 0;
		return true;
	}
    
//}

}

#endif
