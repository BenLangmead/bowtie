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
  $Id: pipe_namer.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_NAMER_H
#define SEQAN_HEADER_PIPE_NAMER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template < typename TCompare >
    struct Namer;

	template < typename TInput, typename TCompare >
    struct Value< Pipe< TInput, Namer<TCompare> > > {
        typedef Pair<
            typename Value<TInput>::Type::T1,
			typename Size<TInput>::Type,
			Compressed
		> Type;
    };


/**
.Spec.Namer:
..cat:Pipelining
..general:Class.Pipe
..summary:Extends the input stream by a second field which names the elements.
..signature:Pipe<TInput, Namer<TCompare> >
..param.TInput:The type of the pipeline module this module reads from.
..param.TCompare:A binary function (see STL's $binary_function$) with result type $int$.
...remarks:Should return $0$ if and only if two elements are equal.
..remarks:The output type is a @Class.Pair@ of input type and size type (i.e. $Pair<Value<TInput>::Type, Size<TInput>::Type>$).
..remarks:The first output field is the original input stream.
..remarks:The second output field is the name. This field begins with 0 and increases by 1 for every distinct element. Two elements gets the same name, if and only if they are equal.
*/

    //////////////////////////////////////////////////////////////////////////////
    // namer class
    template < typename TInput, typename TCompare >
    struct Pipe< TInput, Namer<TCompare> >
    {
		TInput                          &in;
        TCompare                        C;
        typename Value<Pipe>::Type      tmp;
        typename Value<TInput>::Type    last;

/**
.Memfunc.Namer#Pipe:
..class:Spec.Namer
..summary:Constructor
..signature:Pipe<TInput, Namer<TCompare> > (in)
..signature:Pipe<TInput, Namer<TCompare> > (in, comp)
..param.in:Reference to an input pipe.
..param.comp:A $TCompare$ object (copy constructor).
*/
        Pipe(TInput& _in):
            in(_in) {}
        
        Pipe(TInput& _in, const TCompare& CC) :
            in(_in),
            C(CC) {}
        
        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = getValueI1(*in);
            return tmp;
        }

        inline Pipe& operator++() {
            ++in;
            if (!eof(in) && C(last, *in) != 0) {
                #ifdef SEQAN_TEST
                    SEQAN_ASSERT(C(last, *in) < 0);
                #endif
                last = *in;
                ++tmp.i2;
            }
			return *this;
        }

        bool unique() const {
            return tmp.i2 == (length(in) - 1);
        }
    };
    

    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, typename TCompare >
	inline bool control(Pipe< TInput, Namer<TCompare> > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        if (!eof(me.in)) {
            me.last = *me.in;
            me.tmp.i1 = me.last.i1;
        }
        me.tmp.i2 = 0;
		return true;
	}
    
//}

}

#endif
