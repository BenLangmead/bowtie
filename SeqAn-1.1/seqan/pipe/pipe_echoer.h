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
  $Id: pipe_echoer.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_ECHOER_H
#define SEQAN_HEADER_PIPE_ECHOER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    //////////////////////////////////////////////////////////////////////////////
	// some metaprogramming to unrool fixed-size loops
    struct _EchoerFillWorker {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.tmp.i2[I-1] = *(arg.in); ++(arg.in);
        }
    };
    
    struct _EchoerClearWorker {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
			arg.i2[I] = typename Value< typename Value<Arg, 2>::Type >::Type ();
        }
    };
    
    struct _EchoerShiftWorker {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.i2[I] = arg.i2[I-1];
        }
    };


    template < unsigned echoRepeats, bool omitFirst >
    struct Echoer;

    template < typename TInput, unsigned echoRepeats, bool omitFirst >
    struct Value< Pipe< TInput, Echoer< echoRepeats, omitFirst > > > {
        typedef Tuple<typename Value<TInput>::Type, echoRepeats>	EchoType;
        typedef Pair<typename Size<TInput>::Type, EchoType>			Type;
    };


/**
.Spec.Echoer:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs tuples of the $echoRepeats$ last elements of the input stream.
..signature:Pipe<TInput, Echoer<echoRepeats, omitFirst> >
..param.TInput:The type of the pipeline module this module reads from.
..param.echoRepeats:The tuple length.
...remarks:The tuples contain elements $in[i]in[i-1]...in[i-(echoRepeats-1)]$.
..param.omitFirst:Omit half filled tuples.
..param.omitFirst:If $true$, the output stream is $echoRepeats-1$ elements shorter than the input stream.
..param.omitFirst:If $false$, the lengths are identical and the tuple is filled with blanks (default constructed elements) for undefined entries.
..remarks:The output type is a @Class.Tuple@ of input elements and length $echoRepeats$ (i.e. $Tuple<Value<TInput>::Type, echoRepeats>$).
..remarks:The tuples are sequences of the form $in[i]in[i-1]in[i-2]..in[i-echoRepeats+1]$. For $omitFirst=false$ $i$ begins with 0 and for $omitFirst=true$ $i$ begins with $echoRepeats-1$.
*/

    //////////////////////////////////////////////////////////////////////////////
    // echoer class
    template < typename TInput, unsigned echoRepeats, bool omitFirst >
    struct Pipe< TInput, Echoer<echoRepeats, omitFirst> >
    {
        TInput                      &in;
        typename Value<Pipe>::Type	tmp;
        
        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			++in;
            if (eof(in)) return *this;
            LOOP_REVERSE<_EchoerShiftWorker, echoRepeats - 1>::run(this->tmp);
			++tmp.i1;
            tmp.i2[0] = *in;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned echoRepeats, bool omitFirst >
	inline bool control(Pipe< TInput, Echoer< echoRepeats, omitFirst > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.tmp.i1 = 0;
        LOOP<_EchoerClearWorker, echoRepeats - 1>::run(me.tmp);
        if (!eof(me.in)) me.tmp.i2[0] = *me.in;
		return true;
	}
    
    template < typename TInput, unsigned echoRepeats >
    inline bool control(Pipe< TInput, Echoer< echoRepeats, true > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command) || size(me.in) < echoRepeats - 1) return false;
        me.tmp.i1 = 0;
        LOOP_REVERSE<_EchoerFillWorker, echoRepeats - 1>::run(me);
        if (!eof(me.in)) me.tmp.i2[0] = *me.in;
		return true;
    }

    template < typename TInput, unsigned echoRepeats >
    inline Size< Pipe< TInput, Echoer< echoRepeats, true > > >
    length(Pipe< TInput, Echoer< echoRepeats, true > > const &me) {
        return length(me.in) - (echoRepeats - 1);
    }

//}

}

#endif
