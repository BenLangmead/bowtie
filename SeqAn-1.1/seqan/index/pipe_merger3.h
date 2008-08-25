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
  $Id: pipe_merger3.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_MERGER3_H
#define SEQAN_HEADER_INDEX_MERGER3_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct Merger3;

    template < typename TInput0, typename TInput12 >
    struct Value< Pipe< Bundle2< TInput0, TInput12 >, Merger3 > > {
        typedef typename Size<TInput0>::Type Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // merger3 class
    template < typename TInput0, typename TInput12 >
    struct Pipe< Bundle2< TInput0, TInput12 >, Merger3 >
    {
        typedef typename Value<TInput0>::Type   InType0;
        typedef typename Value<TInput12>::Type  InType12;
        typedef typename Size<Pipe>::Type       SizeType;

        Bundle2 <
            TInput0,
            TInput12 >  in;
        SizeType        N, tmp;
        int             minStream;
        bool            twoStreams;

        Pipe(Bundle2< TInput0, TInput12 > _in):
            in(_in) {}
        
        static inline bool less1(const InType0& i0, const InType12& i12)
        { // lexic. order for pairs
            return (i0.i3[0] <  i12.i3[0] ||
                    i0.i3[0] == i12.i3[0] && i0.i2[0] < i12.i2[1]);
        }

        static inline bool less2(const InType0& i0, const InType12& i12)
        { // and triples
            if (i0.i3[0] < i12.i3[0]) return true;
            if (i0.i3[0] > i12.i3[0]) return false;
            if (i0.i3[1] < i12.i3[1]) return true;
            if (i0.i3[1] > i12.i3[1]) return false;
            if (i0.i2[1] < i12.i2[1]) return true;
            return false;
        }

        inline void getMin() {
            if (twoStreams)
                if ((*in.in2).i1 % 3 == 2 ? less1(*in.in1, *in.in2): less2(*in.in1, *in.in2)) {
                    minStream = 0;
                    tmp = N - (*in.in1).i1;
                } else {
                    minStream = 1;
                    tmp = N - (*in.in2).i1;
                }
            else
                if (minStream >= 0)
                    tmp = minStream? N - (*in.in2).i1: N - (*in.in1).i1;
        }

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }
        
        inline Pipe& operator++() {
            if (minStream) {
                #ifdef SEQAN_TEST_SKEW3
                    InType12 a = *in.in2;
                #endif
                ++in.in2;
                if (eof(in.in2)) {
                    minStream = (twoStreams)? 0: -1;
                    twoStreams = false;
                }
                #ifdef SEQAN_TEST_SKEW3
                    else {
                        InType12 b = *in.in2;
                        SEQAN_ASSERT((a.i3[0] < b.i3[0]) || (a.i3[0] == b.i3[0] && a.i3[1] <= b.i3[1]));
                        SEQAN_ASSERT(a.i2[0] < b.i2[0]);
                    }
                #endif
            } else {
                #ifdef SEQAN_TEST_SKEW3
                    InType0 a = *in.in1;
                #endif
                ++in.in1;
                if (eof(in.in1)) {
                    minStream = (twoStreams)? 1: -1;
                    twoStreams = false;
                }
                #ifdef SEQAN_TEST_SKEW3
                    else {
                        InType0 b = *in.in1;
                        SEQAN_ASSERT((a.i3[0] < b.i3[0]) || (a.i3[0] == b.i3[0] && a.i3[1] <= b.i3[1]));
                        SEQAN_ASSERT((a.i3[0] < b.i3[0]) || (a.i3[0] == b.i3[0] && a.i2[0] < b.i2[0]));
                    }
                #endif
            }
            getMin();
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
	inline bool control(Pipe< TInput, Merger3 > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.twoStreams = !eof(me.in.in1);
        me.minStream = eof(me.in.in2)? -1: 1;
        me.N = length(me);
        me.getMin();
		return true;
	}
    
    template < typename TInput >
    inline typename Size< Pipe< TInput, Merger3 > >::Type
    length(Pipe< TInput, Merger3 > const &me) {
        return length(me.in.in1) + length(me.in.in2);
    }

//}

}

#endif
