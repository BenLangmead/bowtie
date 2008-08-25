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
  $Id: pipe_shifter.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_SHIFTER_H
#define SEQAN_HEADER_PIPE_SHIFTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

//    template < int delta, bool omitBlank = false, bool _echoing = (delta < 0) >
//    struct Shifter;
    template < int delta, bool omitBlank = false, bool _echoing = true >
    struct Shifter;


/**
.Spec.Shifter:
..cat:Pipelining
..general:Class.Pipe
..summary:Shifts the input stream by $delta$ elements.
..signature:Pipe<TInput, Shifter<delta, omitBlank> >
..param.TInput:The type of the pipeline module this module reads from.
..param.delta:The shift size. For the output stream holds $out[i]=in[i+delta]$.
...remarks:For $delta>0$ the input stream is cut of at the beginning and for $delta<0$ at the end.
..param.omitBlank:Omit undefined entries.
..param.omitBlank:If $true$, the output stream is $|delta|$ elements shorter than the input stream.
..param.omitBlank:If $false$, the lengths are equal and blanks (default constructed elements) are inserted on the cut-off-opposite side.
...default:$false$
..remarks:The output type equals the input type.
*/

    //////////////////////////////////////////////////////////////////////////////
    // echoer class
    template < typename TInput, int delta, bool omitBlank >
    struct Pipe< TInput, Shifter<delta, omitBlank, true> >
    {
        TInput                      &in;
        typename Size<Pipe>::Type	blankCounter, charCounter;
		typename Value<Pipe>::Type	blank;
        
        Pipe(TInput& _in):
            in(_in),
			blank()	{}

        inline typename Value<Pipe>::Type const & operator*() {
			if (blankCounter)	return blank;
			else				return *in;
        }

        inline Pipe& operator++() {
			if (blankCounter)
				--blankCounter;
			else {
				++in;
				--charCounter;
			}
            return *this;
        }
    };


	template < typename TInput, int delta, bool omitBlank >
    struct Pipe< TInput, Shifter<delta, omitBlank, false> >
    {
        TInput                      &in;
        typename Size<Pipe>::Type	blankCounter, charCounter;
		typename Value<Pipe>::Type	blank;
        
        Pipe(TInput& _in):
            in(_in),
			blank()	{}

        inline typename Value<Pipe>::Type const & operator*() {
			if (charCounter)	return *in;
			else				return blank;
        }

        inline Pipe& operator++() {
			if (charCounter) {
				++in;
				--charCounter;
			} else
				--blankCounter;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, int delta, bool omitBlank >
	inline bool control(Pipe< TInput, Shifter< delta, omitBlank, false > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
		for(typename Size<TInput>::Type i = 0; i < delta && !eof(me.in); ++i)
			++me.in;
		me.blankCounter = (omitBlank)? 0: delta;
		me.charCounter = length(me.in) - delta;
		return !eof(me.in);
	}

    template < typename TInput, int delta, bool omitBlank >
	inline bool control(Pipe< TInput, Shifter< delta, omitBlank, true > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
		me.blankCounter = (omitBlank)? 0: -delta;
		me.charCounter = length(me.in) + delta;
		return true;
	}




    template < typename TInput, int delta, bool _echoing >
    inline Size< Pipe< TInput, Shifter< delta, true, _echoing > > >
    length(Pipe< TInput, Shifter< delta, true, _echoing > > const &me) {
        return length(me.in) - abs(delta);
    }

    template < typename TInput, int delta, bool omitBlank, bool _echoing >
	inline bool control(Pipe< TInput, Shifter< delta, omitBlank, _echoing > > &me, ControlEof const &command) {
		return me.charCounter == 0 && me.blankCounter == 0;
    }

    template < typename TInput, int delta, bool omitBlank, bool _echoing >
	inline bool control(Pipe< TInput, Shifter< delta, omitBlank, _echoing > > &me, ControlEos const &command) {
		return control(me, ControlEof());
    }

//}

}

#endif
