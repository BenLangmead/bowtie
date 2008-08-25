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
  $Id: pump_separator7.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_PUMP_SEPARATOR_H
#define SEQAN_HEADER_INDEX_PUMP_SEPARATOR_H

namespace SEQAN_NAMESPACE_MAIN
{


    template < 
		typename TInput, typename TFunctor,
		typename TOut1, typename TOut2, typename TOut4 
	>
    static void skew7_separate_slices(
		TInput &in, TFunctor const &funcSlice,
		TOut1 &out1, TOut2 &out2, TOut4 &out4)
    {
		beginRead(in);

		resize(out1, funcSlice.n1);
		resize(out2, funcSlice.n2);
		resize(out4, funcSlice.n4);

		beginWrite(out1);
		beginWrite(out2);
		beginWrite(out4);

		typename Value<TInput>::Type i;
		while (!eof(in)) {
			pop(in, i);
			if (i.i1 < funcSlice.n4) {
				push(out4, i);
			} else 
				if (i.i1 < funcSlice.n24) {
					i.i1 -= funcSlice.n4;
					push(out2, i);
				} else {
					i.i1 -= funcSlice.n24;
					push(out1, i);
				}
		}

		endWrite(out4);
		endWrite(out2);
		endWrite(out1);
		endRead(in);
    }
    
//}

}

#endif
