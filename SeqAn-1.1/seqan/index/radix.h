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
  $Id: radix.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_RADIX_H
#define SEQAN_HEADER_RADIX_H

namespace SEQAN_NAMESPACE_MAIN
{

	// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText >
	void radixPass(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K)					// alphabet size
	{
		typedef typename Value<TCountArray>::Type TSize;

		TSize i, sum = 0, n = length(a);
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		
		for (i = 0;  i < n;  i++)										// count occurences
			c[ordValue(getValue(r,getValue(a,i)))]++;
			
		for (i = 0;  i < K;  i++)
		{																// exclusive prefix sums
			TSize t = getValue(c,i);
			c[i] = sum;
			sum += t;
		}
		
		for (i = 0;  i < n;  i++)
		{
			TSize j = getValue(a,i);									// sort
			b[c[ordValue(getValue(r,j))]++] = j;
		}
	}

	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText, 
		typename TShift >
	void radixPass(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K,					// alphabet size
		TShift shift)				// shift value
	{
		typedef typename Value<TCountArray>::Type TSize;

		TSize i, sum = 0, n = length(a), sn = length(r);
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		
		for (i = 0;  i < n;  i++)
		{
			TSize j = getValue(a,i) + shift;							// count occurences
			if (j < sn) c[ordValue(getValue(r,j))]++;
			else        sum++;
		}
		
		for (i = 0;  i < K;  i++)
		{																// exclusive prefix sums
			TSize t = getValue(c,i);
			c[i] = sum;
			sum += t;
		}
		
		for (i = 0, sum = 0;  i < n;  i++)
		{																// sort
			TSize j = getValue(a,i);
			TSize k = j + shift;
			if (k < sn) 
				b[c[ordValue(getValue(r,k))]++] = j;	// On Exception: Make sure, you have resized your sufarray
			else
				b[sum++] = j;							// before calling createSuffixArray(..) to length(text)?
		}
	}

    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText >
	void radixExtend(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K)					// alphabet size
	{
		typedef typename Value<TCountArray>::Type TSize;

		TSize i, sum = 0, n = length(a);
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		
		for (i = 0;  i < n;  i++)										// count occurences
			c[ordValue(getValue(r,getValue(a,i)-1))]++;
			
		for (i = 0;  i < K;  i++) 
		{																// exclusive prefix sums
			TSize t = getValue(c,i);
			c[i] = sum;
			sum += t;
		}
		
		for (i = 0;  i < n;  i++)
		{
			TSize j = getValue(a,i);									// sort
			b[c[ordValue(getValue(r,j-1))]++] = j - 1;
		}
	}

    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText >
	void radixExtendClip(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K)					// alphabet size
	{
		typedef typename Value<TCountArray>::Type TSize;

		TSize i, sum = 0, n = length(a);
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		
		for (i = 0;  i < n;  i++)
		{																// count occurences
	        TSize j = getValue(a,i);
			if (j > 0) c[ordValue(getValue(r,j-1))]++;
		}
		
		for (i = 0;  i < K;  i++)
		{																// exclusive prefix sums
			TSize t = getValue(c,i);
			c[i] = sum;
			sum += t;
		}
		
		for (i = 0;  i < n;  i++)
		{
			TSize j = getValue(a,i);									// sort
			if (j > 0) b[c[ordValue(getValue(r,j-1))]++] = j - 1;
		}
	}

}

#endif
