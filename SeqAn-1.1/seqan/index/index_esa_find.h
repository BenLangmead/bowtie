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
  $Id: index_esa_find.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_ESA_FIND_H
#define SEQAN_HEADER_INDEX_ESA_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ESA finders

/**
.Tag.Index Find Algorithm
..summary:Tag to specify the index search algorithm.
..remarks:These tag can be used to specify the @Function.find@ algorithm 
for @Class.Index@ based substring searches.
..cat:Index

..tag.ESA_FIND_MLR:Binary search with mlr-heuristic.
...remarks:Exact string matching using a suffix array binary search with the mlr-heuristic.

..tag.ESA_FIND_LCPE:Binary search using lcp values.
...remarks:Exact string matching using a suffix array binary search and a lcp-interval tree.

..see:Class.Finder
..see:Spec.Index_ESA
..see:Spec.Index_QGram
*/

	struct _Finder_MLR;		// simple Suffix Array finder with mlr-heuristic
	struct _Finder_LCPE;	// Suffix Array finder using an enhanced LCP-Table

	typedef Tag<_Finder_MLR> const ESA_FIND_MLR;
	typedef Tag<_Finder_LCPE> const ESA_FIND_LCPE;

//____________________________________________________________________________


	template < typename TText, typename TSpec >
	struct DefaultFinder< Index<TText, Index_ESA<TSpec> > > {
        typedef ESA_FIND_MLR Type;	// standard suffix array finder is mlr-heuristic
    };


	//////////////////////////////////////////////////////////////////////////////
	// different layouts of a suffix array or lcp table

	struct SortedList {};			// classical sorted list (suffix array, sorted list, ...)
	struct LeftCompleteTree {};		// flattened search tree root, left child, right child, left child's left child, left child's right child, ...

	template < unsigned BlockSize = 4096 >
	struct BTree {};				// b-tree compacts nodes and its children to blocks of BlockSize

	template < typename TString, typename TSpec >
	class SearchTreeIterator {};

	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//

	template <typename TString>
	class SearchTreeIterator< TString, SortedList >
	{
	public:

		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		inline SearchTreeIterator(TString &string):
			first(begin(string, Standard())),
			count(length(string))
		{
			count2 = count / 2;
			_mid = first;
			goFurther(_mid, count2);
		}
			
        inline const TValue& operator*() const {
			return *_mid;
		}

        inline const TValue* operator->() const {
			return &*_mid;
		}

		inline TSize mid() {
			return count2;
		}

		// descend left
		inline SearchTreeIterator & left()
		{
			count = count2;
			count2 /= 2;
			_mid = first;
			goFurther(_mid, count2);
			return *this;
		}

		// descend right
        inline SearchTreeIterator & right()
		{
			first = ++_mid, count -= count2 + 1;
			count2 = count / 2;
			goFurther(_mid, count2);
			return *this;
		}

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline bool eof() {
            return !count;
        }

		inline operator TIterator & () {
			return _mid;
		}

	private:
		TIterator	first, _mid;
		TSize		count, count2;
	};


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//

	template <typename TString>
	class SearchTreeIterator< TString, LeftCompleteTree > 
	{
	public:
		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		TIterator	it;
		TSize		size;

		inline SearchTreeIterator(TString &string):
			it(begin(string, Standard())),
			size(length(string))
		{
			_left = 0;
			_lSize = 1;
            _iSize = size;
			for(_xSize = 1; _xSize < size; _xSize <<= 1);
            if (!size) _xSize = 0;
		}

		inline SearchTreeIterator():
			it(),
            size(0),
            _xSize(0) {}

        inline const TValue& operator*() const {
			return *it;
		}

        inline const TValue* operator->() const {
			return &*it;
		}

		inline TSize mid() {
			return _xSize >> 1;
		}

		inline TSize leftSize() {
			return mid();
		}

		inline TSize rightSize() {
			return _iSize - mid();
		}

		// descend left
        inline SearchTreeIterator & left()
		{
            if (_xSize <= 1) {
                _xSize = 0;
                return *this;
            }
            _descendLeft();
			_iSize = _xSize;	    // = mid();
            return *this;
        }

		// descend right
        inline SearchTreeIterator & right()
		{
            if (_xSize <= 1) {
                _xSize = 0;
                return *this;
            }
			_iSize -= mid();
            SEQAN_ASSERT(_iSize != 0);    // _xSize/2 is less than _iSize by invariant

            // step down right
            _descendRight();

            // step down left until we have two sons or are a leaf
            while (_iSize <= (_xSize >> 1))
                _descendLeft();

            return *this;
        }

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline SearchTreeIterator& operator--() {
            --it;
            return *this;
        }

        inline SearchTreeIterator operator--(int) {
            SearchTreeIterator before = *this;
            --it;
            return before;
        }

        inline SearchTreeIterator & operator++() {
            ++it;
            return *this;
        }

        inline SearchTreeIterator operator++(int) {
            SearchTreeIterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(SearchTreeIterator const &I) {
            return (_xSize == I._xSize) && (_xSize == 0 || it == I.it);
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool eof() {
            return !_xSize;
        }

		inline operator TIterator & () {
			return it;
		}

	private:
		TSize _left;		// left iterator offset of current interval
		TSize _lSize;		// iterator elements of current level
		TSize _xSize;		// max interval size of current level
		TSize _iSize;		// current interval size

        inline void _descendLeft() 
		{
			goFurther(it, _left + _lSize);
			_left <<= 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }

        inline void _descendRight() 
		{
			goFurther(it, _left + 1 + _lSize);
			_left = (_left << 1) + 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }
    };


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search b-tree like a real tree
	//

	template <typename TString, unsigned BlockSize>
	class SearchTreeIterator< TString, BTree<BlockSize> > 
	{
	public:
		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		enum { BlockHeight = Log2Floor<BlockSize>::VALUE };
		enum { BlockElements = (1 << BlockHeight) - 1 };
		enum { BlockInnerElements = (1 << (BlockHeight - 1)) - 1 };
		enum { BlockLeafs = 1 << (BlockHeight - 1) };

		TIterator	it;
		TSize		size;

		inline SearchTreeIterator(TString &string):
			it(begin(string, Standard())),
			size(length(string))
		{
			//_left = 0;
			//_lSize = 0;
   //         _iSize = size;

			_heightHigh = 1;
			_stepSizeLow = 1;
			for(TSize _xSizeLow = 2; _xSizeLow <= size; _xSizeLow <<= 1) {
				if (_stepSizeLow == BlockLeafs) {
					_stepSizeLow = 1;
					++_heightHigh;
				} else
					_stepSizeLow <<= 1;
			}
			
			_stepSizeLow >>= 1;
			for(_xSizeHigh = 1; _xSizeHigh * BlockSize <= size; _xSizeHigh *= BlockSize);

			_leftLow = (_stepSizeLow << 1) - 1;		// point to the middle
			_leftHigh = 0;
			_lSizeHigh = 1;

			it += _leftLow;

			_updateElements();

			//if (!size) _xSizeLow = 0;
   //         if (!size) _xSizeHigh = 0;
		}

		inline SearchTreeIterator():
			it() {}

        inline const TValue& operator*() const {
			return *it;
		}

		inline const TValue* operator->() const {
			return &*it;
		}

		// descend left
        inline SearchTreeIterator & left() 
		{
            if (_heightHigh == 1 && !_stepSizeLow) {
                _heightHigh = 0;
                return *this;
            }
            _descendLeft();
            return *this;
        }

		// descend right
        inline SearchTreeIterator & right() 
		{
            if (_heightHigh == 1 && !_stepSizeLow) {
				++it;
                _heightHigh = 0;
                return *this;
            }

            // step down right
            _descendRight();

            // step down left until we have two sons or are a leaf
            while (_elements <= _leftLow)
                _descendLeft();

            return *this;
        }

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline SearchTreeIterator& operator--() {
            --it;
            return *this;
        }

        inline SearchTreeIterator operator--(int) {
            SearchTreeIterator before = *this;
            --it;
            return before;
        }

        inline SearchTreeIterator& operator++() {
            ++it;
            return *this;
        }

        inline SearchTreeIterator operator++(int) {
            SearchTreeIterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(SearchTreeIterator const &I) {
			return it == I.it || (eof() && I.eof());
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool eof() {
            return !_heightHigh;
        }

		inline operator TIterator & () {
			return it;
		}

	private:
		TSize _heightHigh;	// height measured in BBlocks
		unsigned _elements;		// elements in current BBlock

		unsigned _leftLow;		// left iterator offset of current interval
		unsigned _stepSizeLow;	// left/right step size of current level

		TSize _leftHigh;		// left BBlock offset of current interval
		TSize _lSizeHigh;	// BBlocks of current level
		TSize _xSizeHigh;	// max BBlocks of current level

		inline void _descendLeft() 
		{
			if (_stepSizeLow) {
				it -= _stepSizeLow;
				_leftLow -= _stepSizeLow;
				_stepSizeLow >>= 1;
			} else
				if (--_heightHigh) {
					_descendBBlock(_leftLow);

					_leftLow = BlockInnerElements;		// point to the middle
					_stepSizeLow = BlockLeafs / 2;		// with step size of half an interval length
				}
        }

		inline void _descendRight() 
		{
			if (_stepSizeLow) {
				it += _stepSizeLow;
				_leftLow += _stepSizeLow;
				_stepSizeLow >>= 1;
			} else
				if (--_heightHigh) {
					_descendBBlock(_leftLow + 1);

					_leftLow = BlockInnerElements;		// point to the middle
					_stepSizeLow = BlockLeafs / 2;		// with step size of half an interval length
				}
        }

		inline void _updateElements() 
		{
			TSize firstElement = (1 + _leftHigh * BlockSize) * _xSizeHigh - 1;
			TSize lastElement = (1 + (_leftHigh + 1) * BlockSize) * _xSizeHigh - 2;

			if (lastElement >= size)
				_elements = (size - firstElement) / _xSizeHigh;
			else
				_elements = BlockElements;
		}

        inline void _descendBBlock(TSize _childIndex) 
		{
			// goFurther to the begin of the current BBlock and further to the destination BBlock
			goFurther(it, BlockSize * (_leftHigh * (BlockSize - 1) + _childIndex + _lSizeHigh) + BlockInnerElements - _leftLow);

			_leftHigh = _leftHigh * BlockSize + _childIndex;
			_xSizeHigh /= BlockSize;
			_lSizeHigh = (size / _xSizeHigh + BlockSize - 1) / BlockSize;

			_updateElements();
        }
    };


	//////////////////////////////////////////////////////////////////////////////
	// substring search with SA table and w/o LCP-table or enhancement
	//

	template <
		typename TText,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_lowerBoundSA(
		TText &text,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query)
	{	// find first element not before query, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;

		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.eof(); )
		{	// divide and conquer, find half that contains answer

			TSuffix		suf = suffix(text, *treeIter);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TDiff		lcp = _min(lcpLower, lcpUpper);

			goFurther(t, lcp);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && *t == *q) {
				++t;
				++q;
				++lcp;
			}
			
            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q)) {
				treeIter.right();
				lcpLower = lcp;
			} else {
				treeIter.left();
				lcpUpper = lcp;
			}
		}
		return treeIter;
	}

	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TText,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_upperBoundSA(
		TText &text,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query)
	{	// find first element that query is before, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;

		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.eof(); )
		{	// divide and conquer, find half that contains answer

			TSuffix		suf = suffix(text, *treeIter);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TDiff		lcp = _min(lcpLower, lcpUpper);

			goFurther(t, lcp);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && *t == *q) {
				++t;
				++q;
				++lcp;
			}
			
            // is text <= query ?
			if (q == qEnd || t == tEnd || !(*q < *t)) {
				treeIter.right();
				lcpLower = lcp;
			} else {
				treeIter.left();
				lcpUpper = lcp;
			}
		}
		return treeIter;
	}

	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TText,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA, Standard>::Type >
	_equalRangeSA(
		TText &text,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query)
	{	// find range equivalent to query, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;

		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.eof(); )
		{	// divide and conquer, check midpoint

			TSuffix		suf = suffix(text, *treeIter);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TDiff		lcp = _min(lcpLower, lcpUpper);

			goFurther(t, lcp);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && *t == *q) {
				++t;
				++q;
				++lcp;
			}
			
            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q))
			{	// range begins above mid, loop
				treeIter.right();
				lcpLower = lcp;
			}
            // is text > query ?
			else if (q != qEnd && (t != tEnd && *q < *t))
			{	// range in first half, loop
				treeIter.left();
				lcpUpper = lcp;
			} else
            // is text == query ?
			{	// range straddles mid, find each end and return
				return Pair<TSAIter> (
					_lowerBoundSA(text, treeIter.leftChild(), query),
					_upperBoundSA(text, treeIter.rightChild(), query)
				);
			}
		}
		return Pair<TSAIter> (treeIter, treeIter);	// empty range
	}


	//////////////////////////////////////////////////////////////////////////////
	// Finder wrappers (return iterators instead of positions)

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	lowerBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	upperBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA const, Standard>::Type >
	equalRangeSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _equalRangeSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}


	//////////////////////////////////////////////////////////////////////////////
	// workarounds for the Visual Studio array problem

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	lowerBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	upperBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA const, Standard>::Type >
	equalRangeSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _equalRangeSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}


	//////////////////////////////////////////////////////////////////////////////
	// wrappers with nice interfaces

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find first element not before query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element not before query, using operator<
		return lowerBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find first element that query is before, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return upperBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA, Standard>::Type > itPair = 
			_equalRangeSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return equalRangeSA(text, sa, query, SortedList());
	}


	//////////////////////////////////////////////////////////////////////////////
	// workarounds for the Visual Studio array problem

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find first element not before query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element not before query, using operator<
		return lowerBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find first element that query is before, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element that query is before, using operator<
		return upperBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA, Standard>::Type > itPair = 
			_equalRangeSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element that query is before, using operator<
		return equalRangeSA(text, sa, query, SortedList());
	}


	//////////////////////////////////////////////////////////////////////////////
	// substring search with enhanced LCP-table
	//

	template <
		typename TText,
		typename TSA,
		typename TLCP,
		typename TSpec,
		typename TQuery,
		typename TDiff_
	>
	inline typename Iterator<TSA, Standard>::Type
	_lowerBoundLCPE(
		TText &text,
		TSA &sa,
		SearchTreeIterator< TLCP, TSpec > treeIter,
		TQuery &query,
		TDiff_ lcpLower,
		TDiff_ lcpUpper)
	{	// find first element not before query, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;
		typedef SearchTreeIterator< TLCP, TSpec >			TLCPTreeIt;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;

		TDiff delta = length(sa) - 1;
		TDiff lcp;
		#ifdef SEQAN_PROFILE_LCPEFIND
			TDiff skippedCompares = 0;	// difference of char compares related to xxx_bound_sa
		#endif

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());
		TSAIter first = begin(sa, Standard());

		// binary search with intervals >= 3 elements
		for (; 1 < delta; )
		{	// divide and conquer, find half that contains answer
			TDiff delta2 = treeIter.leftSize();
			TSAIter mid = first;
			goFurther(mid, delta2);

			if (lcpLower > lcpUpper) 
			{
                TLCPTreeIt leftChild = treeIter;
                leftChild.left();
				TDiff lcpMidLower = *leftChild;

				if (lcpMidLower > lcpLower) 
				{
					// second half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += lcpLower - lcpUpper;
					#endif
					first = mid;
					treeIter.right();
					if ((delta -= delta2) == 1) {
						++first;
						delta = 0;
					}
					continue;
				} 
				else if (lcpMidLower < lcpLower) 
				{
					// first half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += lcpMidLower - lcpUpper;
					#endif
					lcpUpper = lcpMidLower;
					treeIter = leftChild;
					if ((delta = delta2) == 1)
						delta = 0;
					continue;
				}
				lcp = lcpLower;
			} 
			else if (lcpLower < lcpUpper) 
			{
                TLCPTreeIt rightChild = treeIter;
                rightChild.right();
				TDiff lcpMidUpper = *rightChild;

				if (lcpMidUpper > lcpUpper) 
				{
					// first half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += lcpUpper - lcpLower;
					#endif
					treeIter.left();
					if ((delta = delta2) == 1)
						delta = 0;
					continue;
				} 
				else if (lcpMidUpper < lcpUpper) 
				{
					// second half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += lcpMidUpper - lcpLower;
					#endif
					lcpLower = lcpMidUpper;
					first = mid;
                    treeIter = rightChild;
					if ((delta -= delta2) == 1) {
						++first;
						delta = 0;
					}
					continue;
				}
				lcp = lcpUpper;
			} else
				lcp = lcpUpper;


			TSuffix		suf = suffix(text, *mid);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;

			// lcp search changes MIN to MAX here
//			TDiff lcp = _max(lcpLower, lcpUpper);
			#ifdef SEQAN_PROFILE_LCPEFIND
				skippedCompares += lcp - _min(lcpLower, lcpUpper);
			#endif
			goFurther(t, lcp);
			goFurther(q, lcp);
			for(TDiff i = _min(difference(t, tEnd), difference(q, qEnd)); 
				i && *t == *q;
				--i, ++t, ++q, ++lcp);

            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q)) 
			{
				// second half
				lcpLower = lcp;
				first = mid;
				treeIter.right();
				if ((delta -= delta2) == 1) {
					++first;
					delta = 0;
				}
			} else {
				// first half
				lcpUpper = lcp;
				treeIter.left();
				if ((delta = delta2) == 1)
					delta = 0;
			}
		}

		#ifdef SEQAN_PROFILE_LCPEFIND
			SEQAN_PROADD(SEQAN_PROEXTRA3, skippedCompares);
		#endif

        // binary search for intervals of 2 or less elements
        lcp = _min(lcpLower, lcpUpper);
		TQueryIter q = qBegin;
		goFurther(q, lcp);

		while (true) {
			TSuffix		suf = suffix(text, *first);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());

			goFurther(t, lcp);
			for(TDiff i = _min(difference(t, tEnd), difference(q, qEnd)); 
            	i && *t == *q;
            	 --i, ++t, ++q);

            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q)) {
				// second half
				++first;
				if (!delta) return first;
                --delta;
			} else {
				// first half -> end
				return first;
			}
        }
	}

	template <
		typename TText,
		typename TSA,
		typename TLCP,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_lowerBoundLCPE(
		TText &text,
		TSA &sa,
		SearchTreeIterator< TLCP, TSpec > &treeIter,
		TQuery &query)
	{
		return _lowerBoundLCPE(text, sa, treeIter, query, 0, 0);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TText,
		typename TSA,
		typename TLCP,
		typename TSpec,
		typename TQuery,
		typename TDiff_
	>
	inline typename Iterator<TSA, Standard>::Type
	_upperBoundLCPE(
		TText &text,
		TSA &sa,
		SearchTreeIterator< TLCP, TSpec > treeIter,
		TQuery &query,
		TDiff_ lcpLower,
		TDiff_ lcpUpper)
	{	// find first element not before query, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;
		typedef SearchTreeIterator< TLCP, TSpec >			TLCPTreeIt;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;

		TDiff delta = length(sa) - 1;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());
		TSAIter first = begin(sa, Standard());

        // binary search with intervals >= 3 elements
		for (; 1 < delta; )
		{	// divide and conquer, find half that contains answer
			TDiff delta2 = treeIter.leftSize();
			TSAIter mid = first;
			goFurther(mid, delta2);

			if (lcpLower > lcpUpper) 
			{
                TLCPTreeIt leftChild = treeIter;
                leftChild.left();
				TDiff lcpMidLower = *leftChild;

				if (lcpMidLower > lcpLower) 
				{
					// second half
					first = mid;
					treeIter.right();
					if ((delta -= delta2) == 1) {
						++first;
						delta = 0;
					}
					continue;
				} 
				else if (lcpMidLower < lcpLower) 
				{
					// first half
					lcpUpper = lcpMidLower;
					treeIter = leftChild;
					if ((delta = delta2) == 1)
						delta = 0;
					continue;
				}
			} 
			else if (lcpLower < lcpUpper) 
			{
                TLCPTreeIt rightChild = treeIter;
                rightChild.right();
				TDiff lcpMidUpper = *rightChild;

				if (lcpMidUpper > lcpUpper) 
				{
					// first half
					treeIter.left();
					if ((delta = delta2) == 1)
						delta = 0;
					continue;
				}
				else if (lcpMidUpper < lcpUpper) 
				{
					// second half
					lcpLower = lcpMidUpper;
					first = mid;
                    treeIter = rightChild;
					if ((delta -= delta2) == 1) {
						++first;
						delta = 0;
					}
					continue;
				}
			}

			TSuffix		suf = suffix(text, *mid);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;

			// lcp search changes MIN to MAX here
			TDiff lcp = _max(lcpLower, lcpUpper);
			goFurther(t, lcp);
			goFurther(q, lcp);
			TDiff max = _min(difference(t, tEnd), difference(q, qEnd));
            TDiff i = max;
			for(; i && *t == *q; --i, ++t, ++q);
			lcp += max - i;

            // is text <= query ?
			if (q == qEnd || t == tEnd || !(*q < *t)) 
			{
				// second half
				lcpLower = lcp;
				first = mid;
				treeIter.right();
				if ((delta -= delta2) == 1) {
					++first;
					delta = 0;
				}
			} else {
				// first half
				lcpUpper = lcp;
				treeIter.left();
				if ((delta = delta2) == 1)
					delta = 0;
			}
		}

        // binary search for intervals of 2 or less elements
        TDiff lcp = _min(lcpLower, lcpUpper);
		TQueryIter q = qBegin;
		goFurther(q, lcp);

		while (true) {
			TSuffix		suf = suffix(text, *first);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());

			goFurther(t, lcp);
			TDiff i = _min(difference(t, tEnd), difference(q, qEnd));
            for(; i && *t == *q; --i, ++t, ++q);

            // is text <= query ?
			if (q == qEnd || t == tEnd || !(*q < *t)) {
				// second half
				++first;
				if (!delta) return first;
                --delta;
			} else {
				// first half -> end
				return first;
			}
        }
	}

	template <
		typename TText,
		typename TSA,
		typename TLCP,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_upperBoundLCPE(
		TText &text,
		TSA &sa,
		SearchTreeIterator< TLCP, TSpec > &treeIter,
		TQuery &query)
	{
		return _upperBoundLCPE(text, sa, treeIter, query, 0, 0);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TText,
		typename TSA,
		typename TLCP,
		typename TSpec,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA, Standard>::Type >
	_equalRangeLCPE(
		TText &text,
		TSA &sa,
		SearchTreeIterator< TLCP, TSpec > treeIter,
		TQuery &query)
	{	// find first element not before query, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;
		typedef SearchTreeIterator< TLCP, TSpec >			TLCPTreeIt;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;

		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;
        TDiff delta = length(sa) - 1;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());
		TSAIter first = begin(sa, Standard());
		TSAIter last = end(sa, Standard());

        // binary search with intervals >= 3 elements
		for (; 1 < delta; )
		{	// divide and conquer, find half that contains answer
			TDiff delta2 = treeIter.leftSize();
			TSAIter mid = first;
			goFurther(mid, delta2);

			if (lcpLower > lcpUpper) 
			{
                TLCPTreeIt leftChild = treeIter;
                leftChild.left();
				TDiff lcpMidLower = *leftChild;

				if (lcpMidLower > lcpLower) 
				{
					// second half
					first = mid;
					treeIter.right();
					if ((delta -= delta2) == 1) {
						++first;
						delta = 0;
					}
					continue;
				}
				else if (lcpMidLower < lcpLower) 
				{
					// first half
					lcpUpper = lcpMidLower;
					treeIter = leftChild;
					if ((delta = delta2) == 1)
						delta = 0;
					continue;
				}
			} 
			else if (lcpLower < lcpUpper) 
			{
                TLCPTreeIt rightChild = treeIter;
                rightChild.right();
				TDiff lcpMidUpper = *rightChild;

				if (lcpMidUpper > lcpUpper) 
				{
					// first half
					treeIter.left();
					if ((delta = delta2) == 1)
						delta = 0;
					continue;
				} 
				else if (lcpMidUpper < lcpUpper) 
				{
					// second half
					lcpLower = lcpMidUpper;
					first = mid;
                    treeIter = rightChild;
					if ((delta -= delta2) == 1) {
						++first;
						delta = 0;
					}
					continue;
				}
			}

			TSuffix		suf = suffix(text, *mid);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;

			// lcp search changes MIN to MAX here
			TDiff lcp = _max(lcpLower, lcpUpper);
			goFurther(t, lcp);
			goFurther(q, lcp);
			TDiff max = _min(difference(t, tEnd), difference(q, qEnd));
            TDiff i = max;
			for(; i && *t == *q; --i, ++t, ++q);
			lcp += max - i;

            // is text < query ?
			if (q != qEnd && (t == tEnd || !(*q < *t))) {
				// second half
				lcpLower = lcp;
				first = mid;
				treeIter.right();
				if ((delta -= delta2) == 1) {
					++first;
					delta = 0;
				}
			} else
            // is text > query ?
			if (q != qEnd && t != tEnd && (*q < *t)) {
				// first half
				lcpUpper = lcp;
				treeIter.left();
				if ((delta = delta2) == 1)
					delta = 0;
			} 
			else
			{	// range straddles mid, find each end and return
				typename Infix<TSA>::Type leftRange(sa, first, mid);
				TSAIter First2 = _lowerBoundLCPE(
					text, leftRange, treeIter.leftChild(), query, lcpLower, lcp);

				treeIter.right();
				if ((delta -= delta2) == 1) {
					delta = 0;
					++first;
				}

				typename Infix<TSA>::Type rightRange(sa, mid, last);
				TSAIter Last2 = _upperBoundLCPE(
					text, rightRange, treeIter, query, lcp, lcpUpper);

				return Pair<TSAIter> (First2, Last2);
			}
		}

		// range straddles mid, find each end and return
		typename Infix<TSA>::Type midRange(sa, first, last);
		return Pair<TSAIter> (
			_lowerBoundLCPE(text, midRange, treeIter, query, lcpLower, lcpUpper),
			_upperBoundLCPE(text, midRange, treeIter, query, lcpLower, lcpUpper)
		);
	}

	template <
		typename TText,
		typename TSA,
		typename TLCP,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA, Standard>::Type >
	_equalRangeLCPE(
		TText &text,
		TSA &sa,
		TLCP &lcp,
		TQuery &query)
	{
		return _equalRangeLCPE(text, sa, SearchTreeIterator<TLCP, LeftCompleteTree>(lcp), query);
	}


    //////////////////////////////////////////////////////////////////////////////
	// little helpers (not used)
/*
	template < typename LCPFwdIt, typename TSize >
	TSize lcp(LCPFwdIt first, TSize count) {
		if (count > 1) {
			TSize lcp = *first;
			++first;
			count-=2;
			while (count) {
				if (lcp < *first) lcp = *first;
				++first;
				--count;
			}
			return lcp;
		} else
			return 0;
	}
*/

	//////////////////////////////////////////////////////////////////////////////
	// wrappers with nice interfaces

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline typename Position<TLCPE>::Type 
	lowerBoundLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery const &query)
	{	// find first element not before query, using operator<
		return _lowerBoundLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query, 0, 0) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline typename Position<TLCPE>::Type 
	upperBoundLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return upperBoundLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query, 0, 0) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA, Standard>::Type > itPair = 
			_equalRangeLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	lowerBoundLCPEIterator(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery const &query)
	{	// find first element not before query, using operator<
		return _lowerBoundLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query, 0, 0);
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	upperBoundLCPEIterator(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return upperBoundLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query, 0, 0) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA, Standard>::Type >
	equalRangeLCPEIterator(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _equalRangeLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query);
	}

	//////////////////////////////////////////////////////////////////////////////
	// workarounds for the Visual Studio array problem

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline typename Position<TLCPE>::Type 
	lowerBoundLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery *query)
	{	// find first element not before query, using operator<
		return _lowerBoundLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query, 0, 0) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline typename Position<TLCPE>::Type 
	upperBoundLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery *query)
	{	// find first element that query is before, using operator<
		return upperBoundLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query, 0, 0) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA, Standard>::Type > itPair = 
			_equalRangeLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA, Standard>::Type >
	equalRangeLCPEIterator(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _equalRangeLCPE(text, sa, SearchTreeIterator<TLCPE const, LeftCompleteTree>(lcpe), query);
	}

	//////////////////////////////////////////////////////////////////////////////
	// wrappers with nice interfaces
/*
	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TSubText >
	inline Pair<typename Iterator<TSA const>::Type> equalRangeLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TSubText const &subtext)
	{
		return _equalRangeLCPE(
			begin(text), end(text),
			begin(sa), end(sa),
            SearchTreeIterator<typename Iterator<TLCPE const>::Type, LeftCompleteTree>(begin(lcpe), (length(text)>1)?length(text)-1:0),
			begin(subtext), end(subtext));
	}
*/

//////////////////////////////////////////////////////////////////////////////
// _findFirstIndex implementation

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void 
	_findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		ESA_FIND_MLR const)
	{
		Index<TText, TSpec> &index = haystack(finder);
		indexRequire(index, ESA_SA());
		finder.range = equalRangeSAIterator(indexText(index), indexSA(index), pattern);
	}

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void 
	_findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		ESA_FIND_LCPE const)
	{
		Index<TText, TSpec> &index = haystack(finder);
		indexRequire(index, ESA_SA());
		indexRequire(index, ESA_LCPE());
		finder.range = equalRangeLCPEIterator(indexText(index), indexSA(index), indexLCPE(index), pattern);
	}


}
#endif
