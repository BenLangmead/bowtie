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
  $Id: index_esa_algs.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_ESA_ALGS_H
#define SEQAN_HEADER_INDEX_ESA_ALGS_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	//////////////////////////////////////////////////////////////////////////////
	// more sophisticated algorithms on enhanced suffix arrays of 1 sequence
	// (also called virtual suffix trees)
	//////////////////////////////////////////////////////////////////////////////


/**
.Spec.SuperMaxRepeats Iterator:
..cat:Index
..general:Spec.BottomUp Iterator
..summary:Iterator to search for all supermaximal repeats.
..signature:Iterator<TContainer, SuperMaxRepeats>::Type
..signature:Iter<TContainer, VSTree< BottomUp<SuperMaxRepeats> > >
..param.TContainer:Type of an index that can be iterated with a bottom-up iterator.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
*/

	//////////////////////////////////////////////////////////////////////////////
	// super-maximal repeats - suffix tree version
	//////////////////////////////////////////////////////////////////////////////

	template < typename TSTree >
	struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > > {
		typedef PostorderEmptyEdges	Type;
	};

	template < typename TSTree >
	class Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > >:
		public Iter< TSTree, VSTree< BottomUp<> > > 
	{
	public:
		typedef Iter< TSTree, VSTree< BottomUp<> > >	TBase;
		typedef typename Size<TSTree>::Type				TSize;
		typedef typename Value<TSTree>::Type			TValue;
//____________________________________________________________________________

		TSize						minLength;
		typename Set<TValue>::Type	charSet;
//____________________________________________________________________________

		Iter(TSTree &_tree):
			TBase(_tree),
			minLength(1)
		{
			indexRequire(_tree, ESA_ChildTab());
			indexRequire(_tree, ESA_BWT());
			goNext(*this);	// the iterator starts in a suffix, i.e. not a repeat
		}

		Iter(TSTree &_tree, MinimalCtor):
			TBase(_tree, MinimalCtor()) {}

		Iter(TSTree &_tree, TSize _minLength):
			TBase(_tree),
			minLength(_minLength)
		{
			indexRequire(_tree, ESA_ChildTab());
			indexRequire(_tree, ESA_BWT());
			goNext(*this);	// the iterator starts in a suffix, i.e. not a repeat
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			charSet(_origin.charSet) {}
	};

	template < typename TSTree >
	inline void goNext(Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > &it) {
		do {
			goNext(it, PostorderEmptyEdges());
		} while (!atEnd(it) && 
			     !(	childrenAreLeaves(it) && 
					(repLength(it) >= it.minLength) &&
					!isPartiallyLeftExtensible(it, it.charSet)) );
	}
		

/**
.Spec.SuperMaxRepeatsFast Iterator:
..cat:Index
..general:Spec.BottomUp Iterator
..summary:Iterator to search for all supermaximal repeats (for enh. suffix arrays only).
..signature:Iterator<TContainer, SuperMaxRepeatsFast>::Type
..signature:Iter<TContainer, VSTree< BottomUp<SuperMaxRepeatsFast> > >
..param.TContainer:Type of an index based on enhanced suffix array.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
*/

	//////////////////////////////////////////////////////////////////////////////
	// supermaximal repeats - specialized for Enhanced Suffix Arrays
	//////////////////////////////////////////////////////////////////////////////

	template < typename TSTree >
	struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<SuperMaxRepeatsFast> > > > {
		typedef PostorderEmptyEdges	Type;
	};

	template < typename TText, typename TSpec >
	class Iter< Index<TText, Index_ESA<TSpec> >, VSTree< BottomUp<SuperMaxRepeatsFast> > >:
		public Iter< Index<TText, Index_ESA<TSpec> >, VSTree< BottomUp<> > > 
	{
	public:
		typedef Index< TText, Index_ESA<TSpec> >		TIndex;
		typedef Iter< TIndex, VSTree< BottomUp<> > >	TBase;
		typedef typename Size<TIndex>::Type				TSize;
		typedef typename Value<TIndex>::Type			TValue;

		typedef typename Iterator<typename Fibre<TIndex, ESA_LCP>::Type const>::Type	TLCPIter;
//____________________________________________________________________________

		TSize		minLength;
		TLCPIter	lIter, lEnd;	// lcp table iterators (optimization)
		TSize		lValueLast;		// current l-value of interval
		bool 		rising;			// is the left interval border valid
		typename Set<TValue>::Type	charSet;
//____________________________________________________________________________

		Iter(TIndex &_index):
			TBase(_index),
			minLength(1),
			lValueLast(0),
			rising(true)
		{
			this->vDesc.range = Pair<TSize>(0,0);
			indexRequire(_index, ESA_BWT());
			lIter = begin(indexLCP(this->index));
			lEnd  = end(indexLCP(this->index));
			goNext(*this);
		}

		Iter(TIndex &_index, MinimalCtor):
			TBase(_index, MinimalCtor()) {}

		Iter(TIndex &_index, TSize _minLength):
			TBase(_index),
			minLength(_minLength),
			lValueLast(0),
			rising(true)
		{
			this->vDesc.i1 = Pair<TSize>(0,0);
			indexRequire(_index, ESA_BWT());
			lIter = begin(indexLCP(this->index));
			lEnd  = end(indexLCP(this->index));
			goNext(*this);
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			lIter(_origin.lIter),
			lEnd(_origin.lEnd),
			lValueLast(_origin.lValueLast),
			rising(_origin.rising) {}
	};

	template < typename TText, typename TSpec >
	inline void goNext(Iter< Index<TText, Index_ESA<TSpec> >, VSTree< BottomUp<SuperMaxRepeatsFast> > > &it) 
	{
		typedef Index<TText, Index_ESA<TSpec> >		TIndex;
		typename Size<TIndex>::Type					lcp;

		while (it.lIter != it.lEnd) {
			lcp = *it.lIter;

			if (lcp < it.lValueLast) {
				if (it.rising) {
					if (it.lValueLast > it.minLength) {
						_dfsLCP(it) = it.lValueLast;
						++_dfsRange(it).i2;
						++it.lIter;
						it.lValueLast = lcp;
						if (!isPartiallyLeftExtensible(it, it.charSet)) return;
						continue;
					}
					it.rising = false;
				}
			} else
			if (lcp > it.lValueLast) {
				_dfsRange(it).i1 = _dfsRange(it).i2;
				it.rising = true;
			}

			++_dfsRange(it).i2;
			++it.lIter;
			it.lValueLast = lcp;
		}
		_dfsRange(it).i2 = 0;
		return;
	}

	
/**
.Spec.MaxRepeats Iterator:
..cat:Index
..general:Spec.BottomUp Iterator
..summary:Iterator to search for all maximal repeats.
..signature:Iterator<TContainer, MaxRepeats>::Type
..signature:Iter<TContainer, VSTree< BottomUp<MaxRepeats> > >
..param.TContainer:Type of an index that can be iterated with a bottom-up iterator.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
*/

	//////////////////////////////////////////////////////////////////////////////
	// maximal repeats - suffix tree version
	//////////////////////////////////////////////////////////////////////////////

	// contains a list of indices of the same bwt value (fraction)
	template <typename TSize>
	struct _FractionHeader {
		TSize	begin, end;
		TSize	size;
		_FractionHeader() {}
		_FractionHeader(TSize _begin, TSize _end, TSize _size):
			begin(_begin), end(_end), size(_size) {}
	};

	// contains a set of fractions (one for each bwt value) 
	// and a fraction for the undefined bwt value (for the virtual character at position -1)
	template <typename TValue, typename TSize>
	struct _FractionCompound {
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))
		typedef typename Set<TFraction>::Type	TSet;		// c..char, begin/end indices in posList

		TSet			set;
		TFractionHeader	leftmost;

		_FractionCompound():
			leftmost(0,0,0) {}
	};

	// tests existence of left maximal repeats
	// right maximality is implicitly given
	template <typename TValue, typename TSize>
	int _haveMaximalRepeats(
		_FractionCompound<TValue, TSize> const &a,
		_FractionCompound<TValue, TSize> const &b)
	{
		TSize cs = length(a.set), ps = length(b.set);

		if (a.leftmost.size > 0) ++cs;
		if (b.leftmost.size > 0) ++ps;

		if (cs == 0 || ps == 0) return false;
		if (cs  > 1 || ps  > 1) return true;

		if (a.leftmost.size > 0 || b.leftmost.size > 0)
			return true;

		return (keyOf(begin(a.set)) != keyOf(begin(b.set)));
	}


	template <typename TValue, typename TSize>
	int _haveMaximalRepeats(
		_FractionCompound<TValue, TSize> const &a,
		_FractionCompound<TValue, TSize> const &b,
		TValue &equalKey)
	{
		TSize cs = length(a.set), ps = length(b.set);

		if (a.leftmost.size > 0) ++cs;
		if (b.leftmost.size > 0) ++ps;

		if (cs == 0 || ps == 0) return 0;
		if (cs  > 1 || ps  > 1) return 2;	// more than 2

		if (a.leftmost.size > 0 || b.leftmost.size > 0)
			return 2;

		if ((equalKey = keyOf(begin(a.set))) != keyOf(begin(b.set)))
			return 2;
		else
			return 1;
	}


	template < typename TSTree, typename TSpec >
	struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > > {
		typedef PostorderEmptyEdges	Type;
	};

	template < typename TSTree, typename TSpec >
	class Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > >:
		public Iter< TSTree, VSTree< BottomUp<> > > 
	{
	public:
		typedef Iter< TSTree, VSTree< BottomUp<> > >	TBase;
		typedef typename Value<TSTree>::Type			TValue;
		typedef typename Size<TSTree>::Type				TSize;

		typedef _FractionCompound<TValue, TSize>		TFractionCompound;
		typedef String<TFractionCompound, Block<> >		TSetStack;
		typedef String<TSize>							TPositionList;
		
		typedef typename TFractionCompound::TSet		TSet;
		typedef typename Iterator<TSet>::Type			TSetIterator;
		typedef typename Iterator<TSet const>::Type		TConstSetIterator;

		typedef typename TBase::TStackEntry				TStackEntry;

//____________________________________________________________________________

		TSize			minLength;
		TSetStack		setStack;
		TPositionList	posList;	// this list is indexed just as SA is and contains the next entry's index
		bool			canMerge;	// is false, if parent node appears after its first child on stack
//____________________________________________________________________________

		Iter(TSTree &_index):
			TBase(_index, MinimalCtor()),
			minLength(1),
			canMerge(true)
		{
			indexRequire(_index, ESA_SA());
			indexRequire(_index, ESA_LCP());
			indexRequire(_index, ESA_BWT());
			resize(posList, length(_index));

			if (!empty(indexSA(_index))) 
			{
				_dfsOnPush(*this, TStackEntry(0,0));
				goNext(*this);
			}
		}

		Iter(TSTree &_tree, MinimalCtor):
			TBase(_tree, MinimalCtor()) {}

		Iter(TSTree &_index, TSize _minLength):
			TBase(_index, MinimalCtor()),
			minLength(_minLength),
			canMerge(true)
		{
			indexRequire(_index, ESA_SA());
			indexRequire(_index, ESA_LCP());
			indexRequire(_index, ESA_BWT());
			resize(posList, length(_index));

			if (!empty(indexSA(_index))) 
			{
				_dfsOnPush(*this, TStackEntry(0,0));
				goNext(*this);
			}
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			setStack(_origin.setStack),
			posList(_origin.posList),
			canMerge(_origin.canMerge) {}

//____________________________________________________________________________

		inline bool hasRepeats() 
		{
			if (length(setStack) < 2) return false;
			return _haveMaximalRepeats(top(setStack), topPrev(setStack)) > 0;
		}

		inline TSize countRepeats() const
		{
			if (length(setStack) < 2) return 0;

			TFractionCompound const &child  = top(setStack);
			TFractionCompound const &parent = topPrev(setStack);

			TConstSetIterator childFraction	= begin(child.set);
			TConstSetIterator childEnd		= end(child.set);
			TConstSetIterator parentFraction	= begin(parent.set);
			TConstSetIterator parentEnd		= end(parent.set);

			TSize sum = 0;
			for(; childFraction != childEnd; ++childFraction) {
				for(; parentFraction != parentEnd; ++parentFraction) {
					if (keyOf(childFraction) != keyOf(parentFraction))
						sum += objectOf(childFraction).size * objectOf(parentFraction).size;

					sum += child.leftmost.size * objectOf(parentFraction).size;
				}
				sum += objectOf(childFraction).size * parent.leftmost.size;
			}
			sum += child.leftmost.size * parent.leftmost.size;
			return sum;
		}
//____________________________________________________________________________

		inline void _dump() const {
			::std::cerr << "SETSTACK of " << representative(*this) << ":" << ::std::endl;
			typename Iterator<TSetStack const>::Type it = begin(setStack), itEnd = end(setStack);
			while (it != itEnd) {
				TSet const &set = (*it).set;
				typename Iterator<TSet const>::Type sit = begin(set), sitEnd = end(set);

				while (sit != sitEnd) {
					::std::cerr << keyOf(sit) << "::";
					typename TFractionCompound::TFractionHeader head = objectOf(sit);
					TSize i = head.begin;
					while (!_isSizeInval(i)) {
						::std::cerr << i << "  ";
						i = posList[i];
					}
					::std::cerr << ::std::endl;
					++sit;
				}

				::std::cerr << "_________________________" << ::std::endl;
				++it;
			}
		}
	};
//____________________________________________________________________________

	// add bwt partitions of child to parent node
	template < typename TSTree, typename TSpec, typename TValue, typename TSize >
	inline void _fractionMerge(
		Iter<TSTree, VSTree< BottomUp<TSpec> > > &it, 
		_FractionCompound<TValue, TSize> &parent,
		_FractionCompound<TValue, TSize> &child)
	{
		typedef _FractionCompound<TValue, TSize>	TCompound;
		typedef typename TCompound::TFraction		TFraction;
		typedef typename TCompound::TFractionHeader	TFractionHeader;
		typedef typename TCompound::TSet			TSet;
		typedef typename Iterator<TSet>::Type		TSetIterator;

		TSetIterator _end = end(child.set);
		for(TSetIterator i = begin(child.set); i != _end; ++i) {
			if (in(keyOf(i), parent.set)) {	// append child fraction to parent's fraction
				TFractionHeader &parent_header = objectOf(find(keyOf(i), parent.set));
				TFractionHeader const &child_header = objectOf(i);
				it.posList[parent_header.end] = child_header.begin;
				parent_header.end = child_header.end;
				parent_header.size += child_header.size;
			} else
				insert(TFraction(keyOf(i), objectOf(i)), parent.set);	// insert child fraction in parent's set
		}
		if (parent.leftmost.size > 0) {
			if (child.leftmost.size > 0) {
				it.posList[parent.leftmost.end] = child.leftmost.begin;
				parent.leftmost.end = child.leftmost.end;
				parent.leftmost.size += child.leftmost.size;
			}
		} else
			parent.leftmost = child.leftmost;
	}

	// maximal repeat push/leaf handlers of lcp-dfs-traversal
	template < typename TSTree, typename TElement, typename TSpec >
	inline void _dfsOnPush(Iter<TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > &it, TElement const &e) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_dfsOnPush((TBase&)it, e);

		if (it.canMerge)
			push(it.setStack);
/*
		::std::cerr << "PUSH ";
		_dumpHistoryStack(it);
		it._dump();
*/	}

	template < typename TSTree, typename TSpec >
	inline void _dfsOnLeaf(Iter<TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > &it) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_dfsOnLeaf((TBase&)it);

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef typename SAValue<TSTree>::Type	TSAValue;
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;

		push(it.setStack);

		TSTree &index = container(it);

		TSize		gPos = posGlobalize(_dfsRange(it).i1, stringSetLimits(index));
		TSAValue	lPos;
		posLocalize(lPos, _dfsRange(it).i1, stringSetLimits(index));

		if (!posAtFirstLocal(lPos))
			insert(
				TFraction(
					bwtAt(gPos, container(it)),
					TFractionHeader(gPos, gPos, 1)), 
				top(it.setStack).set);
		else
			top(it.setStack).leftmost = TFractionHeader(gPos, gPos, 1);

		_setSizeInval(it.posList[gPos]);
/*
		::std::cerr << "LEAF ";
		_dumpHistoryStack(it);
		it._dump();
*/	}

//____________________________________________________________________________

	template < typename TSTree, typename TSpec >
	inline void goNext(Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > &it) {
		do {
			if (it.canMerge && length(it.setStack) >= 2) {
				_fractionMerge(it, topPrev(it.setStack), top(it.setStack));
				pop(it.setStack);
			}
			goNext(it, PostorderEmptyEdges());
			if (empty(it.history))
				it.canMerge = false;
			else
				it.canMerge = !_dfsReversedOrder(it);
		} while (!eof(it) && !(it.canMerge && (repLength(it) >= it.minLength) && it.hasRepeats()));
	}
//____________________________________________________________________________

	template < typename TSTree, typename TSpec >
	inline typename VertexDescriptor<TSTree>::Type 
	value(Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > const &it) 
	{
		if (empty(it.history))
			return it.vDesc;
		typedef typename VertexDescriptor<TSTree>::Type TDesc;
		return TDesc(top(it.history).i1, it.vDesc.range.i2, 0);
	}
//____________________________________________________________________________

	template < typename TSTree, typename TSpec >
	inline typename Size<TSTree>::Type 
	repLength(Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > const &it) 
	{
		return top(it.history).i2;
	}
//____________________________________________________________________________

///.Function.length.param.object.type:Spec.MaxRepeats Iterator
	template < typename TSTree, typename TSpec >
	inline typename Size<TSTree>::Type 
	length(Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > const &it) {
		return it.countRepeats();
	}
//____________________________________________________________________________

///.Function.begin.param.object.type:Spec.MaxRepeats Iterator
	template < typename TSTree, class TSpec >
	inline typename Iterator< Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > >::Type
	begin(Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > &it) 
	{
		typedef Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > >	TIter;
		return typename Iterator<TIter>::Type(it);
	}
//____________________________________________________________________________

///.Function.end.param.object.type:Spec.MaxRepeats Iterator
	template < typename TSTree, class TSpec >
	inline typename Iterator< Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > >::Type
	end(Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > &it) 
	{
		typedef Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > >	TIter;
		return typename Iterator<TIter>::Type(it, MinimalCtor());
	}




	//////////////////////////////////////////////////////////////////////////////
	// maximal repeat representation (deprecated)
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSTree>
	struct MaxRepeat {
		Iter< TSTree, VSTree<BottomUp<MaxRepeats> > > &it;
	};

	template <typename TSTree>
	struct Value< MaxRepeat<TSTree> > {
		typedef Pair< typename SAValue<TSTree>::Type > Type;
	};

	template <typename TSTree>
	struct Size< MaxRepeat<TSTree> > {
		typedef typename Size<TSTree>::Type Type;
	};


	template <typename TSTree>
	inline typename Size< MaxRepeat<TSTree> >::Type 
	length(MaxRepeat<TSTree> const &repeat) {
		return repeat.it.countRepeats();
	}

/*
	template <typename TSTree>
	inline typename Iterator< MaxRepeat<TSTree> >::Type 
	begin(MaxRepeat<TSTree> &repeat) {
		return Iterator< MaxRepeat<TSTree> >::Type(repeat.it);
	}

	template <typename TSTree>
	inline typename Iterator< MaxRepeat<TSTree> const >::Type 
	begin(MaxRepeat<TSTree> const &repeat) {
		return Iterator< MaxRepeat<TSTree> >::Type(repeat.it);
	}
*/


	template <typename TSTree>
	class Iter< MaxRepeat<TSTree>, MaxRepeatOccurrences > {
	public:

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef	Pair<TSize>						TPair;

		typedef _FractionCompound<TValue, TSize>	TFractionCompound;
		typedef typename TFractionCompound::TSet	TSet;
		typedef typename Iterator<TSet const>::Type	TSetIterator;

		TSize			childPtr, parentPtr;
		TSetIterator	childFraction, childEnd;
		TSetIterator	parentFraction, parentEnd;
		bool			_atEnd;
		TPair			tmp;
		bool			leftmostChild, leftmostParent;
		
		Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const *maxIt;

		inline Iter(Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const &_maxIt):
			maxIt(&_maxIt)
		{
			_init();
		}
		
		inline Iter(Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const &_maxIt, MinimalCtor):
			maxIt(&_maxIt),
			_atEnd(true) {}
		
		inline bool _innerStep() {
			if (_isSizeInval(childPtr = maxIt->posList[childPtr])) {
				if (_isSizeInval(parentPtr = maxIt->posList[parentPtr])) return false;
				childPtr = objectOf(childFraction).begin;
			}
			return true;
		}

		inline void _firstParentFraction() {
			TFractionCompound const &parent = topPrev(maxIt->setStack);

			parentFraction	= begin(parent.set);
			parentEnd		= end(parent.set);

			if (parentFraction != parentEnd) {
				leftmostParent = false;
				parentPtr = objectOf(parentFraction).begin;
			} else {
				leftmostParent = true;
				parentPtr = parent.leftmost.begin;
			}
		}

		inline void _firstChildFraction() {
			TFractionCompound const &child = top(maxIt->setStack);

			childFraction	= begin(child.set);
			childEnd		= end(child.set);

			if (childFraction != childEnd) {
				leftmostChild = false;
				childPtr = objectOf(childFraction).begin;
			} else {
				leftmostChild = true;
				childPtr = child.leftmost.begin;
			}
		}

		inline bool _nextParentFraction() {
			if (leftmostParent)
				return false;

			if (++parentFraction == parentEnd) {
				if (topPrev(maxIt->setStack).leftmost.size > 0) {
					leftmostParent = true;
					parentPtr = topPrev(maxIt->setStack).leftmost.begin;
				} else
					return false;
			} else
				parentPtr = objectOf(parentFraction).begin;

			return true;
		}

		inline bool _nextChildFraction() {
			if (leftmostChild)
				return false;

			if (++childFraction == childEnd) {
				if (top(maxIt->setStack).leftmost.size > 0) {
					leftmostChild = true;
					childPtr = top(maxIt->setStack).leftmost.begin;
				} else
					return false;
			} else
				childPtr = objectOf(childFraction).begin;

			return true;
		}

		inline bool _outerStep() {
			do {
				if (!_nextChildFraction()) {
					_firstChildFraction();
					if (!_nextParentFraction()) {
						_atEnd = true;
						return false;
					}
				}
				if (leftmostChild || leftmostParent) break;
			} while (keyOf(childFraction) == keyOf(parentFraction));		// ignore occurences with equal bwt entries
			return true;
		}

		inline void _init() 
		{
			if (length(maxIt->setStack) < 2) {
				_atEnd = true;
				return;
			}

			_firstChildFraction();
			_firstParentFraction();

			if (!leftmostChild && !leftmostParent &&
				(keyOf(childFraction) == keyOf(parentFraction)))
				_atEnd = !_outerStep();
			else
				_atEnd = false;

			if (!_atEnd) {
				tmp.i1 = saAt(parentPtr, container(*maxIt));
				tmp.i2 = saAt(childPtr, container(*maxIt));
			}
		}
	};

//____________________________________________________________________________

	template < typename TRepeat >
	inline typename Value< Iter<TRepeat, MaxRepeatOccurrences> >::Type &
	value(Iter<TRepeat, MaxRepeatOccurrences> const &it)  {
		return it.tmp;
	}

	template < typename TRepeat >
	inline typename Value< Iter<TRepeat, MaxRepeatOccurrences> >::Type &
	value(Iter<TRepeat, MaxRepeatOccurrences> &it)  {
		return it.tmp;
	}
//____________________________________________________________________________

	template < typename TRepeat >
	inline Iter<TRepeat, MaxRepeatOccurrences> &
	goNext(Iter<TRepeat, MaxRepeatOccurrences> &it)  {
		if (it._innerStep()) {
			it.tmp.i1 = saAt(it.parentPtr, container(*it.maxIt));
			it.tmp.i2 = saAt(it.childPtr, container(*it.maxIt));
			return it;
		}
		if (it._outerStep()) {
			it.tmp.i1 = saAt(it.parentPtr, container(*it.maxIt));
			it.tmp.i2 = saAt(it.childPtr, container(*it.maxIt));
		}
		return it;
	}
//____________________________________________________________________________

	template < typename TRepeat >
	inline Iter<TRepeat, MaxRepeatOccurrences> &
	goBegin(Iter<TRepeat, MaxRepeatOccurrences> &it) 
	{
		it._init;
	}

//____________________________________________________________________________

	template < typename TRepeat >
	inline Iter<TRepeat, MaxRepeatOccurrences> &
	goEnd(Iter<TRepeat, MaxRepeatOccurrences> &it) 
	{
		it._atEnd = true;
	}

	template < typename TRepeat >
	inline bool atEnd(Iter<TRepeat, MaxRepeatOccurrences> const &it) {
		return it._atEnd;
	}

	template < typename TRepeat >
	inline bool atEnd(Iter<TRepeat, MaxRepeatOccurrences> &it) {
		return it._atEnd;
	}

//____________________________________________________________________________

	template < typename TRepeat >
	inline bool
	operator == (
		Iter<TRepeat, MaxRepeatOccurrences> const &itA,
		Iter<TRepeat, MaxRepeatOccurrences> const &itB)
	{
		if (itA._atEnd && itB._atEnd) return true;
		if (itA._atEnd || itB._atEnd) return false;
		return (itA.childPtr == itB.childPtr) && (itA.parentPtr == itB.parentPtr);
	}

	template < typename TRepeat >
	inline bool
	operator != (
		Iter<TRepeat, MaxRepeatOccurrences> const &itA,
		Iter<TRepeat, MaxRepeatOccurrences> const &itB)
	{
		if (itA._atEnd && itB._atEnd) return false;
		if (itA._atEnd || itB._atEnd) return true;
		return (itA.childPtr != itB.childPtr) || (itA.parentPtr != itB.parentPtr);
	}
//____________________________________________________________________________


	template <typename TSTree, typename TSpec>
	struct Size< Iter<TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > > {
		typedef typename Size<TSTree>::Type Type;
	};

	template <typename TSTree>
	struct Size< Iter<MaxRepeat<TSTree>, MaxRepeatOccurrences> > {
		typedef typename Size<TSTree>::Type Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// Iterator wrappers
	//////////////////////////////////////////////////////////////////////////////

	// iterates over all supermaximal repeats
	template <typename TSTree>
	struct Iterator< TSTree, SuperMaxRepeats > {
		typedef Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > Type;
	};
//____________________________________________________________________________

	// iterates over all maximal unique matches
	template <typename TSTree>
	struct Iterator< TSTree, SuperMaxRepeatsFast > {
		typedef Iter< TSTree, VSTree< BottomUp<SuperMaxRepeatsFast> > > Type;
	};
//____________________________________________________________________________

	// iterates over all maximal repeat structures
	template <typename TSTree, typename TSpec>
	struct Iterator< TSTree, _MaxRepeats<TSpec> > {
		typedef Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > Type;
	};
//____________________________________________________________________________

	// iterates over all maximal repeat pairs of a repeat structure

	// Iterator of Iterator<TIndex, MaxRepeats>
	template <typename TSTree, typename TSpec>
	struct Iterator< Iter< TSTree, TSpec >, MaxRepeatOccurrences > {
		typedef Iter< MaxRepeat<TSTree>, MaxRepeatOccurrences > Type;
	};

	template <typename TSTree, typename TSpec>
	struct DefaultIteratorSpec< Iter< TSTree, VSTree< BottomUp<_MaxRepeats<TSpec> > > > > {
		typedef MaxRepeatOccurrences Type;
	};

	// alternative (use MaxRepeat<TIndex> as a wrapper in between)
	template <typename TSTree>
	struct Iterator< MaxRepeat<TSTree>, MaxRepeatOccurrences > {
		typedef Iter <MaxRepeat<TSTree>, MaxRepeatOccurrences > Type;
	};

	template <typename TSTree>
	struct DefaultIteratorSpec< MaxRepeat<TSTree> > {
		typedef MaxRepeatOccurrences Type;
	};

//}

}

#endif
