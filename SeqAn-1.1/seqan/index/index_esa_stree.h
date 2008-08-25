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
  $Id: index_esa_stree.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_ESA_STREE_H
#define SEQAN_HEADER_INDEX_ESA_STREE_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Spec.VSTree Iterator:
..cat:Index
..summary:Abstract iterator for Suffix Trees.
..signature:Iter<TContainer, VSTree<TSpec> >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..param.TSpec:The specialization type.
..remarks:This iterator is a pointer to a node in the Suffix Tree (given by the Enhanced Suffix Array @Spec.Index_ESA@).
Every node can uniquely be mapped to an interval of the Suffix Array containing all suffixes of the node's subtree.
This interval is the @Function.value@ of the iterator.
*/

	template < typename TIndex, typename TSpec >
    struct Value< Iter< TIndex, VSTree<TSpec> > > {
		typedef typename VertexDescriptor<TIndex>::Type Type;
	};
 
	template < typename TIndex, typename TSpec >
	struct Size< Iter< TIndex, VSTree<TSpec> > > {
		typedef typename Size<TIndex>::Type Type;
	};
 
	template < typename TIndex, typename TSpec >
	struct Position< Iter< TIndex, VSTree<TSpec> > > {
		typedef typename Position<TIndex>::Type Type;
	};
 

/**
.Spec.TopDown Iterator:
..cat:Index
..general:Spec.VSTree Iterator
..summary:Iterator for Suffix Trees that can go down and right beginning from the root.
..signature:Iterator<TContainer, TopDown<TSpec> >::Type
..signature:Iter<TContainer, VSTree< TopDown<TSpec> > >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..param.TSpec:The specialization type.
*/


	template < typename TIndex, class TSpec >
	class Iter< TIndex, VSTree< TopDown<TSpec> > > 
	{
	public:

		typedef typename VertexDescriptor<TIndex>::Type	TVertexDesc;
		typedef Iter									iterator;

		TIndex const	&index;		// container of all necessary tables
		TVertexDesc		vDesc;		// current interval in suffix array and
									// right border of parent interval (needed in goRight)

		// pseudo history stack (to go up at most one node)
		TVertexDesc		_parentDesc;

		Iter(TIndex &_index):
			index(_index)
		{
			_indexRequireTopDownIteration(_index);
			goRoot(*this);
		}

		Iter(TIndex &_index, MinimalCtor):
			index(_index),
			vDesc(MinimalCtor()) {}

		Iter(TIndex &_index, TVertexDesc const &_vDesc):
			index(_index),
			vDesc(_vDesc)
		{
			_indexRequireTopDownIteration(_index);
		}

		Iter(Iter const &_origin):
			index(container(_origin)),
			vDesc(value(_origin)) {}
	};


/**
.Spec.TopDownHistory Iterator:
..cat:Index
..general:Spec.TopDown Iterator
..summary:Iterator for Suffix Trees that can go down, right, and up.
..signature:Iterator<TContainer, TopDown< ParentLinks<TSpec> > >::Type
..signature:Iter<TContainer, VSTree< TopDown< ParentLinks<TSpec> > > >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..implements:Concept.Iterator
..param.TSpec:The specialization type.
*/

	template < typename TVSTreeIter >
	struct _HistoryStackEntry;
	
	template < typename TIndex, typename TSpec >
	struct _HistoryStackEntry< Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > >
	{
		typedef Pair<typename Size<TIndex>::Type>	Type;
	};

	template < typename TIndex, class TSpec >
	class Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > >:
		public Iter< TIndex, VSTree< TopDown<> > >
	{
	public:

		typedef Iter< TIndex, VSTree< TopDown<> > >		TBase;
		typedef	typename _HistoryStackEntry<Iter>::Type	TStackEntry;
		typedef String<TStackEntry, Block<> >			TStack;
		typedef Iter									iterator;

		TStack			history;	// contains all previously visited intervals (allows to go up)

		Iter(TIndex &_index):
			TBase(_index) {}

		Iter(TIndex &_index, MinimalCtor):
			TBase(_index, MinimalCtor()) {}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			history(_origin.history) {}
	};


/**
.Spec.BottomUp Iterator:
..cat:Index
..general:Spec.VSTree Iterator
..summary:Iterator for an efficient postorder depth-first search in a suffix tree.
..signature:Iterator<TContainer, BottomUp<TSpec> >::Type
..signature:Iter<TContainer, VSTree< BottomUp<TSpec> > >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..implements:Concept.Iterator
..param.TSpec:The specialization type.
*/

	template < typename TIndex, typename TSpec >
	class Iter< TIndex, VSTree< BottomUp<TSpec> > > 
	{
	public:

		typedef typename VertexDescriptor<TIndex>::Type	TVertexDesc;
		typedef typename Size<TIndex>::Type				TSize;
		typedef	Pair<TSize>								TStackEntry;
		typedef String<TStackEntry, Block<> >			TStack;
		typedef Iter									iterator;

		TIndex	const	&index;			// container of all necessary tables
		TVertexDesc		vDesc;			// current interval in suffix array and
										// right border of parent interval (unused here)
		TSize			lValue;			// current l-value
		TStack			history;		// contains all left borders of current l-intervals (== left borders of history intervals)

		Iter(TIndex &_index):
			index(_index),
			vDesc(MinimalCtor()),
			lValue(0)
		{
			_indexRequireBottomUpIteration(_index);
			goBegin(*this);
		}

		Iter(TIndex &_index, MinimalCtor):
			index(_index),
			vDesc(MinimalCtor()),
			lValue(0) {}

		Iter(Iter const &_origin):
			index(container(_origin)),
			vDesc(value(_origin)),
			lValue(_dfsLCP(_origin)),
			history(_origin.history) {}
	};


	//////////////////////////////////////////////////////////////////////////////
	// Iterator wrappers
	//////////////////////////////////////////////////////////////////////////////

	template <typename TObject, typename TSpec>
	struct Iterator< TObject, BottomUp<TSpec> > {
		typedef Iter< TObject, VSTree< BottomUp<TSpec> > > Type;
	};

	template <typename TObject, typename TSpec>
	struct Iterator< TObject, TopDown<TSpec> > {
		typedef Iter< TObject, VSTree< TopDown<TSpec> > > Type;
	};




	template < typename TIndex, typename TSpec >
	inline void _dumpHistoryStack(Iter<TIndex, VSTree<TSpec> > &it) {
		for(typename Size<TIndex>::Type i = 0; i < length(it.history); ++i)
			::std::cerr << it.history[i] << '\t';
		::std::cerr << value(it) << ::std::endl;
	}

	template <typename TText, typename TSpec>
	inline void
	_dump(Index<TText, Index_ESA<TSpec> > &index)
	{
		::std::cout << "  SA" << ::std::endl;
		for(unsigned i=0; i < length(indexSA(index)); ++i)
			::std::cout << i << ":  " << indexSA(index)[i] << "  " << suffix(indexText(index), indexSA(index)[i]) << ::std::endl;

		::std::cout << ::std::endl << "  LCP" << ::std::endl;
		for(unsigned i=0; i < length(indexLCP(index)); ++i)
			::std::cout << i << ":  " << indexLCP(index)[i] << ::std::endl;

		::std::cout << ::std::endl << "  ChildTab" << ::std::endl;
		for(unsigned i=0; i < length(indexChildTab(index)); ++i)
			::std::cout << i << ":  " << indexChildTab(index)[i] << ::std::endl;

		::std::cout << ::std::endl;
	}


	template < typename TIndex, typename TSpec >
	inline bool _dfsReversedOrder(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it) {
        return lcpAt(_dfsRange(it).i2 - 1, container(it)) > top(it.history).i2;
	}

	// standard push/pop handlers of lcp-dfs-traversal
	template < typename TIndex, typename TSpec, typename TSize >
	inline void _dfsOnPop(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, TSize const) {
        _dfsRange(it).i1 = top(it.history).i1;
		_dfsLCP(it) = top(it.history).i2;
		pop(it.history);
	}

	template < typename TIndex, typename TSpec, typename TElement >
	inline void _dfsOnPush(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, TElement const &e) {
		push(it.history, e);
	}

	template < typename TIndex, typename TSpec >
	inline void _dfsOnLeaf(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it) {
		_setSizeInval(_dfsLCP(it));
	}


	// postorder bottom up iterator (dfs)
	template < typename TIndex, typename TSpec, typename THideEmptyEdges >
	inline void goNextImpl(
		Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, 
		VSTreeIteratorTraits<_Postorder, THideEmptyEdges> const) 
	{
		TIndex const &index = container(it);
		do {
			// postorder dfs via lcp-table
			if (isRoot(it)) {
				_dfsClear(it);
				return;
			}

			if (_dfsRange(it).i2)
			{
				typedef typename Size<TIndex>::Type TSize;
				typedef typename Iter<TIndex, VSTree< BottomUp<TSpec> > >::TStackEntry TStackEntry;
				TStackEntry	_top_ = top(it.history);
				TSize		lcp_i = lcpAt(_dfsRange(it).i2 - 1, index);

				if (lcp_i < _top_.i2) {
					_dfsOnPop(it, lcp_i);
					if (nodePredicate(it)) return;
					else continue;
				}

				if (lcp_i > _top_.i2) {
					_top_.i1 = _dfsRange(it).i1;
					_top_.i2 = lcp_i;
					_dfsOnPush(it, _top_);
				}

	// innerer Knoten:
	// wenn kein Pop, aber Push -> begehe mind. 2. Teilbaum irgendeines Vorfahrs
	// wenn kein Pop, kein Push -> verlasse mind. zweites Kindblatt
	// wenn Pop, aber kein Push -> verlasse Wurzel des mind. 2.Teilbaums
	// wenn Pop und Push        -> verlasse ersten Teilbaum (sieht Vater zum ersten Mal und pusht jenen)

	// wenn nach Pop ein Pop folgen würde	-> Vater ist Top of Stack
	// wenn nach Pop ein Push folgen würde	-> Vater erst beim Push auf Stack (-> zwischenspeichern)
			}

			// last lcp entry (== 0) causes removal of toplevel interval
			if ((_dfsRange(it).i1 = _dfsRange(it).i2++) == length(index)) {
				_dfsOnPop(it, 0);
				_dfsRange(it).i2 = _dfsRange(it).i1;
			} else {
				// skip $ leafs (empty edges)
				if (THideEmptyEdges::VALUE &&
					suffixLength(saAt(_dfsRange(it).i1, index), index) == lcpAt(_dfsRange(it).i1, index))
					continue;

				_dfsOnLeaf(it);
	// Blatt:
	// wenn danach kein Pop, aber Push -> Vater wird erst noch gepusht
	// wenn danach Pop				   -> Vater ist Top des Stack
			}
			if (nodePredicate(it)) return;
		} while (true);
	}

/**
.Function.repLength:
..summary:Returns the length of the substring representing the path from root to $iterator$ node.
..cat:Index
..signature:repLength(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The length of the sequence returned by @Function.representative@
...type:Metafunction.Size|Size type of the underlying index
*/

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type 
	repLength(Iter< TIndex, VSTree<BottomUp<TSpec> > > const &it) 
	{
		typename Size<TIndex>::Type lcp;
		if (!_isSizeInval(lcp = it.lValue))
			return lcp;
		else
			return suffixLength(getOccurrence(it), container(it));
	}


	template < typename TIndex, typename TSize >
	inline typename Size<TIndex>::Type
	repLength(TIndex const &index, VertexESA<TSize> const &vDesc) 
	{
		if (_isLeaf(vDesc)) return suffixLength(saAt(vDesc.range.i1, index), index);
		if (_isRoot(vDesc)) return 0;

		// get l-value of suffix array range
		TSize lval = _getUp(vDesc.range.i2, index);
		if (!(vDesc.range.i1 < lval && lval < vDesc.range.i2))
			lval = _getDown(vDesc.range.i1, index);
		
		// retrieve the longest-common-prefix length of suffices in range
		return lcpAt(lval - 1, index);
	}

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type
	repLength(Iter< TIndex, VSTree<TopDown<TSpec> > > const &it) 
	{
		return repLength(container(it), value(it));
	}

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type
	nodeDepth(Iter< TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > const &it) 
	{
		return length(it.history);
	}

/**
.Function.parentRepLength:
..summary:Returns the length of the substring representing the path from root to $iterator$'s parent node.
..cat:Index
..signature:parentRepLength(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The length of the sequence returned by @Function.parentRepresentative@
...type:Metafunction.Size|Size type of the underlying index
*/

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type
	parentRepLength(Iter< TIndex, VSTree<TopDown<TSpec> > > const &it) 
	{
		return repLength(container(it), nodeUp(it));
	}


/**
.Function.emptyParentEdge:
..summary:Returns $true$ iff the edge label from the $iterator$ node to its parent is empty.
..cat:Index
..signature:bool emptyParentEdge(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
..returns:$true$ if @Function.parentEdgeLength@$ returns 0, otherwise $false$.
...type:Metafunction.Size|Size type of the underlying index
*/

	template < typename TIndex, typename TSpec >
	inline bool
	emptyParentEdge(Iter< TIndex, VSTree<TopDown<TSpec> > > const &it) 
	{
		// the following is more efficient than 
		// return parentEdgeLength(it) == 0;
		TIndex const &index = container(it);
		typename SAValue<TIndex>::Type pos = getOccurrence(it);
		return getSeqOffset(pos, stringSetLimits(index)) + parentRepLength(it)
			== sequenceLength(getSeqNo(pos, stringSetLimits(index)), index);
	}


/**
.Function.lca:
..summary:Returns the last common ancestor of two tree nodes.
..cat:Index
..signature:bool lca(a, b, result)
..param.a:The first node.
...type:Spec.TopDownHistory Iterator
..param.b:The second node.
...type:Spec.TopDownHistory Iterator
..param.result:A reference to the resulting lca node.
...type:Spec.TopDownHistory Iterator
..returns:$false$ if the lca of $a$ and $b$ is the root node, otherwise $true$.
*/

	template < typename TIndex, class TSpec1, class TSpec2 >
	inline bool lca(
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > > &a, 
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > > &b, 
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > > &_lca)
	{
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > >::TStack::const_iterator iA;
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > >::TStack::const_iterator iB;

		typedef typename Size<TIndex>::Type TSize;

		// push current intervals
		push(a.history, value(a).range);
		push(b.history, value(b).range);

		TSize s = min(a.history.size(), b.history.size()), i0 = 0;
		
		while (s) {
			TSize m = s / 2;
			iA = a.history.begin() + i0 + m;
			iB = b.history.begin() + i0 + m;
			if ((*iA).i1 == (*iB).i1 && (*iA).i2 == (*iB).i2) {
				i0 += m + 1;
				s -= m + 1;
			} else
				s = m;
		}

		_lca.history.resize(i0);
		copy(a.history.begin(), a.history.begin() + i0, _lca.history.begin());

		// pop current intervals
		pop(a.history);
		pop(b.history);
		goUp(_lca);

		return i0;
	}

/**
.Function.lcp:
..summary:Returns the length of the longest-common-prefix of two Suffix Tree nodes.
..cat:Index
..signature:lcp(a, b)
..param.a:The first node.
...type:Spec.TopDownHistory Iterator
..param.b:The second node.
...type:Spec.TopDownHistory Iterator
..returns:The lcp-length of $a$ and $b$.
*/

	// return the lcp of a and b by seeking the lca of them
	template < typename TIndex, class TSpec1, class TSpec2 >
	inline typename Size<TIndex>::Type lcp(
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > > &a, 
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > > &b) 
	{
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > >::TStack::const_iterator iA;
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > >::TStack::const_iterator iB;

		typedef typename Size<TIndex>::Type				TSize;
		typedef typename VertexDescriptor<TIndex>::Type	TDesc;

		// push current intervals
		push(a.history, value(a).range);
		push(b.history, value(b).range);

		TSize s = min(a.history.size(), b.history.size()), i0 = 0;
		
		while (s) {
			TSize m = s / 2;
			iA = a.history.begin() + i0 + m;
			iB = b.history.begin() + i0 + m;
			if ((*iA).i1 == (*iB).i1 && (*iA).i2 == (*iB).i2) {
				i0 += m + 1;
				s -= m + 1;
			} else
				s = m;
		}

		TSize _lcp = (i0 > 0)? repLength(container(a), TDesc(a.history[i0 - 1], 0)): 0;

		// pop current intervals
		pop(a.history);
		pop(b.history);

		return _lcp;
	}

///.Function.container.param.iterator.type:Spec.VSTree Iterator

	template < typename TIndex, class TSpec >
	inline TIndex const & container(Iter< TIndex, VSTree<TSpec> > const &it) { 
		return it.index; 
	}

	template < typename TIndex, class TSpec >
	inline TIndex & container(Iter< TIndex, VSTree<TSpec> > &it) { 
		return const_cast<TIndex&>(it.index); 
	}


///.Function.value.param.object.type:Spec.VSTree Iterator

	template < typename TIndex, class TSpec >
	inline typename VertexDescriptor<TIndex>::Type & 
	value(Iter< TIndex, VSTree<TSpec> > &it) { 
		return it.vDesc;
	}

	template < typename TIndex, class TSpec >
	inline typename VertexDescriptor<TIndex>::Type const & 
	value(Iter< TIndex, VSTree<TSpec> > const &it) { 
		return it.vDesc;
	}


/**
.Function.getOccurrence:
..summary:Returns an occurence of the @Function.representative@ substring in the index text.
..cat:Index
..signature:getOccurrence(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:A position where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA Index Fibres.ESA_Text@).
...metafunction:Metafunction.SAValue.<TIndex>
*/

	template < typename TIndex, class TSpec >
	inline typename SAValue<TIndex>::Type 
	getOccurrence(Iter< TIndex, VSTree<TSpec> > const &it) {
		return saAt(value(it).range.i1, container(it));
	}


/**
.Function.countOccurrences:
..summary:Returns the number of occurences of @Function.representative@ in the index text.
..cat:Index
..signature:countOccurrences(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The number of positions where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA Index Fibres.ESA_Text@).
If $iterator$'s container type is $TIndex$ the return type is $Size<TIndex>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type 
	countOccurrences(Iter< TIndex, VSTree<TSpec> > const &it) 
	{
		if (_isSizeInval(value(it).range.i2))
			return length(indexSA(container(it))) - value(it).range.i1;
		else
			return value(it).range.i2 - value(it).range.i1;
	}

/**
.Function.getOccurrences:
..summary:Returns all occurences of the @Function.representative@ substring in the index text.
..cat:Index
..signature:getOccurrences(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:All positions where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA Index Fibres.ESA_Text@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_SA>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, Fibre_SA>::Type const >::Type 
	getOccurrences(Iter< TIndex, VSTree<TSpec> > const &it) 
	{
		if (_isSizeInval(value(it).range.i2))
			return infix(indexSA(container(it)), value(it).range.i1, length(indexSA(container(it))));
		else
			return infix(indexSA(container(it)), value(it).range.i1, value(it).range.i2);
	}

/**
.Function.alignment:
..summary:Returns an alignment of the occurences of the @Function.representative@ substring in the index text.
..cat:Index
..signature:alignment(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:A local alignment corresponding to the seed of the $iterator$.
..remarks:The @Function.representative@ must uniquely occur in every sequence (e.g. in MUMs), 
otherwise the seed returned is one many.
*/

	template < typename TString, typename TSSetSpec, typename TIndexSpec, class TSpec >
	inline Align<TString, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, TSSetSpec>, TIndexSpec >, VSTree<TSpec> > &it) 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;

		Align<TString, ArrayGaps> align;
		TIndex &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurrences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs, Standard()), occEnd = end(occs, Standard());
		while (occ != occEnd) {
			typename Size<TIndex>::Type seqNo = getSeqNo(*occ, stringSetLimits((TIndex const&)index));
			typename Size<TIndex>::Type seqOfs = getSeqOffset(*occ, stringSetLimits((TIndex const&)index));
			setSource(row(align, seqNo), value(indexText(index), seqNo), seqOfs, seqOfs + repLen);
			++occ;
		}
		return align;
	}
/*
	template < typename TString, typename TConcSpec, typename TIndexSpec, class TSpec >
	inline typename Align<TString const, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec >, VSTree<TSpec> > const &it) 
	{
		typedef Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;

		Align<TString const, ArrayGaps> align;
		TIndex const &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurrences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs, Standard()), occEnd = end(occs, Standard());
		while (occ != occEnd) {
			typename Size<TIndex>::Type seqNo = getSeqNo(*occ, stringSetLimits(index));
			typename Size<TIndex>::Type globOfs = posGlobalize(*occ, stringSetLimits(index));
			setSource(row(align, seqNo), concat(indexText(index)), globOfs, globOfs + repLen);
			++occ;
		}
		return align;
	}
*/
	template < typename TString, typename TConcSpec, typename TIndexSpec, class TSpec >
	inline Align<TString, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec >, VSTree<TSpec> > &it) 
	{
		typedef Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;

		Align<TString, ArrayGaps> align;
		TIndex &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurrences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs, Standard()), occEnd = end(occs, Standard());
		while (occ != occEnd) {
			typename Size<TIndex>::Type seqNo = getSeqNo(*occ, stringSetLimits((TIndex const&)index));
			typename Size<TIndex>::Type globOfs = posGlobalize(*occ, stringSetLimits((TIndex const&)index));
			setSource(row(align, seqNo), concat(indexText(index)), globOfs, globOfs + repLen);
			++occ;
		}
		return align;
	}

/**
.Function.getOccurrencesBWT:
..summary:Returns the characters left beside all occurence of the @Function.representative@ substring in the index text.
..cat:Index
..signature:getOccurrencesBWT(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:All positions where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA Index Fibres.ESA_Text@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_BWT>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, ESA_BWT>::Type const >::Type 
	getOccurrencesBWT(Iter< TIndex, VSTree<TSpec> > const &it) 
	{
		if (_isSizeInval(value(it).range.i2))
			return infix(indexBWT(container(it)), value(it).range.i1, length(indexSA(container(it))));
		else
			return infix(indexBWT(container(it)), value(it).range.i1, value(it).range.i2);
	}

/**
.Function.representative:
..summary:Returns a substring representing the path from root to $iterator$ node.
..cat:Index
..signature:representative(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:An @Spec.InfixSegment@ of the raw text of an index (see @Tag.ESA Index Fibres.ESA_Text@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_RawText>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, Fibre_Text>::Type const >::Type 
	representative(Iter< TIndex, VSTree<TSpec> > const &it) 
	{
		return infixWithLength(indexText(container(it)), getOccurrence(it), repLength(it));
	}


/**
.Function.countChildren:
..summary:Count the number of children of a tree node.
..cat:Index
..signature:countChildren(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The number of children of a tree node.
If $iterator$'s container type is $TIndex$, the return type is $Size<TIndex>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type 
	countChildren(Iter< TIndex, VSTree<TSpec> > const &it) 
	{
		if (_isLeaf(it, EmptyEdges())) return 0;

		typedef typename Size<TIndex>::Type TSize;

		TSize i = _getUp(value(it).range.i2, container(it));
		if (!(value(it).range.i1 < i && i < value(it).range.i2))
			i = _getDown(value(it).range.i1, container(it));

		TSize result = (isRoot(it))? 1: 2;
		while (_isNextl(i, container(it))) {
			i = _getNextl(i, container(it));
			++result;
		}
		return result;
	}

	// get the interval of SA of the subtree under the edge beginning with character c
	template < typename TText, class TSpec, typename TValue >
	inline bool 
	_getNodeByChar(
		Iter< Index<TText, Index_ESA<TSpec> >, VSTree<TSpec> > const &it, 
		TValue c, 
		typename VertexDescriptor< Index<TText, Index_ESA<TSpec> > >::Type &childDesc)
	{
		typedef Index<TText, Index_ESA<TSpec> >				TIndex;
		typedef typename Size<TIndex>::Type					TSize;

		if (_isLeaf(it, EmptyEdges())) return false;

		Pair<TSize> child(value(it).range.i1, _getUp(value(it).range.i2, container(it)));
		if (!(value(it).range.i1 < child.i2 && child.i2 < value(it).range.i2))
			child.i2 = _getDown(value(it).range.i1, container(it));

		TSize _lcp = lcpAt(child.i2 - 1, container(it));
		if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
			childDesc.range = child;
			childDesc.parentRight = value(it).range.i2;
			return true;
		}
		child.i1 = child.i2;
		while (_isNextl(child.i2, container(it))) 
		{
			child.i2 = _getNextl(child.i2, container(it));
			if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				childDesc.range = child;
				childDesc.parentRight = value(it).range.i2;
				return true;
			}
			child.i1 = child.i2;
		}

		if (!isRoot(it)) {
			if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				childDesc.range.i1 = child.i1;
				childDesc.range.i2 = childDesc.i2 = value(it).range.i2;
				return true;
			}
		}
		return false;
	}


/**
.Function.nodePredicate:
..summary:If $false$ this node will be skipped during the bottom-up traversal.
..cat:Index
..signature:bool nodePredicate(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:
*/

	template < typename TIndex, class TSpec >
	inline bool
	nodePredicate(Iter<TIndex, TSpec> &)
	{
		return true;
	}


/**
.Function.nodeHullPredicate:
..summary:If $false$ this node and its subtree is concealed.
..cat:Index
..signature:bool nodeHullPredicate(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:
*/

	template < typename TIndex, class TSpec >
	inline bool
	nodeHullPredicate(Iter<TIndex, TSpec> &)
	{
		return true;
	}

//____________________________________________________________________________

	template < typename TText, typename TIndexSpec, class TSpec >
	inline void goRoot(Iter<Index<TText, TIndexSpec>, VSTree<TSpec> > &it) 
	{
		_historyClear(it);
		clear(it);							// start in root node with range (0,infty)
		_setSizeInval(value(it).range.i2);	// infty is equivalent to length(index) and faster to compare
	}

//____________________________________________________________________________

///.Function.begin.param.object.type:Class.Index
	template < typename TText, typename TIndexSpec, class TSpec >
	inline typename Iterator<Index<TText, TIndexSpec>, TSpec >::Type
	begin(Index<TText, TIndexSpec> &index, TSpec const) 
	{
		return typename Iterator<
			Index<TText, TIndexSpec>, 
			TSpec 
		>::Type (index);
	}

///.Function.goBegin.param.iterator.type:Spec.VSTree Iterator
	template < typename TText, typename TIndexSpec, class TSpec >
	inline void goBegin(Iter<Index<TText, TIndexSpec>, VSTree<TSpec> > &it) 
	{
		typedef Iter<Index<TText, TIndexSpec>, VSTree<TSpec> >	TIter;
		typedef typename GetVSTreeIteratorTraits<TIter>::Type	TTraits;
		typedef typename TTraits::HideEmptyEdges				THideEmptyEdges;

		goRoot(it);

		if (TYPECMP<typename TTraits::DFSOrder, _Postorder>::VALUE) {
			while (goDown(it));
			return;
		}

		// if root doesn't suffice predicate, do a dfs-step
		if ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it))
			goNext(it);
	}

	template < typename TText, typename TIndexSpec, class TSpec >
	inline void goBegin(Iter<Index<TText, Index_ESA<TIndexSpec> >, VSTree< BottomUp<TSpec> > > &it) 
	{
		typedef Index<TText, Index_ESA<TIndexSpec> >		TIndex;
		typedef Iter<TIndex, VSTree< BottomUp<TSpec> > >	TIter;
		typedef typename VertexDescriptor<TIndex>::Type		TVertexDesc;
		typedef typename Size<TIndex>::Type					TSize;
		typedef	Pair<TSize>									TStackEntry;

		_dfsClear(it);
		clear(it);							// start in root node with range (0,infty)
		if (!empty(indexSA(container(it)))) {
			_dfsOnPush(it, TStackEntry(0,0));
			goNextImpl(it, typename GetVSTreeIteratorTraits< TIter >::Type());
		}
	}

//____________________________________________________________________________

///.Function.end.param.object.type:Class.Index
	template < typename TText, typename TIndexSpec, class TSpec >
	inline typename Iterator<Index<TText, TIndexSpec>, TSpec >::Type
	end(Index<TText, TIndexSpec> &index, TSpec const) 
	{
		return typename Iterator< 
			Index<TText, TIndexSpec>, 
			TSpec
		>::Type (index, MinimalCtor());
	}

///.Function.goEnd.param.iterator.type:Spec.BottomUp Iterator
///.Function.goEnd.param.iterator.type:Spec.TopDownHistory Iterator
	template < typename TText, typename TIndexSpec, class TSpec >
	inline void goEnd(Iter<Index<TText, Index_ESA<TIndexSpec> >, VSTree<TSpec> > &it) 
	{
		_historyClear(it);
		clear(it);
	}

	template < typename TText, typename TIndexSpec, class TSpec >
	inline void goEnd(Iter<Index<TText, Index_ESA<TIndexSpec> >, VSTree< BottomUp<TSpec> > > &it) 
	{
		_dfsClear(it);
		clear(it);
	}

//____________________________________________________________________________

///.Function.goNext.param.iterator.type:Spec.BottomUp Iterator
///.Function.goNext.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, typename TSpec >
	inline void goNext(Iter<TIndex, VSTree<TSpec> > &it) {
		goNext(it, typename GetVSTreeIteratorTraits< Iter<TIndex, VSTree<TSpec> > >::Type());
	}

	template < typename TIndex, typename TSpec, typename TTraits >
	inline void goNext(Iter<TIndex, VSTree<TSpec> > &it, TTraits const traits) {
		goNextImpl(it, traits);
	}


/**
.Function.goDown:
..summary:Iterates down one edge or a path in a tree.
..cat:Index
..signature:bool goDown(iterator)
..signature:bool goDown(iterator, char)
..signature:bool goDown(iterator, text[, lcp])
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDown Iterator
..param.char:$iterator$ goes down the edge beginning with $char$.
..param.text:$iterator$ goes down the path representing $text$. If $text$ ends within an edge, $iterator$ will point to the child-end of this edge.
..param.lcp:A reference of a size type. When $goDown$ returns, $lcp$ contains the length of the longest-common-prefix of $text$ and a path beginning at the $iterator$ node.
...type:Class.String
...type:Class.Segment
..remarks:$goDown(iterator)$ goes down the leftmost edge in the Suffix Tree, i.e. the edge beginning with the lexicographically smallest character.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
*/

    //////////////////////////////////////////////////////////////////////////////
	// unified history stack access for goDown(..)
       
	template < typename TIndex, class TSpec >
	inline void 
	_historyClear(Iter< TIndex, VSTree<TSpec> > &) {}

	template < typename TIndex, class TSpec >
	inline void 
	_historyClear(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		clear(it.history);
	}
/*
	template < typename TText, class TIndexSpec, class TSpec >
	inline void 
	_historyPush(Iter< Index<TText, Index_ESA<TIndexSpec> >, VSTree<TSpec> > &it) {
		value(it).parentRight = value(it).range.i2;
	}
*/	template < typename TText, class TIndexSpec, class TSpec >
	inline void 
	_historyPush(Iter< Index<TText, Index_ESA<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it) {
		it._parentDesc = value(it);
		value(it).parentRight = value(it).range.i2;
	}
	template < typename TText, class TIndexSpec, class TSpec >
	inline void 
	_historyPush(Iter< Index<TText, Index_ESA<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		value(it).parentRight = value(it).range.i2;
		push(it.history, value(it).range);
	}


	// standard down/right/up handlers of top-down-traversal
	template < typename TIndex, typename TSpec >
	inline void _onGoDown(Iter<TIndex, VSTree< TopDown<TSpec> > > &) {}

	template < typename TIndex, typename TSpec >
	inline void _onGoRight(Iter<TIndex, VSTree< TopDown<TSpec> > > &) {}

	template < typename TIndex, typename TSpec >
	inline void _onGoUp(Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &) {}


    //////////////////////////////////////////////////////////////////////////////
	// goDown

	// go down the leftmost edge (including empty $-edges)
	template < typename TText, class TIndexSpec, class TSpec, typename TDFSOrder >
	inline bool _goDown(
		Iter< Index<TText, Index_ESA<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, False> const)
	{
		typedef Index<TText, Index_ESA<TIndexSpec> >	TIndex;

		if (_isLeaf(it, EmptyEdges())) return false;
		_historyPush(it);

		TIndex const &index = container(it);

		typename Size<TIndex>::Type lval = _getUp(value(it).range.i2, index);
		if (!(value(it).range.i1 < lval && lval < value(it).range.i2))
			lval = _getDown(value(it).range.i1, index);
		value(it).range.i2 = lval;
		return true;
	}

	// go down the leftmost edge (skip empty $-edges)
	template < typename TText, class TIndexSpec, class TSpec, typename TDFSOrder >
	inline bool _goDown(
		Iter< Index<TText, Index_ESA<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, True> const)
	{
		typedef Index<TText, Index_ESA<TIndexSpec> >	TIndex;
		
		if (_isLeaf(it, EmptyEdges())) return false;
		_historyPush(it);

		TIndex const &index = container(it);

		typename Size<TIndex>::Type lval = _getUp(value(it).range.i2, index);
		if (!(value(it).range.i1 < lval && lval < value(it).range.i2))
			lval = _getDown(value(it).range.i1, index);
		value(it).range.i2 = lval;

		typename Size<TIndex>::Type lcp = lcpAt(lval - 1, index);
		//typename typename StringSetLimits<TIndex const>::Type &limits = stringSetLimits(index);
		
		typename SAValue<TIndex>::Type pos = getOccurrence(it);
		if (getSeqOffset(pos, stringSetLimits(index)) + lcp == sequenceLength(getSeqNo(pos, stringSetLimits(index)), index)
			|| !nodeHullPredicate(it)) 
		{
			if (!goRight(it)) {
				_goUp(it);
				return false;
			}
		}
		return true;
	}

	// go down the leftmost edge
	template < typename TIndex, class TSpec >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) {
		if (_goDown(it, typename GetVSTreeIteratorTraits< Iter<TIndex, VSTree< TopDown<TSpec> > > >::Type())) {
			_onGoDown(it);
			return true;
		} else
			return false;
	}


    //////////////////////////////////////////////////////////////////////////////
	// goDown a specific edge (chosen by the first character)

	// go down the edge beginning with c (returns false iff this edge doesn't exists)
	template < typename TIndex, class TSpec, typename TValue >
	inline bool _goDownChar(Iter< TIndex, VSTree< TopDown<TSpec> > > &it, TValue c) 
	{
		typename VertexDescriptor<TIndex>::Type nodeDesc;
		if (_getNodeByChar(it, c, nodeDesc)) {
			_historyPush(it);
			value(it) = nodeDesc;
			return true;
		}
		return false;
	}

	// go down the path corresponding to pattern
	// lcp is the longest prefix of pattern and path
	template < typename TIndex, typename TSpec, typename TString, typename TSize >
	inline bool
	_goDownString(
		Iter< TIndex, VSTree< TopDown<TSpec> > > &node,
		TString const &pattern, 
		TSize &lcp) 
	{
		typedef typename Fibre<TIndex, Fibre_Text>::Type const		TText;
		typedef typename Infix<TText>::Type							TInfix;
		typedef typename Iterator<TInfix, Standard>::Type			IText;
		typedef typename Iterator<TString const, Standard>::Type	IPattern;
		
		IPattern p_begin = begin(pattern, Standard()), p_end = end(pattern, Standard());
		IText t_begin, t_end;

		if (p_begin == p_end) {
			lcp = 0;
			return true;
		}

		TSize parentRepLen = 0;
		// go down the edge beginning with a pattern character
		while (_goDownChar(node, *p_begin))
		{
			TInfix t = representative(node);
			t_begin = begin(t, Standard()) + parentRepLen;
			t_end = end(t, Standard());

			while (t_begin != t_end && p_begin != p_end) 
			{
				// compare each character along the edge
				if (*p_begin != *t_begin) {
					lcp = p_begin - begin(pattern, Standard());
					return false;
				}
				++t_begin;
				++p_begin;
			}

			// was the whole pattern found?
			if (p_begin == p_end) {
				lcp = length(pattern);
				return true;
			}
			parentRepLen = length(t);
		}
		lcp = p_begin - begin(pattern, Standard());
		return false;
	}

	template < typename TIndex, typename TSpec, typename TObject >
	inline bool 
	_goDownObject(
		Iter< TIndex, VSTree< TopDown<TSpec> > > &it, 
		TObject const &obj,
		False)
	{
		return _goDownChar(it, obj);
	}

	template < typename TIndex, typename TSpec, typename TObject >
	inline bool 
	_goDownObject(
		Iter< TIndex, VSTree< TopDown<TSpec> > > &it, 
		TObject const &obj,
		True)
	{
		typename Size<TIndex>::Type dummy;
		return _goDownString(it, obj, dummy);
	}


	// public interface for goDown(it, ...)
	template < typename TIndex, typename TSpec, typename TObject >
	inline bool
	goDown(
		Iter< TIndex, VSTree< TopDown<TSpec> > > &it, 
		TObject const &obj) 
	{
		return _goDownObject(it, obj, typename IsSequence<TObject>::Type());
	}

	template < typename TIndex, typename TSpec, typename TString, typename TSize >
	inline bool 
	goDown(
		Iter< TIndex, VSTree< TopDown<TSpec> > > &it, 
		TString const &pattern,
		TSize &lcp)
	{
		return _goDownString(it, pattern, lcp);
	}

		
/**
.Function.goUp:
..summary:Iterates up one edge to the parent in a tree.
..cat:Index
..signature:goUp(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
*/

	// go up one edge (returns false if in root node)
	// can be used at most once, as no history stack is available
	template < typename TIndex, class TSpec >
	inline bool 
	_goUp(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) 
	{
		if (!isRoot(it)) {
			value(it) = it._parentDesc;
			return true;
		}
		return false;
	}

	// go up one edge (returns false if in root node)
	template < typename TIndex, class TSpec >
	inline bool 
	_goUp(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) 
	{
		if (!empty(it.history)) {
			value(it).range = top(it.history);
			pop(it.history);
			if (!empty(it.history))
				value(it).parentRight = top(it.history).i2;	// copy right boundary of parent's range
			return true;
		}
		return false;
	}

	// go up one edge
	template < typename TIndex, class TSpec >
	inline bool goUp(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		if (_goUp(it)) {
			_onGoUp(it);
			return true;
		} else
			return false;
	}

	// return vertex descriptor of parent's node
	template < typename TIndex, class TSpec >
	inline typename VertexDescriptor<TIndex>::Type
	nodeUp(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > const &it) 
	{
		if (!empty(it.history)) {
			typename Size<TIndex>::Type parentRight = 0;
			if (length(it.history) > 2)
				parentRight = topPrev(it.history).i2;
			return typename VertexDescriptor<TIndex>::Type(top(it.history), parentRight);
		} else
			return value(it);
	}

	// nodeUp adaption for non-history iterators
	// ATTENTION: Do not call nodeUp after a goDown that returned false (or after _goUp)!
	template < typename TIndex, class TSpec >
	inline typename VertexDescriptor<TIndex>::Type const &
	nodeUp(Iter< TIndex, VSTree< TopDown<TSpec> > > const &it) 
	{
		return it._parentDesc;
	}

/**
.Function.goRight:
..summary:Iterates to the next sibling in a tree.
..cat:Index
..signature:goRight(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDown Iterator
*/

	// go right to the lexic. next sibling
	template < typename TText, class TIndexSpec, class TSpec, typename TDFSOrder, typename THideEmptyEdges >
	inline bool _goRight(
		Iter< Index<TText, Index_ESA<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it, 
		VSTreeIteratorTraits<TDFSOrder, THideEmptyEdges> const) 
	{
		typedef Index<TText, Index_ESA<TIndexSpec> > TIndex;

		if (isRoot(it)) return false;		

		typename Size<TIndex>::Type right = value(it).parentRight;
		if (_isSizeInval(right)) right = length(indexSA(container(it)));

		do {
			if (value(it).range.i2 == right)				// not the right-most child?
				return false;

			if (_isNextl(value(it).range.i2, container(it))) 
			{
				value(it).range.i1 = value(it).range.i2;	// go right
				value(it).range.i2 = _getNextl(value(it).range.i2, container(it));
			} else {
				value(it).range.i1 = value(it).range.i2;	// now it is the right-most child
				value(it).range.i2 = value(it).parentRight;
			}

		} while ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it));
		return true;
	}

	// go down the leftmost edge
	template < typename TIndex, class TSpec >
	inline bool goRight(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) {
		if (_goRight(it, typename GetVSTreeIteratorTraits< Iter<TIndex, VSTree< TopDown<TSpec> > > >::Type())) {
			_onGoRight(it);
			return true;
		} else
			return false;
	}

/**
.Function.parentEdgeLength:
..summary:Returns the length of the edge from the $iterator$ node to its parent.
..cat:Index
..signature:parentEdgeLength(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
..returns:The returned value is equal to $length(parentEdgeLabel(iterator))$.
*/

	template < typename TText, class TIndexSpec, class TSpec >
	inline typename Size< Index<TText, Index_ESA<TIndexSpec> > >::Type
	parentEdgeLength(Iter< 
		Index<TText, Index_ESA<TIndexSpec> >, 
		VSTree< TopDown< ParentLinks<TSpec> > > > const &it) 
	{
		return repLength(it) - parentRepLength(it);
	}

/**
.Function.parentEdgeLabel:
..summary:Returns a substring representing the edge from an $iterator$ node to its parent.
..cat:Index
..signature:parentEdgeLabel(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
..returns:An @Spec.InfixSegment@ of the raw text of an index (see @Tag.ESA Index Fibres.ESA_RawText@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_RawText>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, Fibre_Text>::Type const >::Type
	parentEdgeLabel(Iter< TIndex, VSTree< TopDown<TSpec> > > const &it)
	{
		return infixWithLength(
			indexText(container(it)), 
			posAdd(getOccurrence(it), parentRepLength(it)),
			parentEdgeLength(it));
	}

/**
.Function.parentEdgeFirstChar:
..summary:Returns the first character of the edge from an $iterator$ node to its parent.
..cat:Index
..signature:parentEdgeFirstChar(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
..returns:A single character of type $Value<TIndex>::Type$ which is identical to $Value<Fibre<TIndex, ESA_RawText>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Value<TIndex>::Type 
	parentEdgeFirstChar(Iter< TIndex, VSTree<TSpec> > const &it) 
	{
		return infixWithLength(
			indexText(container(it)),
			posAdd(getOccurrence(it), parentRepLength(it)),
			1)[0];
	}


	template < typename TIndex, class TSpec >
	inline void clear(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		value(it) = typename VertexDescriptor<TIndex>::Type(MinimalCtor());
    }

	template < typename TIndex, class TSpec >
	inline void _dfsClear(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		clear(it.history);
    }


    //////////////////////////////////////////////////////////////////////////////
	// dfs traversal for ParentLink iterators

	template < typename TIndex, typename TSpec, typename THideEmptyEdges >
	inline void goNextImpl(
		Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, 
		VSTreeIteratorTraits<_Preorder, THideEmptyEdges> const)
	{
		// preorder dfs
		do {
			if (!goDown(it) && !goRight(it))
				while (goUp(it) && !goRight(it));
			if (isRoot(it)) {
				clear(it);
				return;
			}
		} while (!nodePredicate(it));
	}

	template < typename TIndex, typename TSpec, typename THideEmptyEdges >
	inline void goNextImpl(
		Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, 
		VSTreeIteratorTraits<_Postorder, THideEmptyEdges> const)
	{
		// postorder dfs
		do {
			if (goRight(it))
				while (goDown(it));
			else
				if (!goUp(it)) {
					clear(it);
					return;
				}
		} while (!nodePredicate(it));
	}


    //////////////////////////////////////////////////////////////////////////////
	// boolean functions

	template < typename TIndex, class TSpec >
	inline bool eof(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		return !value(it).range.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool eof(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return !value(it).range.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool empty(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		return !value(it).range.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool empty(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return !value(it).range.i2;
	}

///.Function.atEnd.param.iterator.type:Spec.BottomUp Iterator
///.Function.atEnd.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, class TSpec >
	inline bool atEnd(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		return !value(it).range.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool atEnd(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return !value(it).range.i2;
	}

/**
.Function.isRoot:
..summary:Test whether iterator points to the root node.
..cat:Index
..signature:bool isRoot(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to the root of the tree, otherwise $false$.
*/

	template < typename TIndex, class TSpec >
	inline bool isRoot(Iter<TIndex, VSTree< BottomUp<TSpec> > > const &it) 
	{
		return empty(it.history);
	}

	template < typename TIndex, class TSpec >
	inline bool isRoot(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return _isRoot(value(it));
	}

	template < typename TSize >
	inline bool _isRoot(VertexESA<TSize> const &value) 
	{
		return _isSizeInval(value.range.i2);
	}

/**
.Function.isRightTerminal:
..summary:Test whether iterator points to a suffix.
..cat:Index
..signature:bool isRightTerminal(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to the node representing a suffix, otherwise $false$.
..remarks:Every leaf is also a right terminal (see @Function.Index#isLeaf@), but not vice versa.
*/

	template < typename TIndex, class TSpec >
	inline bool isRightTerminal(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		// do we reach a leaf in a suffix tree with trailing '$'
		typename SAValue<TIndex>::Type pos = getOccurrence(it);
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		return (getSeqOffset(pos, limits) + repLength(it) 
			== sequenceLength(getSeqNo(pos, limits), index));
	}

/**
.Function.isLeftMaximal:
..summary:Test whether the occurences of an iterator's @Function.representative@ mutually differ in the character left of the hits.
..cat:Index
..signature:bool isLeftMaximal(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if there are at least two different characters left of the occurences, otherwise $false$.
..see:Function.getOccurrences
*/

	template < typename TIndex, class TSpec >
	inline bool isLeftMaximal(Iter<TIndex, VSTree<TSpec> > const &it)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type	TOccs;
		typedef typename Infix< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TOccsBWT;
		typedef typename Value< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TValue;

		typedef typename Iterator<TOccs, Standard>::Type	TIter;
		typedef typename Iterator<TOccsBWT, Standard>::Type TIterBWT;
		
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		TOccs occs = getOccurrences(it);
		TOccsBWT bwts = getOccurrencesBWT(it);

		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
		TIterBWT bw = begin(bwts, Standard());

		if (oc == ocEnd) return true;
		if (posAtFirstLocal(*oc, limits)) return true;

		TValue seen = *bw;
		++oc; 
		++bw;
		if (oc == ocEnd) return true;

		do {
			if (posAtFirstLocal(*oc, limits)) return true;
			if (seen != *bw) return true;
			++oc;
			++bw;
		} while (oc != ocEnd);

		return false;
	}

/**
.Function.isPartiallyLeftExtensible:
..summary:Test whether the characters left of the two occurences of @Function.representative@ are equal.
..cat:Index
..signature:bool isPartiallyLeftExtensible(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if there are at least two different characters left of the occurences, otherwise $false$.
..see:Function.getOccurrences
*/

	template < typename TIndex, class TSpec, typename TSet >
	inline bool isPartiallyLeftExtensible(Iter<TIndex, VSTree<TSpec> > const &it, TSet &charSet)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type	TOccs;
		typedef typename Infix< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TOccsBWT;
		typedef typename Value< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TValue;

		typedef typename Iterator<TOccs, Standard>::Type	TIter;
		typedef typename Iterator<TOccsBWT, Standard>::Type TIterBWT;
		
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		clear(charSet);

		TOccs occs = getOccurrences(it);
		TOccsBWT bwts = getOccurrencesBWT(it);

		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
		TIterBWT bw = begin(bwts, Standard());

		while (oc != ocEnd) {
			if (!posAtFirstLocal(*oc, limits)) {
				TValue c = *bw;
				if (in(c, charSet)) return true;
				insert(c, charSet);
			}
			++oc;
			++bw;
		}

		return false;
	}

	template < typename TIndex, class TSpec >
	inline bool isPartiallyLeftExtensible(Iter<TIndex, VSTree<TSpec> > const &it)
	{
		typename Set<typename Value<TIndex>::Type>::Type set;
		return isPartiallyLeftExtensible(it, set);
	}

/**
.Function.isUnique:
..summary:Test whether the @Function.representative@ occurs only once in every sequence.
..cat:Index
..signature:bool isUnique(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if there are at least two different characters left of the occurences, otherwise $false$.
..see:Function.getOccurrences
*/

	template < typename TIndex, class TSpec, typename TSet >
	inline bool isUnique(Iter<TIndex, VSTree<TSpec> > const &it, TSet &set)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;
		typedef typename Size<TIndex>::Type TSize;

		TIndex const &index = container(it);

		clear(set);

		TOccs occs = getOccurrences(it);
		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());

		while (oc != ocEnd) {
			TSize seqNo = getSeqNo(*oc, stringSetLimits(index));
			if (in(seqNo, set)) return false;
			insert(seqNo, set);
			++oc;
		}

		return true;
	}

	template < typename TIndex, class TSpec >
	inline bool isUnique(Iter<TIndex, VSTree<TSpec> > const &it) {
		_VectorSet<
			typename Size<TIndex>::Type,
			Alloc<> 
		> set(countSequences(container(it)));
		return isUnique(it, set);
	}

/**
.Function.getFrequency:
..summary:Returns the number of sequences, which contain the @Function.representative@ as a substring.
..cat:Index
..signature:int getFrequency(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The number of different sequences containing the @Function.representative@.
..see:Function.getOccurrences
*/

	template < typename TIndex, class TSpec, typename TSet >
	inline typename Size<TIndex>::Type
	getFrequency(Iter<TIndex, VSTree<TSpec> > const &it, TSet &set)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;
		typedef typename Size<TIndex>::Type TSize;

		TIndex const &index = container(it);

		clear(set);

		TOccs occs = getOccurrences(it);
		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());

		int counter = 0;
		while (oc != ocEnd) {
			TSize seqNo = getSeqNo(*oc, stringSetLimits(index));
			if (!in(seqNo, set)) {
				++counter;
				insert(seqNo, set);
			}
			++oc;
		}

		return counter;
	}

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type
	getFrequency(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		_VectorSet<
			typename Size<TIndex>::Type,
			Alloc<> 
		> set(countSequences(container(it)));
		return getFrequency(it, set);
	}

/**
.Function.childrenAreLeaves:
..summary:Test whether iterator points to a node with only leaf-children.
..cat:Index
..signature:bool childrenAreLeaves(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to an inner node of the tree, whose children are leaves. Otherwise it is $false$.
*/

	template < typename TIndex, class TSpec >
	inline bool childrenAreLeaves(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return countChildren(it) == countOccurrences(it);
	}

/**
.Function.Index#isLeaf:
..summary:Test whether iterator points to a leaf.
..cat:Index
..signature:bool isLeaf(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to a leaf of the tree, otherwise $false$.
*/

	template < typename TSize >
	inline bool _isLeaf(VertexESA<TSize> const &vDesc)
	{
		// is this a leaf?
		return vDesc.range.i1 + 1 >= vDesc.range.i2;
	}

	// is this a leaf? (including empty $-edges)
	template < typename TIndex, class TSpec, typename TDFSOrder >
	inline bool _isLeaf(
		Iter<TIndex, VSTree<TSpec> > const &it,
		VSTreeIteratorTraits<TDFSOrder, False> const)
	{
		return _isLeaf(value(it));
	}

	// is this a leaf? (hide empty $-edges)
	template < typename TIndex, class TSpec, typename TDFSOrder >
	inline bool _isLeaf(
		Iter<TIndex, VSTree<TSpec> > const &it,
		VSTreeIteratorTraits<TDFSOrder, True> const)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;

		if (_isLeaf(value(it))) return true;

		TIndex const &index = container(it);

		// get representative length (see repLength)
		typename Size<TIndex>::Type lval = _getUp(value(it).range.i2, index);
		if (!(value(it).range.i1 < lval && lval < value(it).range.i2))
			lval = _getDown(value(it).range.i1, index);
		typename Size<TIndex>::Type lcp = lcpAt(lval - 1, index);

		// test suffices in range for empty edges
		TOccs occs = getOccurrences(it);
		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
		for (; oc != ocEnd; ++oc)
			if (getSeqOffset(*oc, stringSetLimits(index)) + lcp != 
				sequenceLength(getSeqNo(*oc, stringSetLimits(index)), index))
				return false;
		return true;
	}

	template < typename TIndex, class TSpec >
	inline bool isLeaf(Iter<TIndex, VSTree<TSpec> > const &it)
	{
		return _isLeaf(it, typename GetVSTreeIteratorTraits< Iter<TIndex, VSTree<TSpec> > >::Type());
	}


	//////////////////////////////////////////////////////////////////////////////
	// (more or less) internal functions for accessing the childtab

	template < typename TSize, typename TIndex >
	inline bool _isNextl(TSize i, TIndex const &index) 
	{
		if (i >= length(index)) return false;
		TSize j = childAt(i, index);
		return (j > i) && lcpAt(j - 1, index) == lcpAt(i - 1, index);
	}

	template < typename TSize, typename TIndex >
	inline bool _isUp(TSize i, TIndex const &index) 
	{
		if (i >= length(index)) return false;
		TSize j = childAt(i, index);
		return (j <= i) && lcpAt(j - 1, index) > lcpAt(i - 1, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getNextl(TSize i, TIndex const &index) 
	{
		return childAt(i, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getUp(TSize i, TIndex const &index) 
	{
		if (!_isSizeInval(i))
			return childAt(i - 1, index);
		else
			return childAt(0, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getDown(TSize i, TIndex const &index) 
	{
		return childAt(i, index);
	}


	//////////////////////////////////////////////////////////////////////////////
	// depth-first search 

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> &
	_dfsRange(Iter< TIndex, VSTree< BottomUp<TSpec> > > &it)
	{
		return value(it).range;
	}

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> const & 
	_dfsRange(Iter< TIndex, VSTree< BottomUp<TSpec> > > const &it) 
	{
		return value(it).range;
	}

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type & _dfsLCP(Iter< TIndex, VSTree< BottomUp<TSpec> > > &it)
	{
		return it.lValue;
	}

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type _dfsLCP(Iter< TIndex, VSTree< BottomUp<TSpec> > > const &it)
	{
		return it.lValue;
	}


}

#endif
