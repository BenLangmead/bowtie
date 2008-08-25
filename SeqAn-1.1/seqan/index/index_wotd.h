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
  $Id: index_wotd.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_WOTD_H
#define SEQAN_HEADER_INDEX_WOTD_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// wotd tree index fibres

	typedef Fibre_Text		Wotd_Text;
	typedef Fibre_RawText	Wotd_RawText;
	typedef Fibre_SA		Wotd_SA;
	typedef Fibre_RawSA		Wotd_RawSA;
	typedef Fibre_Dir		Wotd_Dir;


//////////////////////////////////////////////////////////////////////////////
// wotd tree index

/**
.Spec.Index_Wotd:
..summary:An index based on a lazy suffix tree (see Giegerich et al., "Efficient implementation of lazy suffix trees").
..cat:Index
..general:Class.Index
..signature:Index<TText, Index_Wotd<> >
..param.TText:The text type.
...type:Class.String
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a partially sorted suffix array (see @Tag.Wotd_SA@) and the wotd tree (see @Tag.Wotd_Dir@).
*/

	struct WotdOriginal_;
	typedef Tag<WotdOriginal_> const WotdOriginal;

	template < typename TSpec = void >
	struct Index_Wotd;

/*
	template < typename TObject, typename TSpec >
	struct Fibre< Index<TObject, Index_Wotd<TSpec> >, Fibre_Dir> 
	{
		typedef Index<TObject, Index_Wotd<TSpec> > TIndex;
		typedef String< 
			typename typename Size<TIndex>::Type,
			Alloc<>
		> Type;
	};
*/
	
	template < typename TObject, typename TSpec >
	class Index<TObject, Index_Wotd<TSpec> > {
	public:
		typedef typename Fibre<Index, Wotd_Text>::Type		TText;
		typedef typename Fibre<Index, Wotd_SA>::Type		TSA;
		typedef typename Fibre<Index, Wotd_Dir>::Type		TDir;

		typedef typename Value<Index>::Type					TValue;
		typedef typename Size<TText>::Type					TSize;
		typedef String<TSize, Alloc<> >						TCounter;
		typedef String<typename Value<TSA>::Type, Alloc<> >	TTempSA;
		typedef typename Cargo<Index>::Type					TCargo;

		// 1st word flags
		static TSize const LEAF          = (TSize)1 << (BitsPerValue<TSize>::VALUE - 1); // this node is a leaf
		static TSize const LAST_CHILD    = (TSize)1 << (BitsPerValue<TSize>::VALUE - 2); // this node is the last child
		// 2nd word flag
		static TSize const UNEVALUATED   = (TSize)1 << (BitsPerValue<TSize>::VALUE - 1); // this node is partially evalutated and has no evaluated children
		static TSize const SENTINELS     = (TSize)1 << (BitsPerValue<TSize>::VALUE - 2); // the children of this node have solely $-edges

		static TSize const BITMASK0      = ~(LEAF | LAST_CHILD);
		static TSize const BITMASK1      = ~(UNEVALUATED | SENTINELS);


		Holder<TText>	text;	// underlying text
		TSA				sa;		// suffix array sorted by the first q chars
		TDir			dir;	// bucket directory
		TCargo			cargo;	// user-defined cargo


		TTempSA			tempSA;
		TCounter		tempOcc;
		TCounter		tempBound;

		TSize			sentinelOcc;
		TSize			sentinelBound;

		Index() {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text) {}
	};
/*
    template < typename TText, typename TSpec >
    struct Value< Index<TText, Index_Wotd<TSpec> > > {
		typedef typename Value< typename Fibre< Index<TText, Index_Wotd<TSpec> >, Wotd_RawText >::Type >::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, Index_Wotd<TSpec> > > {
		typedef typename Size< typename Fibre< Index<TText, Index_Wotd<TSpec> >, Wotd_RawText >::Type >::Type Type;
    };
*/


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, Index_Wotd<TSpec> >, Fibre_SA> {
        typedef Default Type;
    };

	template < typename TText, typename TSSetSpec, typename TSpec >
	struct DefaultIndexCreator<Index<StringSet<TText, TSSetSpec>, Index_Wotd<TSpec> >, Fibre_SA> {
        typedef Default Type;
    };


//////////////////////////////////////////////////////////////////////////////


	template <typename TSize>
	struct _VertexWotdOriginal {
		TSize		node;			// position of current node entry in directory
		TSize		parentRepLen;	// representative length of parent node
		TSize		edgeLen;		// length of edge above current node

		_VertexWotdOriginal() {}
		_VertexWotdOriginal(MinimalCtor):
			parentRepLen(0),
			edgeLen(0) 
		{
			_setSizeInval(node);
		}
	};

	template <typename TSize>
	struct _VertexWotdModified {
		TSize		node;			// position of current node entry in directory
		TSize		parentRepLen;	// representative length of parent node
		TSize		edgeLen;		// length of edge above current node
		Pair<TSize> range;			// current SA interval of hits
		TSize		parentRight;	// right boundary of parent node's range (allows to go right)

		_VertexWotdModified() {}
		_VertexWotdModified(MinimalCtor):
			node(0),
			parentRepLen(0),
			edgeLen(0),
			range(0,0),
			parentRight(0) {}
	};

	template < typename TText >
	struct VertexDescriptor< Index<TText, Index_Wotd<WotdOriginal> > > {
		typedef typename Size< Index<TText, Index_Wotd<WotdOriginal> > >::Type TSize;
		typedef _VertexWotdOriginal<TSize> Type;
	};

	template < typename TText, typename TSpec >
	struct VertexDescriptor< Index<TText, Index_Wotd<TSpec> > > {
		typedef typename Size< Index<TText, Index_Wotd<TSpec> > >::Type TSize;
		typedef _VertexWotdModified<TSize> Type;
	};

	template < typename TText, typename TSpec >
	void _indexRequireTopDownIteration(Index<TText, Index_Wotd<TSpec> > &index) 
	{
		indexRequire(index, Wotd_SA());
	}

//////////////////////////////////////////////////////////////////////////////
// history stack functions

	template <typename TSize>
	struct _HistoryStackWotdOriginal
	{
		TSize		node;		// position of current node entry in directory
		TSize		edgeLen;	// length of edge above current node
	};

	template <typename TSize>
	struct _HistoryStackWotdModified
	{
		TSize		node;		// position of current node entry in directory
		TSize		edgeLen;	// length of edge above current node
		Pair<TSize> range;		// current SA interval of hits
	};

	template < typename TText, typename TSpec >
	struct _HistoryStackEntry< Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown< ParentLinks<TSpec> > > > > 
	{
		typedef Index<TText, Index_Wotd<WotdOriginal> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;
		typedef _HistoryStackWotdOriginal<TSize>		Type;
	};

	template < typename TText, typename TIndexSpec, typename TSpec >
	struct _HistoryStackEntry< Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > >
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;
		typedef _HistoryStackWotdModified<TSize>		Type;
	};


	template < typename TText, typename TIndexSpec, typename TSpec >
	inline void 
	_historyPush(Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it) 
	{
		it._parentDesc = value(it);
		value(it).parentRepLen += parentEdgeLength(it);
		value(it).parentRight = value(it).range.i2;
	}

	template < typename TText, typename TSpec >
	inline void 
	_historyPush(Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it) 
	{
		typedef typename Size< Index<TText, Index_Wotd<WotdOriginal> > >::Type TSize;
		TSize edgeLen = parentEdgeLength(it);
		_HistoryStackWotdOriginal<TSize> entry = { value(it).node, edgeLen };
		push(it.history, entry);
		value(it).parentRepLen += edgeLen;
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline void 
	_historyPush(Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it) 
	{
		typedef typename Size< Index<TText, Index_Wotd<TIndexSpec> > >::Type TSize;
		TSize edgeLen = parentEdgeLength(it);
		_HistoryStackWotdModified<TSize> entry = { value(it).node, edgeLen, value(it).range };
		push(it.history, entry);
		value(it).parentRepLen += edgeLen;
		value(it).parentRight = value(it).range.i2;
	}

//////////////////////////////////////////////////////////////////////////////

	template < typename TSize >
	inline bool _isRoot(_VertexWotdOriginal<TSize> const &value) {
		return value.node == 0;
	}

	template < typename TSize >
	inline bool _isRoot(_VertexWotdModified<TSize> const &value) {
		return value.node == 0; 
	}

	// is this a leaf? (including empty $-edges)
	template < typename TText, typename TIndexSpec, typename TSpec, typename TDFSOrder >
	inline bool _isLeaf(
		Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree<TSpec> > const &it,
		VSTreeIteratorTraits<TDFSOrder, False> const)
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> > TIndex;
		TIndex const &index = container(it);
		return dirAt(value(it).node, index) & index.LEAF;
	}

	// is this a leaf? (excluding empty $-edges)
	template < typename TText, typename TIndexSpec, typename TSpec, typename TDFSOrder >
	inline bool _isLeaf(
		Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree<TSpec> > const &it,
		VSTreeIteratorTraits<TDFSOrder, True> const)
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		TIndex const &index = container(it);
		if (dirAt(value(it).node, index) & index.LEAF)
			return true;

		// ensure node evaluation and test for sentinel child edges?
		return _wotdEvaluate(it) & index.SENTINELS;
	}


	// parentEdgeLength - ORIGINAL VERSION
	template < typename TIndex, typename TSize >
	inline typename Size<TIndex>::Type
	parentEdgeLength(TIndex const &index, _VertexWotdOriginal<TSize> &vDesc)
	{
		TSize edgeLen = vDesc.edgeLen;
		if (edgeLen != (TSize)-1)
			return edgeLen;

		TSize pos = vDesc.node;
		TSize w0 = dirAt(pos, index);
		if (w0 & index.LEAF)
			return vDesc.edgeLen = suffixLength(w0 & index.BITMASK0, index);

		TSize w1 = dirAt(pos + 1, index);
		if (w1 & index.UNEVALUATED)
			return vDesc.edgeLen = _bucketLCP(
				infix(indexSA(index), w0 & index.BITMASK0, w1 & index.BITMASK1),
				indexText(index));
		else
			return vDesc.edgeLen = _getNodeLP(index, w1) - (w0 & index.BITMASK0);
	}

	// parentEdgeLength - MODIFIED VERSION
	template < typename TIndex, typename TSize >
	inline typename Size<TIndex>::Type
	parentEdgeLength(TIndex const &index, _VertexWotdModified<TSize> &vDesc)
	{
		TSize edgeLen = vDesc.edgeLen;
		if (edgeLen != (TSize)-1)
			return edgeLen;

		TSize pos = vDesc.node;
		TSize w0 = dirAt(pos, index);
		if (w0 & index.LEAF)
			return vDesc.edgeLen =
				suffixLength(saAt(vDesc.range.i1, index), index) - vDesc.parentRepLen;

		TSize w1 = dirAt(pos + 1, index);
		if (w1 & index.UNEVALUATED)
			if (_isSizeInval(vDesc.range.i2))
				return vDesc.edgeLen = _bucketLCP(
					suffix(indexSA(index), vDesc.range.i1),
					indexText(index),
					vDesc.parentRepLen) - vDesc.parentRepLen;
			else
				return vDesc.edgeLen = _bucketLCP(
					infix(indexSA(index), vDesc.range.i1, vDesc.range.i2),
					indexText(index),
					vDesc.parentRepLen) - vDesc.parentRepLen;
		else
			return (dirAt(w1, index) & index.BITMASK0) - vDesc.parentRepLen;
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline typename Size< Index<TText, Index_Wotd<TIndexSpec> > >::Type
	parentEdgeLength(Iter<
		Index<TText, Index_Wotd<TIndexSpec> >, 
		VSTree< TopDown<TSpec> > > const &it)
	{
		typedef Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > TIter;
		return parentEdgeLength(container(it), value(const_cast<TIter&>(it)));
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline typename Size< Index<TText, Index_Wotd<TIndexSpec> > >::Type
	parentRepLength(Iter<
		Index<TText, Index_Wotd<TIndexSpec> >, 
		VSTree< TopDown<TSpec> > > const &it)
	{
		return value(it).parentRepLen;
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline typename Size< Index<TText, Index_Wotd<TIndexSpec> > >::Type
	parentRepLength(Iter<
		Index<TText, Index_Wotd<TIndexSpec> >, 
		VSTree< TopDown< ParentLinks<TSpec> > > > const &it)
	{
		return value(it).parentRepLen;
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline typename Size< Index<TText, Index_Wotd<TIndexSpec> > >::Type
	repLength(Iter<
		Index<TText, Index_Wotd<TIndexSpec> >, 
		VSTree< TopDown<TSpec> > > const &it) 
	{
		return parentRepLength(it) + parentEdgeLength(it);
	}


	// parentEdgeLabel - ORIGINAL VERSION
	template < typename TText, typename TSpec >
	inline typename Infix< typename Fibre<Index<TText, Index_Wotd<WotdOriginal> >, ESA_RawText>::Type const >::Type 
	parentEdgeLabel(Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > const &it) 
	{
		typedef Index<TText, Index_Wotd<WotdOriginal> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		TIndex const &index = container(it);

		if (isRoot(it))
			return infix(indexRawText(index), 0, 0);
		else {
			TSize occ = _getNodeLP(index, value(it).node);
			return infix(indexRawText(index), occ, occ + parentEdgeLength(it));
		}
	}

	// getOccurrence - ORIGINAL VERSION
	template < typename TText, typename TSpec >
	inline typename SAValue<Index<TText, Index_Wotd<WotdOriginal> > >::Type 
	getOccurrence(Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree<TSpec> > const &it) {
		return _getNodeLP(container(it), value(it).node) - value(it).parentRepLen;
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline bool
	emptyParentEdge(Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it) 
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> > TIndex;

		TIndex const &index = container(it);
		typename SAValue<TIndex>::Type pos = getOccurrence(it);
		return getSeqOffset(pos, stringSetLimits(index)) + value(it).parentRepLen
			== sequenceLength(getSeqNo(pos, stringSetLimits(index)), index);
	}

	// to avoid ambiguity
	template < typename TText, typename TIndexSpec, typename TSpec >
	inline bool
	emptyParentEdge(Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > const &it) 
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> > TIndex;

		TIndex const &index = container(it);
		typename SAValue<TIndex>::Type pos = getOccurrence(it);
		return getSeqOffset(pos, stringSetLimits(index)) + value(it).parentRepLen
			== sequenceLength(getSeqNo(pos, stringSetLimits(index)), index);
	}



	template < typename TText, typename TSpec >
	inline void 
	goRoot(Iter<
		Index<TText, Index_Wotd<WotdOriginal> >, 
		VSTree<TSpec> > &it) 
	{
		_historyClear(it);
		value(it).node = 0;			// start in root node (first entry in dir)
		value(it).parentRepLen = 0;	// parent prefix length is 0
		value(it).edgeLen = 0;		// edge length is 0
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline void 
	goRoot(Iter<
		Index<TText, Index_Wotd<TIndexSpec> >, 
		VSTree<TSpec> > &it) 
	{
		_historyClear(it);
		value(it).range.i1 = 0;		// start in root node with range (0,infty)
		_setSizeInval(value(it).range.i2);	// infty is equivalent to length(index) and faster to compare
		value(it).node = 0;			// start in root node (first entry in dir)
		value(it).parentRepLen = 0;	// parent prefix length is 0
		value(it).edgeLen = 0;		// edge length is 0
	}

	template < typename TText, typename TSpec >
	inline bool atEnd(Iter<Index<TText, Index_Wotd<WotdOriginal> >, VSTree<TSpec> > &it) 
	{
		return _isSizeInval(value(it).node);
	}

	template < typename TText, typename TSpec >
	inline bool atEnd(Iter<Index<TText, Index_Wotd<WotdOriginal> >, VSTree<TSpec> > const &it) 
	{
		return _isSizeInval(value(it).node);
	}

		
	// adjust iterator's right border of SA range
	template < typename TText, typename TSpec >
	inline void
	_adjustRightBorder(
		Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &) 
	{}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline void
	_adjustRightBorder(
		Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it)
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		TIndex	const &index = container(it);
		TSize	pos = value(it).node;

		TSize w0 = dirAt(pos, index);
		if (w0 & index.LEAF)
			value(it).range.i2 = value(it).range.i1 + 1;
		else
			if (w0 & index.LAST_CHILD)
				value(it).range.i2 = value(it).parentRight;
			else {
				w0 = dirAt(pos + 2, index);
				value(it).range.i2 = w0 & index.BITMASK0;
			}
	}

	// go down the leftmost edge (including empty $-edges)
	template < typename TText, typename TSpec, typename TDFSOrder, typename THideEmptyEdges >
	inline bool 
	_goDown(
		Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, THideEmptyEdges> const)
	{
		typedef Index<TText, Index_Wotd<WotdOriginal> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		if (_isLeaf(it, EmptyEdges())) return false;

		TIndex &index = container(it);
		_historyPush(it);

		// ensure node evaluation
		TSize childNode = _wotdEvaluate(it);

		if (THideEmptyEdges::VALUE && (childNode & index.SENTINELS) != 0)
			return false;

		// go down
		value(it).node = childNode & index.BITMASK1;
		value(it).edgeLen = -1;

		// go right if parent edge is empty 
		// or hull predicate is false
		if ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it))
			if (!goRight(it)) {
				_goUp(it);
				return false;
			}

		return true;
	}

	// go down the leftmost edge (excluding empty $-edges)
	template < typename TText, typename TIndexSpec, typename TSpec, typename TDFSOrder, typename THideEmptyEdges >
	inline bool 
	_goDown(
		Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, THideEmptyEdges> const)
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		if (_isLeaf(it, EmptyEdges())) return false;
		TIndex const &index = container(it);

		// ensure node evaluation
		TSize childNode = _wotdEvaluate(it);

		if (THideEmptyEdges::VALUE && (childNode & index.SENTINELS) != 0)
			return false;

		// go down
		_historyPush(it);
		value(it).node = childNode & index.BITMASK1;
		value(it).edgeLen = -1;
		_adjustRightBorder(it);

		// go right if parent edge is empty
		// or hull predicate is false
		if ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it))
			if (!goRight(it)) {
				_goUp(it);
				return false;
			}

		return true;
	}

	// go right to the lexic. next sibling
	template < typename TText, typename TSpec, typename TDFSOrder, typename THideEmptyEdges >
	inline bool 
	_goRight(
		Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, THideEmptyEdges> const) 
	{
		typedef Index<TText, Index_Wotd<WotdOriginal> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		TIndex const &index = container(it);

		do {
			TSize w0 = dirAt(value(it).node, index);
			if (w0 & index.LAST_CHILD)
				return false;

			value(it).node += (w0 & index.LEAF)? 1: 2;
			value(it).edgeLen = -1;

			_adjustRightBorder(it);
		} while ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it));
		return true;
	}

	template < typename TText, typename TIndexSpec, typename TSpec, typename TDFSOrder, typename THideEmptyEdges >
	inline bool 
	_goRight(
		Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, THideEmptyEdges> const) 
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		TIndex const &index = container(it);

		do {
			TSize w0 = dirAt(value(it).node, index);
			if (w0 & index.LAST_CHILD)
				return false;

			value(it).node += (w0 & index.LEAF)? 1: 2;
			value(it).edgeLen = -1;

			value(it).range.i1 = value(it).range.i2;
			_adjustRightBorder(it);
		} while ((THideEmptyEdges::VALUE && emptyParentEdge(it)) || !nodeHullPredicate(it));
		return true;
	}

	// go up one edge (returns false if in root node)
	// can be used at most once, as no history stack is available
	template < typename TText, typename TWotdSpec, typename TSpec >
	inline bool 
	_goUp(Iter< Index<TText, Index_Wotd<TWotdSpec> >, VSTree< TopDown<TSpec> > > &it) 
	{
		if (!isRoot(it)) {
			value(it) = it._parentDesc;
			return true;
		}
		return false;
	}

	// go up one edge (returns false if in root node)
	template < typename TText, typename TSpec >
	inline bool 
	_goUp(Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it) 
	{
		typedef typename Size< Index<TText, Index_Wotd<WotdOriginal> > >::Type TSize;

		if (!empty(it.history)) {
			_HistoryStackWotdOriginal<TSize> const &entry = top(it.history);
			value(it).node = entry.node;
			value(it).parentRepLen -= entry.edgeLen;
			value(it).edgeLen = entry.edgeLen;
			pop(it.history);
			return true;
		}
		return false;
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline bool 
	_goUp(Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > &it) 
	{
		typedef typename Size< Index<TText, Index_Wotd<TIndexSpec> > >::Type TSize;

		if (!empty(it.history)) 
		{
			_HistoryStackWotdModified<TSize> const &entry = top(it.history);
			value(it).node = entry.node;
			value(it).parentRepLen -= entry.edgeLen;
			value(it).edgeLen = entry.edgeLen;
			value(it).range = entry.range;
			pop(it.history);
			if (!empty(it.history))
				value(it).parentRight = top(it.history).range.i2;	// copy right boundary of parent's range
			return true;
		}
		return false;
	}


//////////////////////////////////////////////////////////////////////////////


	// Counting sort - Step 2a: Count characters
	template < typename TBuckets, typename TText >
	inline void
	_wotdCountChars(TBuckets &buckets, TText const &text)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type TTextIterator;

		TTextIterator itText = begin(text, Standard());
		TTextIterator itTextEnd = end(text, Standard());
		for (; itText != itTextEnd; ++itText)
			++buckets[ordValue(*itText)];
	}


	template < typename TBuckets, typename TText, typename TSpec >
	inline void
	_wotdCountChars(TBuckets &buckets, StringSet<TText, TSpec> const &stringSet)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type TTextIterator;

		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TText const &text = value(stringSet, seqNo);
			TTextIterator itText = begin(text, Standard());
			TTextIterator itTextEnd = end(text, Standard());
			for (; itText != itTextEnd; ++itText)
				++buckets[ordValue(*itText)];
		}
	}

	// Counting sort - Step 2b: Count the prefixLen'th characters of suffices
	template < typename TBuckets, typename TText, typename TSA, typename TSize >
	inline typename Size<TText>::Type
	_wotdCountCharsWotdOriginal(
		TBuckets &buckets, 
		TText const &text, 
		TSA &sa,
		TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA const, Standard>::Type	TSAIterator;
		typedef typename Size<TText>::Type						TTextSize;

		TTextIterator itText = begin(text, Standard());
		TSAIterator itSA = begin(sa, Standard());
		TSAIterator itSAEnd = end(sa, Standard());

		TTextSize sentinels = 0;
		TTextSize textLength = length(text);
		for (; itSA != itSAEnd; ++itSA) 
		{
			// add prefix on entries in sa
			TTextSize saValue = (*itSA += prefixLen);
			if (textLength > saValue)
				++buckets[ordValue(*(itText + saValue))];
			else
				if (textLength == saValue) ++sentinels;
		}
		return sentinels;
	}

	template < typename TBuckets, typename TText, typename TSA, typename TSize >
	inline typename Size<TText>::Type
	_wotdCountChars(
		TBuckets &buckets, 
		TText const &text, 
		TSA const &sa,
		TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA const, Standard>::Type	TSAIterator;
		typedef typename Size<TText>::Type						TTextSize;

		TTextIterator itText = begin(text, Standard()) + prefixLen;
		TSAIterator itSA = begin(sa, Standard());
		TSAIterator itSAEnd = end(sa, Standard());

		TTextSize sentinels = 0;
		TTextSize textLength = length(text) - prefixLen;
		for (; itSA != itSAEnd; ++itSA) 
		{
			TTextSize saValue = *itSA;
			if (textLength > saValue)
				++buckets[ordValue(*(itText + saValue))];
			else
				if (textLength == saValue) ++sentinels;
		}
		return sentinels;
	}

	template < 
		typename TBuckets, 
		typename TText, 
		typename TSpec, 
		typename TSA, 
		typename TSize 
	>
	inline typename Size<TText>::Type
	_wotdCountChars(
		TBuckets &buckets, 
		StringSet<TText, TSpec> const &stringSet,
		TSA const &sa,
		TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA const, Standard>::Type	TSAIterator;
		typedef typename Size<TText>::Type						TTextSize;

		TTextIterator itText = TTextIterator();
		TSAIterator itSA = begin(sa, Standard());
		TSAIterator itSAEnd = end(sa, Standard());

		TTextSize sentinels = 0;
		TTextSize textLength = 0;
		unsigned lastSeqSeen = -1;
		Pair<unsigned, TTextSize> lPos;
		for (; itSA != itSAEnd; ++itSA) 
		{
			posLocalize(lPos, *itSA, stringSetLimits(stringSet));
			if (lastSeqSeen != getSeqNo(lPos))
			{
				lastSeqSeen = getSeqNo(lPos);

				// shift textBegin and textLength by prefixLen
				textLength = length(stringSet[lastSeqSeen]) - prefixLen;
				itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
			}
			if (textLength > getSeqOffset(lPos))
				++buckets[ordValue(*(itText + getSeqOffset(lPos)))];
			else
				if (textLength == getSeqOffset(lPos)) ++sentinels;
		}
		return sentinels;
	}


//////////////////////////////////////////////////////////////////////////////


	// Counting sort - Step 3: Cumulative sum
	template < typename TBounds, typename TBuckets, typename TSize >
	inline typename Size<TBuckets>::Type
	_wotdCummulativeSum(TBounds &bounds, TBuckets const &buckets, TSize offset)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TBounds, Standard>::Type		TBoundIterator;
		typedef typename Iterator<TBuckets, Standard>::Type		TBucketsIterator;

		TBucketsIterator it = begin(buckets, Standard());
		TBucketsIterator itEnd = end(buckets, Standard());
		TBoundIterator bit = begin(bounds, Standard());

		typename Value<TBounds>::Type	sum = offset;
		typename Size<TBuckets>::Type	requiredSize = 0;
		typename Value<TBuckets>::Type	diff;

		for (; it != itEnd; ++it, ++bit)
			if ((diff = *it) != 0) {
				requiredSize += (diff > 1)? 2: 1;
				*bit = sum;
				sum += diff;
			}

		return requiredSize;
	}

	template < typename TBounds, typename TBuckets >
	inline typename Size<TBuckets>::Type
	_wotdCummulativeSum(TBounds &bounds, TBuckets const &buckets)
	{
		return _wotdCummulativeSum(bounds, buckets, 0);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.createWotdIndex:
..summary:Builds a q-gram index on a sequence. 
..cat:Index
..signature:createWotdIndex(sa, dir, text)
..param.text:The sequence.
..param.sa:The resulting list in which all q-grams are sorted alphabetically.
..param.dir:The resulting array that indicates at which position in index the corresponding q-grams can be found.
..returns:Index contains the sorted list of qgrams. For each possible q-gram pos contains the first position in index that corresponds to this q-gram. 
*/

	// single sequence
	template < typename TIndex >
	typename Size<TIndex>::Type
	_sortFirstWotdBucket(TIndex &index)
	{
	SEQAN_CHECKPOINT
		typedef typename Fibre<TIndex, Wotd_Text >::Type		TText;
		typedef typename Fibre<TIndex, Wotd_SA >::Type			TSA;
		typedef typename TIndex::TCounter						TCounter;

		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA, Standard>::Type			TSAIterator;
		typedef typename Iterator<TCounter, Standard>::Type		TCntIterator;
		typedef typename Size<TText>::Type						TSize;
		
		TText const &text = indexText(index);
		TCounter &occ = index.tempOcc;
		TCounter &bound = index.tempBound;

		// 1. clear counters and copy SA to temporary SA
		arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);

		// 2. count characters
		_wotdCountChars(occ, text);

		// 3. cumulative sum
		TSize requiredSize = _wotdCummulativeSum(bound, occ);

		// 4. fill suffix array
		{
			TSA &sa = indexSA(index);
			TSAIterator saBeg = begin(sa, Standard());
			TCntIterator boundBeg = begin(bound, Standard());

			TTextIterator itText = begin(text, Standard());
			TTextIterator itTextEnd = end(text, Standard());
			for(TSize i = 0; itText != itTextEnd; ++itText, ++i)
				*(saBeg + (*(boundBeg + ordValue(*itText)))++) = i;
		}
		index.sentinelOcc = 0;
		index.sentinelBound = 0;

		return requiredSize;
	}

	// multiple sequences
	template < typename TText, typename TSpec, typename TIndexSpec >
	typename Size< Index<StringSet<TText, TSpec>, TIndexSpec> >::Type
	_sortFirstWotdBucket(Index<StringSet<TText, TSpec>, TIndexSpec> &index)
	{
	SEQAN_CHECKPOINT
		typedef Index<StringSet<TText, TSpec>, TIndexSpec>		TIndex;
		typedef typename Fibre<TIndex, Wotd_SA >::Type			TSA;
		typedef typename TIndex::TCounter						TCounter;

		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA, Standard>::Type			TSAIterator;
		typedef typename Iterator<TCounter, Standard>::Type		TCntIterator;
		typedef typename Size<TText>::Type						TSize;
		typedef typename Value<TText>::Type						TValue;
		
		StringSet<TText, TSpec> const &stringSet = indexText(index);
		TCounter &occ = index.tempOcc;
		TCounter &bound = index.tempBound;

		// 1. clear counters and copy SA to temporary SA
		arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);

		// 2. count characters
		_wotdCountChars(occ, stringSet);

		// 3. cummulative sum
		TSize requiredSize = _wotdCummulativeSum(bound, occ);

		// 4. fill suffix array
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TSA &sa = indexSA(index);
			TSAIterator saBeg = begin(sa, Standard());
			TCntIterator boundBeg = begin(bound, Standard());

			typename Value<TSA>::Type localPos;
			assignValueI1(localPos, seqNo);
			assignValueI2(localPos, 0);

			TText const &text = value(stringSet, seqNo);
			TTextIterator itText = begin(text, Standard());
			TTextIterator itTextEnd = end(text, Standard());
			for(; itText != itTextEnd; ++itText) 
			{
				*(saBeg + (*(boundBeg + ordValue(*itText)))++) = localPos;
				assignValueI2(localPos, getValueI2(localPos) + 1);
			}
		}
		index.sentinelOcc = 0;
		index.sentinelBound = 0;

		return requiredSize;
	}



	// sort bucket using counting sort
	// (nearly) ORIGINAL VERSION:
	// - brings the bucket with the longest suffix (lpBucket) to front
	// - all other buckets are in lexicographical order
	// - adds prefixLen to all SA entries in SA[left,right)
	template < typename TText, typename TSize >
	TSize _sortWotdBucket(
		Index<TText, Index_Wotd<WotdOriginal> > &index,
		TSize left, 
		TSize right,
		TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef Index<TText, Index_Wotd<WotdOriginal> >			TIndex;
		typedef typename Fibre<TIndex, Wotd_SA >::Type			TSA;
		typedef typename TIndex::TCounter						TCounter;

		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA, Standard>::Type			TSAIterator;
		typedef typename Iterator<TCounter, Standard>::Type		TCntIterator;
		typedef typename Size<TText>::Type						TTextSize;

		TText const &text = indexText(index);
		TCounter const &tempSA = index.tempSA;
		TCounter &occ = index.tempOcc;
		TCounter &bound = index.tempBound;

		// 1. clear counters and copy SA to temporary SA
		arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);

		// 2. count characters
		index.tempSA = infix(indexSA(index), left, right);
		index.sentinelBound = 0;
		index.sentinelOcc =
			_wotdCountCharsWotdOriginal(occ, text, tempSA, prefixLen);

		// 3. cumulative sum
		TSize requiredSize = 0;
		if (index.sentinelOcc != 0)
			requiredSize = (index.sentinelOcc > 1)? 2: 1;

		requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
		index.sentinelBound = left;

		// 4. fill suffix array
		{
			TSA &sa = indexSA(index);
			TSAIterator saBeg = begin(sa, Standard());
			TCntIterator boundBeg = begin(bound, Standard());

			TTextIterator itText = begin(text, Standard());
			TCntIterator itSA = begin(tempSA, Standard());
			TCntIterator itSAEnd = end(tempSA, Standard());
			TTextSize textLength = length(text);
			for(; itSA != itSAEnd; ++itSA)
			{
				TTextSize saValue = *itSA;
				if (textLength > saValue)
					*(saBeg + (*(boundBeg + ordValue(*(itText + saValue))))++) = saValue;
				else
					if (textLength == saValue)
						*(saBeg + index.sentinelBound++) = saValue;
			}
		}

		// 5. move lpBucket to front and add saOffset to hash directory entries
		{
			TSize lpBucket = ordValue(text[tempSA[0]]);
			if (lpBucket != 0) {
				TSize lpBucketOcc = occ[lpBucket];
				TSize lpBucketBound = bound[lpBucket];

				TCntIterator itOcc = begin(occ, Standard()) + lpBucket;
				TCntIterator itBound = begin(bound, Standard()) + lpBucket;
				TCntIterator itBeg = begin(bound, Standard());
				for(; itBound != itBeg; --itBound, --itOcc) {
					*itOcc = *(itOcc - 1);
					*itBound = *(itBound - 1);
				}
				if (index.sentinelOcc != 0) {
					// bring first bucket before sentinel bucket
					*itOcc = index.sentinelOcc;
					*itBound = index.sentinelBound;
					index.sentinelOcc = lpBucketOcc;
					index.sentinelBound = lpBucketBound;
				} else {
					*itOcc = lpBucketOcc;
					*itBound = lpBucketBound;
				}
			} else
				if (index.sentinelOcc != 0) {
					// bring first bucket before sentinel bucket
					TSize swap = index.sentinelOcc;
					index.sentinelOcc = occ[0];
					occ[0] = swap;
					swap = index.sentinelBound;
					index.sentinelBound = bound[0];
					bound[0] = swap;
				}
		}

		return requiredSize;
	}





	// sort bucket using counting sort
	// MODIFIED VERSION:
	// - all buckets are in lexicographical order
	// - SA[left,right) contains real SA entries (the beginning positions of the suffices)

	// single sequence
	template < typename TIndex, typename TSize >
	TSize _sortWotdBucket(
		TIndex &index,
		TSize left, 
		TSize right,
		TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Fibre<TIndex, Wotd_Text >::Type		TText;
		typedef typename Fibre<TIndex, Wotd_SA >::Type			TSA;
		typedef typename TIndex::TCounter						TCounter;

		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA, Standard>::Type			TSAIterator;
		typedef typename Iterator<TCounter, Standard>::Type		TCntIterator;
		typedef typename Size<TText>::Type						TTextSize;

		TText const &text = indexText(index);
		TCounter const &tempSA = index.tempSA;
		TCounter &occ = index.tempOcc;
		TCounter &bound = index.tempBound;

		// 1. clear counters and copy SA to temporary SA
		arrayFill(begin(occ, Standard()), end(occ, Standard()), 0);
		index.tempSA = infix(indexSA(index), left, right);

		// 2. count characters
		index.sentinelBound = 0;
		index.sentinelOcc = _wotdCountChars(occ, text, tempSA, prefixLen);

		// 3. cumulative sum
		TSize requiredSize = 0;
		if (index.sentinelOcc != 0)
			requiredSize = (index.sentinelOcc > 1)? 2: 1;

		requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
		index.sentinelBound = left;

		// 4. fill suffix array
		{
			TSA &sa = indexSA(index);
			TSAIterator saBeg = begin(sa, Standard());
			TCntIterator boundBeg = begin(bound, Standard());

			TTextIterator itText = begin(text, Standard()) + prefixLen;
			TCntIterator itSA = begin(tempSA, Standard());
			TCntIterator itSAEnd = end(tempSA, Standard());
			TTextSize textLength = length(text) - prefixLen;
			for(; itSA != itSAEnd; ++itSA)
			{
				TTextSize saValue = *itSA;
				if (textLength > saValue)
					*(saBeg + (*(boundBeg + ordValue(*(itText + saValue))))++) = saValue;
				else
					if (textLength == saValue)
						*(saBeg + index.sentinelBound++) = saValue;
			}
		}

		return requiredSize;
	}

	// multiple sequences
	template < typename TText, typename TSpec, typename TIndexSpec, typename TSize >
	TSize _sortWotdBucket(
		Index<StringSet<TText, TSpec>, TIndexSpec> &index,
		TSize left, 
		TSize right,
		TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef Index<StringSet<TText, TSpec>, TIndexSpec>			TIndex;
		typedef typename Fibre<TIndex, Wotd_SA >::Type				TSA;
		typedef typename TIndex::TCounter							TCounter;
		typedef typename TIndex::TTempSA							TTempSA;

		typedef typename Iterator<TText const, Standard>::Type		TTextIterator;
		typedef typename Iterator<TSA, Standard>::Type				TSAIterator;
		typedef typename Iterator<TTempSA, Standard>::Type			TTempSAIterator;
		typedef typename Iterator<TCounter, Standard>::Type			TCntIterator;
		typedef typename Size<TText>::Type							TTextSize;
		typedef typename Value<TText>::Type							TValue;

		StringSet<TText, TSpec> const &stringSet = indexText(index);
		TTempSA const &tempSA = index.tempSA;
		TCounter &occ = index.tempOcc;
		TCounter &bound = index.tempBound;

		// 1. clear counters and copy SA to temporary SA
		TCntIterator occBeg = begin(occ, Standard());

		arrayFill(occBeg, end(occ, Standard()), 0);
		index.tempSA = infix(indexSA(index), left, right);

		// 2. count characters
		index.sentinelBound = 0;
		index.sentinelOcc = _wotdCountChars(occ, stringSet, tempSA, prefixLen);

		// 3. cumulative sum
		TSize requiredSize = 0;
		if (index.sentinelOcc != 0)
			requiredSize = (index.sentinelOcc > 1)? 2: 1;

		requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
		index.sentinelBound = left;

		// 4. fill suffix array
		{
			TSA &sa = indexSA(index);
			TSAIterator saBeg = begin(sa, Standard());
			TCntIterator boundBeg = begin(bound, Standard());

			TTextIterator itText = TTextIterator();
			TTempSAIterator itSA = begin(tempSA, Standard());
			TTempSAIterator itSAEnd = end(tempSA, Standard());
			TTextSize textLength = 0;
			unsigned lastSeqSeen = -1;
			Pair<unsigned, TTextSize> lPos;
			for(; itSA != itSAEnd; ++itSA)
			{
				posLocalize(lPos, *itSA, stringSetLimits(index));
				if (lastSeqSeen != getSeqNo(lPos))
				{
					lastSeqSeen = getSeqNo(lPos);

					// shift textBegin and textLength by prefixLen
					textLength = length(stringSet[lastSeqSeen]) - prefixLen;
					itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
				}
				if (textLength > getSeqOffset(lPos))
					*(saBeg + (*(boundBeg + ordValue(*(itText + getSeqOffset(lPos)))))++) = *itSA;
				else
					if (textLength == getSeqOffset(lPos))
						*(saBeg + index.sentinelBound++) = *itSA;
			}
		}

		return requiredSize;
	}





	template < typename TSA, typename TText >
	typename Size<TText>::Type 
	_bucketLCP(TSA const &sa, TText const &text)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA const, Standard>::Type	TSAIterator;
		typedef typename Value<TText>::Type						TValue;
		typedef typename Size<TText>::Type						TTextSize;

		TTextSize prefixLen = 0;
		
		if (length(sa) < 2) return prefixLen;

		TTextIterator itText = begin(text, Standard());
		TSAIterator itSAEnd = end(sa, Standard());
		TTextSize textLength = length(text);

		do {
			TSAIterator itSA = begin(sa, Standard());
			TTextSize sa = *itSA;
			if (textLength == sa) return prefixLen;

			TValue c = *(itText + sa);
			for(++itSA; itSA != itSAEnd; ++itSA) {
				sa = *itSA;
				if (textLength == sa || c != *(itText + sa))
					return prefixLen;
			}
			++prefixLen; --textLength;
			++itText;
		} while (true);
	}

	template < typename TSA, typename TText, typename TSize >
	typename Size<TText>::Type 
	_bucketLCP(TSA const &sa, TText const &text, TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA const, Standard>::Type	TSAIterator;
		typedef typename Value<TText>::Type						TValue;
		typedef typename Size<TText>::Type						TTextSize;

		if (length(sa) < 2) return prefixLen;

		TTextIterator itText = begin(text, Standard()) + prefixLen;
		TSAIterator itSAEnd = end(sa, Standard());
		TTextSize textLength = length(text) - prefixLen;

		do {
			TSAIterator itSA = begin(sa, Standard());
			TTextSize sa = *itSA;
			if (textLength == sa) return prefixLen;

			TValue c = *(itText + sa);
			for(++itSA; itSA != itSAEnd; ++itSA) {
				sa = *itSA;
				if (textLength == sa || c != *(itText + sa))
					return prefixLen;
			}
			++prefixLen; --textLength;
			++itText;
		} while (true);
	}

	template < typename TSA, typename TText, typename TSpec, typename TSize >
	typename Size<TText>::Type 
	_bucketLCP(TSA const &sa, StringSet<TText, TSpec> const &stringSet, TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TTextIterator;
		typedef typename Iterator<TSA const, Standard>::Type	TSAIterator;
		typedef typename Value<TText>::Type						TValue;
		typedef typename Size<TText>::Type						TTextSize;

		if (length(sa) < 2) return prefixLen;

		TSAIterator itSAEnd = end(sa, Standard());
		TTextIterator itText;
		TTextSize textLength;

		Pair<unsigned, TTextSize> lPos;
		do {
			TSAIterator itSA = begin(sa, Standard());
			posLocalize(lPos, *itSA, stringSetLimits(stringSet));

			unsigned lastSeqSeen = getSeqNo(*itSA);
			textLength = length(stringSet[lastSeqSeen]) - prefixLen;
			if (textLength == getSeqOffset(lPos)) return prefixLen;

			itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
			TValue c = *(itText + getSeqOffset(*itSA));
			for(++itSA; itSA != itSAEnd; ++itSA)
			{
				posLocalize(lPos, *itSA, stringSetLimits(stringSet));

				if (lastSeqSeen != getSeqNo(lPos))
				{
					lastSeqSeen = getSeqNo(lPos);

					// shift textBegin and textLength by prefixLen
					textLength = length(stringSet[lastSeqSeen]) - prefixLen;
					itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
				}

				if (textLength == getSeqOffset(lPos) || c != *(itText + getSeqOffset(lPos)))
					return prefixLen;
			}
			++prefixLen; --textLength;
			++itText;
		} while (true);
	}


	template <typename TText, typename TSpec, typename TPos>
	inline TPos
	_getNodeLP(
		Index<TText, Index_Wotd<TSpec> > const &index, 
		TPos pos)
	{
		TPos w0 = dirAt(pos, index);
		if (w0 & index.LEAF)
			return w0 & index.BITMASK0;

		TPos w1 = dirAt(pos + 1, index);
		if (w1 & index.UNEVALUATED)
			return saAt(w0 & index.BITMASK1, index);
		else
			return w0 & index.BITMASK0;
	}

	// store buckets into directory 
	// ORIGINAL VERSION: storing SA entries and topology links in Dir
	template <typename TText, typename TPos>
	inline void 
	_storeWotdChildren(
		Index<TText, Index_Wotd<WotdOriginal> > &index,
		TPos dirOfs)
	{
	SEQAN_CHECKPOINT
		typedef Index<TText, Index_Wotd<WotdOriginal> >		TIndex;
		typedef typename Fibre<TIndex, Wotd_Dir>::Type		TDir;
		typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
		typedef typename TIndex::TCounter					TCounter;
		typedef typename Iterator<TCounter, Standard>::Type	TCntIterator;

		typedef typename Value<TCounter>::Type				TValue;

		TDirIterator itDir = begin(indexDir(index), Standard()) + dirOfs;
		TDirIterator itDirEnd = end(indexDir(index), Standard());
		TDirIterator itPrev = itDirEnd;

		TCntIterator it = begin(index.tempOcc, Standard());
		TCntIterator bit = begin(index.tempBound, Standard());
		TCntIterator itEnd = end(index.tempOcc, Standard());

		TValue occ;
		if (index.sentinelOcc != 0)
			if (index.sentinelOcc > 1) { // occurs on multiseqs
				itPrev = itDir;
				*itDir = index.sentinelBound - index.sentinelOcc;	++itDir;
				*itDir = index.sentinelBound | index.UNEVALUATED;	++itDir;
			} else {
				itPrev = itDir;
				*itDir = saAt(index.sentinelBound - index.sentinelOcc, index) | index.LEAF;
				++itDir;
			}

		for (; it != itEnd; ++it, ++bit)
		{
			if ((occ = *it) == 0) continue;
			if (occ > 1) {
				itPrev = itDir;
				*itDir = *bit - occ; 							++itDir;
				*itDir = *bit | index.UNEVALUATED;				++itDir;
			} else {
				itPrev = itDir;
				*itDir = saAt(*bit - occ, index) | index.LEAF;	++itDir;
			}
		}

		// mark the last child
		if (itPrev != itDirEnd)
			*itPrev |= index.LAST_CHILD;
	}

	// store buckets into directory
	// MODIFIED VERSION: storing SA links and topology links in Dir
	template <typename TText, typename TSpec, typename TSize>
	inline void 
	_storeWotdChildren(
		Index<TText, Index_Wotd<TSpec> > &index,
		TSize dirOfs,
		TSize lcp)
	{
	SEQAN_CHECKPOINT
		typedef Index<TText, Index_Wotd<TSpec> >			TIndex;
		typedef typename Fibre<TIndex, Wotd_Dir>::Type		TDir;
		typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
		typedef typename TIndex::TCounter					TCounter;
		typedef typename Iterator<TCounter, Standard>::Type	TCntIterator;

		typedef typename Value<TCounter>::Type				TValue;

		TDirIterator itDirBegin = begin(indexDir(index), Standard()) + dirOfs;
		TDirIterator itDirEnd = end(indexDir(index), Standard());
		TDirIterator itDir = itDirBegin;
		TDirIterator itPrev = itDirEnd;

		TCntIterator it = begin(index.tempOcc, Standard());
		TCntIterator bit = begin(index.tempBound, Standard());
		TCntIterator itEnd = end(index.tempOcc, Standard());

		TValue occ;
		if (index.sentinelOcc != 0)
			if (index.sentinelOcc > 1) { // occurs on multiseqs
				itPrev = itDir;
				*itDir = index.sentinelBound - index.sentinelOcc;	++itDir;
				*itDir = index.sentinelBound | index.UNEVALUATED;	++itDir;
			} else {
				itPrev = itDir;
				*itDir = (index.sentinelBound - index.sentinelOcc) | index.LEAF;
				++itDir;
			}

		for (; it != itEnd; ++it, ++bit)
		{
			if ((occ = *it) == 0) continue;
			if (occ > 1) {
				itPrev = itDir;
				*itDir = *bit - occ; 				++itDir;
				*itDir = *bit | index.UNEVALUATED;	++itDir;
			} else {
				itPrev = itDir;
				*itDir = (*bit - occ) | index.LEAF;	++itDir;
			}
		}

		// first child gets the mutual lcp value of the children (== parent repLength)
		*itDirBegin = ((*itDirBegin) & ~index.BITMASK0) | lcp;

		// mark the last child
		if (itPrev != itDirEnd)
			*itPrev |= index.LAST_CHILD;
	}


	template < typename TText, typename TSpec >
	inline typename Size< Index<TText, Index_Wotd<WotdOriginal> > >::Type 
	_wotdEvaluate(Iter< Index<TText, Index_Wotd<WotdOriginal> >, VSTree<TSpec> > const &it)
	{
		typedef Index<TText, Index_Wotd<WotdOriginal> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		TIndex &index = const_cast<TIndex&>(container(it));
		TSize pos = value(it).node;
		TSize w1 = dirAt(pos + 1, index);

		// test for evaluation
		if (w1 & index.UNEVALUATED)
		{
			TSize w0 = dirAt(pos, index);
			TSize lp = saAt(w0 & index.BITMASK0, index);
			TSize dst = length(indexDir(index));

			TSize size = _sortWotdBucket(
				index, 
				w0 & index.BITMASK0, 
				w1 & index.BITMASK1, 
				parentEdgeLength(it));

			resize(indexDir(index), dst + size, Generous());
			_storeWotdChildren(index, dst);

			// mark nodes with solely empty child edges
			w1 = dst;
			if (index.sentinelOcc > 1)
				if (size == 2) w1 |= index.SENTINELS;

			assert(!(index.sentinelOcc == 1 && size == 1));

			dirAt(pos, index)     = (w0 & ~index.BITMASK0) | lp;
			dirAt(pos + 1, index) = w1;
		}

		return w1;
	}

	template < typename TText, typename TIndexSpec, typename TSpec >
	inline typename Size< Index<TText, Index_Wotd<TIndexSpec> > >::Type 
	_wotdEvaluate(Iter< Index<TText, Index_Wotd<TIndexSpec> >, VSTree<TSpec> > const &it)
	{
		typedef Index<TText, Index_Wotd<TIndexSpec> >	TIndex;
		typedef typename Size<TIndex>::Type				TSize;

		TIndex &index = const_cast<TIndex&>(container(it));
		TSize pos = value(it).node;
		TSize w1 = dirAt(pos + 1, index);

		// test for evaluation
		if (w1 & index.UNEVALUATED) 
		{
			TSize dst = length(indexDir(index));

			TSize size = _sortWotdBucket(
				index, 
				value(it).range.i1, 
				w1 & index.BITMASK1, 
				repLength(it));

/*			if (globalDumpFlag) {
				::std::cerr << '"' << representative(it) << '"' << ::std::endl;
				_dumpFreq(index);
			}
*/
			resize(indexDir(index), dst + size, Generous());
			_storeWotdChildren(index, dst, repLength(it));

			// mark nodes with solely empty child edges
			w1 = dst;
			if (index.sentinelOcc > 1)
				if (size == 2) w1 |= index.SENTINELS;

			dirAt(pos + 1, index) = w1;
		}

		return w1;
	}


	template <typename TText, typename TSpec>
	inline void
	_dump(Index<TText, Index_Wotd<TSpec> > &index)
	{
		typedef Index<TText, Index_Wotd<TSpec> >			TIndex;
		typedef typename Fibre<TIndex, Wotd_Dir>::Type		TDir;
		typedef typename Value<TDir>::Type					TDirValue;

		::std::cout << "  Dir (wotd)" << ::std::endl;
		for(unsigned i=0; i < length(indexDir(index)); ++i) {
			TDirValue d = indexDir(index)[i];
			::std::cout << i << ":  " << (d & index.BITMASK0);
			if (d & index.LEAF)			::std::cout << "  (Leaf/Uneval)";
			if (d & index.LAST_CHILD)	::std::cout << "  (LastChild/SENTINELS)";
			::std::cout << ::std::endl;
		}

		::std::cout << ::std::endl << "  SA" << ::std::endl;
		for(unsigned i=0; i < length(indexSA(index)); ++i)
			::std::cout << i << ":  " << indexSA(index)[i] << "  " << suffix(indexText(index), indexSA(index)[i]) << ::std::endl;

		::std::cout << ::std::endl;

	}

//////////////////////////////////////////////////////////////////////////////
// _goDownChar

	template < typename TText, class TSpec, typename TValue >
	inline bool _goDownChar(
		Iter<Index<TText, Index_Wotd<WotdOriginal> >, VSTree< TopDown<TSpec> > > &it,
		TValue c)
	{
		typedef Index<TText, Index_Wotd<TSpec> >	TIndex;
		typedef typename Value<TIndex>::Type		TIndexValue;

		bool sorted = false;
		if (!goDown(it)) return false;
		do {
			if (parentEdgeLength(it) != 0) {
				TIndexValue edgeChar = parentEdgeLabel(it)[0];
				if (edgeChar == c) return true;		// the edge is found
				if (sorted && edgeChar > c) break;	// too far (except the first one,
			}										// child edges are sorted)
			sorted = true;
		} while (goRight(it));
		_goUp(it);
		return false;
	}

	template < typename TText, class TIndexSpec, class TSpec, typename TValue >
	inline bool _goDownChar(
		Iter<Index<TText, Index_Wotd<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
		TValue c)
	{
		typedef Index<TText, Index_Wotd<TSpec> >	TIndex;
		typedef typename Value<TIndex>::Type		TIndexValue;

		if (!goDown(it)) return false;
		do {
			if (parentEdgeLength(it) != 0) {
				TIndexValue edgeChar = parentEdgeLabel(it)[0];
				if (edgeChar == c) return true;	// the edge is found
				if (edgeChar > c) break;		// too far (child edges are sorted)
			}
		} while (goRight(it));
		_goUp(it);
		return false;
	}

/*
	template < typename TText, typename TSpec, typename TValue >
	inline bool
	_getNodeByChar(
		Iter< Index<TText, Index_Wotd<TSpec> >, VSTree<TSpec> > const &it,
		TValue c, 
		typename VertexDescriptor< Index<TText, Index_Wotd<TSpec> > >::Type &childDesc)
	{
		typedef Index<TText, Index_Wotd<TSpec> >			TIndex;
		typedef typename Fibre<TIndex, Wotd_Dir>::Type		TDir;
		typedef typename Fibre<TIndex, Wotd_SA>::Type		TSA;

		typedef typename Value<TSA>::Type					TSAValue;
		typedef typename Value<TDir>::Type					TDirValue;

		typename Size<TIndex>::Type len = length(index);
		typename VertexDescriptor<TIndex>::Type	desc;

		TSAValue pos = _firstSuffixOfBucket(index, value(it).node);
		while (pos == len || value < (c = textAt(index, pos + value.i2))) {
			value.node += (dirAt(value.node, index) & index.LEAF)? 1: 2;
			pos = _firstSuffixOfBucket(index, value.node);
		}
		assert(pos != len);

		return c == value;
	}
*/

//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation 

	template <typename TText, typename TSpec>
	inline void _wotdCreateFirstLevel(Index<TText, Index_Wotd<TSpec> > &index)
	{
		typedef Index<TText, Index_Wotd<TSpec> >	TIndex;
		typedef typename Value<TIndex>::Type		TValue;
		typedef typename Size<TIndex>::Type			TSize;

		resize(indexSA(index), length(indexRawText(index)));
		resize(index.tempOcc, ValueSize<TValue>::VALUE + 1);
		resize(index.tempBound, ValueSize<TValue>::VALUE + 1);

		TSize size = _sortFirstWotdBucket(index);
		if (size > 0) 
		{
			resize(indexDir(index), size + 2, Generous());
			_storeWotdChildren(index, 2, 0);

			// mark nodes with solely empty child edges
			TSize w1 = 2;
			if (index.sentinelOcc > 1)
				if (size == 2) w1 |= index.SENTINELS;

			dirAt(0, index) = 0 | index.LAST_CHILD;
			dirAt(1, index) = w1;

		} else {
			resize(indexDir(index), 1);
			dirAt(0, index) = 0 | index.LAST_CHILD | index.LEAF;
		}
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, Index_Wotd<TSpec> > &index, Wotd_SA const, Default const)
	{
		_wotdCreateFirstLevel(index);
		return true;
	}

}

#endif //#ifndef SEQAN_HEADER_...
