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
  $Id: index_esa_base.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_ESA_BASE_H
#define SEQAN_HEADER_INDEX_ESA_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	// dfs order
	struct _Preorder;
	struct _Postorder;

	template <typename TDFSOrder = _Postorder, typename THideEmptyEdges = True>
	struct VSTreeIteratorTraits {
		typedef TDFSOrder DFSOrder;
		typedef THideEmptyEdges HideEmptyEdges;
	};

/**
.Tag.Preorder:
..summary:Preorder depth-first search.
..cat:Index
..signature:Preorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a preorder fashion (visit the node before its children).
..see:Tag.Postorder
*/

/**
.Tag.Postorder:
..summary:Postorder depth-first search.
..cat:Index
..signature:Postorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a postorder fashion (visit the node after its children).
..see:Tag.Preorder
*/

/**
.Tag.PreorderEmptyEdges:
..summary:Preorder depth-first search in a suffix tree with leaves for every suffix.
..cat:Index
..signature:PreorderEmptyEdges
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a preorder fashion (visit the node before its children).
Empty edges are traversed also, i.e. for every suffix there is a leaf node representing it.
..see:Tag.PostorderEmptyEdges
..see:Tag.Preorder
*/

/**
.Tag.PostorderEmptyEdges:
..summary:Postorder depth-first search in a suffix tree with leaves for every suffix.
..cat:Index
..signature:PostorderEmptyEdges
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a postorder fashion (visit the node after its children).
Empty edges are traversed also, i.e. for every suffix there is a leaf node representing it.
..see:Tag.Postorder
*/

/**
.Tag.EmptyEdges:
..summary:Consider a suffix tree with leaves for every suffix.
..cat:Index
..signature:EmptyEdges
..remarks:When given as a second parameter in @Function.goDown@, empty edges are traversed also, i.e. for every suffix there is a leaf node representing it.
..see:Tag.HideEmptyEdges
*/

/**
.Tag.HideEmptyEdges:
..summary:Consider a suffix tree with no empty edges (default behaviour).
..cat:Index
..signature:HideEmptyEdges
..remarks:When given as a second parameter in @Function.goDown@, only non-empty edges are traversed.
..see:Tag.EmptyEdges
*/

	// predefined iterator traits
	struct Preorder:			VSTreeIteratorTraits<_Preorder,  True> {};
	struct Postorder:			VSTreeIteratorTraits<_Postorder, True> {};
	struct PreorderEmptyEdges:	VSTreeIteratorTraits<_Preorder,  False> {};	// also iterate over
	struct PostorderEmptyEdges:	VSTreeIteratorTraits<_Postorder, False> {};	// empty edges (with $-label)

	// traits for TopDown iterators (w/o ParentLinks) for which postorder/preorder is ignored
	struct HideEmptyEdges:		VSTreeIteratorTraits<_Postorder, True> {};
	struct EmptyEdges:			VSTreeIteratorTraits<_Postorder, False> {};	// empty edges (with $-label)
	
	// MultiMEMs are more specialized MaxRepeats
	template <typename TSpec = void>
	struct _MaxRepeats;	// base class
	struct _MultiMEMs;	// subclass tag



	// virtual suffix tree iterators
	template <typename TSpec = void>
	struct VSTree;

		// top down traversal iterators
		template <typename TSpec = Preorder>
		struct TopDown;						// starts in the suffix tree root and can go down and go right

			// allows an top-down iterator to go up
			template < typename TSpec = Preorder >
			struct ParentLinks {};			// .. can also go up

		// bottom up traversal iterators
		template <typename TSpec = Postorder>
		struct BottomUp;					// starts in the first node of a depth-first-search and can go next

			struct	SuperMaxRepeats;					// maximal repeat and not part of a longer repeat
			struct	SuperMaxRepeatsFast;
			struct	MUMs;								// Maximal Unique Match (unique in every sequence)

			typedef _MaxRepeats<void>		MaxRepeats;	// maximal repeat
			struct	MaxRepeatOccurrences;
			typedef _MaxRepeats<_MultiMEMs> MultiMEMs;	// Multiple Maximal Exact Match
			struct	MultiMEMOccurences;					// i.e. maximal match over different sequences


/**
.Metafunction.GetVSTreeIteratorTraits:
..cat:Index
..summary:Default behaviour of @Function.goNext@ when no second parameter is given.
..signature:GetVSTreeIteratorTraits<TIterator>::Type
..param.TIterator:A @Spec.VSTree Iterator@.
..returns:$Tag.Postorder$ by default and $Tag.Preorder$ if $TIterator$ is $VSTree<TopDown<ParentLinks<> > >$ or $VSTree<TopDown<ParentLinks<Preorder> > >$.
*/

	template <typename TIterator>
	struct GetVSTreeIteratorTraits:
		DeepestSpec<TIterator> {};

//////////////////////////////////////////////////////////////////////////////

	template <typename TSize>
	struct VertexESA {
		Pair<TSize> range;			// current SA interval of hits (unique node identifier)
		TSize		parentRight;	// right boundary of parent node's range (allows to go right)

		VertexESA() {}

		VertexESA(MinimalCtor):
			range(0,0),
			parentRight(0) {}

		VertexESA(TSize otherRangeLeft, TSize otherRangeRight, TSize otherParentRight):
			range(Pair<TSize>(otherRangeLeft, otherRangeRight)),
			parentRight(otherParentRight) {}

		VertexESA(Pair<TSize> const &otherRange, TSize otherParentRight):
			range(otherRange),
			parentRight(otherParentRight) {}

		VertexESA(VertexESA const &other):
			range(other.range),
			parentRight(other.parentRight) {}
	};

	template < typename TText, typename TSpec >
	struct VertexDescriptor< Index<TText, Index_ESA<TSpec> > > {
		typedef typename Size< Index<TText, Index_ESA<TSpec> > >::Type TSize;
		typedef VertexESA<TSize> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	struct ArrayGaps;

	template <typename TSource, typename TSpec>
	class Align;


//////////////////////////////////////////////////////////////////////////////
// ESA fibres

/**
.Tag.ESA Index Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of an @Spec.Index_ESA.ESA@ index.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of an Enhanced Suffix Array based @Spec.Index_ESA.Index@.
..cat:Index

..tag.ESA_Text:The original text the index should be based on.

..tag.ESA_RawText:The raw text the index is really based on.
...remarks:$ESA_Text$ and $ESA_RawText$ fibres are equal by default.
They differ if the index text is a set of strings. Then, raw text is the concatenation of all strings in this set.

..tag.ESA_SA:The suffix array.
...remarks:The suffix array contains the indices of all suffices of $ESA_RawText$ in lexicographical order.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of the @Metafunction.SAValue@ of $TIndex$.

..tag.ESA_LCP:The lcp table.
...remarks:The lcp table contains the lcp-value of two adjacent suffices in the suffix array $ESA_SA$.
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.ESA_ChildTab:The child table.
...remarks:The child table contains structural information of the suffix tree (see Abhouelda et al.).
...remarks:@Metafunction.Fibre@ returns a @Class.String@ over the alphabet of a size type.

..tag.ESA_BWT:The Burrows-Wheeler table.
...remarks:The Burrows-Wheeler table contains the Burrows-Wheeler transformation of $ESA_RawText$.
The entries are the characters left of the corresponding suffix in the suffix array $ESA_SA$.
...remarks:@Metafunction.Fibre@ returns the same type for $ESA_RawText$ and for $ESA_BWT$.

..see:Metafunction.Fibre
..see:Function.getFibre
..see:Spec.Index_ESA
*/

///.Metafunction.Fibre.param.TSpec.type:Tag.ESA Index Fibres

	typedef Fibre_Text		ESA_Text;
	typedef Fibre_RawText	ESA_RawText;
	typedef Fibre_SA		ESA_SA;
	typedef Fibre_RawSA		ESA_RawSA;
	typedef Fibre_SAE		ESA_SAE;
	typedef Fibre_LCP		ESA_LCP;
	typedef Fibre_LCPE		ESA_LCPE;
	typedef Fibre_ChildTab	ESA_ChildTab;
	typedef Fibre_BWT		ESA_BWT;


//////////////////////////////////////////////////////////////////////////////
// ESA index

/**
.Spec.Index_ESA:
..summary:An index based on an enhanced suffix array.
..cat:Index
..general:Class.Index
..signature:Index<TText, Index_ESA<> >
..param.TText:The text type.
...type:Class.String
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array (see @Tag.ESA Index Fibres.ESA_SA@), a lcp table (see @Tag.ESA Index Fibres.ESA_LCP@), etc.
..remarks:This index can be accessed as a Suffix Tree using the @Spec.VSTree Iterator@ classes.
*/

/*
	already defined in index_base.h

	template <typename TSpec = void>
	struct Index_ESA;
*/

	template < typename TText, typename TSpec >
	class Index<TText, Index_ESA<TSpec> > {
	public:
		Holder<typename Fibre<Index, ESA_Text>::Type>	text;
		typename Fibre<Index, ESA_SA>::Type				sa;			// suffix array 
		typename Fibre<Index, ESA_LCP>::Type			lcp;		// longest-common-prefix table
		typename Fibre<Index, ESA_LCPE>::Type			lcpe;		// extended lcp table
		typename Fibre<Index, ESA_ChildTab>::Type		childtab;	// child table (tree topology)
		typename Fibre<Index, ESA_BWT>::Type			bwt;		// burrows-wheeler table
		typename Cargo<Index>::Type						cargo;		// user-defined cargo

		Index() {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text) {}
	};

//////////////////////////////////////////////////////////////////////////////

	template < typename TText, typename TSpec >
	void _indexRequireTopDownIteration(Index<TText, Index_ESA<TSpec> > &index) 
	{
		indexRequire(index, ESA_SA());
		indexRequire(index, ESA_LCP());
		indexRequire(index, ESA_ChildTab());
	}

	template < typename TText, typename TSpec >
	void _indexRequireBottomUpIteration(Index<TText, Index_ESA<TSpec> > &index) 
	{
		indexRequire(index, ESA_SA());
		indexRequire(index, ESA_LCP());
	}

//////////////////////////////////////////////////////////////////////////////
///.Function.clear.param.object.type:Class.Index

	template <typename TText, typename TSpec>
	inline void clear(Index<TText, Index_ESA<TSpec> > &index) {
		clear(getFibre(index, ESA_SA()));
		clear(getFibre(index, ESA_LCP()));
		clear(getFibre(index, ESA_LCPE()));
		clear(getFibre(index, ESA_ChildTab()));
		clear(getFibre(index, ESA_BWT()));
	}


//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!open(getFibre(index, ESA_Text()), toCString(name), openMode)) && 
			(!open(getFibre(index, ESA_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	open(getFibre(index, ESA_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".lcp");	open(getFibre(index, ESA_LCP()), toCString(name), openMode);
		name = fileName;	append(name, ".child");	open(getFibre(index, ESA_ChildTab()), toCString(name), openMode);
		name = fileName;	append(name, ".bwt");	open(getFibre(index, ESA_BWT()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName) 
	{
		return open(index, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!save(getFibre(index, ESA_Text()), toCString(name), openMode)) && 
			(!save(getFibre(index, ESA_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	save(getFibre(index, ESA_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".lcp");	save(getFibre(index, ESA_LCP()), toCString(name), openMode);
		name = fileName;	append(name, ".child");	save(getFibre(index, ESA_ChildTab()), toCString(name), openMode);
		name = fileName;	append(name, ".bwt");	save(getFibre(index, ESA_BWT()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName)
	{
		return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif
