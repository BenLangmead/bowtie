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
  $Id: index_dfi.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_DFI_H
#define SEQAN_HEADER_INDEX_DFI_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// simple struct to store frequency vectors

	struct _DFIEntry
	{
		unsigned						lastSeqSeen;
		String<unsigned, Array<2> >		freq;
//		String<unsigned>				freq;
	};


//////////////////////////////////////////////////////////////////////////////
// constant frequency predicate

	template <bool RESULT>
	struct _DFIPredDefault 
	{
        inline bool operator()(_DFIEntry const &) const {
			return RESULT;
		}
    };


//////////////////////////////////////////////////////////////////////////////
// DFI - The Deferred Frequency Index

/**
.Spec.Index_DFI:
..summary:The Deferred Frequency Index (see Weese and Schulz, "Efficient string mining under constraints via the
deferred frequency index").
..cat:Index
..general:Class.Index
..signature:Index<TText, Index_Wotd< WotdDFI<TPredHull, TPred> > >
..param.TText:The text type.
...type:Class.String
..param.TPred:An arbitrary frequency predicate
..param.TPredHull:A monotonic hull of $TPred$
..remarks:This index is based on a lazy suffix tree (see @Spec.Index_Wotd@).
All $TPredHull$ sufficing nodes can be iterated using a @Spec.TopDown Iterator@.
To iterate the exact solution set of $TPred$, use a $Spec.TopDownHistory Iterator$ of this index.
*/

	template < 
		typename TPredHull = _DFIPredDefault<true>,
		typename TPred = _DFIPredDefault<true>
	>
	struct WotdDFI;

	template < 
		typename TObject,
		typename TPredHull,
		typename TPred
	>
	class Index<TObject, Index_Wotd< WotdDFI<TPredHull, TPred> > >:
		public Index<TObject, Index_Wotd<> > 
	{
	public:
		
		typedef Index<TObject, Index_Wotd<> >	TBase;

		// extending base class
		typedef typename TBase::TText	TText;
		typedef typename TBase::TValue	TValue;
		typedef typename TBase::TSize	TSize;
		
		using TBase::LEAF;
		using TBase::LAST_CHILD;
		using TBase::UNEVALUATED;
		using TBase::SENTINELS;

		// frequency strings for DFI
		typedef _DFIEntry						TDFIEntry;
		typedef String<
			TDFIEntry, 
			Array<ValueSize<TValue>::VALUE> >	TDFIEntries;
		typedef String<unsigned, Array<3> >		TDFIDatasets;

		// 1st word flags
		static TSize const DFI_PRED_HULL = (TSize)1 << (BitsPerValue<TSize>::VALUE - 3); // this node fulfills all monotonic frequency predicates (e.g. min_freq)
		static TSize const DFI_PRED      = (TSize)1 << (BitsPerValue<TSize>::VALUE - 4); // this node fulfills all frequency predicates (e.g. emerging, minmax_freq)

		static TSize const BITMASK0      = ~(LEAF | LAST_CHILD | DFI_PRED_HULL | DFI_PRED);
		static TSize const BITMASK1      = ~(UNEVALUATED | SENTINELS);

		TDFIEntries		entry;
		TDFIDatasets	ds;

		TPredHull		predHull;
		TPred			pred;
		
		Index() {}

		template <typename _TText>
		Index(_TText &_text):
			TBase(_text) {}

		template <typename _TText>
		Index(_TText &_text, TPredHull const &_predHull):
			TBase(_text),
			predHull(_predHull) {}

		template <typename _TText>
		Index(_TText &_text, TPredHull const &_predHull, TPred const &_pred):
			TBase(_text),
			predHull(_predHull),
			pred(_pred) {}
	};

	
	template < 
		typename TText, 
		typename TPredHull,
		typename TPred,
		typename TSpec
	>
	inline bool nodePredicate(
		Iter<Index<TText, Index_Wotd< WotdDFI<TPredHull, TPred> > >, TSpec> const &it)
	{
		typedef Index<TText, Index_Wotd< WotdDFI<TPredHull, TPred> > > TIndex;
		return dirAt(value(it).node, container(it)) & TIndex::DFI_PRED;
	}

	template < 
		typename TText, 
		typename TPredHull,
		typename TPred,
		typename TSpec
	>
	inline bool nodeHullPredicate(
		Iter<Index<TText, Index_Wotd< WotdDFI<TPredHull, TPred> > >, TSpec> const &it)
	{
		typedef Index<TText, Index_Wotd< WotdDFI<TPredHull, TPred> > > TIndex;
		return dirAt(value(it).node, container(it)) & TIndex::DFI_PRED_HULL;
	}


//////////////////////////////////////////////////////////////////////////////
// we modify counting sort used in the wotd-index
// to count the frequencies in passing

	template < typename TText, typename TSpec, typename TPredHull, typename TPred >
	typename Size< Index<StringSet<TText, TSpec>, Index_Wotd<WotdDFI<TPredHull, TPred> > > >::Type
	_sortFirstWotdBucket(Index<StringSet<TText, TSpec>, Index_Wotd<WotdDFI<TPredHull, TPred> > > &index)
	{
	SEQAN_CHECKPOINT
		typedef Index<StringSet<TText, TSpec>, Index_Wotd<WotdDFI<TPredHull, TPred> > >	TIndex;
		typedef typename Fibre<TIndex, Wotd_SA >::Type			TSA;
		typedef typename TIndex::TCounter						TCounter;

		typedef typename TIndex::TDFIEntry						TDFIEntry;
		typedef typename TIndex::TDFIDatasets					TDFIDatasets;
		typedef typename Iterator<TDFIDatasets, Standard>::Type	TDFIDatasetsIterator;

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

		for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i) {
			TDFIEntry &entry = index.entry[i];
			entry.lastSeqSeen = -1;
			for(unsigned  j = 0; j < length(entry.freq); ++j)
				entry.freq[j] = 0;
		}

		// 2. count characters
		_wotdCountChars(occ, stringSet);

		// 3. cummulative sum
		TSize requiredSize = _wotdCummulativeSum(bound, occ);

		// 4. fill suffix array
		unsigned dsNo = 0;
		TDFIDatasetsIterator currentDS = begin(index.ds, Standard()) + 1;
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			// search for surrounding dataset
			while (seqNo >= *currentDS) {
				++dsNo;
				++currentDS;
			}

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
				unsigned ord = ordValue(*itText);
				TDFIEntry &entry = index.entry[ord];
				// new sequence is seen for <ord> character
				// -> increment frequency of current dataset
				if (entry.lastSeqSeen != seqNo) 
				{
					entry.lastSeqSeen = seqNo;
					++entry.freq[dsNo];
				}

				*(saBeg + (*(boundBeg + ord))++) = localPos;
				assignValueI2(localPos, getValueI2(localPos) + 1);
			}
		}
		index.sentinelOcc = 0;
		index.sentinelBound = 0;

		return requiredSize;
	}

	// sort bucket using radixsort
	// - all buckets are in lexicographical order
	// - SA[left,right) contains real SA entries (the beginning positions of the suffices)
	template < typename TText, typename TSpec, typename TPredHull, typename TPred, typename TSize >
	TSize _sortWotdBucket(
		Index<StringSet<TText, TSpec>, Index_Wotd<WotdDFI<TPredHull, TPred> > > &index,
		TSize left, 
		TSize right,
		TSize prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef Index<StringSet<TText, TSpec>, Index_Wotd<WotdDFI<TPredHull, TPred> > >	TIndex;
		typedef typename Fibre<TIndex, Wotd_SA >::Type				TSA;
		typedef typename TIndex::TCounter							TCounter;
		typedef typename TIndex::TTempSA							TTempSA;
		typedef typename TIndex::TDFIEntry							TDFIEntry;
		typedef typename TIndex::TDFIDatasets						TDFIDatasets;
		typedef typename Iterator<TDFIDatasets, Standard>::Type		TDFIDatasetsIterator;

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

		for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i) {
			TDFIEntry &entry = index.entry[i];
			entry.lastSeqSeen = -1;
			for(unsigned  j = 0; j < length(entry.freq); ++j)
				entry.freq[j] = 0;
		}

		index.sentinelOcc = 0;
		index.sentinelBound = 0;

		// 2. count characters
		{
			TDFIDatasetsIterator currentDS = begin(index.ds, Standard()) + 1;
			TTextIterator itText = TTextIterator();
			TTempSAIterator itSA = begin(tempSA, Standard());
			TTempSAIterator itSAEnd = end(tempSA, Standard());
			TTextSize textLength = 0;
			unsigned lastSeqSeen = -1;
			unsigned dsNo = 0;
			Pair<unsigned, TTextSize> lPos;
			for (; itSA != itSAEnd; ++itSA) 
			{
				posLocalize(lPos, *itSA, stringSetLimits(index));
				if (lastSeqSeen != getSeqNo(lPos))
				{
					lastSeqSeen = getSeqNo(lPos);

					// search for surrounding dataset
					while (lastSeqSeen >= *currentDS) {
						++dsNo;
						++currentDS;
					}

					// shift textBegin and textLength by prefixLen
					textLength = length(stringSet[lastSeqSeen]) - prefixLen;
					itText = begin(stringSet[lastSeqSeen], Standard()) + prefixLen;
				}
				if (textLength > getSeqOffset(lPos)) {
					unsigned ord = ordValue(*(itText + getSeqOffset(lPos)));
					TDFIEntry &entry = index.entry[ord];
					// new sequence is seen for <ord> character
					// -> increment frequency of current dataset
					if (entry.lastSeqSeen != lastSeqSeen) 
					{
						entry.lastSeqSeen = lastSeqSeen;
						++entry.freq[dsNo];
					}
					++*(occBeg + ord);
				} else
					if (textLength == getSeqOffset(lPos)) ++index.sentinelOcc;
			}
		}

		// 3. cumulative sum
		TSize requiredSize = 0;
		if (index.sentinelOcc != 0)
			requiredSize = (index.sentinelOcc > 1)? 2: 1;

		requiredSize += _wotdCummulativeSum(bound, occ, left + index.sentinelOcc);
		index.sentinelBound = left;
/*
		::std::cout << "$=" << index.sentinelOcc<<"@"<<index.sentinelBound << "\t";
		for(int i=0; i<length(occ);++i)
			if (occ[i])
				::std::cout << i << "=" << occ[i]<<"@"<<bound[i] << "\t";
*/
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


	// store buckets into directory
	// storing SA links and topology links in Dir
	template <typename TText, typename TPredHull, typename TPred, typename TSize>
	inline void 
	_storeWotdChildren(
		Index<TText, Index_Wotd<WotdDFI<TPredHull, TPred> > > &index,
		TSize dirOfs,
		TSize lcp)
	{
	SEQAN_CHECKPOINT
		typedef Index<TText, Index_Wotd<WotdDFI<TPredHull, TPred> > >	TIndex;

		typedef typename Fibre<TIndex, Wotd_Dir>::Type		TDir;
		typedef typename TIndex::TCounter					TCounter;
		typedef typename TIndex::TDFIEntries				TEntries;

		typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
		typedef typename Iterator<TCounter, Standard>::Type	TCntIterator;
		typedef typename Iterator<TEntries, Standard>::Type	TEntriesIterator;

		typedef typename Value<TCounter>::Type				TCntValue;
		typedef typename Value<TDir>::Type					TDirValue;

		TDirIterator itDirBegin = begin(indexDir(index), Standard()) + dirOfs;
		TDirIterator itDirEnd = end(indexDir(index), Standard());
		TDirIterator itDir = itDirBegin;
		TDirIterator itPrev = itDirEnd;

		TCntIterator it = begin(index.tempOcc, Standard());
		TCntIterator bit = begin(index.tempBound, Standard());
		TCntIterator itEnd = end(index.tempOcc, Standard());
		TEntriesIterator itEntry = begin(index.entry, Standard());

		TCntValue occ;
		if (index.sentinelOcc != 0)
		{
			TDirValue orMask = (index.predHull(*itEntry))? index.DFI_PRED_HULL: 0;
			if (index.pred(*itEntry)) orMask |= index.DFI_PRED;

			if (index.sentinelOcc > 1) { // occurs on multiseqs
				itPrev = itDir;
				*itDir = (index.sentinelBound - index.sentinelOcc) | orMask;	++itDir;
				*itDir = index.sentinelBound | index.UNEVALUATED;				++itDir;
			} else {
				itPrev = itDir;
				*itDir = (index.sentinelBound - index.sentinelOcc) | index.LEAF | orMask;
				++itDir;
			}
		}

		for (; it != itEnd; ++it, ++bit, ++itEntry)
		{
			if ((occ = *it) == 0) continue;

			TDirValue orMask = (index.predHull(*itEntry))? index.DFI_PRED_HULL: 0;
			if (index.pred(*itEntry)) orMask |= index.DFI_PRED;

			if (occ > 1) {
				itPrev = itDir;
				*itDir = (*bit - occ) | orMask;					++itDir;
				*itDir = *bit | index.UNEVALUATED;				++itDir;
			} else {
				itPrev = itDir;
				*itDir = (*bit - occ) | index.LEAF | orMask;	++itDir;
			}
		}

		// first child gets the mutual lcp value of the children (== parent repLength)
		*itDirBegin = ((*itDirBegin) & ~index.BITMASK0) | lcp;

		// mark the last child
		if (itPrev != itDirEnd)
			*itPrev |= index.LAST_CHILD;
	}

//////////////////////////////////////////////////////////////////////////////
// debug output

	template <typename TText, typename TPredHull, typename TPred>
	inline void
	_dump(Index<TText, Index_Wotd< WotdDFI<TPredHull, TPred> > > &index)
	{
		typedef Index_Wotd< WotdDFI<TPredHull, TPred> >		TSpec;
		typedef Index<TText, TSpec>							TIndex;
		typedef typename Fibre<TIndex, Wotd_Dir>::Type		TDir;
		typedef typename Value<TDir>::Type					TDirValue;

		::std::cout << "  Dir (wotd/DFI)" << ::std::endl;
		for(unsigned i=0; i < length(indexDir(index)); ++i) {
			TDirValue d = indexDir(index)[i];
			::std::cout << i << ":  " << (d & index.BITMASK0);
			if (d & index.LEAF)				::std::cout << "  (Leaf/Uneval)";
			if (d & index.LAST_CHILD)		::std::cout << "  (LastChild/SENTINELS)";
			if (d & index.DFI_PRED_HULL)	::std::cout << "  (PRED_HULL)";
			if (d & index.DFI_PRED)			::std::cout << "  (PRED)";
			::std::cout << ::std::endl;
		}

		::std::cout << ::std::endl << "  SA" << ::std::endl;
		for(unsigned i=0; i < length(indexSA(index)); ++i)
			::std::cout << i << ":  " << indexSA(index)[i] << "  " << suffix(indexText(index), indexSA(index)[i]) << ::std::endl;

		::std::cout << ::std::endl;
	}

	template <typename TText, typename TPredHull, typename TPred>
	inline void
	_dumpFreq(Index<TText, Index_Wotd< WotdDFI<TPredHull, TPred> > > &index)
	{
		typedef WotdDFI<TPredHull, TPred> TSpec;
		typedef Index<TText, Index_Wotd<TSpec> >				TIndex;
		typedef typename Value<TIndex>::Type					TValue;

		for(unsigned i=0; i<length(index.tempOcc); ++i)
			if (index.tempOcc[i] != 0) {
				::std::cout << "  Freq[" << (TValue)i << "] = (";
				for(unsigned d=0; d<length(index.entry[i].freq); ++d) {
					if (d>0) ::std::cout << ",";
					::std::cout << index.entry[i].freq[d];
				}
				::std::cout << ")" << ::std::endl;
			}
	}

//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation 

	template <typename TText, typename TPredHull, typename TPred>
	inline bool indexCreate(Index<TText, Index_Wotd<WotdDFI<TPredHull, TPred> > > &index, Wotd_SA const, Default const)
	{
		typedef Index<TText, Index_Wotd<WotdDFI<TPredHull, TPred> > >	TIndex;
		typedef typename Value<TIndex>::Type							TValue;
		typedef typename TIndex::TBase									TBase;

		resize(index.entry, (unsigned) ValueSize<TValue>::VALUE);
		if (empty(index.ds)) {
			resize(index.ds, 2);
			index.ds[0] = 0;
			index.ds[1] = countSequences(index);
		}
		for(unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
			resize(index.entry[i].freq, length(index.ds) - 1);

		_wotdCreateFirstLevel(index);
		return true;
	}
}

#endif //#ifndef SEQAN_HEADER_...
