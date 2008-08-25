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
  $Id: repeat_base.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_REPEAT_BASE_H
#define SEQAN_HEADER_REPEAT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template <typename TPos, typename TPeriod>
	struct Repeat {
		TPos		beginPosition;
		TPos		endPosition;
		TPeriod		period;
	};

	template <typename TPos, typename TPeriod>
	struct Value< Repeat<TPos, TPeriod> > {
		typedef TPos Type;
	};

	template <typename TPos, typename TPeriod>
	struct Size< Repeat<TPos, TPeriod> > {
		typedef TPeriod Type;
	};



	template <typename TSize>
	struct RepeatFinderParams {
		TSize minRepeatLen;
		TSize maxPeriod;
	};

	// custom TSpec for our customized wotd-Index
	struct TRepeatFinder;

	template <typename TText>
	struct Cargo<Index<TText, Index_Wotd<TRepeatFinder> > > 
	{
		typedef Index<TText, Index_Wotd<TRepeatFinder> >	TIndex;
		typedef typename Size<TIndex>::Type					TSize;
		typedef RepeatFinderParams<TSize>					Type;
	};


	// node predicate
	template <typename TText, typename TSpec>
	bool nodePredicate(Iter<Index<TText, Index_Wotd<TRepeatFinder> >, TSpec> &it) 
	{
//		return countOccurrences(it) * nodeDepth(it) >= cargo(container(it)).minRepeatLen;
		return countOccurrences(it) * repLength(it) >= cargo(container(it)).minRepeatLen;
	}

	// monotonic hull
	template <typename TText, typename TSpec>
	bool nodeHullPredicate(Iter<Index<TText, Index_Wotd<TRepeatFinder> >, TSpec> &it) 
	{
//		return nodeDepth(it) <= cargo(container(it)).maxPeriod;
		return repLength(it) <= cargo(container(it)).maxPeriod;
	}

	template <typename TPos>
	struct _RepeatLess : public ::std::binary_function<TPos, TPos, bool>
	{
		// key less
		inline bool operator() (TPos const &a, TPos const &b) {
			return posLess(a, b);
		}
	};

	// main function
	template <typename TRepeatStore, typename TText, typename TRepeatSize, typename TPeriodSize>
	void findRepeats(TRepeatStore &repString, TText const &text, TRepeatSize minRepeatLen, TPeriodSize maxPeriod) 
	{
		typedef Index<TText, Index_Wotd<TRepeatFinder> >					TIndex;
		typedef typename Size<TIndex>::Type									TSize;
		typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type	TNodeIterator;
		typedef typename Fibre<TIndex, Fibre_SA>::Type const				TSA;
		typedef typename Infix<TSA>::Type									TOccString;
		typedef typename Iterator<TOccString>::Type							TOccIterator;

		typedef typename Value<TRepeatStore>::Type							TRepeat;
		typedef typename Value<TOccString>::Type							TOcc;

		typedef ::std::map<TOcc,TRepeat,_RepeatLess<TOcc> >					TRepeatList;

		if (maxPeriod < 1) return;
		if (maxPeriod == 1) 
		{
			findRepeats(repString, text, minRepeatLen);
			return;
		}

		TIndex		index(text);
		TRepeatList list;

		// set repeat finder parameters
		cargo(index).minRepeatLen = minRepeatLen;
		cargo(index).maxPeriod = maxPeriod;

		TNodeIterator nodeIt(index);
		TOccIterator itA, itB, itRepBegin, itEnd;
		TRepeat rep;
		for (; !atEnd(nodeIt); goNext(nodeIt))
		{
			if (isRoot(nodeIt)) continue;

			// get occurrences
			TOccString occ = getOccurrences(nodeIt);
			itA = begin(occ, Standard());
			itEnd = end(occ, Standard());
			itRepBegin = itB = itA;

			TSize repLen = repLength(nodeIt);		// representative length
			if ((TSize)minRepeatLen <= repLen) continue;

			TSize diff, period = 0;					// period of current repeat
			TSize repeatLen = 0;					// overall length of current repeat
			TSize minLen = minRepeatLen - repLen;	// minimum repeat length minus length of representative

			for (++itB; itB != itEnd; ++itB)
			{
				diff = posSub(*itB, *itA);
				if (diff != period || getSeqNo(*itA) != getSeqNo(*itB))
				{
					// is the repeat long enough?
					if (repeatLen >= minLen)
						// is the repeat self overlapping or connected?
						if (parentRepLength(nodeIt) < period && period <= repLen)
						{
							// insert repeat
							rep.beginPosition = *itRepBegin;
							rep.endPosition = posAdd(*itA, period);
							rep.period = period;
//							::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
							list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
						}
					itRepBegin = itA;
					period = diff;
					repeatLen = 0;
				}
				repeatLen += period;
				itA = itB;
			}

			// is the last repeat long enough?
			if (repeatLen >= minLen)
				// is the repeat self overlapping or connected?
				if (parentRepLength(nodeIt) < period && period <= repLen)
				{
					// insert repeat
					rep.beginPosition = *itRepBegin;
					rep.endPosition = posAdd(*itA, period);
					rep.period = period;
//					::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
					list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
				}
		}

		// copy low-complex regions to result string
		resize(repString, list.size());
		typename TRepeatList::const_iterator lit = list.begin();
		typename TRepeatList::const_iterator litEnd = list.end();
		for (TSize i = 0; lit != litEnd; ++lit, ++i)
			repString[i] = (*lit).second;
	}

	// period-1 optimization
	template <typename TRepeatStore, typename TString, typename TSpec, typename TRepeatSize>
	void findRepeats(TRepeatStore &repString, StringSet<TString, TSpec> const &text, TRepeatSize minRepeatLen) 
	{
		typedef typename Value<TRepeatStore>::Type	TRepeat;
		typedef typename Iterator<TString>::Type	TIterator;
		typedef typename Value<TString>::Type		TValue;
		typedef typename Size<TString>::Type		TSize;

		TRepeat rep;
		rep.period = 1;
		clear(repString);

		for( unsigned i = 0; i < length(text); ++i)
		{
			TIterator it = begin(text[i], Standard());
			TIterator itEnd = end(text[i], Standard());
			if (it == itEnd) continue;

			TValue last = *it;
			TSize repLeft = 0;
			TSize repRight = 1;
			rep.beginPosition.i1 = i;
			rep.endPosition.i1 = i;

			for (++it; it != itEnd; ++it, ++repRight) 
			{
				if (last != *it)
				{
					if ((TRepeatSize)(repRight - repLeft) > minRepeatLen)
					{
						// insert repeat
						rep.beginPosition.i2 = repLeft;
						rep.endPosition.i2 = repRight;
//						::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
						appendValue(repString, rep);
					}
					repLeft = repRight;
					last = *it;
				}
			}
			if ((TRepeatSize)(repRight - repLeft) > minRepeatLen)
			{
				// insert repeat
				rep.beginPosition.i2 = repLeft;
				rep.endPosition.i2 = repRight;
//				::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
				appendValue(repString, rep);
			}
		}
	}

}	// namespace seqan

#endif
