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
  $Id: find_pex.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_PEX_H
#define SEQAN_HEADER_FIND_PEX_H

// uncomment this for verbose debug output
//#define SEQAN_DEBUG_PEX

namespace SEQAN_NAMESPACE_MAIN 
{

struct Hierarchical;			
struct NonHierarchical;

template <typename TVerification, typename TMultiFinder = AhoCorasick>
struct Pex;

typedef Pex<Hierarchical,AhoCorasick>      PexHierarchical;
typedef Pex<NonHierarchical,AhoCorasick>   PexNonHierarchical;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.PexMultiFinder:
..summary:Determines the multiple exact string matching algorithm used by the Pex algorithm.
..signature:PexMultiFinder< Pattern<TNeedle,Pex<TVerification,TMultiFinder> > >::Type
..param.TMultiFinder:The specification for the multiple exact string matching algorithm that should be used with the Pex algorithm.
..returns.param.Type:Pattern type of the multiple exact string matching algorithm for the specified Pattern.
..see:Spec.Pex
..remarks: For a description of Pattern usage see @Class.Pattern@.
..remarks: Overload this Metafunction if you want to use something else for verification then $Pattern<String<Segment<TNeedle> > , TMultiFinder>$.
*/

template<typename T>
struct PexMultiFinder;

template<typename TNeedle,typename TVerification, typename TMultiFinder>
struct PexMultiFinder< Pattern<TNeedle, Pex<TVerification , TMultiFinder > > >
{
  typedef Pattern<String<Segment<TNeedle> > , TMultiFinder> Type;
};
 
//////////////////////////////////////////////////////////////////////////////

template<typename TPosition,typename TScore,typename TVerifier,typename TNeedle>
struct _PexRange{
  TPosition start,end;
  TScore error;
  TVerifier verifier;
  
  _PexRange()
  {}
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Pex:
..summary: Provides a fast approximate string matching filter that splits the needle into several pieces that are searched with a multiple exact string matching algorithm and later verified.
..general:Class.Pattern
..cat:Searching
..signature:Pattern<TNeedle, Pex<TVerification,TMultiFinder> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TVerification: Determines if the hierarchical verification proposed by Navarro and Baeza-Yates is used or not.
..param.TMultiFinder: Specifies the algorithm for the multiple exact string matching algorithm
...type:Spec.AhoCorasick
..remarks: There are two defaults available $PexHierarchical$ and $PexNonHierarchical$ (e.g. $Pattern<String<char> ,PexHierarchical> $) that both use the @Spec.AhoCorasick@ algorithm for the multiple exact string matching.
*/
///.Class.Pattern.param.TSpec.type:Spec.Pex

/**
.Spec.Hierarchical:
..summary: By using this Specialization the hierarchical verification is enabled.
..general:Spec.Pex
..cat:Searching
..signature:Pattern<TNeedle, Pex<Hierarchical,TMultiFinder> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TMultiFinder: Specifies the algorithm for the multiple exact string matching algorithm
*/
///.Spec.Pex.param.TVerification.type:Spec.Hierarchical
/**
.Spec.NonHierarchical:
..summary: By using this Specialization the hierarchical verification is disabled.
..general:Spec.Pex
..cat:Searching
..signature:Pattern<TNeedle, Pex<NonHierarchical,TMultiFinder> >
..param.TNeedle:The needle type.
...type:Class.String
..param.TMultiFinder: Specifies the algorithm for the multiple exact string matching algorithm
*/
///.Spec.Pex.param.TVerification.type:Spec.NonHierarchical

template <typename TNeedle, typename TVerification, typename TMultiFinder>
class Pattern<TNeedle, Pex<TVerification,TMultiFinder > > 
{
 public:
   typedef typename Position<TNeedle>::Type TPosition;
   typedef unsigned TScore;
   typedef Pattern<TNeedle,MyersUkkonen> TVerifier;
   typedef typename PexMultiFinder< 
                       Pattern<TNeedle, Pex<TVerification,TMultiFinder > > 
                                  >::Type TMFinder; 
  
   // the maximal accepted error
   TScore limit;
   // reference to the needle
   Holder<TNeedle> data_needle;
   // pattern object for the multi pattern search
   TMFinder multiPattern;
   // needles for the multi pattern search
   String<Segment<TNeedle> >  splitted_needles;
   
   // data store for the verification tree respectively the splitted needle
   ::std::map<unsigned, _PexRange<TPosition,TScore,TVerifier,TNeedle> > range_table;
   // map leafs of the tree to parts of the needle
   ::std::map<unsigned, unsigned> leaf_map;

   // store the infixes for the verifiers
   String<Segment<TNeedle> > segment_store;
  
   // track position where the last occurence was found
   unsigned lastFPos;
   unsigned lastFNdl;
   
   // indicator to track if we already found an occurence
   bool findNext,patternNeedsInit; 

   unsigned needleLength;

   Pattern() {}

   template <typename TNeedle2>
   Pattern(TNeedle2 const & ndl)
     : limit(1)
   {
     setHost(*this, ndl);
   }

   template <typename TNeedle2>
   Pattern(TNeedle2 const & ndl, int _limit = -1)
     : limit(- _limit)
   {
SEQAN_CHECKPOINT
     setHost(*this, ndl);
   }

   ~Pattern() {
     SEQAN_CHECKPOINT
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2, typename TVerification, typename TMultiFinder>
void setHost (Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me, TNeedle2 const & needle) 
{
  // initialisation of the find-tree etc. will be done when patternInit
  // is called to assure that we already know the scoreLimit
  me.data_needle = needle;
  me.needleLength = length(needle);
  me.findNext = false;
  me.patternNeedsInit = true;
}

template <typename TNeedle, typename TNeedle2, typename TVerification, typename TMultiFinder>
void setHost (Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me, TNeedle2 & needle)
{
  setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TVerification, typename TMultiFinder>
inline typename Host<Pattern<TNeedle, Pex<TVerification,TMultiFinder > > const>::Type & 
host(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me)
{
SEQAN_CHECKPOINT
  return value(me.data_needle);
}

template <typename TNeedle, typename TVerification, typename TMultiFinder>
inline typename Host<Pattern<TNeedle, Pex<TVerification,TMultiFinder > > const>::Type & 
host(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > const & me)
{
SEQAN_CHECKPOINT
  return value(me.data_needle);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TMultiFinder>
int _getRoot(Pattern<TNeedle, Pex<NonHierarchical, TMultiFinder > > & me) 
{
SEQAN_CHECKPOINT
  return length(me.splitted_needles);
}

template <typename TNeedle, typename TMultiFinder>
int _getRoot(Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > &) 
{
SEQAN_CHECKPOINT
  return 1;
}

//////////////////////////////////////////////////////////////////////////////
///.Function.getScore.param.pattern.type:Spec.Pex

template <typename TNeedle, typename TVerification, typename TMultiFinder>
int getScore(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me) 
{
SEQAN_CHECKPOINT
  return getScore(me.range_table[_getRoot(me)].verifier);
}

//////////////////////////////////////////////////////////////////////////////
///.Function.scoreLimit.param.pattern.type:Spec.Pex

template <typename TNeedle, typename TVerification, typename TMultiFinder>
inline int 
scoreLimit(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > const & me)
{
SEQAN_CHECKPOINT
  return - (int) me.limit;
}


//////////////////////////////////////////////////////////////////////////////
///.Function.setScoreLimit.param.pattern.type:Spec.Pex

template <typename TNeedle, typename TScoreValue,typename TVerification, typename TMultiFinder>
inline void 
setScoreLimit(Pattern<TNeedle, Pex<TVerification,TMultiFinder > > & me, 
			  TScoreValue _limit)
{
SEQAN_CHECKPOINT
  me.patternNeedsInit = true;
  me.limit = (- _limit);
}

//////////////////////////////////////////////////////////////////////////////
//   PexNonHierarchical -- functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TFinder, typename TMultiFinder>
void _patternInit(Pattern<TNeedle, Pex<NonHierarchical, TMultiFinder > > &me, TFinder &)
{
SEQAN_CHECKPOINT
  typedef typename Position<TNeedle>::Type TPosition;
  typedef unsigned TScore;
  typedef Pattern<TNeedle,MyersUkkonen> TVerifier;

  // split pattern
  unsigned k = me.limit + 1;
  unsigned seg_len = me.needleLength / k; //::std::floor(me.needleLength/k); //Wer hat den floor-Quatsch hier programmiert?
  
  // 
  clear(me.splitted_needles);
  clear(me.range_table);
  clear(me.segment_store);
  unsigned s = 0;
  unsigned c = 0;
  unsigned i = 0;
  while(s < me.needleLength)
  { 
    _PexRange<TPosition,TScore,TVerifier,TNeedle> pr;
    pr.start = s;
    pr.end = (c == me.limit ? me.needleLength : s + seg_len);
    pr.error = 0;

    insert(me.range_table,i,pr);
    appendValue(me.splitted_needles,infix(value(me.data_needle),pr.start,pr.end));
    s += (c == me.limit ? me.needleLength : seg_len);
    ++c;
    ++i;
  }

  me.lastFPos = 0;
  me.lastFNdl = 0;

  // insert complete needle in range table to use the verifier
  appendValue(me.segment_store,infix(value(me.data_needle),0,me.needleLength));
  _PexRange<TPosition,TScore,TVerifier,TNeedle> pr;
  pr.start = 0;
  pr.end = me.needleLength;
  pr.error = me.limit;
  setHost(pr.verifier,me.segment_store[0]);
  setScoreLimit(pr.verifier, - static_cast<int>(me.limit));  
  insert(me.range_table,length(me.splitted_needles),pr);
  
  // init multipattern finder
  setHost(me.multiPattern,me.splitted_needles);
  
  me.patternNeedsInit = false;
  me.findNext = false;

#ifdef SEQAN_DEBUG_PEX
  ::std::cout << " -------------------------------------------------  " << ::std::endl;
  ::std::cout << "                   PATTERN INIT                     " << ::std::endl;
  ::std::cout << "Needle:   " << value(me.data_needle) << ::std::endl;
  ::std::cout << "|Needle|: " << me.needleLength << ::std::endl;
  ::std::cout << "seg_len:  " << seg_len << ::std::endl;
  ::std::cout << "limit:    " << me.limit << ::std::endl;
  ::std::cout << "k:        " << k << ::std::endl;
  ::std::cout << "computed following needles for multipattern search: " << ::std::endl;
  for(unsigned i = 0;i < length(me.splitted_needles);++i)  ::std::cout << me.splitted_needles[i] << ::std::endl;
  ::std::cout << " -------------------------------------------------  " << ::std::endl;
#endif  

}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TMultiFinder>
inline bool find (TFinder & finder, Pattern<TNeedle, Pex<NonHierarchical, TMultiFinder > > & me)
{
SEQAN_CHECKPOINT

  typedef typename Host<TFinder>::Type    THost;
  typedef Segment<THost>                  THostSegment;
  typedef Finder<THostSegment>            THSFinder;
  TFinder mf(finder);
  unsigned startPos;

  if (empty(finder))
  {
    _finderSetNonEmpty(finder);
  }
  if(me.patternNeedsInit)
  {
     _patternInit(me, finder);
  }

  if(me.findNext){
    startPos = position(finder);
    int start = me.lastFPos - me.range_table[me.lastFNdl].start - me.limit;
    int end   = me.lastFPos + (me.needleLength - me.range_table[me.lastFNdl].start) + me.limit;
    
    // adjust start and end if they point over the edges of host(finder)
    start = (start < 0 ? 0 : start);
    end = (end > static_cast<int>(length(host(finder))) ? length(host(finder)) : end);

    THostSegment s(infix(host(finder),start,end));
    THSFinder f(s);

    while(find(f,me.range_table[_getRoot(me)].verifier))
    {
      unsigned nP = start + position(f);
      if(nP > startPos){
	// compute new position
	unsigned offset = nP - position(finder);
	finder += offset;
	me.findNext = true;
	return true;
      }
    }
    // reset mf finder to old position
    unsigned mf_offset = position(finder) - me.lastFPos;
    mf -= mf_offset;
  }
  me.findNext = false;
  startPos = position(finder);

  while(find(mf,me.multiPattern))
  {
    int s = position(mf) - me.range_table[position(me.multiPattern)].start - me.limit;
    int e   = position(mf) + (me.needleLength - me.range_table[position(me.multiPattern)].start) + me.limit;

    // adjust start and end if they point over the edges of host(finder)
    s = (s < 0 ? 0 : s);
    e = (e > static_cast<int>(length(host(finder))) ? length(host(finder)) : e);

    THostSegment i(infix(host(mf),s,e));
    THSFinder f(i);
    while(find(f,me.range_table[_getRoot(me)].verifier))
    {
      unsigned nP = s + position(f);
      if(nP > startPos){
	// compute new position
	unsigned offset = nP - position(finder);
	finder += offset;
	me.lastFPos = position(mf);
	me.lastFNdl = position(me.multiPattern);
	me.findNext = true;
	return true;
      }
    }
  }
  // set finder to end position
  unsigned t = length(host(finder))- position(finder);
  finder += t;

  return false;
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//   PexHierarchical -- functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TMultiFinder>
void _createTree(Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > &me, unsigned start, unsigned end,
		 unsigned k, unsigned parent, unsigned direction ,unsigned idx, unsigned plen)
{
  //create tree like proposed in Navarro & Raffinot
  // direction == 0 .. choose left child in the tree
  // direction == 1 .. choose right child in the tree

#ifdef SEQAN_DEBUG_PEX
  ::std::cout << "called _createTree:" << ::std::endl;
  ::std::cout << "  start: " << start << ::std::endl;
  ::std::cout << "  end  : " << end << ::std::endl;
  ::std::cout << "  seq  : " << infix(value(me.data_needle),start,end + 1) << ::std::endl;
  ::std::cout << "  k    : " << k << ::std::endl;
  ::std::cout << "  paren: " << parent << ::std::endl;
  ::std::cout << "  direc: " << direction << ::std::endl;
  ::std::cout << "  idx  : " << idx << ::std::endl;
  ::std::cout << "  plen : " << plen << ::std::endl;
  ::std::cout << " ----------------------------- " << ::std::endl;
#endif
  typedef typename Position<TNeedle>::Type TPosition;
  typedef unsigned TScore;
  typedef Pattern<TNeedle,MyersUkkonen> TVerifier; 

  _PexRange<TPosition,TScore,TVerifier,TNeedle> pr;
  pr.start = start;
  pr.end = end;
  pr.error = k;

  appendValue(me.segment_store,infix(value(me.data_needle),pr.start,pr.end + 1));
  setScoreLimit(pr.verifier, - static_cast<int>(pr.error));
  setHost(pr.verifier, me.segment_store[length(me.segment_store) - 1]);
  
  unsigned left = k/2 + 1; //::std::ceil(static_cast<double>(k + 1)/2);
  unsigned cur_idx = (parent << 1) + direction;

  // insert pr into the tree
  insert(me.range_table,cur_idx,pr);
  
  if(k == 0){
    appendValue(me.splitted_needles,infix(value(me.data_needle),pr.start,pr.end + 1));
#ifdef SEQAN_DEBUG_PEX
    ::std::cout << "inserted : " << me.splitted_needles[length(me.splitted_needles) - 1] << " into splitted needles" << ::std::endl;
    ::std::cout << "assign to leaf_map " << length(me.splitted_needles) - 1 << " value " << cur_idx << ::std::endl;
    ::std::cout << " ----------------------------- " << ::std::endl;
#endif
    me.leaf_map[length(me.splitted_needles) - 1] = cur_idx;
  }else{
    // recusivly create the rest of the tree
//    _createTree(me, start, start + left * plen - 1, ::std::floor(static_cast<double>(left * k)/ static_cast<double>(k + 1)),cur_idx,0,idx,plen);
    _createTree(me, start, start + left * plen - 1, left * k / (k+1),cur_idx,0,idx,plen);
//    _createTree(me,  start + left * plen, end, ::std::floor(static_cast<double>((k + 1 - left)*k)/ static_cast<double>(k + 1)),cur_idx,1,idx + left,plen);
    _createTree(me,  start + left * plen, end, (k + 1 - left)*k / (k+1),cur_idx,1,idx + left,plen);
  }
}

template <typename TNeedle, typename TFinder, typename TMultiFinder>
void _patternInit(Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > &me, TFinder &)
{
SEQAN_CHECKPOINT
  typedef typename Position<TNeedle>::Type TPosition;
  typedef unsigned TScore;
  typedef Pattern<TNeedle,MyersUkkonen> TVerifier;

  unsigned k = me.limit + 1;
  unsigned plen = me.needleLength / k; //::std::floor(static_cast<double>(me.needleLength)/static_cast<double>(k));
  
  // reset
  clear(me.splitted_needles);
  clear(me.range_table);
  clear(me.leaf_map);
  clear(me.segment_store);

  // build the verification tree
  _createTree(me, 0, me.needleLength - 1,me.limit, 0, 1 , 0, plen);

  me.lastFPos = 0;
  me.lastFNdl = 0;
  setHost(me.multiPattern,me.splitted_needles);
  me.patternNeedsInit = false;
  me.findNext = false;

#ifdef SEQAN_DEBUG_PEX
  ::std::cout << " -------------------------------------------------  " << ::std::endl;
  ::std::cout << "                   PATTERN INIT                     " << ::std::endl;
  ::std::cout << "Needle:   " << value(me.data_needle) << ::std::endl;
  ::std::cout << "|Needle|: " << me.needleLength << ::std::endl;
  ::std::cout << "limit:    " << me.limit << ::std::endl;
  ::std::cout << "k:        " << k << ::std::endl;
  ::std::cout << "computed following needles for multipattern search: " << ::std::endl;
  for(unsigned i = 0;i < length(me.splitted_needles);++i)  ::std::cout << me.splitted_needles[i] << ::std::endl;
  ::std::cout << " -------------------------------------------------  " << ::std::endl;
#endif  
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFinder, typename TNeedle, typename TMultiFinder>
inline bool find (TFinder & finder, Pattern<TNeedle, Pex<Hierarchical, TMultiFinder > > & me)
{
SEQAN_CHECKPOINT

  typedef typename Host<TFinder>::Type    THost;
  typedef Segment<THost>                  THostSegment;
  typedef Finder<THostSegment>            THSFinder;
  TFinder mf(finder);
  unsigned startPos;

  if (empty(finder))
  {
    _finderSetNonEmpty(finder);
  }
  if(me.patternNeedsInit)
  {
     _patternInit(me, finder);
  }
  if(me.findNext){
    // we found an occurence
    startPos = position(finder);
    unsigned pnode = _getRoot(me); // use root 
    unsigned in = me.range_table[me.leaf_map[me.lastFNdl]].start;
    
    int p1 = me.lastFPos - (in - me.range_table[pnode].start) - me.range_table[pnode].error;
    int p2 = me.lastFPos + (me.range_table[pnode].end - in + 1) + me.range_table[pnode].error;

    // adjust start and end if they point over the edges of host(finder)
    p1 = (p1 < 0 ? 0 : p1);
    p2 = (p2 > static_cast<int>(length(host(finder))) ? length(host(finder)) : p2);
    THostSegment i(infix(host(mf),p1,p2));
    THSFinder f(i);

    while(find(f,me.range_table[pnode].verifier))
    {
      unsigned nP = p1 + position(f);
      if(nP > startPos)
	{
	  // compute new position
	  unsigned offset = nP - position(finder);
	  finder += offset;
	  me.findNext = true;
	  return true;
	}      
    }
    // reset mf finder to old position
    unsigned mf_offset = position(finder) - me.lastFPos;
    mf -= mf_offset;
  }
  me.findNext = false;
  startPos = position(finder);

  while(find(mf,me.multiPattern))
  {
    // get found leaf
    unsigned node = me.leaf_map[position(me.multiPattern)];
    unsigned in = me.range_table[node].start;
    node = node >> 1;
    bool cand = true;

    while( cand && node != 1) // stop when reaching root
    {
      int p1 = position(mf) - (in - me.range_table[node].start) - me.range_table[node].error;
      int p2 = position(mf) + (me.range_table[node].end - in + 1) + me.range_table[node].error;

      // adjust start and end if they point over the edges of host(finder)
      p1 = (p1 < 0 ? 0 : p1);
      p2 = (p2 > static_cast<int>(length(host(finder))) ? length(host(finder)) : p2);
      THostSegment i(infix(host(mf),p1,p2));
      THSFinder f(i);
      cand = find(f,me.range_table[node].verifier);
      node = node >> 1;
    }
    // if we verfied till here .. verify the complete pattern
    if(cand){
      // we found an occurence
      node = _getRoot(me); // use root 
      int p1 = position(mf) - (in - me.range_table[node].start) - me.range_table[node].error;
      int p2 = position(mf) + (me.range_table[node].end - in + 1) + me.range_table[node].error;

      // adjust start and end if they point over the edges of host(finder)
      p1 = (p1 < 0 ? 0 : p1);
      p2 = (p2 > static_cast<int>(length(host(finder))) ? length(host(finder)) : p2);
      THostSegment i(infix(host(mf),p1,p2));
      THSFinder f(i);
      while(find(f,me.range_table[node].verifier))
      {
	unsigned nP = p1 + position(f);
	if(nP > startPos)
	{
	  // compute new position
	  unsigned offset = nP - position(finder);
	  finder += offset;
	  me.lastFPos = position(mf);
	  me.lastFNdl = position(me.multiPattern);
	  me.findNext = true;
	  return true;
	}      
      }
    }
  }
  // nothing more to find -> set finder to end position
  unsigned t = length(host(finder))- position(finder);
  finder += t;

  return false;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_..

