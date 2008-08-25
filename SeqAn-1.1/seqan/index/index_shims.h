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
  $Id: index_shims.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_SHIMS_H
#define SEQAN_HEADER_INDEX_SHIMS_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// Suffix Array creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	// build suffix array with an external pipeling algorithm (skew3, skew7, ...)
	template < 
		typename TSA, 
		typename TObject, 
		typename TAlgSpec >
	void createSuffixArrayExt(
		TSA &suffixArray,
		TObject const &text,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename _MakeUnsigned< typename Value<TObject>::Type >::Type TUValue;

        // specialization
		typedef Pipe< TObject, Source<> >				src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
		typedef Pipe< unsigner_t, TAlgSpec >	        creator_t;

		// instantiation and processing
		src_t		src(text);
        unsigner_t  unsigner(src);
		creator_t	creator(unsigner);

		suffixArray << creator;
		#ifdef SEQAN_TEST_INDEX
			isSuffixArray(suffixArray, text);
		#endif
	}


	// build suffix array (external) for mutliple sequences
	template < 
		typename TSA, 
		typename TString, 
		typename TSpec,
		typename TAlgSpec >
	void createSuffixArrayExt(
		TSA &suffixArray,
		StringSet<TString, TSpec> const &stringSet,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type			TConcat;
        typedef typename _MakeUnsigned< typename Value<TConcat>::Type >::Type	TUValue;
		typedef Multi<
			TAlgSpec, 
			typename Value<TSA>::Type, 
			typename StringSetLimits<StringSet<TString, TSpec> >::Type >		MultiConstrSpec;

        // specialization
		typedef Pipe< TConcat, Source<> >				src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
		typedef Pipe< unsigner_t, MultiConstrSpec >	    creator_t;

		// instantiation and processing
		src_t		src(concat(stringSet));
        unsigner_t  unsigner(src);
		creator_t	creator(unsigner, stringSetLimits(stringSet));

		suffixArray << creator;
		#ifdef SEQAN_TEST_INDEX
			//isSuffixArray(suffixArray, stringSet);
		#endif
	}


/**
.Function.createSuffixArray:
..summary:Creates a suffix array from a given text.
..cat:Index
..signature:createSuffixArray(suffixArray, text[, algo_tag])
..param.suffixArray:A reference to the resulting suffix array.
..param.text:A given text.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $suffixArray$ must be at least $length(text)$ before calling this function.
*/

    template < typename TSA,
               typename TText,
			   typename TAlgSpec >
    inline void createSuffixArray(
		TSA &SA,
		TText const &s,
		TAlgSpec const &alg)
	{
	SEQAN_CHECKPOINT
		// -> call internal memory algorithm with an extended interface (+ alphabet size, max_depth)
		if (BitsPerValue< typename Value<TText>::Type >::VALUE > 16)
			createSuffixArray(SA, s, alg, length(s), 0);
		else
			createSuffixArray(SA, s, alg, ValueSize< typename Value<TText>::Type >::VALUE, 0);
	}

	template < 
		typename TSA, 
		typename TValue, 
		typename TConfig,
		typename TAlgSpec >
	inline void createSuffixArray(
		TSA &SA,
		String< TValue, External<TConfig> > const &s,
		TAlgSpec const alg)
	{
	SEQAN_CHECKPOINT
		// -> explicitly create SA using external memory
        createSuffixArrayExt(SA, s, alg);
	}

	template < 
		typename TSA, 
		typename TValue, 
		typename TSpec,
		typename TSSetSpec,
		typename TAlgSpec >
	inline void createSuffixArray(
		TSA &SA,
		StringSet< String<TValue, TSpec>, TSSetSpec > const &s,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
        createSuffixArrayExt(SA, s, Skew7());
	}


//____________________________________________________________________________



	//////////////////////////////////////////////////////////////////////////////
	// LCP Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	template < 
        typename TLCPTable,
		typename TObject, 
        typename TSA,
		typename TAlgSpec >
	void createLCPTableExt(
		TLCPTable &LCP,
		TObject const &text,
		TSA const &suffixArray,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
		// specialization
		typedef Pipe< TObject, Source<> >							srcText_t;
		typedef Pipe< TSA, Source<> >   							srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, TAlgSpec >		creator_t;

		// instantiation and processing
		srcText_t	srcText(text);
		srcSA_t		srcSA(suffixArray);
		creator_t	creator(bundle2(srcText, srcSA));

		LCP << creator;
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, suffixArray, text);
		#endif
	}


	// build lcp table (external) for mutliple sequences
	template < 
        typename TLCPTable,
		typename TString,
		typename TSpec,
        typename TSA,
		typename TAlgSpec >
	void createLCPTableExt(
		TLCPTable &LCP,
		StringSet<TString, TSpec> const &stringSet,
		TSA const &suffixArray,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type TConcat;
		typedef Multi<
			TAlgSpec, 
			typename Value<TSA>::Type, 
			typename StringSetLimits<StringSet<TString, TSpec> >::Type > MultiConstrSpec;

		// specialization
		typedef Pipe< TConcat, Source<> >								srcText_t;
		typedef Pipe< TSA, Source<> >   								srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, MultiConstrSpec >	creator_t;

		// instantiation and processing
		srcText_t	srcText(concat(stringSet));
		srcSA_t		srcSA(suffixArray);
		creator_t	creator(bundle2(srcText, srcSA), stringSetLimits(stringSet));

		LCP << creator;
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, suffixArray, text);
		#endif
	}


/**
.Function.createLCPTable:
..summary:Creates a lcp table from a given text and suffix array.
..cat:Index
..signature:createLCPTable(lcp, text, suffixArray[, algo_tag])
..param.lcp:A reference to the resulting lcp table.
..param.text:A given text.
..param.suffixArray:The suffix array of $text$.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $lcp$ must be at least $length(text)$ before calling this function.
*/

	template < 
        typename TLCPTable,
        typename TValue,
		typename TConfig,
		typename TSA >
	void createLCPTable(
		TLCPTable &LCP,
		String< TValue, External<TConfig> > const &s,
		TSA const &SA,
		Kasai const)
	{
	SEQAN_CHECKPOINT
        createLCPTableExt(LCP, s, SA, Kasai());
	}

	template < 
        typename TLCPTable,
        typename TValue,
		typename TConfig,
		typename TSSetSpec,
		typename TSA >
	void createLCPTable(
		TLCPTable &LCP,
		StringSet< String<TValue, External<TConfig> >, TSSetSpec > const &s,
		TSA const &SA,
		Kasai const)
	{
        createLCPTableExt(LCP, s, SA, Kasai());
	}


//____________________________________________________________________________



	//////////////////////////////////////////////////////////////////////////////
	// Enhanced LCP Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

    // build enhanced LCP table with an external pipelining algorithm (ext kasai, ...)
	// and a dynamic programming tree construction alg
	// (in contrast to the LCP table the enhanced LCP table contains the lcp-values 
	// of suffix intervals used in the binary search)
	template < 
		typename TValue, 
        typename TSpec,
		typename TObject, 
        typename TSA,
		typename TAlgSpec >
	void createLCPETableExt(
		String< TValue, TSpec > &LCPE,
		TObject const &text,
		TSA const &suffixArray,
		TAlgSpec const)
	{
	SEQAN_CHECKPOINT
		typedef typename Concatenator<TObject>::Type				TConcat;

		// specialization
		typedef Pipe< TConcat, Source<> >						    srcText_t;
		typedef Pipe< TSA, Source<> >					    	    srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, TAlgSpec >		creator_t;

		// instantiation and processing
		srcText_t	srcText(concat(text));
		srcSA_t		srcSA(suffixArray);
		creator_t	creator(bundle2(srcText, srcSA));

		#ifdef SEQAN_TEST_INDEX
			isLCPTable(creator, suffixArray, text);
		#endif
		createLCPBinTree(LCPE, creator);
	}

    // build enhanced LCP table with an lcp algorithm
	// and a dynamic programming tree construction alg
    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA,
		typename TAlgSpec >
    void createLCPETable(
		String< TValue, TSpec > &LCPE,
		TText const &s,
		TSA const &SA,
		TAlgSpec const alg)
	{
	SEQAN_CHECKPOINT
        //TSA LCP;
        //resize(LCP, length(s), Exact());
		// we use LCPE[n-lcpSize..n-1] as a temporary buffer instead of allocating one
		typename Suffix<String< TValue, TSpec > >::Type LCP = suffix(LCPE, length(LCPE) - length(s));

		createLCPTable(LCP, s, SA, alg);
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, SA, s);
		#endif
        createLCPBinTree(LCPE, LCP);
    }

    template <
        typename TValue,
        typename TConfig,
        typename TText,
        typename TSA,
		typename TAlgSpec >
    void createLCPETable(
		String< TValue, External<TConfig> > &LCPE,
		TText const &s,
		TSA const &SA,
		TAlgSpec const alg)
	{
	SEQAN_CHECKPOINT
        createLCPETableExt(LCPE, s, SA, alg);
    }

    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA>
    inline void createLCPETable(
		String< TValue, TSpec > &LCPE,
		TText &s,
		TSA &SA)
	{
	SEQAN_CHECKPOINT
		CreateLCPETable(LCPE, s, SA, Kasai());
    }


//____________________________________________________________________________



	//////////////////////////////////////////////////////////////////////////////
	// Burrows-Wheeler-Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	template < typename TBWT,
               typename TText,
			   typename TSA >
    void createBWTableExt(
		TBWT &bwt,
		TText const &s,
		TSA const &SA)
	{
	SEQAN_CHECKPOINT
		// specialization
		typedef Pipe< TText, Source<> >						srcText_t;
		typedef Pipe< TSA, Source<> >   					srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, BWT >	creator_t;

		// instantiation and processing
		srcText_t	srcText(s);
		srcSA_t		srcSA(SA);
		creator_t	creator(bundle2(srcText, srcSA));

		bwt << creator;
	}

/**
.Function.createBWTable:
..summary:Creates a Burrows-Wheeler table from a given text and suffix array.
..cat:Index
..signature:createBWTable(bwt, text, suffixArray[, algo_tag])
..param.bwt:A reference to the resulting Burrows-Wheeler table.
..param.text:A given text.
..param.suffixArray:The suffix array of $text$.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $bwt$ must be at least $length(text)$ before calling this function.
*/

	// default
	template < typename TBWT, typename TText, typename TSA, typename _TTextRandom >
    inline void _createBWTableWrapper(TBWT &bwt, TText const &s, TSA const &sa,		_TTextRandom const)
	{
	SEQAN_CHECKPOINT
		createBWTableExt(bwt, concat(s), sa);
	}

	// text supports fast random access
	template < typename TBWT, typename TText, typename TSA >
    inline void _createBWTableWrapper(TBWT &bwt, TText const &s, TSA const &sa,		True const)
	{
	SEQAN_CHECKPOINT
		createBWTableInt(bwt, concat(s), sa);
	}

	template < typename TBWT, typename TText, typename TSA >
    inline void createBWTable(TBWT &bwt, TText const &s, TSA const &sa)
	{
	SEQAN_CHECKPOINT
		_createBWTableWrapper(bwt, s, sa, typename AllowsFastRandomAccess<TText>::Type());
	}


//////////////////////////////////////////////////////////////////////////////

	template <typename TOccValue>
	struct _SAValueLess:
		public ::std::less<TOccValue> {};

	template <typename T1, typename T2, typename TCompression>
	struct _SAValueLess< Pair<T1,T2,TCompression> >:
		public ::std::binary_function< Pair<T1,T2,TCompression>, Pair<T1,T2,TCompression>, bool> 
	{
		inline bool operator()(Pair<T1,T2,TCompression> const &a, Pair<T1,T2,TCompression> const &b) const {
			return	getValueI1(a) < getValueI1(b) ||
					getValueI1(a) == getValueI1(b) && getValueI2(a) < getValueI2(b);
		}
	};

/**
.Function.orderOccurrences:
..summary:Sorts a string of occurrences.
..cat:Index
..signature:orderOccurrences(occstring)
..param.occstring:String of occurrences.
...remarks:Contains suffix array values returned by @Function.getOccurrences@.
..remarks:The occurrences are sorted by increasing positions.
..see:Function.getOccurrences
..see:Metafunction.SAValue
*/
	template <typename TValue, typename TSpec>
	inline void orderOccurrences(String<TValue, TSpec> &occString)
	{
	SEQAN_CHECKPOINT
		::std::sort(begin(occString, Standard()), end(occString, Standard()), _SAValueLess<TValue>());
	}


//////////////////////////////////////////////////////////////////////////////
// fibre creators

/**
.Function.indexCreate:
..summary:Creates a specific @Metafunction.Fibre@.
..cat:Index
..signature:indexCreate(index, fibre_tag[, algo_tag])
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA Index Fibres.ESA_SA@).
..param.algo_tag:A tag that identifies the algorithm which is used to create the fibre.
...default:The result of @Metafunction.DefaultIndexCreator@.
..returns:A $bool$ which is $true$ on a successful creation.
..remarks:$indexCreate$ calls the fibre corresponding $createXXX(..)$ function (e.g. @Function.createSuffixArray@).
*/

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, Fibre_SA, TSpecAlg const alg) {
	SEQAN_CHECKPOINT
		resize(indexSA(index), length(indexRawText(index)), Exact());
		createSuffixArray(indexSA(index), indexText(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, Fibre_LCP, TSpecAlg const alg) {
	SEQAN_CHECKPOINT
		resize(indexLCP(index), length(indexRawText(index)), Exact());
		createLCPTable(indexLCP(index), indexText(index), indexSA(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, Fibre_LCPE, TSpecAlg const alg) {
	SEQAN_CHECKPOINT
	//TODO: separate LCP from LCPE (for now LCPE = LCP + extra)
//		resize(indexLCP(index), length(indexRawText(index)), Exact());
//		createLCPETable(indexLCPE(index), indexRawText(index), indexSA(index), alg);
		return false;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, TSpec> &index, Fibre_BWT, BWT const) {
	SEQAN_CHECKPOINT
		resize(indexBWT(index), length(indexRawText(index)), Exact());
		createBWTable(indexBWT(index), indexText(index), indexRawSA(index));
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, TSpec> &index, Fibre_ChildTab, ChildTab const) {
	SEQAN_CHECKPOINT
		resize(indexChildTab(index), length(indexRawText(index)), Exact());
		createChildTable(indexChildTab(index), indexLCP(index));
		return true;
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexCreate(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
	SEQAN_CHECKPOINT
		return indexCreate(index, fibre, typename DefaultIndexCreator<Index<TText, TSpec>, Tag<TFibre> const>::Type());
	}


//////////////////////////////////////////////////////////////////////////////
// automatic fibre creation

/**
.Function.indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..cat:Index
..signature:indexSupplied(index, fibre_tag)
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA Index Fibres.ESA_SA@).
..returns:A $bool$ which is $true$, iff the fibre is present.
*/

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSupplied(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
	SEQAN_CHECKPOINT
		return !empty(getFibre(index, fibre));
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRequire:
..summary:On-demand creation of a specific @Metafunction.Fibre@.
..cat:Index
..signature:indexRequire(index, fibre_tag)
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA Index Fibres.ESA_SA@).
..returns:A $bool$ which is $true$ on a successful creation.
..remarks:If the fibre already exists (@Function.indexSupplied@ is true) then $indexRequire$ does nothing.
If the fibre doesn't exist then @Function.indexCreate@ is called to create it.
*/

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexRequire(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
	SEQAN_CHECKPOINT
		if (indexSupplied(index, fibre)) return true;				// if the table doesn't exist,
		if (!indexSolveDependencies(index, fibre)) return false;	// fulfill requirements
		return indexCreate(index, fibre);							// and create table
	}


//////////////////////////////////////////////////////////////////////////////
// index cargo interface

	template <typename TText, typename TSpec>
	inline typename Reference< typename Cargo<Index<TText, TSpec> >::Type >::Type
	cargo(Index<TText, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return me.cargo;
	}

	template <typename TText, typename TSpec>
	inline typename Reference< typename Cargo<Index<TText, TSpec> const>::Type >::Type
	cargo(Index<TText, TSpec> const & me)
	{
	SEQAN_CHECKPOINT
		return me.cargo;
	}

//////////////////////////////////////////////////////////////////////////////
// solve dependencies

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, TSpec> &, Tag<TFibre> const) {
	SEQAN_CHECKPOINT
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Fibre_LCP) {
	SEQAN_CHECKPOINT
		return indexRequire(index, Fibre_SA());
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Fibre_LCPE) {
	SEQAN_CHECKPOINT
		return indexRequire(index, Fibre_LCP());
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Fibre_ChildTab) {
	SEQAN_CHECKPOINT
		return indexRequire(index, Fibre_LCP());
	}

	template <typename TText, typename TSpec>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Fibre_BWT) {
	SEQAN_CHECKPOINT
		return indexRequire(index, Fibre_SA());
	}


//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TValue, typename TSpec >
	inline bool open(String<TValue, TSpec> &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		String<TValue, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(string, extString);
		return true;
	}
	template < typename TValue, typename TSpec >
	inline bool open(String<TValue, TSpec> &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return open(string, fileName, OPEN_RDONLY);
	}

	template < typename THost, typename TSpec >
	inline bool open(Segment<THost, TSpec> &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(string, extString);
		return true;
	}
	template < typename THost, typename TSpec >
	inline bool open(Segment<THost, TSpec> &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return open(string, fileName, OPEN_RDONLY);
	}

	template < typename TValue, typename TSpec, typename TSSSpec >
	inline bool open(StringSet<String<TValue, TSpec>, TSSSpec> &multi, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		char id[12]; // 2^32 has 10 decimal digits + 1 (0x00)
		unsigned i = 0;
		clear(multi);
		while (true) {
			sprintf(id, ".%u", i);
			String<char> name;
			name = fileName;	append(name, id);
			{
				String<TValue, External< ExternalConfigLarge<> > > extString;
				if (!open(extString, toCString(name), openMode | OPEN_QUIET)) break;
				resize(multi, i + 1);
				assign(value(multi, i), extString);
			}
			++i;
		}
		return i > 1;
	}
	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool open(StringSet<String<TValue, TSpec>, TSSSpec> &multi, const char *fileName) {
	SEQAN_CHECKPOINT
		return open(multi, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TValue, typename TSpec >
	inline bool save(String<TValue, TSpec> const &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
// 
//		if (length(string) == 0) return true;
		String<TValue, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(extString, string);
		return true;
	}
	template < typename TValue, typename TSpec >
	inline bool save(String<TValue, TSpec> &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

	template < typename THost, typename TSpec >
	inline bool save(Segment<THost, TSpec> const &string, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		if (length(string) == 0) return true;
		String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(extString, string);
		return true;
	}
	template < typename THost, typename TSpec >
	inline bool save(Segment<THost, TSpec> const &string, const char *fileName) {
	SEQAN_CHECKPOINT
		return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool save(StringSet<String<TValue, TSpec>, TSSSpec> const &multi, const char *fileName, int openMode) {
	SEQAN_CHECKPOINT
		if (length(multi) == 0) return true;
		char id[12]; // 2^32 has 10 decimal digits + 2 ('.' and 0x00)
		for(unsigned i = 0; i < length(multi); ++i) {
			sprintf(id, ".%u", i);
			String<char> name;
			name = fileName;	
			append(name, &(id[0]));
			if (!save(getValue(multi, i), toCString(name), openMode))
				return false;
		}
		return true;
	}
	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool save(StringSet<String<TValue, TSpec>, TSSSpec> const &multi, const char *fileName) {
	SEQAN_CHECKPOINT
		return save(multi, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif
