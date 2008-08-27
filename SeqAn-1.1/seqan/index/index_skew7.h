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
  $Id: index_skew7.h,v 1.2 2008/08/27 00:27:57 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_SKEW7_H
#define SEQAN_HEADER_INDEX_SKEW7_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct Skew7 {};


    //////////////////////////////////////////////////////////////////////////////
    // external Skew7 algorithm
    //////////////////////////////////////////////////////////////////////////////

	template <typename T>
	struct _SkewDC<7, T> {
		static const unsigned VALUE[];
	};

	template <typename T>
	const unsigned _SkewDC<7, T>::VALUE[] = { 3,   1, 2, 4 };


	// *** COMPARATORS & MAPS ***
        
    template <typename InType, typename Result = int>
    struct skew7_ncomp : public ::std::binary_function<InType,InType,Result> {
        inline Result operator()(const InType &a, const InType &b) const
        {
			typedef typename InType::T1 SizeType;
            typedef typename InType::T2 Septet;
            const typename Septet::T *sa = a.i2.i;
            const typename Septet::T *sb = b.i2.i;

            SizeType n = Septet::size;
            if (a.i1 < n) n = a.i1;
            if (b.i1 < n) n = b.i1;
            for(SizeType i = 0; i < n; i++, ++sa, ++sb) {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }
            if (n < Septet::size) {
                return (a.i1 < b.i1)? -1 : 1;
            } else
                return 0;
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, typename Result>
    struct skew7_ncomp< Pair<T1, Tuple<TValue, 7, Compressed>, Compressed >, Result > :
        public ::std::binary_function<
            Pair<T1, Tuple<TValue, 7, Compressed>, Compressed >,
            Pair<T1, Tuple<TValue, 7, Compressed>, Compressed >,
            Result> {       
        inline Result operator()(
            const Pair<T1, Tuple<TValue, 7, Compressed>, Compressed > &a,
            const Pair<T1, Tuple<TValue, 7, Compressed>, Compressed > &b) const
        {
            if (a.i2 < b.i2) return -1;
            if (a.i2 > b.i2) return 1;
            if (a.i1 < 7 || b.i1 < 7) 
                return (a.i1 < b.i1)? -1 : 1;
            return 0;
        }
    };

    template <typename InType, typename Result = typename InType::T1>
    struct skew7_nmap_linear : public ::std::unary_function<InType,Result> {
        Result BN4, BN;
        skew7_nmap_linear(Result _BN):BN4(_BN+1),BN(_BN) { }
        inline Result operator()(const InType& x) const
		{ Result i = x.i1; return (i%7 == 4)? BN4-(i-(i/7)*4): BN-(i-(i/7)*4); }
    };

    template <typename InType, typename Result = typename InType::T1>
    struct skew7_nmap_sliced : public ::std::unary_function<InType,Result> {
        Result off[5];
        skew7_nmap_sliced(Result _BN)
        { 
			off[0] = 0;
			off[1] = _BN - 1; 
			off[2] = (2*_BN)/3 - 1; 
			off[3] = 0;
			off[4] = _BN/3 - 1; 
		}
        inline Result operator()(const InType& x) const
        { return off[x.i1 % 7] - x.i1/7; }
    };


    template <typename InType, typename Result = InType>
    struct skew7_unslicer_func : public ::std::unary_function<InType,Result> {
        Result o1, o2, o4, n4, n24;
        skew7_unslicer_func(Result N):
            o1(N - (N + 6) % 7),
            o2(N - (N + 5) % 7),
            o4(N - (N + 3) % 7),
            n4((N + 3) / 7),
            n24((N + 5) / 7 + n4) { }
        
        inline Result operator()(const InType& x) const
        { return (x < n4)  ? o4 -  x        * 7:
                 (x < n24) ? o2 - (x - n4)  * 7: 
                             o1 - (x - n24) * 7; }
    };

    template <typename InType, typename Result = typename InType::T2::T>
    struct skew7_nmap_extended : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const
        { return x.i2[0]; }
    };

    template <typename InType, const int EXT_LENGTH, typename Result = int>
    struct skew7_extend_comp : public ::std::binary_function<InType,InType,Result> {
        inline Result operator()(const InType &a, const InType &b) const
        {
            for(unsigned int i = 0; i < EXT_LENGTH; i++) {
                if (a.i3[i] == b.i3[i]) continue;
                return (a.i3[i] <  b.i3[i])? -1 : 1;
            }
            return (a.i2[0] < b.i2[0])? -1 : 1;
        }
    };

    // optimized for bitvectors
    template <typename T1, typename T2, typename T, const int _size, const int EXT_LENGTH, typename Result>
    struct skew7_extend_comp< Triple<T1,T2,Tuple<T,_size,Compressed>, Compressed>, EXT_LENGTH, Result> :
        public ::std::binary_function<
            Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed>,
            Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed>,
            Result> 
    {
        inline Result operator()(
            const Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> &a,
            const Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> &b) const
        {
            if (a.i3 < b.i3) return -1;
            if (a.i3 > b.i3) return 1;
            return (a.i2[0] < b.i2[0])? -1 : 1;
        }
    };

	template < typename TInput >
    struct Value< Pipe< TInput, Skew7 > > {
        typedef typename Size<TInput>::Type Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // skew7 class
    template < typename TInput >
    struct Pipe< TInput, Skew7 >
    {
  
        // *** SPECIALIZATION ***

        // use compression if lessorequal 16 different values per char
        typedef typename IF< 
            (BitsPerValue<_TypeOf(TInput)>::VALUE > 0) && 
            (BitsPerValue<_TypeOf(TInput)>::VALUE <= 4), 
            Compressed, 
            void>::Type compress;
//        typedef void compress;

        // use skew7 for recursion (more I/O-efficient)
        typedef Skew7 recurseSpec;

        // step 1
		typedef Pipe< TInput, Sampler<7, compress> >  TSamplerDC7;          
                                        typedef skew7_ncomp<_TypeOf(TSamplerDC7)> ncomp_t;
        typedef Pool< _TypeOf(TSamplerDC7), SorterSpec< SorterConfigSize<ncomp_t, _TSizeOf(TSamplerDC7) > > > TSortTuples;
		typedef Pipe< TSortTuples, Namer<ncomp_t> > TNamer;
                                        typedef skew7_nmap_sliced<_TypeOf(TNamer)> nmap_sliced_t;
                                        typedef skew7_nmap_linear<_TypeOf(TNamer)> nmap_linear_t;
        typedef Pool< _TypeOf(TNamer), MapperSpec< MapperConfigSize< nmap_sliced_t, _TSizeOf(TNamer) > > > TNames_Sliced;

        // unique names - shortcut
        typedef Pool< _TypeOf(TNames_Sliced), MapperSpec< MapperConfigSize< nmap_linear_t, _TSizeOf(TNames_Sliced) > > > TNames_Linear_Unique;

        // non-unique names
        typedef Pipe< TNames_Sliced, Filter< filterI2<_TypeOf(TNames_Sliced)> > > TFilter;

			// recursion
			typedef Pipe< TFilter, recurseSpec > TRecurse;
										typedef skew7_unslicer_func<_TypeOf(TRecurse)> unslicer_func_t;
			typedef Pipe< TRecurse, Filter<unslicer_func_t> > TUnslicer;
			typedef Pipe< TUnslicer, Counter > TRenamer;

			// no recursion inMemory shortcut
			typedef Pipe< TFilter, LarssonSadakane > TInMem;
			typedef Pipe< TInMem, Filter<unslicer_func_t> > TUnslicerInMem;
			typedef Pipe< TUnslicerInMem, Counter > TRenamerInMem;

        typedef Pool< _TypeOf(TRenamer), MapperSpec< MapperConfigSize< nmap_linear_t, _TSizeOf(TRenamer) > > > TNames_Linear;
        
        // step 2
        typedef Pipe< Bundle2< TInput, TNames_Linear >, Extender7<compress> > TExtender;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out0),3> extend0_comp_t;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out6),2> extend6_comp_t;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out5),1> extend5_comp_t;
                                        typedef skew7_extend_comp<_TypeOf(typename TExtender::Out3),1> extend3_comp_t;
        typedef Pool< _TypeOf(typename TExtender::Out0), SorterSpec< SorterConfigSize< extend0_comp_t, _TSizeOf(typename TExtender::Out0) > > > TSorterS0;
        typedef Pool< _TypeOf(typename TExtender::Out6), SorterSpec< SorterConfigSize< extend6_comp_t, _TSizeOf(typename TExtender::Out6) > > > TSorterS6;
        typedef Pool< _TypeOf(typename TExtender::Out5), SorterSpec< SorterConfigSize< extend5_comp_t, _TSizeOf(typename TExtender::Out5) > > > TSorterS5;
        typedef Pool< _TypeOf(typename TExtender::Out3), SorterSpec< SorterConfigSize< extend3_comp_t, _TSizeOf(typename TExtender::Out3) > > > TSorterS3;

        // step 3
                                        typedef skew7_nmap_extended<_TypeOf(typename TExtender::Out124)> nmap_extended_t;
		typedef Pool< _TypeOf(typename TExtender::Out124), MapperSpec< MapperConfigSize< nmap_extended_t, _TSizeOf(typename TExtender::Out124) > > > TSorterS124;
        typedef Pipe< Bundle5< TSorterS0, TSorterS3, TSorterS5, TSorterS6, TSorterS124 >, Merger7 > TMerger;

        TSorterS0   sortedS0;
        TSorterS3   sortedS3;
        TSorterS5   sortedS5;
        TSorterS6   sortedS6;
        TSorterS124 sortedS124;
        TMerger     in;
            
        Pipe():
			in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124)) {}

		Pipe(TInput& _textIn):
			in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124))
		{
			process(_textIn);
		}
        
	    template < typename _TInput >
        bool process(_TInput &textIn) {

            SEQAN_PROADD(SEQAN_PRODEPTH, 1);
            SEQAN_PROMARK("Rekursionsabstieg");
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "enter level " << SEQAN_PROVAL(SEQAN_PRODEPTH) << " compression: ";
				::std::cerr << TYPECMP<compress, Compressed>::VALUE << " " << BitsPerValue<_TypeOf(TInput)>::VALUE << ::std::endl;
            #endif
            {


            // *** INSTANTIATION ***

            // step 1
            TSamplerDC7                 sampler(textIn);
            TSortTuples                 sorter;
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "  sort names (" << length(sampler)<< ")" << ::std::endl;
            #endif
            sorter << sampler;
            SEQAN_PROMARK("Sorter (2) - 7-lets sortieren");

            TNamer                      namer(sorter);
            nmap_sliced_t               map_sliced(length(namer));
            nmap_linear_t               map_linear(length(namer));
            TNames_Sliced               names_sliced(map_sliced);
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "  slice names" << ::std::endl;
            #endif
            names_sliced << namer;

            if (namer.unique() || empty(names_sliced)) {
                // unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - s124 konstruieren");
                TNames_Linear_Unique        names_linear(map_linear);

                #ifdef SEQAN_DEBUG_INDEX
                    ::std::cerr << "  make names linear" << ::std::endl;
                #endif
                names_linear << names_sliced;
				clear(names_sliced);
                SEQAN_PROMARK("Mapper (10) - ISA124 konstruieren");

                // step 2
                skew7_extend(textIn, names_linear, sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);

            } else {
                // non-unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - s124 konstruieren");

                TFilter                     filter(names_sliced);
				TNames_Linear               names_linear(map_linear);

				if (length(filter) > 128*1024*1024) 
				{
					// recursion
					TRecurse                    recurse(filter);

					#ifdef SEQAN_TEST_SKEW7
					{
						String<typename Value<TFilter>::Type, Alloc<> > _text;
						_text << filter;
						SEQAN_DO(isSuffixArray(recurse, _text));
					}
					#endif

					clear(filter);
					unslicer_func_t             func(length(textIn));
					TUnslicer                   unslicer(recurse, func);
					TRenamer                    renamer(unslicer);

					#ifdef SEQAN_DEBUG_INDEX
						::std::cerr << "  rename names" << ::std::endl;
					#endif

					names_linear << renamer;
				}
				else
				{
					TInMem						inMem(filter);

					clear(filter);
					unslicer_func_t				func(length(textIn));
					TUnslicerInMem              unslicer(inMem, func);
					TRenamerInMem               renamer(unslicer);

					#ifdef SEQAN_DEBUG_INDEX
						::std::cerr << "  rename names" << ::std::endl;
					#endif

					names_linear << renamer;
				} 

                SEQAN_PROMARK("Mapper (10) - ISA124 konstruieren");
               
                // step 2
                #ifdef SEQAN_DEBUG_INDEX
                    ::std::cerr << "  prepare merge" << ::std::endl;
                #endif
                skew7_extend(textIn, names_linear, sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);
                SEQAN_PROMARK("Mapper (12), Sorter (13-16) - SA124, SA3, SA5, SA6, SA0 verschmelzen");
            }
            
            // step 3
            // ... is done on-demand by merger
            }
            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "left level " << SEQAN_PROVAL(SEQAN_PRODEPTH) << ::std::endl;
            #endif
            SEQAN_PROMARK("Rekursionsaufstieg");
            SEQAN_PROSUB(SEQAN_PRODEPTH, 1);

            return true;
        }

        inline typename Value<Pipe>::Type const operator*() {
            return *in;
        }
        
        inline Pipe& operator++() {
            ++in;
            return *this;
        }        
    };

    // not sure which interface is more intuitive, we support both
    // you can call "skew << pipe" or "skew_t skew(pipe); skew.process()"
    // for the first we would need no _in member
	template < typename TInput, typename TObject >
    inline bool operator<<(Pipe< TInput, Skew7 > &me, TObject &textIn) {
        return me.process(textIn);
    }

	template < 
		typename TSA, 
		typename TValue, 
		typename TConfig >
	inline void createSuffixArray(
		TSA &SA,
		String< TValue, External<TConfig> > &s,
		Skew7 const &spec,
		unsigned K,
        unsigned maxdepth)
	{
        createSuffixArrayExt(SA, s, spec);
	}


    //////////////////////////////////////////////////////////////////////////////
    // internal Skew7 algorithm
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
	// typedefs and helpers

	// compares n characters and in case of equality the names a2 and b2 (no clipping)
    template <typename TTextIter, typename TSize> inline
    bool _leqSkew7(TTextIter a1, TSize a2,   TTextIter b1, TSize b2,   TSize n)
    { // lexic. order for n-tupels
        for (; n != 0; --n, ++a1, ++b1) {
            if (lexLess(*a1, *b1)) return true;
            if (lexLess(*b1, *a1)) return false;
        }
        return (a2 <= b2);
    }

	// compares at most the last n characters (a) with b (clipping)
    template <typename TTextIter, typename TSize> inline
    bool _leqSkew7(TTextIter a,   TTextIter b,   TSize n)
    { // lexic. order for n-tupels
        for (; n != 0; --n, ++a, ++b) {
            if (lexLess(*a, *b)) return true;
            if (lexLess(*b, *a)) return false;
        }
        return true;	// a is shorter than b
    }

	// compares two suffixes of residue classes a and b
    template <typename TTextIter, typename TSize, typename TString> inline
    bool _leqSkew7(unsigned a, unsigned b,   TTextIter spos[], const TSize tpos[], const bool islast[], const TString &s124, const long adjust[7][7])
    {
        TTextIter sa = spos[a];
        TTextIter sb = spos[b];
		TSize shft = _SkewShift<7>::VALUE[a][b];
        if (sa > sb) {
            if ((a != 0) && (a < shft) && islast[a]) // do we need to clip?
                return _leqSkew7 (sa,   sb,   a);
        } else {
            if ((b != 0) && (b < shft) && islast[b]) // do we need to clip?
                return !_leqSkew7 (sb,   sa,   b);
        }
        return _leqSkew7 (sa, s124[tpos[a] + adjust[a][b]],   sb, s124[tpos[b] + adjust[b][a]],   shft);
    }


	//////////////////////////////////////////////////////////////////////////////
	// Alternative Skew Implementation
	// find the suffix array SA of s[0..n-1] in {0..K}^n
	//
	// the following algorithm divides the suffixes in seven residue classes
	// that results in a more space and time efficient division
	// 
	// difference cover is {3,5,6} and corresponds to {1,2,4}'
	//
	// * no trailing 0's required
	// * no dummy triples in special cases

    template < typename TSA,
               typename TText >
    void createSuffixArray(
		TSA &SA,
		TText &s,
		Skew7 const &,
		unsigned K,
        unsigned maxdepth,
		unsigned depth)
    {
		typedef typename Value<TSA>::Type TSize;
		typedef typename Value<TText>::Type TValue;
		typedef typename Iterator<TSA, Standard>::Type TSAIter;
		typedef typename Iterator<TText const, Standard>::Type TValueIter;

		SEQAN_ASSERT(IsContiguous<TText>::VALUE);
		SEQAN_ASSERT(IsContiguous<TSA>::VALUE);

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cerr << "WARNING: TSize is more than 32 bit long (Skew3). This is probably not what you want." << ::std::endl;
        #endif

		TSize n = length(s);
        if (n < 1) return;

		TSize _n[7];
		TSize _o[7];
	    
		_n[0] = n/7;
		_o[0] = n%7;
		TSize j = n + 6;
		for(int i = 1; i < 7; ++i, --j) {
			_n[i] = j/7;
			_o[i] = j%7;
		}

		TSize _n24  = _n[2]+_n[4];
		TSize _n124 = _n[1]+_n24;

        SEQAN_PROSET(SEQAN_PRODEPTH, depth);
        SEQAN_PROSET(SEQAN_PROEXTRA1, K);
        SEQAN_PROMARK("Rekursionsabstieg");
        #ifdef SEQAN_DEBUG_INDEX
			::std::cerr << "enter level " << depth << " (" << n << ")" << ::std::endl;
        #endif

        String<TSize, Alloc<> > s124;
        resize(s124, _n124, Exact());
		// we use SA[n-n124..n-1] as a temporary buffer instead of allocating one
		typename Suffix<TSA>::Type SA124 = suffix(SA, n - _n124);


		// generate positions of mod 3, mod 5 and mod 6 suffixes
		{
			TSize j = 0;
			if (_n[2] > _n[4]) s124[j++] = _o[2];
			if (_n[1] > _n[4]) s124[j++] = _o[1];

			for (TSize i=_o[4];  i < n;  i+=7) {
				s124[j++] = i;
				s124[j++] = i + 2;
				s124[j++] = i + 3;
			}
		}


		// lsb radix sort the mod 3, mod 5 and mod 6 7-tupels
		{
			String<TSize, Alloc<> > cnt;
			resize(cnt, K, Exact());	// counter array

			radixPass(SA124, s124,  s, cnt, K, 6);
			radixPass(s124,  SA124, s, cnt, K, 5);
			radixPass(SA124, s124,  s, cnt, K, 4);
			radixPass(s124,  SA124, s, cnt, K, 3);
			radixPass(SA124, s124,  s, cnt, K, 2);
			radixPass(s124,  SA124, s, cnt, K, 1);
			radixPass(SA124, s124,  s, cnt, K);
		}
        SEQAN_PROMARK("7-lets sortiert");

		// find lexicographic names of 7-tupel
		TSize name = 0;
		{
			TSize ofs[7] = {0, _n24, _n[4], 0, 0, 0, 0};
			bool differ = true;
			//TValue c0 = TValue(), c1 = TValue(), c2 = TValue(), c3 = TValue(), c4 = TValue(), c5 = TValue(), c6 = TValue();
			TValue c0 = 0, c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0;
			for (TSize i = 0, clip = n - 6, l;  i < _n124;  i++) {
				if ((l = SA124[i]) < clip) {
					if (differ || s[l] != c0 || s[l+1] != c1 || s[l+2] != c2 || s[l+3] != c3 ||
												s[l+4] != c4 || s[l+5] != c5 || s[l+6] != c6) {
						name++;  c0 = s[l];  c1 = s[l+1];  c2 = s[l+2];  c3 = s[l+3];
											 c4 = s[l+4];  c5 = s[l+5];  c6 = s[l+6];
						differ = false;
					}
				} else {
					name++;
					differ = true;  // the last 6 7-tupels always differ from the rest
				}
				s124[(l/7) + ofs[(n-l) % 7]] = name - 1;   // select a third
			}
		}
        SEQAN_PROMARK("s12 konstruiert");

		// recurse if names are not yet unique
		if (name < _n124) {
            if (depth != maxdepth)
            {
			    createSuffixArray(SA124, s124, Skew7(), name, maxdepth, depth + 1);
			    #ifdef SEQAN_TEST_SKEW7
				    SEQAN_DO(isSuffixArray(SA124, s124));
			    #endif
            }
			// store unique names in s124 using the suffix array
			for (TSize i = 0;  i < _n124;  i++) s124[SA124[i]] = i;
		} else // generate the suffix array of s124 directly
			for (TSize i = 0;  i < _n124;  i++) SA124[s124[i]] = i;


		// use SA[0...n3-1] and SA[n3...n3+n5-1] as a temporary buffers instead of allocating some
		// and allocate SA0, SA3, SA5 and SA6

		{
			typename Infix<TSA>::Type s3 = infix(SA, 0, _n[3]), s5 = infix(SA, _n[3], _n[3] + _n[5]);
			String<TSize, Alloc<> > SA0, SA3, SA5, SA6;
			resize(SA0, _n[0], Exact());
			resize(SA3, _n[3], Exact());
			resize(SA5, _n[5], Exact());
			resize(SA6, _n[6], Exact());

			// stably sort the mod 5 and mod 3 suffixes from SA124 by their first character
			{
				for (TSize i=0, j3=0, j5=0, l;  i < _n124;  i++) {
					l = SA124[i];
					if (l < _n[4]) {
						if ((l = _o[4] + (7 * l)) > 0)
							s5[j5++] = l - 1;
					} else if (l < _n24) {
						if ((l = _o[2] + (7 * (l - _n[4]))) > 0)
							s3[j3++] = l - 1;
					}
				}

				{
					String<TSize, Alloc<> > cnt;
					resize(cnt, K, Exact());	// counter array

					radixPass(SA3, s3, s, cnt, K);
                    SEQAN_PROMARK("SA3 konstruiert");

					radixPass(SA5, s5, s, cnt, K);
                    SEQAN_PROMARK("SA5 konstruiert");
	    
					// stably sort the mod 6 suffixes from SA5 by their first character

					if (_n[5] == _n[6]) radixExtend    (SA6, SA5, s, cnt, K);
					else                radixExtendClip(SA6, SA5, s, cnt, K);
                    SEQAN_PROMARK("SA6 konstruiert");

					// stably sort the mod 0 suffixes from SA6 by their first character

					if (_n[6] == _n[0]) radixExtend    (SA0, SA6, s, cnt, K);
					else                radixExtendClip(SA0, SA6, s, cnt, K);	    
                    SEQAN_PROMARK("SA0 konstruiert");
				}
			}

			// MULTIWAY MERGE all SA_ streams
			{
				// a helper matrix to lex-name-compare every combination of suffixes
				long adjust[7][7] = // we use long instead of TSize because we need a sign
					//      0               1              2             3             4              5               6
				   {{0             , _n124-_n[0]   , _n24-_n[0]  , _n124-_n[0] , _n[4]-_n[0] , _n[4]-_n[0]   , _n24-_n[0]    },  // 0
					{1-_n[1]       , 0             , 0           , 1-_n[1]     , 0           , 1-_n[1]-_n[2] , 1-_n[1]-_n[2] },  // 1*
					{1-_n[2]       , 0             , 0           , _n[1]       , 0           , _n[1]         , 1-_n[2]       },  // 2*
					{1+_n[4]-_n[3] , 1+_n[4]-_n[3] , _n24-_n[3]  , 0           , _n124-_n[3] , _n24-_n[3]    , _n124-_n[3]   },  // 3
					{_n[1]+_n[2]   , 0             , 0           ,_n[2]        , 0           , _n[1]+_n[2]   , _n[2]         },  // 4*
					{_n24-_n[5]    , _n124-_n[5]   , _n[4]-_n[5] , _n[4]-_n[5] , _n24-_n[5]  , 0             , _n124-_n[5]   },  // 5
					{_n124-_n[6]   , _n24-_n[6]    , _n124-_n[6] , _n[4]-_n[6] , _n[4]-_n[6] , _n24-_n[6]    , 0             }}; // 6

				TSAIter pos[7] = {begin(SA0)      , begin(SA124)      , begin(SA124)      , 
								  begin(SA3)      , begin(SA124)      , begin(SA5)        , begin(SA6)      };
				TSAIter max[7] = {begin(SA0)+_n[0], begin(SA124)+_n124, begin(SA124)+_n124,
								  begin(SA3)+_n[3], begin(SA124)+_n124, begin(SA5)+_n[5]  , begin(SA6)+_n[6]};
				TValueIter spos[7];
				TSize tpos[7];
				bool islast[7];

				int a, b, rank[5];
				int fill = 0;
				TSize k = 0;

				#define SEQAN_GET_ISKEW7(ii) (ii < _n[4] ? (ii * 7) + _o[4] : (ii < _n24 ? ((ii - _n[4]) * 7) + _o[2] : ((ii - _n24) * 7) + _o[1]))
				#define SEQAN_GET_ASKEW7(ii) (ii < _n[4] ? 4 : (ii < _n24 ? 2 : 1))

				// fill the stream ranking list
				for(int i = 0; i < 7; ++i) {
					if (!_n[i]) continue;
					if (i == 2 || i == 4) continue; // insert only the least suffix of SA124

					if (i == 1) {

						TSize ii  = *(pos[1]);
						a         = SEQAN_GET_ASKEW7(ii);
						TSize j	  = SEQAN_GET_ISKEW7(ii);

						tpos[a]   = ii;
						spos[a]   = begin(s) + j;
						islast[a] = (j + 7 >= n);

					} else {

						a = i;
						TSize j   = *(pos[a]);

						tpos[a]   = j / 7;
						spos[a]   = begin(s) + j;
						islast[a] = (j + 7 >= n);

					}

					// get the rank of stream a's suffix
					int j;
					for(j = 0;  j < fill;  ++j) {
						b = rank[j];
						if (_leqSkew7 (a,  b,   spos, tpos, islast, s124, adjust)) break;
					}

					// insert the suffix
					for(int i = fill; i > j; --i)
						rank[i] = rank[i-1];
					rank[j] = a;
					fill++;
				}

				// main merge loop
				// in order to find the least suffix in every step we use a stream ranking list
				// so we only need to keep up the ordering and thus rank[0] is always the least suffix
				while (fill > 1) {

					// add the least suffix to SA and get the next of the corresponding stream
					a = rank[0];
					SA[k++] = spos[a] - begin(s);
					if (a == 1 || a == 2 || a == 4)
						pos[4] = pos[2] = ++pos[1];
					else
						++pos[a];

					if (pos[a] < max[a]) {

						// set corresponding spos, tpos, islast values and adapt a if necessary
						if (a == 1 || a == 2 || a == 4) {

							TSize ii  = *(pos[1]);
							a		  = SEQAN_GET_ASKEW7(ii);
							TSize j   = SEQAN_GET_ISKEW7(ii);

							tpos[a]   = ii;
							spos[a]   = begin(s) + j;
							islast[a] = (j + 7 >= n);

						} else {

							TSize j   = *(pos[a]);

							tpos[a]   = j / 7;
							spos[a]   = begin(s) + j;
							islast[a] = (j + 7 >= n);

						}

						// get the rank of stream a's suffix

						// linear search
						int right;
						for(right = 1;  right < fill;  right++) {
							b = rank[right];
							if (_leqSkew7 (a,  b,   spos, tpos, islast, s124, adjust)) break;
						}
/*
						// binary search
						int left = 0;
						int right = fill;
						while (left + 1 != right) {
							int middle = (left + right) >> 2;
							if (leq<TValue, TSize> (a,  rank[middle],   spos, tpos, islast, s124, adjust))
								right = middle;
							else
								left = middle;
						}*/

						// remove the least suffix ...
						for(int i = 1; i < right; ++i)
							rank[i-1] = rank[i];

						// ... and insert the new one
						rank[right-1] = a;

					} else {
						// only remove the least suffix
						fill--;
						for(int i = 0; i < fill; ++i)
							rank[i] = rank[i+1];
					}
				}

				// only one (or less) stream left to fill SA with
				a = rank[0];
				if (a == 1 || a == 2 || a == 4)
					for (;  k < n;  ++k) { TSize ii = *(pos[1]++); SA[k] = SEQAN_GET_ISKEW7(ii); }
				else
					for (;  k < n;  ++k) SA[k] = *(pos[a]++);
			}	    
		}
        SEQAN_PROMARK("SA124, SA3, SA5, SA6, SA0 verschmolzen");

        #ifdef SEQAN_DEBUG_INDEX
            ::std::cerr << "left level " << depth << ::std::endl;
        #endif
        SEQAN_PROMARK("Rekursionsaufstieg");
        SEQAN_PROSUB(SEQAN_PRODEPTH, 1);
	}

    template < typename TSA,
               typename TText >
    inline void createSuffixArray(
		TSA &SA,
		TText &s,
		Skew7 const &alg,
		unsigned K,
        unsigned maxdepth)
	{
		createSuffixArray(SA, s, alg, K, maxdepth, 1);
	}

    // creates suffix array sorted by the first maxLCP chars of suffixes
    template < typename TSA,
               typename TText,
               typename TSize >
    inline void createSuffixArrayPart(
		TSA &SA,
		TText &s,
		Skew7 const &_dummy,
        TSize maxLCP,
        unsigned K = ValueSize< typename Value<TText>::Type >::VALUE)
    {
        unsigned depth = 0;
        for(TSize i = 1; i < maxLCP; i*=7) ++depth;
        createSuffixArray(SA, s, _dummy, K, depth);
    }

//}

}

#endif
