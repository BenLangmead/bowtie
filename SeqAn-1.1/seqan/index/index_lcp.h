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
  $Id: index_lcp.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_LCP_H
#define SEQAN_HEADER_INDEX_LCP_H

namespace SEQAN_NAMESPACE_MAIN
{
	
//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct Kasai {};
	struct KasaiOriginal {};	// original, but more space-consuming algorithm


    //////////////////////////////////////////////////////////////////////////////
    // external LCP algorithm (modified Kasai et al. for pipelining)
    //////////////////////////////////////////////////////////////////////////////


    template < typename TTextInput, typename TSuffixArrayInput >
    struct Value< Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Kasai > > {
        typedef typename Size<TTextInput>::Type Type;
    };

	template <typename InType, typename Result = typename InType::T2::T>
	struct map_inverse : public ::std::unary_function<InType,Result> {
        inline Result operator()(const InType& x) const
        { return x.i2[0]; }
    };

	//////////////////////////////////////////////////////////////////////////////
    // lcp class
    template < typename TTextInput, typename TSuffixArrayInput >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Kasai >
    {
        // *** SPECIALIZATION ***

        typedef Pipe< TSuffixArrayInput, Echoer<2,false> > TEchoer;
                                        typedef map_inverse<_TypeOf(TEchoer)> map_inverse_t;
		                                typedef typename Size<TTextInput>::Type	TSize;
		typedef Pool< _TypeOf(TEchoer), MapperSpec< MapperConfigSize< map_inverse_t, TSize> > > TInverter;
		                                typedef Pair<TSize> TCoreType;
		typedef Pool< TCoreType, MapperSpec< MapperConfigSize< filterI1<TCoreType>, TSize > > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<TCoreType> > > TFilter;

		TTextInput				*textIn;
        TSuffixArrayInput		*suffixArrayIn;
        TLinearMapper           mapper;
		TFilter					in;
        const LcpConfig			conf;
        
        Pipe():
            in(mapper) {}

        Pipe(const LcpConfig &_conf):
            in(mapper),
			conf(_conf) {}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn):
            textIn(&_bundleIn.in1),
			suffixArrayIn(&_bundleIn.in2),
            in(mapper)
		{
			process();
		}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn, LcpConfig const &_conf):
            textIn(&_bundleIn.in1),
			suffixArrayIn(&_bundleIn.in2),
            in(mapper),
			conf(_conf)
		{
			process();
		}
        
        inline void process() {
            process(*textIn, *suffixArrayIn);
        }

		template < typename _TTextInput, typename _TSuffixArrayInput >
        bool process(_TTextInput &textIn, _TSuffixArrayInput &suffixArrayIn) {

            // *** INSTANTIATION ***
			
			TEchoer						echoer(suffixArrayIn);
			TInverter					inverter(echoer);

            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "  invert suffix array" << ::std::endl;
            #endif
			inverter << echoer;
			SEQAN_PROMARK("Suffix-Array invertiert");

			lcp_process(textIn, inverter, mapper);
            return true;
        }

        inline typename Value<Pipe>::Type const operator*() const {
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
	template < typename TInput, typename _TTextInput, typename _TSuffixArrayInput >
    inline bool operator<<(Pipe< TInput, Kasai > &me, Bundle2< _TTextInput, _TSuffixArrayInput > const &bundleIn) {
 	    return me.process(bundleIn.in1, bundleIn.in2);
    }



    //////////////////////////////////////////////////////////////////////////////
    // external LCP algorithm (optimized for multiple sequences)
    //////////////////////////////////////////////////////////////////////////////


    template < typename TTextInput, typename TSuffixArrayInput, typename TPair, typename TLimitsString >
    struct Value< Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Multi<Kasai, TPair, TLimitsString> > > {
        typedef typename Size<TTextInput>::Type Type;
    };

	template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
	struct map_inverse_multi : public ::std::unary_function<InType,Result> {
		TLimitsString const &limits;
		map_inverse_multi(TLimitsString const &_limits) : limits(_limits) {}
        inline Result operator()(const InType& x) const
        {
			return posGlobalize(x.i2[0], limits);
		}
    };

    //////////////////////////////////////////////////////////////////////////////
    // lcp class
    template < typename TTextInput, typename TSuffixArrayInput, typename TPair, typename TLimitsString >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Multi<Kasai, TPair, TLimitsString> >
    {
        // *** SPECIALIZATION ***

        typedef Pipe< TSuffixArrayInput, Echoer<2,false> > TEchoer;
                                        typedef map_inverse_multi<_TypeOf(TEchoer), TLimitsString, _TSizeOf(TEchoer)> map_inverse_t;
		                                typedef typename Size<TTextInput>::Type	TSize;
		typedef Pool< _TypeOf(TEchoer), MapperSpec< MapperConfigSize< map_inverse_t, TSize> > > TInverter;
		                                typedef Pair<TSize> TCoreType;
		typedef Pool< TCoreType, MapperSpec< MapperConfigSize< filterI1<TCoreType>, TSize > > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<TCoreType> > > TFilter;

        TLinearMapper           mapper;
		TFilter					in;
		TLimitsString const		&limits;
        const LcpConfig			conf;
        
        Pipe(TLimitsString const &_limits):
            in(mapper),
			limits(_limits)	{}

        Pipe(TLimitsString const &_limits, LcpConfig const &_conf):
            in(mapper),
			limits(_limits),
			conf(_conf) {}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn, TLimitsString const &_limits):
            in(mapper),
			limits(_limits)
		{
			process(_bundleIn.in1, _bundleIn.in2);
		}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn, TLimitsString const &_limits, LcpConfig const &_conf):
            in(mapper),
			limits(_limits),
			conf(_conf)
		{
			process(_bundleIn.in1, _bundleIn.in2);
		}
        
		template < typename _TTextInput, typename _TSuffixArrayInput >
        bool process(_TTextInput &textIn, _TSuffixArrayInput &suffixArrayIn) {

            // *** INSTANTIATION ***
			
			TEchoer						echoer(suffixArrayIn);
			map_inverse_t				map_inverse(limits);
			TInverter					inverter(echoer, map_inverse);

            #ifdef SEQAN_DEBUG_INDEX
                ::std::cerr << "  invert suffix array" << ::std::endl;
            #endif
			inverter << echoer;
			SEQAN_PROMARK("Suffix-Array invertiert");

			lcp_process_multi(textIn, limits, inverter, mapper);
            return true;
        }

        inline typename Value<Pipe>::Type const operator*() const {
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
	template < typename TInput, typename _TTextInput, typename _TSuffixArrayInput, typename TPair, typename TLimitsString >
    inline bool operator<<(Pipe< TInput, Multi<Kasai, TPair, TLimitsString> > &me, Bundle2< _TTextInput, _TSuffixArrayInput > const &bundleIn) {
 	    return me.process(bundleIn.in1, bundleIn.in2);
    }

	
	
	//////////////////////////////////////////////////////////////////////////////
    // internal Kasai algorithm
    //////////////////////////////////////////////////////////////////////////////


    template < typename TLCPTable,
               typename TText,
               typename TSA >
    void createLCPTable(
		TLCPTable &LCP,
		TText const &s,
		TSA const &SA,
		KasaiOriginal const)
	{
		typedef typename Value<TSA>::Type TSize;

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cerr << "WARNING: TSize size is greater 4 (Kasai)" << ::std::endl;
        #endif

		TSize n = length(s);
        if (n < 2) return;

        #ifdef SEQAN_DEBUG_INDEX
            TSize lcpMax = 0, lcpAvrg = 0, lcpNumer = 0, sigma = 1;	// for lcpMax, lcpMean, |Sigma|
        #endif

        String<TSize, Alloc<> > ISA;
		resize(ISA, n, Exact());

		for(TSize i = 0; i < n; ++i)
			ISA[SA[i]] = i;
		
		SEQAN_PROMARK("Suffix-Array invertiert");

		typename Iterator<TText const>::Type Ibegin = begin(s);
        typename Iterator<TText const>::Type I = Ibegin, J;
        for(TSize i = 0, h = 0, j, isa; i < n; ++i) {
			if ((isa = ISA[i])) {
				J = Ibegin + h + (j = SA[isa - 1]);
                for(TSize hMax = _min(n - i, n - j); h < hMax && *I == *J; ++I, ++J, ++h);
				LCP[isa - 1] = h;
                #ifdef SEQAN_DEBUG_INDEX
                    if ((lcpNumer += h) > n) {
                        lcpNumer -= n;
                        ++lcpAvrg;
                    }
                    if (lcpMax < h) lcpMax = h;
					if (!h) ++sigma;
                #endif
			}
			if (h) --h;
            else ++I;
        }
		LCP[n - 1] = 0;
        #ifdef SEQAN_DEBUG_INDEX
            ::std::cerr << "  n: " << n;
            ::std::cerr << "  lcpMax: " << lcpMax;
            ::std::cerr << "  lcpAvrg: " << (TSize)(lcpAvrg + (lcpNumer + n/2) / n);
            ::std::cerr << "  sigma: " << sigma << ::std::endl;
        #endif
	}

	// HINT:
	// In contrast to the upper functions 
	// createLCPTableInPlace expects the lcp table to be of size n
    template < typename TLCPTable,
               typename TText,
               typename TSA >
    void createLCPTable(
		TLCPTable &LCP,
		TText const &s,
		TSA const &SA,
		Kasai const)
	{
		typedef typename Value<TSA>::Type TSize;

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cerr << "WARNING: TSize size is greater 4 (Kasai)" << ::std::endl;
        #endif

		TSize n = length(s);
        if (n < 2) return;

        #ifdef SEQAN_DEBUG_INDEX
            TSize lcpMax = 0, lcpAvrg = 0, lcpNumer = 0, sigma = 1;	// for lcpMax, lcpMean, |Sigma|
        #endif

		TSize mark = ~(~0u>>1);
		TSize mask =   ~0u>>1;

		for(TSize i = 0; i < n; ++i)
			LCP[SA[i]] = i;
		
		SEQAN_PROMARK("Suffix-Array invertiert");
        #ifdef SEQAN_DEBUG_INDEX
			::std::cerr << "Suffix-Array invertiert" << ::std::endl;
		#endif

		typename Iterator<TText const>::Type Ibegin = begin(s);
        typename Iterator<TText const>::Type I = Ibegin, J;
        for(TSize i = 0, h = 0, j, isa; i < n; ++i) {
			if ((isa = LCP[i] + 1) < n) {
				J = Ibegin + h + (j = SA[isa]);
                for(TSize hMax = _min(n - i, n - j); h < hMax && *I == *J; ++I, ++J, ++h);
				LCP[i] = h | mark;
                #ifdef SEQAN_DEBUG_INDEX
                    if ((lcpNumer += h) > n) {
                        lcpNumer -= n;
                        ++lcpAvrg;
                    }
                    if (lcpMax < h) lcpMax = h;
					if (!h) ++sigma;
                #endif
			}
			if (h) --h;
            else ++I;
        }
		LCP[SA[n - 1]] = mark;
		SEQAN_PROMARK("permutierte LCP-Tabelle erzeugt");
        #ifdef SEQAN_DEBUG_INDEX
			::std::cerr << "permutierte LCP-Tabelle erzeugt" << ::std::endl;
		#endif
        for(TSize i = 0, j, tmp; i < n; ++i)
			if (LCP[i] & mark) {
				j = i;
				tmp = LCP[j];
				while (SA[j] != i) {
					LCP[j] = LCP[SA[j]] & mask;
					j = SA[j];
				}
				LCP[j] = tmp & mask;
			}
        #ifdef SEQAN_DEBUG_INDEX
			::std::cerr << "LCP-Tabelle erzeugt" << ::std::endl;
		#endif

        #ifdef SEQAN_DEBUG_INDEX
            ::std::cerr << "  n: " << n;
            ::std::cerr << "  lcpMax: " << lcpMax;
            ::std::cerr << "  lcpAvrg: " << (TSize)(lcpAvrg + (lcpNumer + n/2) / n);
            ::std::cerr << "  sigma: " << sigma << ::std::endl;
        #endif
	}


	// Kasai in-place for multiple sequences
    template < typename TLCPTable,
			   typename TString,
			   typename TSpec,
               typename TSA >
    void createLCPTable(
		TLCPTable &LCP,
		StringSet<TString, TSpec> const &sset,
		TSA const &SA,
		Kasai const)
	{
		typedef typename Concatenator<StringSet<TString, TSpec> const>::Type	TText;
		typedef typename StringSetLimits<StringSet<TString, TSpec> const>::Type	TLimitsString;
		typedef typename Value<TSA>::Type										TPair;
		typedef _PairDecrementer<TPair, TLimitsString>							TDecrementer;
		typedef typename Value<TLCPTable>::Type									TSize;

		#ifdef SEQAN_DEBUG_INDEX
			if (sizeof(TSize) > 4)
				::std::cerr << "WARNING: TSize size is greater 4 (Kasai)" << ::std::endl;
        #endif

		TText &s = concat(sset);
		TSize n = length(s);
        if (n < 2) return;

        #ifdef SEQAN_DEBUG_INDEX
            TSize lcpMax = 0, lcpAvrg = 0, lcpNumer = 0, sigma = 1;	// for lcpMax, lcpMean, |Sigma|
        #endif

		TSize mark = ~(~0u>>1);
		TSize mask =   ~0u>>1;

		TLimitsString const &limits = stringSetLimits(sset);
		{
			typename Iterator<TSA const>::Type itSA = begin(SA);
			for(TSize i = 0; i < n; ++i, ++itSA)
				LCP[posGlobalize(*itSA, limits)] = i;
		}
		
		SEQAN_PROMARK("Suffix-Array invertiert");
        #ifdef SEQAN_DEBUG_INDEX
			::std::cerr << "Suffix-Array invertiert" << ::std::endl;
		#endif

		typename Iterator<TText const>::Type Ibegin = begin(s);
        typename Iterator<TText const>::Type I = Ibegin, J;

		TDecrementer dec(limits);

		for(TSize i = 0, h = 0, j, isa; i < n; ++i, --dec) {
			if ((isa = LCP[i] + 1) < n) {
				j = posGlobalize(SA[isa], limits);
				J = Ibegin + h + j;

                for(TSize hMax = _min(getValueI2((TPair)dec), n - j); h < hMax && *I == *J; ++I, ++J, ++h);
				LCP[i] = h | mark;
                #ifdef SEQAN_DEBUG_INDEX
                    if ((lcpNumer += h) > n) {
                        lcpNumer -= n;
                        ++lcpAvrg;
                    }
                    if (lcpMax < h) lcpMax = h;
					if (!h) ++sigma;
                #endif
			}
			if (h) 
				--h;
			else 
				++I;
        }

		LCP[posGlobalize(SA[n - 1], limits)] = mark;

		SEQAN_PROMARK("permutierte LCP-Tabelle erzeugt");
        #ifdef SEQAN_DEBUG_INDEX
			::std::cerr << "permutierte LCP-Tabelle erzeugt" << ::std::endl;
		#endif
        for(TSize sa_j, i = 0, j, tmp; i < n; ++i)
			if (LCP[i] & mark) {
				j = i;
				tmp = LCP[j];
				sa_j = posGlobalize(SA[j], limits);
				while (sa_j != i) {
					LCP[j] = LCP[sa_j] & mask;
					j = sa_j;
					sa_j = posGlobalize(SA[j], limits);
				}
				LCP[j] = tmp & mask;
			}
        #ifdef SEQAN_DEBUG_INDEX
			::std::cerr << "LCP-Tabelle erzeugt" << ::std::endl;
		#endif

        #ifdef SEQAN_DEBUG_INDEX
            ::std::cerr << "  n: " << n;
            ::std::cerr << "  lcpMax: " << lcpMax;
            ::std::cerr << "  lcpAvrg: " << (TSize)(lcpAvrg + (lcpNumer + n/2) / n);
            ::std::cerr << "  sigma: " << sigma << ::std::endl;
        #endif
	}


//}

}

#endif
