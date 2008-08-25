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
  $Id: pipe_merger7.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_MERGER7_H
#define SEQAN_HEADER_INDEX_MERGER7_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	template <int I, typename T = void>
	struct _SkewShift;

	template <int I, typename T = void>
	struct _SkewNIndx;


	template <typename T>
	struct _SkewShift<7, T> {
		static const unsigned VALUE[7][7];
	};

	template <typename T>
	struct _SkewNIndx<7, T> {
		static const unsigned VALUE[7][7];
	};


	template <typename T>
	const unsigned _SkewShift<7, T>::VALUE[7][7] =
								  {{0,6,5,6,3,3,5},
                                   {6,0,0,6,0,4,4},
                                   {5,0,0,1,0,1,5},
                                   {6,6,1,0,2,1,2},
                                   {3,0,0,2,0,3,2},
                                   {3,4,1,1,3,0,4},
                                   {5,4,5,2,2,4,0}};


	template <typename T>
	const unsigned _SkewNIndx<7, T>::VALUE[7][7] =
								  {{0,0,0,0,0,1,2},
                                   {0,0,0,0,1,1,2},
                                   {0,1,1,1,1,2,2},
                                   {0,0,1,1,1,1,2},
                                   {0,0,1,2,2,2,2},
                                   {0,0,0,1,2,2,2},
                                   {0,0,0,0,1,2,2}};


    template <typename TValue>
    struct SkewDCStream {
        TValue      i;
        unsigned    stream;
    };
    
    template <typename TValue>
	std::ostream& operator<<(std::ostream &out, const SkewDCStream<TValue> &s) {
		out << "< " << s.i.i1 << " , [ ";
		for(unsigned i = 0; i < 3; ++i)
			out << s.i.i2[i] << " ";
		out << "] , [ ";
		for(unsigned i = 0; i < 6; ++i)
			out << s.i.i3[i] << " ";
		out << "] , " << s.stream << " >";
		return out;
	}


	// greater-operator for two SkewDCStreams
    template <typename TValue>
    struct CompareSkewDCStream :
        public ::std::binary_function < SkewDCStream<TValue>,
                                      SkewDCStream<TValue>,
                                      bool >
    {
		template <typename TSize>
		inline static bool crossBoarderCompare(TSize const a, TSize const b) {
			return a < b;
		}

		template <typename T1, typename T2, typename TCompression>
		inline static bool crossBoarderCompare(
			Pair<T1, T2, TCompression> const &a,
			Pair<T1, T2, TCompression> const &b)
		{
			return (getValueI1(a) >  getValueI1(b)) ||
				  ((getValueI1(a) == getValueI1(b)) && (getValueI2(a) < getValueI2(b)));
		}

        inline bool operator()(const SkewDCStream<TValue> &a,
			                   const SkewDCStream<TValue> &b) const 
        {
            typedef typename TValue::T2::T SizeType;
            int shft = _SkewShift<7>::VALUE[a.stream][b.stream];
            for(int i = 0; i < shft; ++i) {
                if (a.i.i3[i] < b.i.i3[i]) return false;
                if (a.i.i3[i] > b.i.i3[i]) return true;
            }
            SizeType na = a.i.i2[_SkewNIndx<7>::VALUE[a.stream][shft]];
            SizeType nb = b.i.i2[_SkewNIndx<7>::VALUE[b.stream][shft]];
            if (na < nb) return false;
            if (na > nb) return true;

			// we get here, only if a septet crosses the boarder of
			// 1) the single text (a/b.i.i1 is an ordinal number)
			// 2) a sequence in a multiple sequence text (a/b.i.i1 is a Pair)
			return crossBoarderCompare(getValueI1(a.i), getValueI1(b.i));
        }
    };


	// greater-operator for two SkewDCStreams (optimized for bit-compressed character tuples)
    template <typename T1, typename T2, typename T, unsigned _size>
    struct CompareSkewDCStream< Triple<T1,T2,Tuple<T,_size,Compressed>, Compressed> > :
        public ::std::binary_function < SkewDCStream<Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> >,
                                      SkewDCStream<Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> >,
                                      bool >
    {
		template <typename TSize>
		inline static bool crossBoarderCompare(TSize const a, TSize const b) {
			return a < b;
		}

		template <typename _T1, typename _T2, typename _TCompression>
		inline static bool crossBoarderCompare(
			Pair<_T1, _T2, _TCompression> const &a,
			Pair<_T1, _T2, _TCompression> const &b)
		{
			return (getValueI1(a) >  getValueI1(b)) ||
				  ((getValueI1(a) == getValueI1(b)) && (getValueI2(a) < getValueI2(b)));
		}

        typedef Tuple<T,_size,Compressed> T3;
        inline bool operator()(const SkewDCStream<Triple<T1,T2,T3,Compressed> > &a,
			                   const SkewDCStream<Triple<T1,T2,T3,Compressed> > &b) const 
        {
            typedef typename T2::T SizeType;
            int shft = _SkewShift<7>::VALUE[a.stream][b.stream];
            typename T3::CT mask = ~((1 << ((_size - shft) * T3::bitSize)) - 1);

            if ((a.i.i3.i & mask) < (b.i.i3.i & mask)) return false;
            if ((a.i.i3.i & mask) > (b.i.i3.i & mask)) return true;
            SizeType na = a.i.i2[_SkewNIndx<7>::VALUE[a.stream][shft]];
            SizeType nb = b.i.i2[_SkewNIndx<7>::VALUE[b.stream][shft]];
            if (na < nb) return false;
            if (na > nb) return true;

			// we get here, only if a septet crosses the boarder of
			// 1) the single text (a/b.i.i1 is an ordinal number)
			// 2) a sequence in a multiple sequence text (a/b.i.i1 is a Pair)
			return crossBoarderCompare(getValueI1(a.i), getValueI1(b.i));
		}
    };




	struct Merger7;

    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124 >
    struct Value< Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7 > > {
		typedef typename Value<TInput0>::Type::T1 Type;
    };


	template <typename TLimitsString>
	struct Merger7Multi;

    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124, typename TLimitsString >
    struct Value< Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7Multi<TLimitsString> > > {
		typedef typename Value<TInput0>::Type::T1 Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // merger7 class
    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124 >
    struct Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7 >
    {
        typedef typename Value<TInput0>::Type       InType0;
        typedef typename Value<TInput3>::Type       InType3;
        typedef typename Value<TInput5>::Type       InType5;
        typedef typename Value<TInput6>::Type       InType6;
        typedef typename Value<TInput124>::Type     InType124;
        typedef typename Size<Pipe>::Type           SizeType;

        typedef typename InType0::T3::T             Type;
                
        typedef SkewDCStream<InType0>               TSkewDCStream;

		TSkewDCStream								inValue[7];
		unsigned									rank[5];
		unsigned									first;
        CompareSkewDCStream<InType0>				streamGreater;

        Bundle5 <
            TInput0,
            TInput3,
            TInput5,
            TInput6,
            TInput124 > in;
        SizeType        N;
        
        Pipe(Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 > _in):
            in(_in) {}

        template <typename T1, typename T2, typename T3>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,T3,Compressed> const &src) {
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            memcpy(&dst.i.i3, &src.i3, sizeof(T3));
        }

        template <typename T1, typename T2, typename T, unsigned _size>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> const &src) 
		{
			typedef typename InType0::T3 CharTuple;
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            dst.i.i3 = src.i3;
            dst.i.i3 <<= CharTuple::size - _size;
        }

        inline typename Value<Pipe>::Type const operator*() {
            return getValueI1(inValue[rank[first]].i);
        }

		inline void insertStream(unsigned stream) {
            switch (stream) {
                case 0:
					if (eof(in.in1)) {
						++first;
						return;
					}
					inValue[0].i.i1 = N - (*in.in1).i1;
					_copy(inValue[0], *in.in1);
					++in.in1;
                    break;

                case 3:
					if (eof(in.in2)) {
						++first;
						return;
					}
					inValue[3].i.i1 = N - (*in.in2).i1;
					_copy(inValue[3], *in.in2);
					++in.in2;
                    break;

                case 5:
					if (eof(in.in3)) {
						++first;
						return;
					}
					inValue[5].i.i1 = N - (*in.in3).i1;
					_copy(inValue[5], *in.in3);
					++in.in3;
                    break;

                case 6:
					if (eof(in.in4)) {
						++first;
						return;
					}
					inValue[6].i.i1 = N - (*in.in4).i1;
					_copy(inValue[6], *in.in4);
					++in.in4;
					break;

				default:	// case 1, 2, or 4
					if (eof(in.in5)) {
						++first;
						return;
					}
					// calculate residue class from the suffix length
					inValue[1].stream = (*in.in5).i1 % 7;
					inValue[1].i.i1 = N - (*in.in5).i1;
					_copy(inValue[1], *in.in5);
					++in.in5;
            }

			// linear search
			int right;
			for(right = first + 1;  right < 5;  ++right)
				if (!streamGreater(inValue[stream],inValue[rank[right]])) break;

			// remove the least suffix ...
			for(int i = first + 1;  i < right;  ++i)
				rank[i-1] = rank[i];

			// ... and insert the new one
			rank[right-1] = stream;
		}
        
        Pipe& operator++() {
			insertStream(rank[first]);
            return *this;
        }

        void fill() {
			inValue[0].stream = 0;
			inValue[3].stream = 3;
			inValue[5].stream = 5;
			inValue[6].stream = 6;

			first = 4;	insertStream(0);
			--first;    insertStream(3);
			--first;    insertStream(5);
			--first;    insertStream(6);
			--first;	insertStream(1);
        }
	};
    

    //////////////////////////////////////////////////////////////////////////////
    // merger7 class for multiple sequences
    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124, typename TLimitsString >
    struct Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7Multi<TLimitsString> >
    {
        typedef typename Value<TInput0>::Type       InType0;
        typedef typename Value<TInput3>::Type       InType3;
        typedef typename Value<TInput5>::Type       InType5;
        typedef typename Value<TInput6>::Type       InType6;
        typedef typename Value<TInput124>::Type     InType124;
        typedef typename Size<Pipe>::Type           SizeType;

        typedef typename InType0::T3::T             Type;
        
        typedef SkewDCStream<InType0>               TSkewDCStream;

		TSkewDCStream								inValue[7];
		unsigned									rank[5];
		unsigned									first;
        CompareSkewDCStream<InType0>				streamGreater;

        Bundle5 <
            TInput0,
            TInput3,
            TInput5,
            TInput6,
            TInput124 >		in;
		TLimitsString const &limits;
        
        Pipe(Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 > _in, TLimitsString const &_limits):
            in(_in),
			limits(_limits) {}

        template <typename T1, typename T2, typename T3>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,T3,Compressed> const &src) {
            memcpy(&dst.i.i1, &src.i1, sizeof(T1));
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            memcpy(&dst.i.i3, &src.i3, sizeof(T3));
        }

        template <typename T1, typename T2, typename T, unsigned _size>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> const &src) 
		{
			typedef typename InType0::T3 CharTuple;
            memcpy(&dst.i.i1, &src.i1, sizeof(T1));
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            dst.i.i3 = src.i3;
            dst.i.i3 <<= CharTuple::size - _size;
        }

        inline typename Value<Pipe>::Type const operator*() {
            return getValueI1(inValue[rank[first]].i);
        }

		inline void insertStream(unsigned stream) {
            switch (stream) {
                case 0:
					if (eof(in.in1)) {
						++first;
						return;
					}
					_copy(inValue[0], *in.in1);
					++in.in1;
                    break;

                case 3:
					if (eof(in.in2)) {
						++first;
						return;
					}
					_copy(inValue[3], *in.in2);
					++in.in2;
                    break;

                case 5:
					if (eof(in.in3)) {
						++first;
						return;
					}
					_copy(inValue[5], *in.in3);
					++in.in3;
                    break;

                case 6:
					if (eof(in.in4)) {
						++first;
						return;
					}
					_copy(inValue[6], *in.in4);
					++in.in4;
					break;

				default:	// case 1, 2, or 4
					if (eof(in.in5)) {
						++first;
						return;
					}
					// calculate residue class from the suffix length
					typename Value<InType124,1>::Type in_in5_i1 = getValueI1(*in.in5);
					unsigned seqNo = getValueI1(in_in5_i1);
					inValue[1].stream = ((limits[seqNo + 1] - limits[seqNo]) - getValueI2(in_in5_i1)) % 7;

					_copy(inValue[1], *in.in5);
					++in.in5;
            }

			// linear search
			int right;
			for(right = first + 1;  right < 5;  ++right)
				if (!streamGreater(inValue[stream],inValue[rank[right]])) break;

			// remove the least suffix ...
			for(int i = first + 1;  i < right;  ++i)
				rank[i-1] = rank[i];

			// ... and insert the new one
			rank[right-1] = stream;
		}
        
        Pipe& operator++() {
			insertStream(rank[first]);
            return *this;
        }

        void fill() {
			inValue[0].stream = 0;
			inValue[3].stream = 3;
			inValue[5].stream = 5;
			inValue[6].stream = 6;

			first = 4;	insertStream(0);
			--first;    insertStream(3);
			--first;    insertStream(5);
			--first;    insertStream(6);
			--first;	insertStream(1);
        }
	};
    

    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
	inline bool control(Pipe< TInput, Merger7 > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.N = length(me);
        me.fill();
		return true;
	}
    
    template < typename TInput >
	inline bool control(Pipe< TInput, Merger7 > &me, ControlEof const &) {
		return me.first == 5;
    }

    template < typename TInput >
	inline bool control(Pipe< TInput, Merger7 > &me, ControlEos const &) {
		return control(me, ControlEof());
    }

    template < typename TInput >
    inline typename Size< Pipe< TInput, Merger7 > >::Type
    length(Pipe< TInput, Merger7 > const &me) {
        return length(me.in.in1) +
               length(me.in.in2) + 
               length(me.in.in3) +
               length(me.in.in4) +
               length(me.in.in5);
    }




    template < typename TInput, typename TLimitsString >
	inline bool control(Pipe< TInput, Merger7Multi<TLimitsString> > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.fill();
		return true;
	}
    
    template < typename TInput, typename TLimitsString >
	inline bool control(Pipe< TInput, Merger7Multi<TLimitsString> > &me, ControlEof const &) {
		return me.first == 5;
    }

    template < typename TInput, typename TLimitsString >
	inline bool control(Pipe< TInput, Merger7Multi<TLimitsString> > &me, ControlEos const &) {
		return control(me, ControlEof());
    }

    template < typename TInput, typename TLimitsString >
    inline typename Size< Pipe< TInput, Merger7Multi<TLimitsString> > >::Type
    length(Pipe< TInput, Merger7Multi<TLimitsString> > const &me) {
        return length(me.in.in1) +
               length(me.in.in2) + 
               length(me.in.in3) +
               length(me.in.in4) +
               length(me.in.in5);
    }

//}

}

#endif
