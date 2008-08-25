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
  $Id: pipe_sampler.h,v 1.1 2008/08/25 16:20:02 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_SAMPLER_H
#define SEQAN_HEADER_PIPE_SAMPLER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	template <int I, typename T = void>
	struct _SkewDC;

//////////////////////////////////////////////////////////////////////////////

	template < unsigned m, typename TCompression = void >
	struct Sampler;

    template < typename TInput, unsigned m, typename TCompression >
    struct Value< Pipe< TInput, Sampler<m, TCompression> > > {
        typedef Tuple<typename Value<TInput>::Type, m, TCompression>	mTuple;
        typedef Pair<typename Size<TInput>::Type, mTuple, Compressed>	Type;
    };

//////////////////////////////////////////////////////////////////////////////

    template < typename TInput, unsigned m, typename TCompression, typename TPair, typename TLimitsString >
    struct Value< Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> > > {
        typedef Tuple<typename Value<TInput>::Type, m, TCompression>	mTuple;
        typedef Pair<TPair, mTuple, Compressed>							Type;
    };

//////////////////////////////////////////////////////////////////////////////


/**
.Spec.Sampler:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs m-tuples beginning at a position of difference cover DC.
..signature:Pipe<TInput, Sampler<m, DC[, TCompression]> >
..param.TInput:The type of the pipeline module this module reads from.
..param.m:The tuple size.
..param.DC:A set of non-negative integers less than $m$.
..param.DC:$DC[0]$ contains the size of the set and $DC[1..DC[0]]$ contains the distinct and ordered elements.
..param.TCompression:Enable/Disable compression.
..param.TCompression:If $void$, no compression is used.
..param.TCompression:If $Compressed$, bit-compressed @Class.Tuple@s are used.
...default:void.
..example:The set ${1,2,4}$ is represented by $int DC[] = { 3, 1, 2, 4 }$.
..remarks:The output type is a @Class.Pair@ of size type and @Class.Tuple@ of input elements and length m (i.e. $Pair<Size<TInput>::Type, Tuple<Value<TInput>::Type, m, TCompression> >$).
..remarks:The first output field contains the beginning position of the m-tuple in the second field.
The m-tuples are substrings of the input stream beginning at positions $i$, with $i mod m$ is element of the set DC.
*/

    //////////////////////////////////////////////////////////////////////////////
    // sampler class
    template < typename TInput, unsigned m, typename TCompression >
    struct Pipe< TInput, Sampler<m, TCompression> >
    {
        typedef typename Value<Pipe>::Type  OutType;
        typedef typename Size<Pipe>::Type   SizeType;

		TInput		&in;
        bool        filter[m];
        SizeType    idx, _size, _rest;
        unsigned    idxMod;
        OutType     tmp1, tmp2;
        OutType     *outRef, *tmpRef;
        bool        last;
        
        Pipe(TInput& _in):
            in(_in),
            outRef(&tmp1),
            tmpRef(&tmp2) {}
        
        inline void prepare() {
            memset<sizeof(filter), 0>(filter);
			for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
                filter[_SkewDC<m>::VALUE[i]] = true;

            idx = length(in);
            idxMod = idx % m;

            while (!filter[idxMod] && !eof(in)) {
                ++in;
                if (idxMod == 0) idxMod = m;
                --idxMod; --idx;
            }
            _rest = length(*this);
            fill();
            swap();
        }
        
        inline void fill() {
            unsigned i;
            for(i = 0; i < m && !eof(in); ++i, ++in)
                tmpRef->i2.i[i] = *in;
            last = eof(in);
            for(; i < m; ++i)
                tmpRef->i2.i[i] = 0;
            tmpRef->i1 = idx;
        }
        
        inline void rotate(unsigned r) {
            for(unsigned i = 0; i < m; ++i, ++r) {
                if (r == m) r = 0;
                tmpRef->i2.i[i] = outRef->i2.i[r];
            }
        }
        
        inline void swap() {
            OutType *newOutRef = tmpRef;
            tmpRef = outRef;
            outRef = newOutRef;
        }
        
        inline OutType const& operator*() {
            return *outRef;
        }
        
        Pipe& operator++() {
            unsigned skipped = 0;
			if (--_rest) {
				if (!last)
					do {
						outRef->i2.i[skipped++] = *in;
						++in;
						if (idxMod == 0) idxMod = m;
						--idxMod; --idx;
						if (eof(in)) {
							last = true;
							while (!filter[idxMod]) {
								outRef->i2.i[skipped++] = 0;
								if (idxMod == 0) idxMod = m;
								--idxMod; --idx;
							};
							break;
						}
					} while (!filter[idxMod]);
				else
					do {
						outRef->i2.i[skipped++] = 0;
						if (idxMod == 0) idxMod = m;
						--idxMod; --idx;
					} while (!filter[idxMod]);
			}
            rotate(skipped);
            tmpRef->i1 = idx;
            swap();
            return *this;
        }        
    };


    //////////////////////////////////////////////////////////////////////////////
    // sampler class (uses bit compression)
    template < typename TInput, unsigned m >
    struct Pipe< TInput, Sampler<m, Compressed> >
    {
        typedef typename Value<Pipe>::Type  OutType;
        typedef typename Size<Pipe>::Type   SizeType;
        typedef typename OutType::T2        TTuple;

		TInput		&in;
        bool        filter[m];
        SizeType    _size, _rest;
        unsigned    idxMod;
        OutType     tmp;
        bool        last;
        
        Pipe(TInput& _in):
            in(_in) {}
        
        inline void prepare() {
            memset<sizeof(filter), 0>(filter);
            for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
                filter[_SkewDC<m>::VALUE[i]] = true;

            tmp.i1 = length(in);
            idxMod = tmp.i1 % m;

            while (!filter[idxMod] && !eof(in)) {
                ++in;
                if (idxMod == 0) idxMod = m;
                --idxMod; --tmp.i1;
            }
            _rest = length(*this);
            fill();
        }
        
        inline void fill() {
            unsigned i;
            clear(tmp.i2);
            for(i = 0; i < m && !eof(in); ++i, ++in) {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
            }
            last = eof(in);
            tmp.i2 <<= (m - i);
        }
        
        inline OutType const& operator*() {
            return tmp;
        }
        
        Pipe& operator++() {
			if (--_rest) {
				if (!last)
					do {
						tmp.i2 <<= 1;
						tmp.i2 |= *in;
						++in;
						if (idxMod == 0) idxMod = m;
						--idxMod; --tmp.i1;
						if (eof(in)) {
							last = true;
							while (!filter[idxMod]) {
								tmp.i2 <<= 1;
								if (idxMod == 0) idxMod = m;
								--idxMod; --tmp.i1;
							};
							break;
						}
					} while (!filter[idxMod]);
				else
					do {
						tmp.i2 <<= 1;
						if (idxMod == 0) idxMod = m;
						--idxMod; --tmp.i1;
					} while (!filter[idxMod]);
			}
            return *this;
        }        
    };



    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned m, typename TCompression >
	inline bool control(Pipe< TInput, Sampler<m, TCompression> > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.prepare();
		return true;
    }

    template < typename TInput, unsigned m, typename TCompression >
	inline bool control(Pipe< TInput, Sampler<m, TCompression> > &me, ControlEof const &command) {
		return me._rest == 0;
    }

    template < typename TInput, unsigned m, typename TCompression >
    inline typename Size< Pipe< TInput, Sampler<m, TCompression> > >::Type
	length(Pipe< TInput, Sampler<m, TCompression> > const &me) {
        typename Size< Pipe< TInput, Sampler<m> > >::Type _size = 0, n = length(me.in);
        for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
            if (_SkewDC<m>::VALUE[i])
                _size += (n + m - _SkewDC<m>::VALUE[i]) / m;
            else
                _size += n / m;
        return _size;
    }



//////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////
    // sampler class
    template < typename TInput, unsigned m, typename TCompression, typename TPair, typename TLimitsString >
    struct Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> >
    {
        typedef typename Value<Pipe>::Type  OutType;
        typedef typename Size<Pipe>::Type   SizeType;

		typedef _PairDecrementer<TPair, TLimitsString, m>	Decrementer;

		TInput		&in;
        bool		filter[m];
        Decrementer	localPos;
		SizeType	_size, _rest;
        OutType		tmp1, tmp2;
        OutType		*outRef, *tmpRef;
        bool		last;

		TLimitsString const &limits;
        
        Pipe(TInput& _in, TLimitsString const &_limits):
            in(_in),
            outRef(&tmp1),
            tmpRef(&tmp2),
			limits(_limits) {}
       
        inline void prepare() {
            memset<sizeof(filter), 0>(filter);
			for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
                filter[_SkewDC<m>::VALUE[i]] = true;

			setHost(localPos, limits);

			while (!filter[localPos.residue] && !eof(in)) {
                ++in;
				--localPos;
            }
            _rest = length(*this);
            fill();
            swap();
        }
        
        inline void fill() {
            unsigned i;
            for(i = 0; i < m && !eof(in); ++i, ++in)
                tmpRef->i2.i[i] = *in;
            last = eof(in);
            for(; i < m; ++i)
                tmpRef->i2.i[i] = 0;
            tmpRef->i1 = localPos;
        }
        
        inline void rotate(unsigned r) {
            for(unsigned i = 0; i < m; ++i, ++r) {
                if (r == m) r = 0;
                tmpRef->i2.i[i] = outRef->i2.i[r];
            }
        }
        
        inline void swap() {
            OutType *newOutRef = tmpRef;
            tmpRef = outRef;
            outRef = newOutRef;
        }
        
        inline OutType const& operator*() {
            return *outRef;
        }
        
        Pipe& operator++() {
            unsigned skipped = 0;
			if (--_rest) {
				if (!last)
					do {
						outRef->i2.i[skipped++] = *in;
						++in;
						--localPos;
						if (eof(in)) {
							last = true;
							while (!filter[localPos.residue]) {
								outRef->i2.i[skipped++] = 0;
								--localPos;
							};
							break;
						}
					} while (!filter[localPos.residue]);
				else
					do {
						outRef->i2.i[skipped++] = 0;
						--localPos;
					} while (!filter[localPos.residue]);
				rotate(skipped);
				tmpRef->i1 = localPos;
				swap();
			}
            return *this;
        }        
    };


    //////////////////////////////////////////////////////////////////////////////
    // sampler class (uses bit compression)
	template < typename TInput, unsigned m, typename TPair, typename TLimitsString >
    struct Pipe< TInput, Multi<Sampler<m, Compressed>, TPair, TLimitsString> >
    {
        typedef typename Value<Pipe>::Type  OutType;
        typedef typename Size<Pipe>::Type   SizeType;
        typedef typename OutType::T2        TTuple;

		typedef _PairDecrementer<TPair, TLimitsString, m>	Decrementer;

		TInput		&in;
        bool		filter[m];
        SizeType	_size, _rest;
        Decrementer	localPos;
        OutType		tmp;
        bool		last;

		TLimitsString const &limits;
        
        Pipe(TInput& _in, TLimitsString const &_limits):
            in(_in),
			limits(_limits) {}
        
        inline void prepare() {
            memset<sizeof(filter), 0>(filter);
            for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
                filter[_SkewDC<m>::VALUE[i]] = true;

			setHost(localPos, limits);

			while (!filter[localPos.residue] && !eof(in)) {
                ++in;
				--localPos;
            }
            _rest = length(*this);
            fill();
        }
        
        inline void fill() {
            unsigned i;
            clear(tmp.i2);
            for(i = 0; i < m && !eof(in); ++i, ++in) {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
            }
            last = eof(in);
            tmp.i2 <<= (m - i);
			tmp.i1 = localPos;
        }
        
        inline OutType const& operator*() {
            return tmp;
        }
        
        Pipe& operator++() {
			if (--_rest) {
				if (!last)
					do {
						tmp.i2 <<= 1;
						tmp.i2 |= *in;
						++in;
						--localPos;
						if (eof(in)) {
							last = true;
							while (!filter[localPos.residue]) {
								tmp.i2 <<= 1;
								--localPos;
							};
							break;
						}
					} while (!filter[localPos.residue]);
				else
					do {
						tmp.i2 <<= 1;
						--localPos;
					} while (!filter[localPos.residue]);
			}
            tmp.i1 = localPos;
            return *this;
        }        
    };



    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
	template < typename TInput, unsigned m, typename TCompression, typename TPair, typename TLimitsString >
	inline bool control(Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.prepare();
		return true;
    }

	template < typename TInput, unsigned m, typename TCompression, typename TPair, typename TLimitsString >
	inline bool control(Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> > &me, ControlEof const &command) {
		return me._rest == 0;
    }

	template < typename TInput, unsigned m, typename TCompression, typename TPair, typename TLimitsString >
	inline bool control(Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> > &me, ControlEos const &command) {
		return control(me, ControlEof());
    }

	template < typename TInput, unsigned m, typename TCompression, typename TPair, typename TLimitsString >
    inline typename Size< Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> > >::Type
	length(Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> > const &me)
	{
		typedef typename Size< Pipe< TInput, Multi<Sampler<m, TCompression>, TPair, TLimitsString> > >::Type TSize;
		typename Iterator<TLimitsString const, Standard>::Type it = begin(me.limits), itEnd = end(me.limits);
		
		if (it == itEnd) return 0;

		TSize sum = 0;
		TSize size;
		TSize old = *it; ++it;

		while (it != itEnd) {
			size = *it - old;
			old = *it;
			
			for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
				if (_SkewDC<m>::VALUE[i])
					sum += (size + m - _SkewDC<m>::VALUE[i]) / m;
				else
					sum += size / m;

			++it;
		}

        return sum;
    }
//}

}

#endif
