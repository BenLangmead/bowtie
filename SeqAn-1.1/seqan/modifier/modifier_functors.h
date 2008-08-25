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
  $Id: modifier_functors.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_FUNCTORS_H
#define SEQAN_HEADER_MODIFIER_FUNCTORS_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// text transformation
	//////////////////////////////////////////////////////////////////////////////

/**
.Class.FunctorUpcase:
..cat:Modifier
..summary:Functor that returns the upper case character to a given character.
..signature:FunctorUpcase<TValue>
..param.TValue:The input value type.
..remarks:This Functor is a derivation of the STL unary function.
*/

    template <typename InType, typename Result = InType>
	struct FunctorUpcase : public ::std::unary_function<InType,Result> 
	{
        inline Result operator()(InType x) const {
			if (('a' <= x) && (x <= 'z')) return (x + ('A' - 'a'));
			return x; 
		}
    };

/**
.Class.FunctorLowcase:
..cat:Modifier
..summary:Functor that returns the lower case character to a given character.
..signature:FunctorLowcase<TValue>
..param.TValue:The input value type.
..remarks:This Functor is a derivation of the STL unary function.
*/

    template <typename InType, typename Result = InType>
    struct FunctorLowcase : public ::std::unary_function<InType,Result> 
	{
        inline Result operator()(InType x) const {
			if (('A' <= x) && (x <= 'Z')) return (x + ('a' - 'A'));
			return x; 
		}
    };


	//////////////////////////////////////////////////////////////////////////////
	// alphabet transformation
	//////////////////////////////////////////////////////////////////////////////

/**
.Class.FunctorConvert:
..cat:Modifier
..summary:Functor that converts a $TInValue$ type to a $TOutValue$ type character.
..signature:FunctorConvert<TInValue, TOutValue>
..param.TInValue:The input value type.
..param.TOutValue:The output value type.
..remarks:This Functor is a derivation of the STL unary function.
*/

    template <typename InType, typename OutType>
    struct FunctorConvert : public ::std::unary_function<InType,OutType> 
	{
        inline OutType operator()(InType x) const {
			return x; 
		}
    };

	//////////////////////////////////////////////////////////////////////////////
	// DNA complement
	//////////////////////////////////////////////////////////////////////////////

/**
.Class.FunctorComplement:
..cat:Modifier
..summary:Functor that returns the complement nucleotide to a given nucleotide.
..signature:FunctorComplement<TValue>
..param.TValue:The input value type.
...type:Spec.Dna
...type:Spec.Dna5
..remarks:This Functor is a derivation of the STL unary function.
*/

    template <typename TValue>
    struct FunctorComplement;
	
//____________________________________________________________________________

	template <typename T = void>
	struct _Translate_Table_Dna5_2_Dna5Complement
	{
		static char const VALUE[5];
	};
	template <typename T>
	char const _Translate_Table_Dna5_2_Dna5Complement<T>::VALUE[5] = {'T', 'G', 'C', 'A', 'N'};

//____________________________________________________________________________


    template <>
    struct FunctorComplement<Dna> : public ::std::unary_function<Dna,Dna> 
	{
        inline Dna operator()(Dna x) const {
			return _Translate_Table_Dna5_2_Dna5Complement<>::VALUE[x.value]; 
		}
    };

    template <>
    struct FunctorComplement<Dna5> : public ::std::unary_function<Dna5,Dna5> 
	{
        inline Dna operator()(Dna5 x) const {
			return _Translate_Table_Dna5_2_Dna5Complement<>::VALUE[x.value]; 
		}
    };


}

#endif
