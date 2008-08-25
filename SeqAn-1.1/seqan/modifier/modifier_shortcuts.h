 /*==========================================================================
                SeqAn - The Library for object Analysis
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
  $Id: modifier_shortcuts.h,v 1.1 2008/08/25 16:20:07 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_SHORTCUTS_H
#define SEQAN_HEADER_MODIFIER_SHORTCUTS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

typedef ModView< FunctorComplement<Dna> >	ModComplementDna;
typedef ModView< FunctorComplement<Dna5> >	ModComplementDna5;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >		DnaStringComplement;
typedef ModifiedString<Dna5String, ModView< FunctorComplement<Dna5> > >		Dna5StringComplement;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModReverse>		DnaStringReverse;
typedef ModifiedString<Dna5String, ModReverse>		Dna5StringReverse;

//////////////////////////////////////////////////////////////////////////////
/*
typedef ModifiedString<DnaStringReverse, ModComplementDna>		DnaStringReverseComplement;
typedef ModifiedString<Dna5StringReverse, ModComplementDna5>	Dna5StringReverseComplement;
*/

typedef ModifiedString<
			ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >, 
			ModReverse
		>	DnaStringReverseComplement;

typedef ModifiedString<
			ModifiedString<	Dna5String, ModView< FunctorComplement<Dna5> > >, 
			ModReverse
		>	Dna5StringReverseComplement;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template < typename TSequence >
inline void complementInPlace(TSequence & sequence) 
{
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
} 

template < typename TSequence, typename TSpec >
inline void complementInPlace(StringSet<TSequence, TSpec> & stringSet)
{
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		complementInPlace(stringSet[seqNo]);
}

template < typename TSequence >
inline void reverseComplementInPlace(TSequence & sequence) 
{
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
	reverseInPlace(sequence);
} 

template < typename TSequence, typename TSpec >
inline void reverseComplementInPlace(StringSet<TSequence, TSpec> & stringSet)
{
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		reverseComplementInPlace(stringSet[seqNo]);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Shortcut.ModComplementDna:
..cat:Modifier
..summary:Modifier specialization type for the complement of @Spec.Dna@ alphabet sequences.
..signature:DnaStringComplement
..shortcutfor:Spec.ModView
...signature:ModView< FunctorComplement<Dna> >
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.ModComplementDna5:
..cat:Modifier
..summary:Modifier specialization type for the complement of @Spec.Dna5@ alphabet sequences.
..signature:Dna5StringComplement
..shortcutfor:Spec.ModView
...signature:ModView< FunctorComplement<Dna5> >
..see:Spec.ModView
..see:Class.FunctorComplement
*/


/**
.Shortcut.DnaStringComplement:
..cat:Modifier
..summary:Modifier for the complement of a @Shortcut.DnaString@.
..signature:DnaStringComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >
..see:Shortcut.DnaString
..see:Class.ModifiedString
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.Dna5StringComplement:
..cat:Modifier
..summary:Modifier for the complement of a @Shortcut.Dna5String@.
..signature:Dna5StringComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<Dna5String, ModView< FunctorComplement<Dna5> > >
..see:Shortcut.Dna5String
..see:Class.ModifiedString
..see:Spec.ModView
..see:Class.FunctorComplement
*/


/**
.Shortcut.DnaStringReverse:
..cat:Modifier
..summary:Modifier for the reverse of a @Shortcut.DnaString@.
..signature:DnaStringReverse
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<DnaString, ModReverse>
..see:Shortcut.DnaString
..see:Class.ModifiedString
..see:Spec.ModReverse
*/

/**
.Shortcut.Dna5StringReverse:
..cat:Modifier
..summary:Modifier for the reverse of a @Shortcut.Dna5String@.
..signature:Dna5StringReverse
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<Dna5String, ModReverse>
..see:Shortcut.Dna5String
..see:Class.ModifiedString
..see:Spec.ModReverse
*/


/**
.Shortcut.DnaStringReverseComplement:
..cat:Modifier
..summary:Modifier for the reverse complement of a @Shortcut.DnaString@.
..signature:DnaStringReverseComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >, ModReverseComplement>
..see:Shortcut.DnaString
..see:Class.ModifiedString
..see:Spec.ModReverse
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.Dna5StringReverseComplement:
..cat:Modifier
..summary:Modifier for the reverse complement of a @Shortcut.Dna5String@.
..signature:Dna5StringReverseComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<ModifiedString<Dna5String, ModView< FunctorComplement<Dna> > >, ModReverseComplement>
..see:Shortcut.Dna5String
..see:Class.ModifiedString
..see:Spec.ModReverse
..see:Spec.ModView
..see:Class.FunctorComplement
*/

}

#endif
