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
  $Id: index_pizzachili_find.h,v 1.1 2008/08/25 16:20:05 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_PIZZACHILI_FIND_H
#define SEQAN_HEADER_INDEX_PIZZACHILI_FIND_H

namespace SEQAN_NAMESPACE_MAIN {

struct _PizzaChiliFinder;

/**
.Tag.Index Find Algorithm
..cat:Searching
..tag.PizzaChiliFinder:Finds an occurrence in a @Spec.Pizza & Chili Index@ index.
...remarks:The actual algorithm used for searching depends on the @Tag.Pizza & Chili Index Tags@ used.
..see:Spec.Pizza & Chili Index
*/

typedef Tag<_PizzaChiliFinder> const PizzaChiliFinder;

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
struct DefaultFinder<Index<TText, PizzaChili<TSpec> > > {
    typedef PizzaChiliFinder Type;
};

//////////////////////////////////////////////////////////////////////////////

// We define a special suffix array fibre for our finder. This fibre is used to
// represent the search results and our search results aren't actually part of
// a suffix array but rather a C array of ulong values.
template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> >, Tag<_Fibre_SA> const> {
    typedef typename impl::ulong_t* Type;
};

//////////////////////////////////////////////////////////////////////////////

// The following specialization of the finder is necessary only because we have
// to free the data from the C API in the destructor.
template <typename TText, typename TSpec>
class Finder<Index<TText, PizzaChili<TSpec> >, PizzaChiliFinder>
    : public Finder<Index<TText, PizzaChili<TSpec> >, Default>
{
    typedef Finder<Index<TText, PizzaChili<TSpec> >, Default> TBase;

    typedef typename TBase::TIndex TIndex;
    typedef typename TBase::TSA TSA;
    typedef typename TBase::TIterator TIterator;

public:

    using TBase::index;
    using TBase::range;
    using TBase::data_iterator;

    Finder() : TBase() { }

    Finder(TIndex& index) : TBase(index) { }

    Finder(TIndex const& index) : TBase(index) { }

    ~Finder() {
SEQAN_CHECKPOINT
        if (range.i1 != 0)
            ::std::free(range.i1);
    }
};

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TPattern>
    inline uchar_t* getPizzaChiliString(TPattern const& pattern) {
SEQAN_CHECKPOINT
        typedef
            typename _RemoveConst<
                typename Value<TPattern>::Type
            >::Type alph_t;

        // This const_cast is safe because 'cstr' is only read.
        // On the other hand, this cast is necessary to prevent a copy from
        // being made: this would result in an invalid pointer to local memory.
        String<alph_t, CStyle> const cstr = *const_cast<TPattern*>(&pattern);
        // This const_cast is safe because the return value is only read.
        return
            reinterpret_cast<uchar_t*>(
                const_cast<alph_t*>(static_cast<alph_t const*>(cstr))
            );
    }

    inline uchar_t* getPizzaChiliString(char const* pattern) {
SEQAN_CHECKPOINT
        return reinterpret_cast<uchar_t*>(const_cast<char*>(pattern));
    }
} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec, typename TSpecFinder, typename TPattern>
inline void _findFirstIndex(
    Finder<Index<TText, PizzaChili<TSpec> >, TSpecFinder>& finder,
    TPattern const& pattern,
    PizzaChiliFinder const
) {
SEQAN_CHECKPOINT
    typedef Index<TText, PizzaChili<TSpec> > TIndex;
    typedef typename PizzaChiliCodeProvider<TSpec>::Type TCodeProvider;

    if (finder.range.i1 != 0)
        ::std::free(finder.range.i1);

    TIndex& index = haystack(finder);
    indexRequire(index, PizzaChili_Compressed());

    impl::uchar_t* c_pattern =
        impl::getPizzaChiliString(pattern);
    impl::ulong_t patternlength = length(pattern);

    impl::ulong_t numocc;
    impl::ulong_t* occ;

    impl::error_t e =
        TCodeProvider::locate(
            index.index_handle,
            c_pattern,
            patternlength,
            &occ,
            &numocc
        );

    if (e != 0)
        SEQAN_ABORT(TCodeProvider::error_index(e));

    finder.range.i1 = occ;
    finder.range.i2 = occ + numocc;
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_FIND_H
