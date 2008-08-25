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
  $Id: pipe_source.h,v 1.1 2008/08/25 16:20:02 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_SOURCE_H
#define SEQAN_HEADER_PIPE_SOURCE_H

namespace SEQAN_NAMESPACE_MAIN
{

/*
	template < typename TValue >
    struct Size< Pipe< TValue*, Source<> > > {
		typedef ptrdiff_t Type;
	};
*/

/**
.Spec.Source:
..cat:Pipelining
..general:Class.Pipe
..summary:Pipelining adaptor for arbitrary containers or iterators.
..signature:Pipe<TInput, Source<> >
..param.TInput:The type of container or iterator this module reads from.
*/


	//////////////////////////////////////////////////////////////////////////////
	// source class
	template < typename TInput, typename TSpec >
	struct Pipe< TInput, Source<TSpec> >
    {
		TInput const &in;
		typename Iterator<TInput const, Rooted>::Type cur;

		Pipe(typename _RemoveConst<TInput>::Type &_cont):
			in(_cont) {}

		Pipe(TInput const &_cont):
			in(_cont) {}

        inline typename Value<TInput>::Type const & operator*() const {
            return *cur;
        }
    
        inline Pipe& operator++() {
            ++cur;
            return *this;
        }
    };

    template < typename TInput, typename TSpec, typename TIteratorSpec >
	struct Iterator< Pipe< TInput, Source<TSpec> >, TIteratorSpec> {
		typedef typename Iterator<TInput const, TIteratorSpec>::Type Type;
	};

    template < typename TInput, typename TSpec, typename TTag >
	inline typename Iterator< Pipe< TInput, Source<TSpec> >, Tag<TTag> const >::Type
	begin(Pipe< TInput, Source<TSpec> > &pipe, Tag<TTag> const) {
		return begin(pipe.in);
	}

    template < typename TInput, typename TSpec, typename TTag >
	inline typename Iterator< Pipe< TInput, Source<TSpec> >, Tag<TTag> const >::Type
	end(Pipe< TInput, Source<TSpec> > &pipe, Tag<TTag> const) {
		return end(pipe.in);
	}



	//////////////////////////////////////////////////////////////////////////////
	// simple buffer adaption for fast random accessable containers
    template < typename TContainer >
	struct ContainerBuffer 
	{
		TContainer *cont;

        ContainerBuffer(): cont(NULL) {}
		ContainerBuffer(TContainer &_cont): cont(&_cont) {}

		template <typename TSize>
		inline typename Reference<TContainer>::Type 
		operator[] (TSize i) {
			return (*cont)[i];
		}
    };

    template < typename TContainer, typename TIteratorSpec >
	struct Iterator< ContainerBuffer<TContainer>, TIteratorSpec>:
		Iterator<TContainer, TIteratorSpec> {};

    template < typename TContainer  >
	struct Value< ContainerBuffer<TContainer> >:
		Value<TContainer> {};

    template < typename TContainer >
	struct Size< ContainerBuffer<TContainer> >:
		Size<TContainer> {};


	template < typename TContainer >
    inline typename Size< ContainerBuffer<TContainer> >::Type
    length(ContainerBuffer<TContainer> const &me) {
		if (me.cont)
			return length(*me.cont);
		else
			return 0;
    }


	//////////////////////////////////////////////////////////////////////////////
	// simple buffer adaption for fast random accessable iterators
    template < typename TIterator >
	struct IteratorBuffer 
	{
		TIterator						begin;
		typename Size<TIterator>::Type	size;

        IteratorBuffer(): begin(NULL), size(0) {}
		template <typename TSize>
		IteratorBuffer(TIterator &_begin, TSize _size): begin(&_begin), size(_size) {}

		template <typename TSize>
		inline typename Reference<TIterator>::Type 
		operator[] (TSize i) {
			return *(begin + i);
		}
    };

    template < typename TIterator, typename TIteratorSpec >
	struct Iterator< IteratorBuffer<TIterator>, TIteratorSpec >:
		Iterator<TIterator, TIteratorSpec> {};

    template < typename TIterator  >
	struct Value< IteratorBuffer<TIterator> >:
		Value<TIterator> {};

    template < typename TIterator >
	struct Size< IteratorBuffer<TIterator> >:
		Size<TIterator> {};


	template < typename TIterator >
    inline typename Size< IteratorBuffer<TIterator> >::Type
    length(IteratorBuffer<TIterator> const &me) {
		return me.size;
    }


    //////////////////////////////////////////////////////////////////////////////
	// handler that emulates a buffer by using a container buffer
    struct _SourceNonCachingSpec;
	typedef Tag<_SourceNonCachingSpec> SourceNonCachingSpec;

	template < typename TPipe >
	struct BufferHandler< TPipe, SourceNonCachingSpec >
    {
        typedef typename Source<TPipe>::Type	TSource;
        typedef ContainerBuffer<TSource const>	TBuffer;

		TPipe	&pipe;

		BufferHandler(TPipe &_pipe):
			pipe(_pipe) {}

		template <typename TSize>
		BufferHandler(TPipe &_pipe, TSize):
			pipe(_pipe) {}

        inline TBuffer first() {
            return TBuffer(pipe.in);
        }

        inline TBuffer next() {
            return TBuffer();
        }

        inline void process() {}
        inline void end() {}
        inline void cancel() {}
    };

	template < typename TPipe >
    struct Value< BufferHandler< TPipe, SourceNonCachingSpec > > {
        typedef ContainerBuffer< typename Source<TPipe>::Type const > Type;
    };


    //////////////////////////////////////////////////////////////////////////////
	// caching buffer handler
    struct _SourceCachingSpec;
	typedef Tag<_SourceCachingSpec> SourceCachingSpec;

	template < typename TPipe >
	struct BufferHandler< TPipe, SourceCachingSpec >
    {
        typedef typename Value<TPipe>::Type			TValue;
        typedef typename Size<TPipe>::Type			TSize;

        typedef SimpleBuffer<TValue>						TBuffer;
		typedef IPipeIterator<TPipe>						ISource;
		typedef typename Iterator<TBuffer, Standard>::Type	ITarget;

		TPipe		&pipe;
		size_t		bufferSize;
        TSize		rest;
        TBuffer		buffer;
		ISource		source;

		BufferHandler(TPipe &_pipe, size_t requestedSize):
			pipe(_pipe),
			bufferSize(requestedSize),
            rest(0) {}

        inline TBuffer& first() {
            rest = length(pipe);
			allocPage(buffer, _min(bufferSize, rest), *this);
			source = ISource(pipe);
			for(ITarget target = buffer.begin; target != buffer.end; ++target) {
				*target = *source;
				++source;
			}
            if (!(rest -= size(buffer))) source = ISource();
			return buffer;
        }

        inline TBuffer& next() {
			resize(buffer, _min(bufferSize, rest));
			ITarget _end = buffer.begin + size(buffer);
			for(ITarget target = buffer.begin; target != _end; ++target) {
				*target = *source;
				++source;
			}
            if (!(rest -= size(buffer))) source = ISource();
			return buffer;
        }

        inline void process() {}
        inline void end() { cancel(); }
        inline void cancel() { source = ISource(); freePage(buffer, *this); }
    };

	template < typename TPipe >
    struct Value< BufferHandler< TPipe, SourceCachingSpec > > {
        typedef SimpleBuffer< typename Value<TPipe>::Type > Type;
    };

    //////////////////////////////////////////////////////////////////////////////
	// buffer handler optimized for external string sources
    struct _ExtStringSourceCachingSpec;
	typedef Tag<_ExtStringSourceCachingSpec> ExtStringSourceCachingSpec;

	template < typename TSequence, typename TSpec >
	struct BufferHandler< Pipe<TSequence, TSpec>, ExtStringSourceCachingSpec >
    {
		typedef typename Value<TSequence>::Type		TValue;
        typedef typename Size<TSequence>::Type		TSize;
        typedef SimpleBuffer<TValue>				TBuffer;
        typedef Pipe<TSequence, TSpec>				TPipe;

		typedef typename Iterator<TSequence const, Standard>::Type	ISource;
		typedef typename Iterator<TBuffer, Standard>::Type			ITarget;

		TPipe		&pipe;
		unsigned	bufferSize;
        TSize		rest;
        TBuffer		buffer;
		ISource		source;

		BufferHandler(TPipe &_pipe, unsigned requestedSize):
			pipe(_pipe),
			bufferSize(requestedSize),
            rest(0) {}

        inline TBuffer& first() {
            rest = length(pipe.in);
			allocPage(buffer, _min(bufferSize, rest), *this);
			source = begin(pipe.in);
			for(ITarget target = buffer.begin; target != buffer.end; ++target) {
				*target = *source;
				++source;
			}
            if (!(rest -= size(buffer))) source = ISource();
			return buffer;
        }

        inline TBuffer& next() {
			resize(buffer, _min(bufferSize, rest));
			ITarget _end = buffer.begin + size(buffer);
			for(ITarget target = buffer.begin; target != _end; ++target) {
				*target = *source;
				++source;
			}
            if (!(rest -= size(buffer))) source = ISource();
			return buffer;
        }

        inline void process() {}
        inline void end() { cancel(); }
        inline void cancel() { source = ISource(); freePage(buffer, *this); }
    };

	template < typename TSequence, typename TSpec >
    struct Value< BufferHandler< Pipe<TSequence, TSpec>, ExtStringSourceCachingSpec > > {
		typedef SimpleBuffer< typename Value<TSequence>::Type > Type;
    };

    //////////////////////////////////////////////////////////////////////////////
    // global functions

    // choose the most efficient buffer handler
    template < typename TInput, typename TSpec >
    struct BufReadHandler< Pipe<TInput, TSpec> > {
		typedef			
			typename IF< 
				AllowsFastRandomAccess<TInput>::VALUE,
				BufferHandler< Pipe<TInput, TSpec>, SourceNonCachingSpec>,
//				BufferHandler< Pipe<TInput, TSpec>, SourceCachingSpec>
				BufferHandler< Pipe<TInput, TSpec>, ExtStringSourceCachingSpec>
			>::Type Type;
    };


    template < typename TValue, typename TConfig, typename TSpec >
    struct BufReadHandler< Pipe< String<TValue, External<TConfig> >, TSpec > > {
		typedef BufferHandler< 
			Pipe< String<TValue, External<TConfig> > const, TSpec >, 
			ExtStringSourceCachingSpec 
		> Type;
    };

    template < typename TValue, typename TConfig, typename TSpec >
    struct BufReadHandler< Pipe< String<TValue, External<TConfig> > const, TSpec > > {
		typedef BufferHandler< 
			Pipe< String<TValue, External<TConfig> > const, TSpec >, 
			ExtStringSourceCachingSpec 
		> Type;
    };


	template < typename TInput, typename TSpec, typename TCommand >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, TCommand const &) {
        return true;
    }

	template < typename TInput, typename TSpec >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, ControlBeginRead const &) {
		me.cur = begin(me.in, Rooted());
		return true;
	}
	
	template < typename TInput, typename TSpec >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, ControlEndRead const &) {
        me.cur = end(me.in, Rooted());
		return true;
	}
	
	template < typename TInput, typename TSpec >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, ControlEof const &) {
		return atEnd(me.cur);
	}

	template < typename TInput, typename TSpec >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, ControlEos const &) {
		return atEndOfSequence(me.cur);
	}

    template < typename TInput, typename TSpec >
    inline typename Size< Pipe< TInput, Source<TSpec> > >::Type
	length(Pipe< TInput, Source<TSpec> > const &me) {
        return length(me.in);
    }

}

#endif
