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
  $Id: pool_base.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_POOL_BASE_H
#define SEQAN_HEADER_POOL_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Spec.PoolConfigSize:
..cat:Pipelining
..general:Spec.PoolSpec
..summary:Configuration of Pool.
..signature:PoolConfig<TSize, TFile>
..param.TSize:The Pool's size type.
..param.TFile:The underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..see:Spec.PoolConfig
*/

	template < typename TSize,
		       typename TFile = File<> >						// default file type
    struct PoolConfigSize {
		typedef TSize SizeType;
        typedef TFile File;
    };

/**
.Spec.PoolConfig:
..cat:Pipelining
..general:Spec.PoolSpec
..summary:Configuration of Pool.
..signature:PoolConfig<TFile>
..param.TFile:The underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:Using this configuration spec., the Pool's size type is $Size<TFile>::Type$. To use a custom size type @Spec.PoolConfigSize@ should be used.
..see:Spec.PoolConfigSize
*/

	template < typename TFile = File<> >						// default file type
    struct PoolConfig {
		typedef typename Size<TFile>::Type SizeType;
        typedef TFile File;
    };

/**
.Spec.PoolSpec:
..cat:Pipelining
..general:Class.Pool
..summary:Stores/Retrieves all elements to/from disk.
..signature:Pool<TValue, PoolSpec<TConfig> >
..param.TValue:The value type, that is the type of the stream elements.
..param.TConfig:Configuration Spec. Defines destination function, size type, and file type.
...type:Spec.PoolConfig
...type:Spec.PoolConfigSize
..remarks:The Pool's input/output type is $TValue$ and the size type is determined by the $TConfig$.
*/

    template < typename TConfig = PoolConfig<> >
    struct PoolSpec {
        typedef TConfig Config;
    };

/**
.Class.Pool:
..cat:Pipelining
..summary:Pools are push- and pop-passive pipeline modules.
..signature:Pool<TValue, TSpec>
..param.TValue:The value type, that is the type of the stream elements.
..param.TSpec:The specializing type.
...default:PoolSpec<>, see @Spec.PoolSpec@.
..remarks:Use @Metafunction.Value@ to get the output type of a given Pipe (returns $Value<TInput>::Type$ by default).
..remarks:Use @Metafunction.Size@ to get the size type of a given Pipe (returns $Size<TInput>::Type$ by default).
*/

    template < typename TValue,
               typename TSpec = PoolSpec<> >
    struct Pool;


    struct PoolParameters
    {
/*
        enum { DefaultMemBufferSize     = 384 * 1024*1024,			// max memory config
               DefaultPageSize          = 32 * 1024*1024,
               DefaultBucketBufferSize  = 64 * 1024*1024,
               DefaultReadAheadBuffers  = 2,
               DefaultWriteBackBuffers  = 2,
               DefaultWriteBackBuckets  = 16,
               DefaultAbsoluteSizes     = true };
*/


        enum { DefaultMemBufferSize     = 128 * 1024*1024,		// normal memory config
               DefaultPageSize          = 32 * 1024*1024,
               DefaultBucketBufferSize  = 64 * 1024*1024,
               DefaultReadAheadBuffers  = 2,
               DefaultWriteBackBuffers  = 2,
               DefaultWriteBackBuckets  = 16,
               DefaultAbsoluteSizes     = true };


/*
        enum { DefaultMemBufferSize     = 0*8192,//64 * 1024*1024,	// low memory config
               DefaultPageSize          = 2 * 1024*1024,
               DefaultBucketBufferSize  = 6 * 1024*1024,
               DefaultReadAheadBuffers  = 4,
               DefaultWriteBackBuffers  = 4,
               DefaultWriteBackBuckets  = 16,
               DefaultAbsoluteSizes     = true };
*/

        unsigned memBufferSize;
        unsigned pageSize;
        unsigned bucketBufferSize;
        unsigned readAheadBuffers;
        unsigned writeBackBuffers;
        unsigned writeBackBuckets;
        bool     absoluteSizes;     // when false, sizes are measured in units of TValue
                                    // when true, sizes are measured in bytes

        PoolParameters():
            memBufferSize(DefaultMemBufferSize),
            pageSize(DefaultPageSize),
            bucketBufferSize(DefaultBucketBufferSize),
            readAheadBuffers(DefaultReadAheadBuffers),
            writeBackBuffers(DefaultWriteBackBuffers),
            writeBackBuckets(DefaultWriteBackBuckets),
            absoluteSizes(DefaultAbsoluteSizes) {}

        template < typename TValue >
        void absolutize(unsigned aligning, TValue *) {
            if (!absoluteSizes) return;
            memBufferSize = memBufferSize / sizeof(TValue);
            bucketBufferSize = bucketBufferSize / sizeof(TValue);
            pageSize = pageSize / sizeof(TValue);
            pageSize = (pageSize / aligning) * aligning;
            if (!pageSize) pageSize += aligning;
        }
    };


        
    template < typename TBuffer, typename THandler >
    inline TBuffer& processBuffer(TBuffer &h, THandler &) {
        return h;
    }


    //////////////////////////////////////////////////////////////////////////////
	// handler that manages a simple memory buffer
    struct _MemorySpec;
	typedef Tag<_MemorySpec> MemorySpec;

    template < typename TPool >
    struct BufferHandler< TPool, MemorySpec >
    {
		typedef TPool							Pool;
        typedef typename TPool::Type            Type;

        typedef SimpleBuffer<Type>				Buffer;
        typedef Buffer&							BufferRef;

		TPool	&pool;
        Buffer	empty;

		BufferHandler(TPool &_pool):
			pool(_pool) {}

		BufferHandler(TPool &_pool, unsigned):
			pool(_pool) {}

        inline BufferRef first() {
            return pool.memBuffer;
        }

        inline BufferRef next() {
            return empty;
        }

        inline void process() {
            processBuffer(pool.memBuffer, *this);
        }

        inline void end() {}
        inline void cancel() {}
    };


    //////////////////////////////////////////////////////////////////////////////
	// generic block based asynchronous read handler
    struct _ReadFileSpec;
	typedef Tag<_ReadFileSpec> ReadFileSpec;

    template < typename TPool >
    struct BufferHandler< TPool, ReadFileSpec >
    {
		typedef TPool								Pool;
        typedef typename TPool::Type                Type;
        typedef typename TPool::File                File;

        typedef SimpleBuffer<Type>					Buffer;
        typedef PageFrame<Type, File, Dynamic<> >   PageFrame;
		typedef PageFrame&                          PageFrameRef;
        typedef Buffer&								BufferRef;
        typedef PageChain<PageFrame>	            PageChain;

        TPool       &pool;
        PageChain   chain;
        unsigned    pageSize;
        unsigned    readPageNo, _pages;
        Buffer	    empty;

        BufferHandler(TPool &_pool):
            pool(_pool),
            chain(_min(_pool.readAheadBuffers, _pool.pages())),
			pageSize(_pool.pageSize) {}

        BufferHandler(TPool &_pool, unsigned _requestedBufferSize, unsigned _readAheadBuffers = 1):
            pool(_pool),
//            pageSize(alignSize(_min(_pool.size(), _requestedBufferSize), _pool.pageSize)),
            chain(_min(_readAheadBuffers, _pool.pages(pageSize = alignSize(_min(_pool.size(), _requestedBufferSize), _pool.pageSize)))) 
        {
			#ifdef SEQAN_HEADER_PIPE_DEBUG
				::std::cerr << "___BufferHandler___" << ::std::endl;
				::std::cerr << "pagesize: " << pageSize << ::std::endl;
				::std::cerr << "readaheadbuffers: " << chain.maxFrames << ::std::endl;
				::std::cerr << "pages: " << pool.pages(pageSize) << ::std::endl;
			#endif
        }

        ~BufferHandler() {
            end();
        }
        
        inline BufferRef first() {
            _pages = pool.pages(pageSize);
            if (!_pages) return empty;

            // enqueue reading of the <readAheadBuffers> first blocks
            readPageNo = 0;
            PageFrame *p = chain.first;
            while (p) {
                p->pageNo = readPageNo++;
                _read(*p);
                p = p->next;
            }

            // retrieve the very first
            waitFor(*chain.first);
            return processBuffer(*chain.first, *this);
        }

        inline BufferRef next() {
            // step one buffer ahead
            chain.getReadyPage();

            // read ahead
            chain.last->pageNo = readPageNo++;
            _read(*chain.last);
            
            // retrieve the next buffer in order
            waitFor(*chain.first);
            return processBuffer(*chain.first, *this);
        }

        inline void end() {
            cancel();
        }

        inline void cancel() {
            PageFrame *p = chain.first;
            while (p) {
				::seqan::cancel(*p, pool.file);
                freePage(*p, pool.file);
                p = p->next;
            }
        }

        inline void process() {}

    private:
		bool _error() {
//            ::std::cerr << "Error in BufWriteFileHandler::_read " << pool.file.error() << ::std::endl;
            return true;
        }

        inline bool _read(PageFrame &pf) {
            if (pf.pageNo < _pages) {
                // alloc if empty
                if (!pf.begin)
                    allocPage(pf, pageSize, pool.file);

                // set buffer size according to read size
                resize(pf, pool.dataSize(pf.pageNo, pageSize));

                // read asynchronously (if possible) from disk
                return readPage(pf, pool.file) || _error();
            } else {
                // free if allocated
                freePage(pf, pool.file);
                return false;
            }
        }
    };


	//////////////////////////////////////////////////////////////////////////////
	// generic block based asynchronous write handler
    struct _WriteFileSpec;
	typedef Tag<_WriteFileSpec> WriteFileSpec;

    template < typename TPool >
    struct BufferHandler< TPool, WriteFileSpec >
    {
		typedef TPool								Pool;
        typedef typename TPool::Type                Type;
        typedef typename TPool::File                File;

        typedef SimpleBuffer<Type>					Buffer;
        typedef PageFrame<Type, File, Dynamic<> >   PageFrame;
		typedef PageFrame&                          PageFrameRef;
        typedef Buffer&								BufferRef;
        typedef PageChain<PageFrame>	            PageChain;

        TPool       &pool;
        PageChain   chain;
        unsigned    pageSize;
        unsigned    writePageNo, _pages;
        Buffer	    empty;

        BufferHandler(TPool &_pool):
            pool(_pool),
            chain(_min(_pool.writeBackBuffers, _pool.pages())),
			pageSize(_pool.pageSize) {}

        BufferHandler(TPool &_pool, unsigned _requestedBufferSize, unsigned _writeBackBuffers = 1):
            pool(_pool),
            chain(_min(_writeBackBuffers, _pool.pages(pageSize))),
			pageSize(alignSize(_min(_pool.size(), _requestedBufferSize), _pool.pageSize)) {}

        ~BufferHandler() {
            cancel();
        }
        
    public:

        inline BufferRef first() {
            _pages = pool.pages(pageSize);
            if (!_pages) return empty;

            // get a ready page frame
            PageFrameRef pf = *chain.getReadyPage();

            if (!pf.begin)
                allocPage(pf, pageSize, pool.file);

            writePageNo = 0;
            resize(pf, pool.dataSize(pf.pageNo = writePageNo++, pageSize));
            return pf;
        }

        inline BufferRef next() {
            // write previously provided buffer to disk
            _write(processBuffer(*chain.last, *this));

            // step one buffer ahead
            PageFrameRef pf = *chain.getReadyPage();

            if (!pf.begin)
                allocPage(pf, pageSize, pool.file);

            resize(pf, pool.dataSize(pf.pageNo = writePageNo++, pageSize));
            return pf;
        }

        inline void end() {
            // write previously provided buffer to disk
            _write(processBuffer(*chain.last, *this));
            flush();
        }

        inline void flush() {
            PageFrame *p = chain.first;
            while (p) {
                waitFor(*p);
                freePage(*p, pool.file);
                p = p->next;
            }
			::seqan::flush(pool.file);
        }

        inline void cancel() {
            PageFrame *p = chain.first;
            while (p) {
                ::seqan::cancel(*p, pool.file);
                freePage(*p, pool.file);
                p = p->next;
            }
        }

        inline void process() {}

    private:

        bool _error() {
//            ::std::cerr << "Error in BufWriteFileHandler::_write " << pool.file.error() << ::std::endl;
            return true;
        }

        bool _write(PageFrame &pf) {
            if (pf.pageNo < _pages) {
                // write asynchronously (if possible) to disk
                return writePage(pf, pool.file) || _error();
            } else {
                // free if allocated
                freePage(pf, pool.file);
                return false;
            }
        }
    };

	//////////////////////////////////////////////////////////////////////////////
	// generic buffered multiplex handler
    struct _MultiplexSpec;
	typedef Tag<_MultiplexSpec> MultiplexSpec;

	template < typename TBufferHandler1, typename TBufferHandler2 >
    struct BufferHandler< Bundle2< TBufferHandler1, TBufferHandler2 >, MultiplexSpec >
	{
        typedef typename TBufferHandler1::Type   Type;
        typedef typename TBufferHandler1::Buffer Buffer;

		TBufferHandler1 *handler1;
		TBufferHandler2 *handler2;

		template <typename TPool>
		BufferHandler(TPool &_pool) {
			if (_pool.memBuffer.begin || _pool._size == 0) {	// decision to choose handler1 or handler2
				handler1 = new TBufferHandler1(_pool);
				handler2 = NULL;
			} else {
				handler1 = NULL;
				handler2 = new TBufferHandler2(_pool);
			}
		}

		template <typename TPool>
		BufferHandler(TPool &_pool, unsigned _requestedBufferSize) {
			if (_pool.memBuffer.begin || _pool._size == 0) {	// decision to choose handler1 or handler2
				handler1 = new TBufferHandler1(_pool, _requestedBufferSize);
				handler2 = NULL;
			} else {
				handler1 = NULL;
				handler2 = new TBufferHandler2(_pool, _requestedBufferSize);
			}
		}

		~BufferHandler() {
			delete handler1;
			delete handler2;
		}
		
        inline Buffer first() {
			if (handler1)	return handler1->first();
			else			return handler2->first();
		}

		inline Buffer next() {
			if (handler1)	return handler1->next();
			else			return handler2->next();
		}

		inline void end() {
			if (handler1)	return handler1->end();
			else			return handler2->end();
		}

		inline void process() {
			if (handler1)	return handler1->process();
			else			return handler2->process();
		}

		inline void cancel() {
			if (handler1)	return handler1->cancel();
			else			return handler2->cancel();
		}
	};


	template < typename THandler1, typename THandler2 >
    struct Handler< Bundle2< THandler1, THandler2 >, MultiplexSpec >
	{
        typedef typename THandler1::Type            Type;

        THandler1 *handler1;
		THandler2 *handler2;

		template <typename TPool>
		Handler(TPool &_pool) {
			if (_pool.memBuffer.begin || _pool._size == 0) {	// decision to choose handler1 or handler2
				handler1 = new THandler1(_pool);
				handler2 = NULL;
			} else {
				handler1 = NULL;
				handler2 = new THandler2(_pool);
			}
		}

		~Handler() {
			delete handler1;
			delete handler2;
		}
		
        inline bool begin() {
			if (handler1)	return handler1->begin();
			else			return handler2->begin();
		}

        inline Type const & front() const {
            if (handler1)	return handler1->front();
			else			return handler2->front();
        }

        inline void pop() {
            if (handler1)	handler1->pop();
			else			handler2->pop();
        }

        inline void pop(Type &_Ref) {
            if (handler1)	handler1->pop(_Ref);
			else			handler2->pop(_Ref);
        }

        inline void push(Type const & _Val) {
            if (handler1)	handler1->push(_Val);
			else			handler2->push(_Val);
        }

        inline bool eof() const {
            if (handler1)	return handler1->eof();
			else			return handler2->eof();
        }

		inline void end() {
			if (handler1)	return handler1->end();
			else			return handler2->end();
		}

		inline void process() {
			if (handler1)	return handler1->process();
			else			return handler2->process();
		}
	};



	//////////////////////////////////////////////////////////////////////////////
	// character and buffer based handler definitions
    template < typename TValue, typename TSpec >
    struct BufReadHandler< Pool< TValue, TSpec > >
    {
        typedef BufferHandler< Bundle2<
			BufferHandler< Pool< TValue, TSpec >, MemorySpec >, 
			BufferHandler< Pool< TValue, TSpec >, ReadFileSpec >
		>, MultiplexSpec > Type;
    };

    template < typename TValue, typename TSpec >
    struct BufWriteHandler< Pool< TValue, TSpec > >
    {
        typedef BufferHandler< Bundle2<
			BufferHandler< Pool< TValue, TSpec >, MemorySpec >, 
			BufferHandler< Pool< TValue, TSpec >, WriteFileSpec >
		>, MultiplexSpec > Type;
    };


    template < typename TValue, typename TPoolSpec, typename TSpec >
    struct Value< BufferHandler< Pool< TValue, TPoolSpec >, TSpec > > {
        typedef SimpleBuffer< TValue > Type;
    };

    template < typename TPool >
    struct HandlerArgs {
		typedef Nothing Type;
	};


    //////////////////////////////////////////////////////////////////////////////
	// base class of all pool classes like Pool, Mapper, Sorter

    template < typename TValue,
               typename TSpec >
    struct Pool
    {
		typedef TValue	                    Type;
        typedef TSpec                       Spec;
        typedef typename TSpec::Config      Config;
        typedef typename Config::File       File;
        typedef typename Config::SizeType   SizeType;
        typedef SimpleBuffer<TValue>        Buffer;

		// public handlers to read simultanously buffer- or character-wise
        typedef typename ReadHandler<Pool>::Type    ReadHandler;
        typedef typename WriteHandler<Pool>::Type   WriteHandler;
		typedef typename HandlerArgs<Pool>::Type	HandlerArgs;

		File				file;
        bool                _temporary, _ownFile;
		SizeType			_size;
        unsigned            _pages;
        unsigned            pageSize;
        unsigned            bucketBufferSize;
        unsigned            readAheadBuffers;
        unsigned            writeBackBuffers;
        unsigned            writeBackBuckets;

        Buffer				memBuffer;
        unsigned            memBufferSize;
        HandlerArgs         handlerArgs;

		bool				_partiallyFilled;		// the pool is partially filled (it contains undefined values)
		TValue				undefinedValue;			// value to represent undefined (unwritten) entries

    protected:
        unsigned            _lastPageNo;
        unsigned            _lastPageSize;

        int                 listeners;
       
        ReadHandler         *reader;
        WriteHandler        *writer;


    public:

        //////////////////////////////////////////////////////////////////////////////
		// public iterator types

		typedef SizeType				size_type;

        Pool(const PoolParameters &_conf = PoolParameters()):
            file(NULL)
        {
			_init(_conf);
            _setSize(0);
        }
        
		Pool(HandlerArgs const &args, PoolParameters const &_conf = PoolParameters()):
            file(NULL),
			handlerArgs(args)
        {
			_init(_conf);
            _setSize(0);
        }
        
        template < typename TInput, typename TPipeSpec >
        Pool(Pipe<TInput, TPipeSpec> &, const PoolParameters &_conf = PoolParameters()):
            file(NULL)
        {
			_init(_conf);
            _setSize(0);
        }
        
        template < typename TInput, typename TPipeSpec >
		Pool(Pipe<TInput, TPipeSpec> &, HandlerArgs const &args, PoolParameters const &_conf = PoolParameters()):
            file(NULL),
			handlerArgs(args)
        {
			_init(_conf);
            _setSize(0);
        }
        
        Pool(File &_file, const PoolParameters &_conf = PoolParameters()):
            file(_file)
        {
			_init(_conf);
            _ownFile = false;
            _temporary = false;
            memBufferSize = 0;
			_setSize(::seqan::size(file) / sizeof(TValue));
        }
        
        Pool(const char *fileName, const PoolParameters &_conf = PoolParameters())
        {
			_init(_conf);
            _temporary = false;
            memBufferSize = 0;
            if (_ownFile = open(file, fileName))
                _setSize(::seqan::size(file) / sizeof(TValue));
            else
                _setSize(0);
        }
        
        ~Pool() {
            endRead();
            endWrite();
            if (_temporary) 
                clear();
            else
                if (_ownFile) close(file);
        }
        
        inline void clear() {
            resize(0);
        }

		inline size_type size() const {
			return _size;
		}

        // this is not a real resize
        void resize(size_type _newSize) {
			typedef typename Size<File>::Type TFSize;
            if (_newSize == _size) return;

            _freeHandlers();	// if you forgot to call endRead/endWrite we have no trouble

			if (_temporary && _ownFile) {
				if (_size != 0) {
					if (memBuffer.begin)
						freePage(memBuffer, *this);
                    else {
						close(file);
                        SEQAN_PROSUB(SEQAN_PROIOVOLUME, (_proFloat)((TFSize)_size * (TFSize)sizeof(TValue)));
                    }
				}

				if (_newSize != 0) {
					if (_newSize <= (size_type)memBufferSize)
						allocPage(memBuffer, _newSize, *this);
                    else {
						openTemp(file);
                        SEQAN_PROADD(SEQAN_PROIOVOLUME, (_proFloat)((TFSize)_newSize * (TFSize)sizeof(TValue)));
                    }
				}
			}

            _setSize(_newSize);
        }


        //////////////////////////////////////////////////////////////////////////////
		// auto-disposal interface (deprecated)

        inline void addListener() {
            if (!listeners) return;
            ++listeners;
        }

        inline void delListener() {
            if (!listeners) return;
            if (--listeners == 0) {
                clear();
            }
        }


        //////////////////////////////////////////////////////////////////////////////
		// queue interface

        inline Type const & front() const {
            return reader->front();
        }

		inline Type const & operator*() const {
			return reader->front();
		}

        inline void pop() {
			reader->pop();
        }

        inline Pool& operator++() {
			reader->pop();
			return *this;
        }

        inline void pop(Type &_Ref) {
            reader->pop(_Ref);
        }

        inline void push(Type const &_Val) {
            writer->push(_Val);
        }

        inline bool eof() {
            if (reader) return reader->eof();
            if (writer) return writer->eof();
            return true;
        }


        //////////////////////////////////////////////////////////////////////////////
		// flow control

        bool beginWrite() {
            _freeHandlers(); // if you forgot to call endRead/endWrite we have no trouble
			return (
                (writer = new WriteHandler(*this)) && 
                writer->begin());
        }

        bool endWrite() {
			if (writer) {
				writer->end();
				writer->process();
			}
            delete writer;
            writer = NULL;
            return true;
        }

        bool beginRead() {
            _freeHandlers(); // if you forgot to call endRead/endWrite we have no trouble
			return (
                (reader = new ReadHandler(*this)) &&
                reader->begin());
        }

        bool endRead() {
            if (reader) reader->end();
            delete reader;
            reader = NULL;
            delListener();
            return true;
        }


        inline unsigned pages() const {
            return _pages;
        }

        inline unsigned pages(unsigned __pageSize) const {
            return enclosingBlocks(_size, (unsigned)__pageSize);
        }

        // used by buffer handlers ...
        inline unsigned dataSize(unsigned __pageNo) const {
            return (__pageNo != _lastPageNo)? pageSize: _lastPageSize;
        }

        // used by buffer handlers with variable PageSize ...
        inline unsigned dataSize(unsigned __pageNo, unsigned __pageSize) const {
            return (__pageNo != _size / __pageSize)? __pageSize: _size % __pageSize;
        }

    protected:

        inline void _freeHandlers() {
            if (reader) reader->end();
            if (writer) writer->end();
            delete reader;
            delete writer;
            reader = NULL;
            writer = NULL;
        }

        inline void _setSize(size_type _newSize) {
            _size = _newSize;
            _pages = enclosingBlocks(_size, (unsigned)pageSize);
            _lastPageNo = _size / pageSize;
            _lastPageSize = _size % pageSize;
        }

		void _init(PoolParameters _conf = PoolParameters()) {
            _conf.absolutize(sectorSize(file), (Type*)NULL);
            memBufferSize    = _conf.memBufferSize;
            pageSize		 = _conf.pageSize;
            bucketBufferSize = _conf.bucketBufferSize;
            readAheadBuffers = _conf.readAheadBuffers;
            writeBackBuffers = _conf.writeBackBuffers;
            writeBackBuckets = _conf.writeBackBuffers;
            listeners = 0;
            reader = NULL;
            writer = NULL;
            _ownFile = true;
            _temporary = true;
		}

    };

    template < typename TValue, typename TSpec, typename TIteratorSpec>
	struct Iterator< Pool< TValue, TSpec > const, TIteratorSpec> {
		typedef IPipeIterator< Pool< TValue, TSpec > > Type;
	};

    template < typename TValue, typename TSpec, typename TIteratorSpec>
	struct Iterator< Pool< TValue, TSpec >, TIteratorSpec> {
		typedef OPipeIterator< Pool< TValue, TSpec > > Type;
	};

    template < typename TValue, typename TSpec >
	OPipeIterator< Pool< TValue, TSpec > >
	begin(Pool< TValue, TSpec > &pool) {
		return OPipeIterator< Pool< TValue, TSpec > >(pool);
	}

    template < typename TValue, typename TSpec >
	OPipeIterator< Pool< TValue, TSpec > >
	end(Pool< TValue, TSpec > &pool) {
		return OPipeIterator< Pool< TValue, TSpec > >();
	}



    //////////////////////////////////////////////////////////////////////////////
	// a pool is a queue-like container for a large amount of data
    // in contrast to a queue you can read the whole content more than once
    // but can't pop directly after a push
    // instead, access to a pool looks like the following:
    //
    //  1. resize(<new size>)
    //  2. beginWrite()
    //  3. push(a), push(b), push(c) ...
    //  4. endWrite()
    // (5. do something else)
    //  6. beginRead()
    //  7. pop(a), pop(b), pop(c) ...
    //  8. endRead()
    // (9. clear() - if you want to save memory and refill the pool later)
    //
    // you can repeat steps 2-4 and 6-8 independently as often as you like

    //template < typename TValue, typename TFile = File<> >
    //struct Pool: public Pool< TValue, PoolSpec< PoolConfig< TFile > > >
    //{
    //    typedef TValue	                    Type;
    //    typedef TFile                       File;
    //    typedef typename Size<TFile>::Type  SizeType;
    //};


	// seqan namespace traits
    template < typename TValue, typename TSpec >
    struct Value< Pool< TValue, TSpec > > {
        typedef TValue Type;
    };

    template < typename TValue, typename TSpec >
    struct Size< Pool< TValue, TSpec > > {
        typedef typename Pool< TValue, TSpec >::SizeType Type;
    };

    template < typename TValue, typename TSpec >
    struct Position< Pool< TValue, TSpec > > {
        typedef typename Size<Pool< TValue, TSpec > >::Type Type;
    };

    template < typename TValue, typename TSpec >
    struct Difference< Pool< TValue, TSpec > > {
		typedef typename _MakeSigned<typename Size<Pool< TValue, TSpec > >::Type>::Type Type;
    };


///.Function.clear.param.object.type:Class.Pool

    template < typename TValue, typename TSpec >
    inline void clear(Pool<TValue, TSpec> &me)
    {
        return me.clear();
    }

	// deprecated
    template < typename TValue, typename TSpec >
    inline typename Size< Pool<TValue, TSpec> >::Type
    size(Pool<TValue, TSpec> const &me)
    {
        return me.size();
    }

///.Function.length.param.object.type:Class.Pool

	template < typename TValue, typename TSpec >
    inline typename Size< Pool<TValue, TSpec> >::Type
    length(Pool<TValue, TSpec> const &me)
    {
        return me.size();
    }

///.Function.resize.param.object.type:Class.Pool

	template < typename TValue, typename TSpec, typename TSize >
    inline TSize resize(Pool<TValue, TSpec> &me, TSize new_length)
    {
	    me.resize(new_length);
        return me.size();
    }

///.Function.Pipelining#front.param.object.type:Class.Pool

    template < typename TValue, typename TSpec >
	inline typename Value< Pool<TValue, TSpec> >::Type const & front(Pool<TValue, TSpec> &me) {
        return me.front();
    }

///.Function.pop.param.object.type:Class.Pool

    template < typename TValue, typename TSpec >
    inline void pop(Pool<TValue, TSpec> &me) {
        me.pop();
    }

    template < typename TValue, typename TSpec >
    inline void pop(Pool<TValue, TSpec> &me, TValue &_Ref) {
        me.pop(_Ref);
    }

/**
.Function.push:
..cat:Pipelining
..summary:Appends an item at the end of an input stream.
..signature:push(object, val)
..param.object:A push-passive pipeline module.
...type:Class.Pool
..param.val:Item to be pushed.
..remarks:@Function.push@ can only be called within a write process surrounded by @Function.beginWrite@ and @Function.endWrite@.
*/

    template < typename TValue, typename TSpec >
    inline void push(Pool<TValue, TSpec> &me, TValue const &_Val) {
        me.push(_Val);
    }

    template < typename TValue, typename TSpec >
    ::std::ostream& operator<<(::std::ostream &out, Pool<TValue, TSpec> &p) {
        beginRead(p);
        while (!eof(p)) {
		    out << front(p) << ::std::endl;
            pop(p);
        }
        endRead(p);
		return out;
	}


	// the pipe interface of pool classes
	//namespace SEQAN_NAMESPACE_PIPELINING
	//{
		//template < typename TValue, typename TSpec >
		//struct Value< Pool< TValue, TSpec > >
		//{
		//	typedef TValue Type;
		//};

		//template < typename TValue, typename TSpec >
		//struct Size< Pool< TValue, TSpec > >
		//{
		//	typedef typename Size< TFile >::Type Type;
		//};

		template < typename TValue, typename TSpec >
	    inline bool control(Pool< TValue, TSpec > &me, ControlEof const &) {
		    return me.eof();
	    }
    	
		template < typename TValue, typename TSpec >
	    inline bool control(Pool< TValue, TSpec > &me, ControlEos const &) {
		    return me.eof();
	    }
    	
		template < typename TValue, typename TSpec >
	    inline bool control(Pool< TValue, TSpec > &me, ControlClear const &) {
		    me.clear();
		    return true;
	    }
    	
		template < typename TValue, typename TSpec >
	    inline bool control(Pool< TValue, TSpec > &me, ControlBeginRead const &) {
		    return me.beginRead();
	    }
    	
		template < typename TValue, typename TSpec >
	    inline bool control(Pool< TValue, TSpec > &me, ControlEndRead const &) {
		    return me.endRead();
	    }
    	
/**
.Function.beginWrite:
..cat:Pipelining
..summary:Initiates a write process.
..signature:beginWrite(object)
..param.object:A push-passive pipeline module.
...type:Class.Pool
..returns:A $bool$ which is $true$ on success.
..remarks:$beginWrite$ prepares a @Class.Pool@ for succeeding writes.
..remarks:A write process must be terminated with @Function.endWrite@. Nested write processes are not allowed.
..see:Function.endWrite
*/

		template < typename TValue, typename TSpec >
	    inline bool beginWrite(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
		    return me.beginWrite();
        }

/**
.Function.endWrite:
..cat:Pipelining
..summary:Terminates a write process.
..signature:endWrite(object)
..param.object:A push-passive pipeline module.
...type:Class.Pool
..returns:A $bool$ which is $true$ on success.
..remarks:$endWrite$ closes the input stream and frees resources possibly allocated by @Function.beginWrite@.
..see:Function.beginWrite
*/

		template < typename TValue, typename TSpec >
	    inline bool endWrite(Pool< TValue, TSpec > &me) {
		    return me.endWrite();
        }

///.Function.atEnd.param.iterator.type:Class.Pool

		template < typename TValue, typename TSpec >
	    inline bool eof(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
            return control(me, ControlEof());
        }

		template < typename TValue, typename TSpec >
	    inline bool beginRead(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
            return control(me, ControlBeginRead());
        }

		template < typename TValue, typename TSpec >
	    inline bool endRead(Pool< TValue, TSpec > &me) {
SEQAN_CHECKPOINT
            return control(me, ControlEndRead());
        }

	//}

    // pipe/pool -> pool
    template < typename TValue,
               typename TSpec,
               typename TSource >
    inline bool append(Pool<TValue, TSpec> &dest, TSource &src) {
        typename Size<TSource>::Type leftToRead = length(src);
        if (!beginRead(src)) return false;
        while (leftToRead) {
            push(dest, *src);
            ++src;
            --leftToRead;
        }
        endRead(src);
        return true;
    }

    // string -> pool
    template < typename TValue,
               typename TSpec,
               typename TStringSpec >
    inline bool append(Pool<TValue, TSpec> &dest, String<TValue, TStringSpec> &src) {
        typedef typename Iterator< String<TValue, TStringSpec> const, Standard >::Type TIter;
        TIter _cur = begin(src, Standard());
		TIter _end = end(src, Standard());
        while (_cur != _end) {
            push(dest, *_cur);
            ++_cur;
        }
        return true;
    }

///.Function.assign.param.target.type:Class.Pool

    template < typename TValue,
               typename TSpec,
               typename TSource >
    inline bool assign(Pool<TValue, TSpec> &dest, TSource &src) {
        typename Size<TSource>::Type _size = length(src);
        resize(dest, _size);
        return beginWrite(dest) && append(dest, src) && endWrite(dest);
    }

    template < typename TValue,
               typename TSpec,
               typename TSource >
    inline bool operator<<(Pool<TValue, TSpec> &dest, TSource &src) {
        return assign(dest, src);
    }



///.Function.assign.param.source.type:Class.Pool

    // pool -> string
    template < typename TValue1,
			   typename TStringSpec,
			   typename TValue2,
			   typename TSpec >
    inline bool assign(String<TValue1, TStringSpec> &dest, Pool<TValue2, TSpec> &src) {
        typedef typename Iterator< String<TValue1, TStringSpec>, Standard >::Type TIter;
        typename Size< String<TValue1, TStringSpec> >::Type _size = length(src);
        resize(dest, _size);
        if (!beginRead(src)) return false;
        TIter _cur = begin(dest, Standard());
		TIter _end = end(dest, Standard());
        while (_cur != _end) {
            *_cur = *src;
            ++_cur;
            ++src;
        }
        endRead(src);
        return true;
    }

    template < typename TValue1,
			   typename TStringSpec,
			   typename TValue2,
			   typename TSpec >
    inline bool operator<<(String<TValue1, TStringSpec> &dest, Pool<TValue2, TSpec> &src) {
        return assign(dest, src);
    }

}

#endif
