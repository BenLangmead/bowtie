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
  $Id: pool_sorter.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_POOL_SORTER_H
#define SEQAN_HEADER_POOL_SORTER_H

namespace SEQAN_NAMESPACE_MAIN
{

    template < typename TValue, typename Compare >
    struct MergeStreamComparer : public ::std::binary_function < PageBucket<TValue>,
														       PageBucket<TValue>,
														       bool>
    {
        Compare C;
        MergeStreamComparer(Compare &CC): C(CC) { }
        inline bool operator()(const PageBucket<TValue> &a,
			                   const PageBucket<TValue> &b) const {
            return C(static_cast<const TValue&>(*(a.cur)),
                     static_cast<const TValue&>(*(b.cur))) > 0;
        }
    };

    template < typename TCompare >
	struct AdaptorCompare2Less : 
		public ::std::binary_function <
			typename TCompare::first_argument_type, 
			typename TCompare::second_argument_type, 
			bool >
    {
        TCompare C;
        AdaptorCompare2Less(TCompare &CC): C(CC) { }
        inline bool operator() (
			typename TCompare::first_argument_type const &a, 
			typename TCompare::second_argument_type const &b) const 
		{
            return C(a, b) < 0;
        }
    };

/**
.Spec.SorterConfigSize:
..cat:Pipelining
..general:Spec.SorterSpec
..summary:Configuration of Sorter.
..signature:SorterConfigSize<TCompare, TSize, TFile>
..param.TCompare:The compare function (see STL's $binary_function$).
...remarks:Let $comp$ be an object of type $TCompare$. $comp(a,b)$ should return a value less, equal, or greater than $0$ if $a<b$, $a==b$, or $a>b$.
..param.TSize:The Sorter's size type.
..param.TFile:The underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..see:Spec.SorterConfig
*/

	template < typename TCompare,
		       typename TSize,
			   typename TFile = File<> >
    struct SorterConfigSize {
        typedef TCompare	Compare;
		typedef TSize		SizeType;
        typedef TFile		File;
    };



/**
.Spec.SorterConfig:
..cat:Pipelining
..general:Spec.SorterSpec
..summary:Configuration of Sorter.
..signature:SorterConfig<TCompare, TFile>
..param.TCompare:The compare function (see STL's $binary_function$).
...remarks:Let $comp$ be an object of type $TCompare$. $comp(a,b)$ should return a value less, equal, or greater than $0$ if $a<b$, $a==b$, or $a>b$.
..param.TFile:The underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:Using this configuration spec., the Sorter's size type is $Size<TFile>::Type$. To use a custom size type @Spec.SorterConfigSize@ should be used.
..see:Spec.SorterConfigSize
*/

    template < typename TCompare, typename TFile = File<> >
    struct SorterConfig {
        typedef TCompare					Compare;
		typedef typename Size<TFile>::Type	SizeType;
        typedef TFile						File;
    };

/**
.Spec.SorterSpec:
..cat:Pipelining
..general:Class.Pool
..summary:Sorts all elements using a custom compare function.
..signature:Pool<TValue, SorterSpec<TConfig> >
..param.TValue:The value type, that is the type of the stream elements.
..param.TConfig:Configuration Spec. Defines compare function, size type, and file type.
...type:Spec.SorterConfig
...type:Spec.SorterConfigSize
..remarks:The Pool's input/output type is $TValue$ and the size type is determined by the $TConfig$.
..note:If there exists an explicit function mapping input elements to their destined positions in the output stream, @Spec.MapperSpec@ should be preferred.
*/

    template < typename TConfig >
    struct SorterSpec {
        typedef TConfig Config;
    };

    template < typename TValue,
			   typename TConfig >
    struct HandlerArgs< Pool< TValue, SorterSpec<TConfig> > >
    {
        typedef typename TConfig::Compare Type;
    };


    //////////////////////////////////////////////////////////////////////////////
	// cache bucket based synchronous multiway merge
    struct _ReadSorterSpec;
	typedef Tag<_ReadSorterSpec> ReadSorterSpec;

    template < typename TPool >
    struct Handler< TPool, ReadSorterSpec >
    {
        typedef typename TPool::Type                    Type;
        typedef typename TPool::Buffer					Buffer;
        typedef typename TPool::Spec::Config::Compare   Compare;
        typedef PageBucketExtended<Type>                PageBucket;

        typedef MergeStreamComparer<Type, Compare>      StreamComparer;
/*        typedef ::std::priority_queue <
            PageBucket,
            ::std::vector<PageBucket>,
            MergeStreamComparer<Type, Compare> >	    PQueue;
*/
        typedef PriorityType<
			PageBucket, 
			MergeStreamComparer<Type, Compare> >	    PQueue;

		TPool   &pool;
        Buffer	bucketBuffer;
        PQueue  pqueue;

        Handler(TPool &_pool):
            pool(_pool),
            pqueue(StreamComparer(_pool.handlerArgs)) { }

        ~Handler() {
            cancel();
        }

		struct insertBucket : public ::std::unary_function<PageBucket,void> {
			Handler &me;
			insertBucket(Handler &_me): me(_me) {}

			inline void operator() (PageBucket &pb) const {
                pb.pageNo = length(me.pqueue);
                readBucket(pb, pb.pageNo, me.pool.pageSize, me.pool.dataSize(pb.pageNo), me.pool.file);
                push(me.pqueue, pb);
			}
		};

        bool begin() {
            // 1. initially fill priority queue
//    		pqueue.reserve(pool.pages);
            equiDistantDistribution(
                bucketBuffer, pool.bucketBufferSize, *this,
                pool._size, pool.pageSize,
                insertBucket(*this));
            return true;
        }

        inline Type const & front() const {
			return *(top(pqueue).cur);
        }

        inline void pop(Type &_Ref) {
            PageBucket &pb = top(pqueue);
            _Ref = *pb.cur;
            if (++pb.cur == pb.end)
                // bucket is empty, we have to fetch the next bucket
				if (!readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo), pool.file)) {
					::seqan::pop(pqueue);
					return;
				}
			adjustTop(pqueue);
        }

        inline void pop() {
            PageBucket &pb = top(pqueue);
            if (++pb.cur == pb.end)
                // bucket is empty, we have to fetch the next bucket
				if (!readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo), pool.file)) {
					::seqan::pop(pqueue);
					return;
				}
			adjustTop(pqueue);
        }

		inline bool eof() const {
			return empty(pqueue);
		}

        inline void end() {
            cancel();
        }
        
        void cancel()
        {
            clear(pqueue);
            freePage(bucketBuffer, *this);
        }

        inline void process() {}

    };


    template < typename TPool >
    struct BufferHandler< TPool, ReadSorterSpec >
    {
        typedef typename TPool::Type                Type;
        typedef typename TPool::Buffer				Buffer;
        typedef typename TPool::Spec::Compare       Compare;

        typedef Buffer&								BufferRef;
		typedef PageBucketExtended<Type>            PageBucket;

        typedef MergeStreamComparer<Type, Compare>  StreamComparer;
        typedef ::std::priority_queue <
            PageBucket,
            ::std::vector<PageBucket>,
            MergeStreamComparer<Type, Compare> >	PQueue;

        TPool       &pool;
        unsigned    mergeBufferSize;
        Buffer	    bucketBuffer, mergeBuffer, tmpBuffer;
        PQueue	    pqueue;

        BufferHandler(TPool &_pool):
            pool(_pool),
            mergeBufferSize(_pool.pageSize),
            pqueue(StreamComparer(_pool.handlerData.compare)) { }

        BufferHandler(TPool &_pool, unsigned _requestedBufferSize):
            pool(_pool),
            mergeBufferSize(_min(_pool.size(), _requestedBufferSize)),
            pqueue(StreamComparer(_pool.handlerData.compare)) { }
        
        ~BufferHandler() {
            cancel();
        }

		struct insertBucket : public ::std::unary_function<PageBucket,void> {
			BufferHandler &me;
			insertBucket(BufferHandler &_me): me(_me) {}

			inline void operator() (PageBucket &pb) const {
                pb.pageNo = length(pqueue);
                readBucket(pb, pb.pageNo, me.pool.pageSize, me.pool.dataSize(pb.pageNo), pool.file);
                push(pqueue, pb);
			}
		};

        inline BufferRef first() {
            // 1. initially fill priority queue
//    		pqueue.reserve(pool.pages);
            equiDistantDistribution(
                bucketBuffer, pool.bucketBufferSize, *this,
                pool._size, pool.pageSize,
                insertBucket(*this));
			allocPage(mergeBuffer, mergeBufferSize, *this);
			return merge();
        }

		inline BufferRef next() {
			return merge();
		}

        inline void end() {
            cancel();
        }
        
        void cancel()
        {
            clear(pqueue);
			freePage(mergeBuffer, *this);
            freePage(bucketBuffer, *this);
        }

        inline void process() {}

    private:

        BufferRef merge()
        {
            // 2. merge streams into mergeBuffer
            
			typename PQueue::size_type pqsize = length(pqueue);
			if (!pqsize) {
				resize(tmpBuffer, 0);
				return tmpBuffer;
			}

            if (pqsize == 1)
			{
                // only one stream left => retrieve what's left in stream es
                PageBucket &pb = pqueue.top();

                // flush stream
				if (pb.cur != pb.end) {
					tmpBuffer.begin = pb.cur;
					tmpBuffer.end = pb.end;
					pb.cur = pb.end;
					if (pb.pageOfs == pool.dataSize(pb.pageNo))
						pop(pqueue);
					return tmpBuffer;
				}

                // read directly from disk
				pb.begin = mergeBuffer.begin;
				pb.end = pb.cur = pb.begin + pageSize(mergeBuffer);

				resize(mergeBuffer, readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo)));
				if (pb.pageOfs == pool.dataSize(pb.pageNo))
					pop(pqueue);
			}
            else
            {
                for(Type *cur = mergeBuffer.begin; cur != mergeBuffer.end; ++cur) {
                    PageBucket &pb = pqueue.top();
                    *cur = *pb.cur;
					if (++pb.cur == pb.end) {
						// bucket is empty, we have to fetch the next bucket
						if (!readBucket(pb, pb.pageNo, pool.pageSize, pool.dataSize(pb.pageNo), pool.file)) {
							pop(pqueue);
							// queue contains only one stream
							// => we return what we have merged
							if (--pqsize == 1) {
								resize(cur - mergeBuffer.begin + 1);
								return mergeBuffer;
							}
						}
					}
					adjustTop(pqueue);
                }
				resize(mergeBuffer, pageSize(mergeBuffer));
            }
            
            return mergeBuffer;
        }
	};

	template < typename TValue,
			   typename TConfig >
    inline PageFrame< TValue, typename TConfig::File, Dynamic<> > & processBuffer(
        PageFrame< TValue, typename TConfig::File, Dynamic<> > &buf,
        BufferHandler< Pool< TValue, SorterSpec<TConfig> >, WriteFileSpec > &me)
    {
		sort(buf.begin, buf.end, AdaptorCompare2Less<typename TConfig::Compare>(me.pool.handlerArgs));
		return buf;
    }


	//////////////////////////////////////////////////////////////////////////////
	// character and buffer based handler definitions
    template < typename TValue,
               typename TConfig >
    inline SimpleBuffer<TValue> & processBuffer(
        SimpleBuffer< TValue > &buf,
        BufferHandler< Pool< TValue, SorterSpec<TConfig> >, MemorySpec > &me)
    {
        sort(buf.begin, buf.end, AdaptorCompare2Less<typename TConfig::Compare>(me.pool.handlerArgs));
		return buf;
    }

    template < typename TValue,
			   typename TConfig >
    struct BufReadHandler< Pool< TValue, SorterSpec<TConfig> > >
    {
        typedef BufferHandler< Bundle2<
			BufferHandler< Pool< TValue, SorterSpec<TConfig> >, MemorySpec >,
			BufferHandler< Pool< TValue, SorterSpec<TConfig> >, ReadSorterSpec >
		>, MultiplexSpec > Type;
    };

    template < typename TValue,
			   typename TConfig >
    struct ReadHandler< Pool< TValue, SorterSpec<TConfig> > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler < Pool< TValue, SorterSpec<TConfig> >, MemorySpec >, AdapterSpec >,
			Handler<				 Pool< TValue, SorterSpec<TConfig> >, ReadSorterSpec >
		>, MultiplexSpec > Type;
    };


    //template < typename TValue,
    //           typename TCompare,
			 //  typename TFile = File<> >
    //struct Sorter: public Pool< TValue, SorterSpec< SorterConfig<TCompare, TFile> > >
    //{
    //    typedef Pool< TValue, SorterSpec< SorterConfig<TCompare, TFile> > > Base;

    //    typedef TValue	                    Type;
    //    typedef TFile                       File;
    //    typedef typename Size<TFile>::Type  SizeType;

    //    Sorter(TCompare const & C)
    //    {
    //        Base::handlerArgs = C;
    //    }

    //    template < typename TInput, typename TSpec >
    //    Sorter(Pipe<TInput, TSpec> &src):
    //        Base(src) {}

    //    template < typename TInput, typename TSpec >
    //    Sorter(Pipe<TInput, TSpec> &src, TCompare const & C):
    //        Base(src)
    //    {
    //        Base::handlerArgs = C;
    //    }
    //};

 //   template < typename TValue,
 //              typename TCompare,
	//		   typename TFile >
 //   ::std::ostream& operator<<(::std::ostream &out, Sorter<TValue, TCompare, TFile> &p) {
 //       beginRead(p);
 //       while (!eof(p)) {
	//	    out << p.front() << ::std::endl;
 //           p.pop();
 //       }
 //       endRead(p);
	//	return out;
	//}

}

#endif
