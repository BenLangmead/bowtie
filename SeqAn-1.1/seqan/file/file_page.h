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
  $Id: file_page.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_PAGE_H
#define SEQAN_HEADER_FILE_PAGE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


	//////////////////////////////////////////////////////////////////////////////
	// base class for memory buffers

	template < typename TValue >
	struct SimpleBuffer {
        typedef TValue      Type;
        typedef size_t		SizeType;

		typedef	TValue&		TypeRef;
		typedef TValue*     TypePtr;
		typedef TValue*     Iterator;

		Iterator            begin;      // the beginning of the buffer
        Iterator            end;        // end of valid data
        SizeType            pageSize;   // size of allocated memory

        SimpleBuffer():
            begin(NULL),
            end(NULL) {}

        SimpleBuffer(TypePtr _begin, TypePtr _end):
            begin(_begin),
            end(_end) {}

        SimpleBuffer(TypePtr _begin, SizeType _size):
            begin(_begin),
            end(_begin + _size) {}

        SimpleBuffer(SizeType _pageSize):
            begin(NULL),
            end(NULL),
            pageSize(_pageSize) {}

        inline Type& operator[](SizeType i) { return begin[i]; }
        inline Type const & operator[](SizeType i) const { return begin[i]; }
	};

    template < typename TValue >
    struct Iterator< SimpleBuffer<TValue>, Standard > {
        typedef TValue* Type;
    };

    template < typename TValue >
    struct Iterator< SimpleBuffer<TValue> const, Standard > {
        typedef TValue const * Type;
    };

    template < typename TValue >
    struct Value< SimpleBuffer<TValue> >
    {
        typedef TValue Type;
    };

    template < typename TValue >
    struct Size< SimpleBuffer<TValue> >
    {
		typedef typename SimpleBuffer<TValue>::SizeType Type;
    };

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    pageSize(SimpleBuffer<TValue> &me) {
        return me.pageSize;
    }

    template < typename TValue, typename TSize >
    inline void setPageSize(SimpleBuffer<TValue> &me, TSize size) {
        me.pageSize = size;
    }

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    size(SimpleBuffer<TValue> const &me) {
        return me.end - me.begin;
    }

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    length(SimpleBuffer<TValue> const &me) {
        return me.end - me.begin;
    }

    template < typename TValue, typename TSize >
    inline void resize(SimpleBuffer<TValue> &me, TSize size) {
        me.end = me.begin + size;
    }

    template < typename TValue, typename TSize, typename T >
	inline void allocPage(SimpleBuffer<TValue> &pf, TSize size, T const & me) {
        setPageSize(pf, size);
        allocate(me, pf.begin, pageSize(pf));
        resize(pf, size);
	}

	template < typename TValue, typename T > inline
	void freePage(SimpleBuffer<TValue> &pf, T const & me) {
		deallocate(me, pf.begin, pageSize(pf));
		pf.begin = NULL;
        resize(pf, 0);
        setPageSize(pf, 0);
	}

	template < typename TValue >
	inline TValue* begin(SimpleBuffer<TValue> &pf, Standard) {
		return pf.begin;
	}

	template < typename TValue >
	inline TValue const * begin(SimpleBuffer<TValue> const &pf, Standard) {
		return pf.begin;
	}

	template < typename TValue >
	inline TValue * end(SimpleBuffer<TValue> &pf, Standard) {
		return pf.end;
	}

	template < typename TValue >
	inline TValue const * end(SimpleBuffer<TValue> const &pf, Standard) {
		return pf.end;
	}


    //////////////////////////////////////////////////////////////////////////////
	// a bucket is a structure to represent a small window of a page
    // used by algorithms which need a global view of all pages (merge sort, mapper)

	template < typename TValue >
    struct PageBucket {
        unsigned    pageOfs;                // begin of bucket window with relation to page begin
        TValue  	*begin, *cur, *end;     // begin/end of buckets memory buffer and a pointer
    };

    template < typename TValue >
    struct PageBucketExtended : public PageBucket< TValue > {
		int     	pageNo;		            // related page (needed by merger sort)
    };

	template < typename TValue >
    ::std::ostream& operator<<(::std::ostream &out, const PageBucketExtended<TValue> &pb) {
        for(TValue *cur = pb.begin; cur != pb.end; cur++)
            out << *cur << " ";
        return out;
    }


    template < typename TValue, typename TFile, typename TSpec >
    struct PageFrame {};    


	//////////////////////////////////////////////////////////////////////////////
	// page frame of dynamic size

    template < typename TSpec = void >
    struct Dynamic;

    // forward declaration
    template < typename TPageFrame >
    struct PageChain;

    template < typename TValue, typename TFile >
	struct PageFrame< TValue, TFile, Dynamic<> >: public SimpleBuffer< TValue >
	{
        typedef TValue                          Type;
        typedef unsigned                        SizeType;
		typedef TFile							File;

		typedef	TValue&							TypeRef;
		typedef TValue*           	            TypePtr;
        typedef PageChain<PageFrame>		    PageChain;
		typedef SimpleBuffer<TValue>	        Base;
        typedef typename aRequest<TFile>::Type  aRequest;

		bool			dirty;		// data needs to be written to disk before freeing
		unsigned   		pageNo;		// maps frames to pages (reverse vector mapper)
        aRequest        request;    // request structure of the async io process

        enum Status		{ READY, READING, WRITING };
		Status status;

        PageFrame       *next;      // next buffer in a chained list
//        PageChain       *chain;     // related chain 

        PageFrame(/*BufChain *_chain = NULL*/):
			Base(),
            dirty(false),
            pageNo(-1),
			status(READY),
//            chain(_chain) 
			next(NULL) {}
    };
    

	//////////////////////////////////////////////////////////////////////////////
	// page frame of static size

    template < unsigned _PageSize >
    struct Fixed;

    typedef ::std::list<int>		PageLRUList;    // least recently usage list
	typedef PageLRUList::iterator	PageLRUEntry;

    template < typename TValue,
               typename TFile,
               unsigned _PageSize >
	struct PageFrame<TValue, TFile, Fixed<_PageSize> >
	{
        typedef TValue                          Type;
        typedef unsigned                        SizeType;
		typedef TFile							File;
		enum { PageSize = _PageSize };

		typedef	TValue&							TypeRef;
		typedef VolatilePtr<TValue>	            TypePtr;
        typedef typename aRequest<TFile>::Type	aRequest;

		bool			dirty;		// data needs to be written to disk before freeing
		int     		pageNo;		// maps frames to pages (reverse vector mapper)
		TypePtr			begin;	    // start address of page memory
        aRequest        request;    // request structure of the async io process

		enum Status		{ READY, READING, WRITING };
		enum DataStatus	{ ON_DISK = -1, UNINITIALIZED = -2 };
		enum Priority	{ NORMAL_LEVEL = 0, PREFETCH_LEVEL = 1, ITERATOR_LEVEL = 2, PERMANENT_LEVEL = 3 };

		Status status;
		DataStatus dataStatus;

		PageLRUEntry	lruEntry;   // priority based lru
        Priority        priority;

		PageFrame():
            dirty(false),
            pageNo(-1),
			begin(NULL),
			status(READY),
            priority(NORMAL_LEVEL) {}

        inline Type& operator[](SizeType i) { return begin[i]; }
        inline Type const & operator[](SizeType i) const { return begin[i]; }
	};

    template < typename TValue, typename TFile, typename TSpec, typename TSize >
    inline void resize(PageFrame<TValue, TFile, Dynamic<TSpec> > &me, TSize size) {
        me.end = me.begin + size;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    size(PageFrame<TValue, TFile, Fixed<_PageSize> > &/*me*/) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    length(PageFrame<TValue, TFile, Fixed<_PageSize> > const &/*me*/) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    pageSize(PageFrame<TValue, TFile, Fixed<_PageSize> > &/*me*/) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize, typename TSize >
    inline void resize(PageFrame<TValue, TFile, Fixed<_PageSize> > &/*me*/, TSize /*size*/) {}


	//////////////////////////////////////////////////////////////////////////////
	// various page frame methods

	template < typename TValue, typename TFile, typename TSpec >
    ::std::ostream& operator<<(::std::ostream &out, const PageFrame<TValue, TFile, TSpec > &pf) {
        out << "PageFrame @ " << pf.pageNo;
        if (pf.dirty)
            out << " DIRTY";
        else
            out << " CLEAN";

        switch (pf.status) {
			case PageFrame<TValue, TFile, TSpec >::READY:
                out << " READY";
                break;
			case PageFrame<TValue, TFile, TSpec >::READING:
                out << " READING";
                break;
			case PageFrame<TValue, TFile, TSpec >::WRITING:
                out << " WRITING";
        }

        if (pf.dataStatus == pf.ON_DISK)
            out << " ON_DISK";
        else
            out << " UNITIALIZED";

        out << " Prio:" << pf.priority;
        out << " Buffer:" << (TValue*)pf.begin;

        return out;
	}

	template < typename TValue, typename TFile, typename TSpec, typename T > inline
	void allocPage(PageFrame<TValue, TFile, TSpec> &pf, T const & me) {
		TValue* tmp = NULL;
		allocate(me, tmp, pageSize(pf));
		pf.begin = tmp;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "allocPage: " << ::std::hex << tmp << ::std::dec << ::std::endl;
		#endif
	}

	template < typename TValue, typename TFile, typename TSpec, typename T > inline
	void freePage(PageFrame<TValue, TFile, TSpec> &pf, T const & me) {
		#ifdef SEQAN_VVERBOSE
			if ((TValue*)pf.begin)
				::std::cerr << "freePage:  " << ::std::hex << (TValue*)pf.begin << ::std::dec << ::std::endl;
		#endif
        nukeCopies(pf.begin);
		deallocate(me, (TValue*)pf.begin, pageSize(pf));
		pf.begin = NULL;
        resize(pf, 0);
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool readPage(int pageNo, PageFrame<TValue, TFile, TSpec> &pf, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READING;
//        resize(pf, pageSize(pf));
		return areadAt(file, (TValue*)pf.begin, size(pf), (pos_t)pageNo * (pos_t)pageSize(pf), pf.request);
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writePage(PageFrame<TValue, TFile, TSpec> &pf, int pageNo, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.status = pf.WRITING;
//        resize(pf, pageSize(pf));
		return awriteAt(file, (TValue*)pf.begin, size(pf), (pos_t)pageNo * (pos_t)pageSize(pf), pf.request);
	}

	template < typename TValue, typename TFile, typename TSpec, typename TSize> inline
    bool readLastPage(int pageNo, PageFrame<TValue, TFile, TSpec> &pf, TFile &file, TSize size) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return readAt(file, (TValue*)pf.begin, size, (pos_t)pageNo * (pos_t)pageSize(pf));
	}

	template < typename TValue, typename TFile, typename TSpec, typename TSize > inline
	bool writeLastPage(PageFrame<TValue, TFile, TSpec> &pf, int pageNo, TFile &file, TSize size) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return writeAt(file, (TValue*)pf.begin, size, (pos_t)pageNo * (pos_t)pageSize(pf));
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool waitFor(PageFrame<TValue, TFile, TSpec> &pf) {
		if ((pf.status != pf.READY) && waitFor(pf.request)) {
			pf.status = pf.READY;
			pf.dirty = false;
            return true;
		}
        return false;
	}

	template < typename TValue, typename TFile, typename TSpec, typename TTime > inline
	bool waitFor(PageFrame<TValue, TFile, TSpec> &pf, TTime timeOut) {
		if ((pf.status != pf.READY) && waitFor(pf.request, timeOut)) {
			pf.status = pf.READY;
			pf.dirty = false;
			return true;
		}
        return false;
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool cancel(PageFrame<TValue, TFile, TSpec> &pf, TFile &file) {
        waitFor(pf, 0);
		if (pf.status != pf.READY) {
            if (!cancel(file, pf.request)) return false;
            pf.status = pf.READY;
        }
        return true;
	}

	//////////////////////////////////////////////////////////////////////////////
	// page based read/write methods used by Pool classes

    template < typename TValue, typename TFile, typename TSpec > inline
	bool readPage(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, TFile &file) {
        if (size(pf) == pageSize(pf))
            return readPage(pf.pageNo, pf, file);
        else
            return readLastPage(pf.pageNo, pf, file, size(pf));
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writePage(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, TFile &file) {
        if (size(pf) == pageSize(pf))
            return writePage(pf, pf.pageNo, file);
        else
            return writeLastPage(pf, pf.pageNo, file, size(pf));
	}

	template < typename TValue, typename TFile > inline
	unsigned readBucket(PageBucket<TValue> &b, int pageNo, unsigned pageSize, unsigned dataSize, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
        unsigned readSize = _min(dataSize - b.pageOfs, (unsigned)(b.end - b.begin));
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readBucket:  " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)pageNo * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << readSize << ::std::endl;
		#endif
        if (readSize && readAt(file, b.begin, readSize, (pos_t)pageNo * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += readSize;
            b.cur = b.begin;
            b.end = b.begin + readSize;
            return readSize;
        } else
            return 0;
	}

	template < typename TValue, typename TFile > inline
	bool writeBucket(PageBucket<TValue> &b, int pageNo, unsigned pageSize, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)pageNo * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << b.cur - b.begin << ::std::endl;
		#endif
        if ((b.cur == b.begin) || writeAt(file, b.begin, b.cur - b.begin, (pos_t)pageNo * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += b.cur - b.begin;
            b.cur = b.begin;
            return true;
        } else
            return false;
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writeBucket(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, unsigned &pageOfs, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << pf.begin;
			::std::cerr << " from page " << ::std::dec << pf.pageNo << " at " << (pos_t)pf.pageNo * (pos_t)pageSize(pf) + pageOfs;
			::std::cerr << " size " << size(pf) << ::std::endl;
		#endif
        if (pf.end == pf.begin) return true;
        if (awriteAt(file, pf.begin, size(pf), (pos_t)pf.pageNo * (pos_t)pageSize(pf) + pageOfs, pf.request)) {
            pf.status = pf.WRITING;
            pageOfs += size(pf);
            return true;
        } else
            return false;
	}



    template < typename TPageFrame >
    struct PageChain {
        typedef typename TPageFrame::Type		Type;
        typedef typename TPageFrame::SizeType	SizeType;

		typedef	Type&							TypeRef;
		typedef Type*							TypePtr;
		typedef TPageFrame						PageFrame;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef TPageFrame * 		iterator;
		typedef TPageFrame const *	const_iterator;
		typedef TPageFrame      	value_type;
		typedef TPageFrame &     	reference;
		typedef TPageFrame const &  const_reference;
		typedef TPageFrame * 		pointer;
		typedef TPageFrame const *  const_pointer;
		typedef unsigned       		size_type;
		typedef int                 difference_type;

        pointer             first, last;
/*        pointer             firstFree, lastFree;
        Semaphore           freeBuffers;
        Mutex               dataLock;*/

        unsigned            frames, maxFrames;
        
        PageChain(unsigned _maxFrames = UINT_MAX):
            first(NULL),
            last(NULL),
            frames(0),
            maxFrames(_maxFrames)
//            freeBuffers(_buffers)
        {
            for(unsigned i = 0; i < _maxFrames; ++i)
                pushBack();
        }
        
        ~PageChain()
        {
            while (first)
                popFront();
        }
        
        inline reference operator[](int k) {
            pointer p = first;
            while (k) {
                p = p->next;
                --k;
            }
            return *p;
        }
        
        inline const_reference operator[](int k) const {
            pointer p = first;
            while (k) {
                p = p->next;
                --k;
            }
            return *p;
        }

        inline pointer getReadyPage() {
            if (!first || (frames < maxFrames && waitFor(*first, 0)))
                return pushBack();
            else {
                waitFor(*first);
                return firstToEnd();
            }
        }

        template < typename TFile >
        inline void cancelAll(TFile &file) {
            pointer p = first;
            while (p) {
                cancel(*p, file);
                p = p->next;
            }
        }

        inline void waitForAll() {
            pointer p = first;
            while (p) {
                waitFor(*p);
                p = p->next;
            }
        }

/*        BufHeader* getFreeBuffer() {
            freeBuffers.lock();
            dataLock.lock();
            BufHeader *h = firstFree;
            if (!(firstFree = firstFree->next)) lastFree = NULL;
            dataLock.unlock();
            return h;
        }
        
        void addFreeBuffer(BufHeader *h) {
            dataLock.lock();
            if (lastFree) {
                lastFree->next = h;
                lastFree = h;
            } else {
                firstFree = h;
                lastFree = h;
            }
            dataLock.unlock();
            freeBuffers.unlock();
        }*/
        
    private:

        inline pointer firstToEnd() {
            last->next = first;
            last = first;
            first = first->next;
            last->next = NULL;
            return last;
        }

        inline pointer pushBack() {
            pointer p = new PageFrame();
            if (p) {
                if (last)
                    last->next = p;
                else
                    first = p;
                last = p;
                ++frames;
            }
            return p;
        }

        inline pointer popFront() {
            pointer p = first;
            if (p) {
                first = first->next;
                if (!first) last = NULL;
                --frames;
                delete p;
            }
            return p;
        }
    };


/*    
    template < typename ST >
    struct BufArrayConfig
    {
        typedef ST SizeType;

        unsigned int    disks;
        unsigned int    buffersPerDisk;
        SizeType        bufferSize;         // prefetch/cache buffer size
        
        BufArrayConfig():
            disks(1),
            buffersPerDisk(3),
            bufferSize(2 * 1024 * 1024) { }
    };
    
    template < typename T, typename ST, typename _File >
    struct BufArray
    {
        typedef T Type;
        typedef ST SizeType;
        typedef BufChain<T,ST,_File> BufChain;
        typedef BufArrayConfig<ST> Config;
        
        Config		conf;
        BufChain**	disk;
        
        BufArray(const Config &_conf):
            conf(_conf)
        {
            disk = new BufChain*[conf.disks];
            for(unsigned int i = 0; i < conf.disks; i++)
                disk[i] = new BufChain(conf.buffersPerDisk, conf.bufferSize, i);
            BufChain *prev = disk[conf.disks - 1];
            for(unsigned int i = 0; i < conf.disks; i++)
                prev = prev->next = disk[i];
        }

        ~BufArray() {
            for(unsigned int i = 0; i < conf.disks; i++)
                delete disk[i];
            delete[] disk;
        }

        BufChain& operator[](const int k) {
            return *disk[k];
        }

        const BufChain& operator[](const int k) const {
            return *disk[k];
        }
    };*/



	//////////////////////////////////////////////////////////////////////////////
	// page container with lru mechanism
	// the lru strategy uses different priorities
	// the page with the least priority is used first
	// 0..random access pages
	// 1..forward iterator pages
	// 2..quasi permanent pages

    template < typename TPageFrame,
               unsigned _Frames,
               unsigned _PriorityLevels = TPageFrame::PERMANENT_LEVEL + 1 >
	struct PageContainer : public ::std::vector<TPageFrame>
	{
        typedef typename TPageFrame::Type			Type;
        typedef typename TPageFrame::SizeType		SizeType;
		typedef typename ::std::vector<TPageFrame>	Base;

		typedef	Type&							TypeRef;
		typedef Type*							TypePtr;
		typedef TPageFrame						PageFrame;

		enum { PriorityLevels = _PriorityLevels };

		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef typename Base::iterator			iterator;
		typedef typename Base::const_iterator	const_iterator;
		typedef TPageFrame      				value_type;
		typedef TPageFrame &     				reference;
		typedef TPageFrame const &				const_reference;
		typedef TPageFrame * 					pointer;
		typedef TPageFrame const *				const_pointer;
		typedef unsigned       					size_type;
		typedef int								difference_type;
		
        PageLRUList lruList[_PriorityLevels];

		PageContainer():
			Base(_Frames)
		{
            for(size_type i = 0; i < _Frames; ++i)
			    (*this)[i].lruEntry = lruList[0].insert(lruList[0].end(), i);
        }

		inline void reserve(unsigned _Count) {
			Base::reserve(_Count);
		}

		void resize(unsigned _Count) {
			unsigned _Size = Base::size();
			if (_Size < _Count)
                for(unsigned i = _Size; i < _Count; ++i)
                    push_back();
			else 
				if (_Size > _Count)
					for(unsigned i = _Count; i < _Size; ++i)
						pop_back();
		}

		inline void push_back() {
			Base::push_back(PageFrame());
			Base::back().lruEntry = lruList[0].insert(lruList[0].end(), Base::size() - 1);
		}

		inline void erase(int frameNo) {
			lruList[(*this)[frameNo].priority].erase((*this)[frameNo].lruEntry);
            Base::erase(Base::begin() + frameNo);
		}

        inline void rename(int frameNo) {
            *((*this)[frameNo].lruEntry) = frameNo;
        }

		inline void pop_back() {
			lruList[Base::back().priority].erase(Base::back().lruEntry);
			Base::pop_back();
		}


		//////////////////////////////////////////////////////////////////////////////
		// lru strategy interface

		inline void upgrade(const PageFrame &pf) {
			lruList[pf.priority].splice(lruList[pf.priority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].begin();
		}

		inline void downgrade(const PageFrame &pf) {
			lruList[pf.priority].splice(lruList[pf.priority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].end() - 1;
		}

		inline void upgrade(PageFrame &pf, int newPriority) {
			lruList[newPriority].splice(lruList[newPriority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].begin();
			pf.priority = static_cast<typename PageFrame::Priority> (newPriority);
		}

		inline void downgrade(PageFrame &pf, int newPriority) {
			lruList[newPriority].splice(lruList[newPriority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].end() - 1;
			pf.priority = static_cast<typename PageFrame::Priority> (newPriority);
		}

		inline void _dump() {
			for(unsigned i = 0; i < _PriorityLevels; ++i) {
                ::std::cerr << "|";
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
                    PageFrame &pf = (*this)[*I];
                    ::std::cerr << pf.pageNo;
                    if (pf.dirty) ::std::cerr << "*";
                    else          ::std::cerr << " ";
                    if (pf.status == PageFrame::READY) ::std::cerr << "  ";
                    else                               ::std::cerr << ". ";
				};
            }
            ::std::cerr << ::std::endl;
		}

        // Function is a functor which is called with a PageFrame object,
        // that is dirty or not READY (in an IO transfer)
		template <class Function>
		inline int mru(Function _Func, unsigned maxLevel = _PriorityLevels - 1) {
			for(unsigned i = 0; i <= maxLevel; ++i) {
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
					PageFrame& pf =(*this)[*I];
					if (pf.status == PageFrame::READY && !pf.dirty)
						return *I;
					else
						if (_Func(pf)) return *I;
				};
            }
			#ifdef SEQAN_VVERBOSE
				::std::cerr << "ALL PAGES DIRTY OR IN USE (try to use const iterators) :-(" << ::std::endl;
			#endif
			return -1;
		}

		inline int mruDirty() {
            for(unsigned i = 0; i < _PriorityLevels; ++i)
                if (!lruList[i].empty())
                    return lruList[i].back();
			return -1;
		}
	};




    template < typename TValue, typename TSize, typename T, class Function >
    inline bool equiDistantDistribution(
        SimpleBuffer<TValue> &_clusterBuffer, unsigned _bufferSize, T const &me,
        TSize _size, unsigned _pageSize,
        Function const &_Func)
    {
        unsigned _pages         = enclosingBlocks(_size, (unsigned)_pageSize);
        if (!_pages) {
			::std::cerr << "equiDistantDistribution: _pages is null!" << ::std::endl;
            return false;
        }

        if (_bufferSize < _pages) {
			::std::cerr << "equiDistantDistribution: clusterBufferSize is too small -> raised to " << _pages << ::std::endl;
            _bufferSize = _pages;
        }

        unsigned lastPageSize   = _size % _pageSize;
        unsigned pages          = _pages;

        if ((TSize)_bufferSize > _size)
            _bufferSize = _size;

        allocPage(_clusterBuffer, _bufferSize, me);
        PageBucketExtended<TValue> pb;
        pb.begin = _clusterBuffer.begin;

        unsigned clusterSize = _bufferSize / pages;
        if (lastPageSize > 0 && clusterSize >= lastPageSize) {
            // last page bucket would get more memory than page would need
            // --> exclude from equi-size distribution
            if (--pages) {
                _bufferSize -= lastPageSize;
                clusterSize = _bufferSize / pages;
            }
        }

        if (pages) {
            unsigned remainder = _bufferSize % pages;
            for(unsigned i = 0, numerator = 0; i < pages; ++i) {
                pb.end = pb.begin + clusterSize;
                if ((numerator += remainder) >= pages) {    // simple bresenham for distribution
                    numerator -= pages;
                    ++pb.end;
                }
                pb.cur = pb.begin;
                pb.pageOfs = 0;
			    _Func(pb);
                pb.begin = pb.end;
            }
        }

        if (pages < _pages) {
            pb.end = pb.begin + lastPageSize;
            pb.cur = pb.begin;
            pb.pageOfs = 0;
			_Func(pb);
        }

        return true;
    }

    template < typename TValue, typename TSize, typename T, class Function >
    inline unsigned equiDistantAlignedDistribution(
        SimpleBuffer<TValue> &_clusterBuffer, unsigned aligning, unsigned _bufferSize, T const &me,
        TSize _size, unsigned _pageSize,
        Function const &_Func)
    {
        unsigned _pages         = enclosingBlocks(_size, (unsigned)_pageSize);
        if (!_pages) {
			::std::cerr << "equiDistantDistribution: _pages is null!" << ::std::endl;
            return 0;
        }

        if (_bufferSize < _pages) {
			::std::cerr << "equiDistantAlignedDistribution: clusterBufferSize is too small -> raised to " << _pages << ::std::endl;
            _bufferSize = _pages;
        }

        unsigned lastPageSize   = _size % _pageSize;
        unsigned pages          = _pages;

        if ((TSize)_bufferSize > _size)
            _bufferSize = _size;

        unsigned clusterSize = _bufferSize / pages;
        unsigned aclusterSize = (clusterSize / aligning) * aligning;
        if (clusterSize - aclusterSize > aligning / 2)
            aclusterSize += aligning;

		if (aclusterSize != 0) {

			if (lastPageSize > 0 && aclusterSize > lastPageSize) {
				// last page bucket would get more memory than page would need
				// --> exclude from equi-size distribution
				--pages;
				allocPage(_clusterBuffer, aclusterSize * pages + lastPageSize, me);
			} else
				allocPage(_clusterBuffer, aclusterSize * pages, me);

			PageBucketExtended<TValue> pb;
			pb.begin = _clusterBuffer.begin;

			if (pages) {
				for(unsigned i = 0; i < pages; ++i) {
					pb.end = pb.begin + aclusterSize;
					pb.cur = pb.begin;
					pb.pageOfs = 0;
					_Func(pb);
					pb.begin = pb.end;
				}
			}

			if (pages < _pages) {
				pb.end = pb.begin + lastPageSize;
				pb.cur = pb.begin;
				pb.pageOfs = 0;
				_Func(pb);
			}
		}

        return aclusterSize;
    }

}

#endif
