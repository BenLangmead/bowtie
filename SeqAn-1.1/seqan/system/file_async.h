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
  $Id: file_async.h,v 1.1 2008/08/25 16:20:08 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_FILE_ASYNC_H
#define SEQAN_HEADER_FILE_ASYNC_H

namespace SEQAN_NAMESPACE_MAIN
{

 
	template <typename TSpec /* = void */>
	struct Async;


#ifdef PLATFORM_WINDOWS


    static DWORD _transferedBytes;  // for reporting

	template <typename TSpec>
	class File<Async<TSpec> >
    {
    public:

        typedef LONGLONG    FilePtr;
        typedef ULONGLONG   SizeType;
        typedef DWORD       _SizeType;
        typedef HANDLE      Handle;

		Handle              hFile, hFileAsync;
        bool                noBuffering;

        File():
            hFile(INVALID_HANDLE_VALUE) {}

        File(void *): // to be compatible with the FILE*(NULL) constructor
            hFile(INVALID_HANDLE_VALUE) {}

        bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) {
			SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
            noBuffering = getExtraFlags(openMode | OPEN_ASYNC) & (FILE_FLAG_NO_BUFFERING | FILE_FLAG_OVERLAPPED);
            hFileAsync = CreateFile(fileName,
                                getFileAccess(openMode | OPEN_ASYNC),
                                FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL,
                                getCreationFlags(openMode | OPEN_ASYNC),
                                getExtraFlags(openMode | OPEN_ASYNC),
                                NULL);

            if (hFileAsync == INVALID_HANDLE_VALUE) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Open failed on file " << fileName << ". (ErrNo=" << GetLastError() << ")" << ::std::endl;
                return false;
            }
            #ifdef SEQAN_VERBOSE
				if (!(openMode & OPEN_QUIET))
	                ::std::cerr << "file opened asynchronously " << fileName << " handle " << ::std::hex << hFileAsync << ::std::dec << ::std::endl;
            #endif

            if (noBuffering) {
                hFile = CreateFile(fileName,                // in this case io must be sector aligned
                                getFileAccess(openMode),    // so we open a second file, for unaligned access
                                FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL,
                                OPEN_EXISTING,
                                getExtraFlags(openMode & ~OPEN_ASYNC),
                                NULL);
                if (hFile == INVALID_HANDLE_VALUE) {
					if (!(openMode & OPEN_QUIET))
	                	::std::cerr << "Open failed on secondary file " << fileName << ". (ErrNo=" << GetLastError() << ")" << ::std::endl;
                    return false;
                }
	            #ifdef SEQAN_VERBOSE
					if (!(openMode & OPEN_QUIET))
	                	::std::cerr << "async file opened  " << fileName << " handle " << ::std::hex << hFile << ::std::dec << ::std::endl;
                #endif
            } else
                hFile = hFileAsync;

            return true;
        }

        bool openTemp(int openMode = DefaultOpenTempMode<File>::VALUE) {
            char szTempName[MAX_PATH];
#ifdef SEQAN_DEFAULT_TMPDIR
            static const char szTempPath[MAX_PATH] = SEQAN_DEFAULT_TMPDIR;
#else
            char szTempPath[MAX_PATH];
            if (!GetTempPath(MAX_PATH, szTempPath)) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Couldn't get a temporary path name. (ErrNo=" << GetLastError() << ")" << ::std::endl;
                return false;
            }
#endif
            if (!GetTempFileName(szTempPath, "GNDX", 0, szTempName)) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Couldn't get a temporary file name. (ErrNo=" << GetLastError() << ")" << ::std::endl;
                return false;
            }
            return open(szTempName, openMode | OPEN_TEMPORARY);
        }

        inline bool close() {
            BOOL result = true;
            #ifdef SEQAN_VERBOSE
                ::std::cerr << "files closed handles " << ::std::hex << hFileAsync << " and " << hFile << ::std::dec << ::std::endl;
            #endif
            if (hFile != hFileAsync)
                result &= CloseHandle(hFileAsync);
            result &= CloseHandle(hFile);
            hFileAsync = INVALID_HANDLE_VALUE;
            hFile = INVALID_HANDLE_VALUE;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return result;
        }

        inline bool read(void *memPtr, _SizeType count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    bool result = ReadFile(hFile, memPtr, count, &_transferedBytes, NULL);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
        }

        inline bool write(void const *memPtr, _SizeType count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    bool result = WriteFile(hFile, memPtr, count, &_transferedBytes, NULL);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
        }

		inline FilePtr seek(FilePtr _pos, DWORD origin = FILE_BEGIN) {
//          LARGE_INTEGER li = _pos;
//			return SetFilePointer(hFileAsync, li.LowPart, &li.HighPart, MoveMethod);
            LARGE_INTEGER new_pos, pos;
            pos.QuadPart = _pos;
            SetFilePointerEx(hFile, pos, &new_pos, origin);
//            position = new_pos.QuadPart;
            return new_pos.QuadPart;
		}

		inline FilePtr tell() {
			return seek(0, FILE_CURRENT);
        }

		inline FilePtr size() const {
            LARGE_INTEGER result;
            DWORD dwError, high;
            result.LowPart = GetFileSize(hFile, &high);
            result.HighPart = high;
            if (result.LowPart == INVALID_FILE_SIZE && (dwError = GetLastError()) != NO_ERROR) {
				::std::cerr << "Couldn't get file size. (ErrNo=" << dwError << ")" << ::std::endl;
                return 0;
            }
            return result.QuadPart;
        }

        inline bool setEOF() const {
            return SetEndOfFile(hFile);
        }

		inline static DWORD error() {
			return GetLastError();
		}

        operator bool () const {
            return (hFile != INVALID_HANDLE_VALUE) && (hFileAsync != INVALID_HANDLE_VALUE);
        }

    protected:

        DWORD getFileAccess(int openMode) {
            switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    return GENERIC_READ;
                case OPEN_WRONLY:
                    return GENERIC_WRITE;
                case OPEN_RDWR:
                    return GENERIC_READ | GENERIC_WRITE;
				default:
					return 0;
            }
        }

        DWORD getCreationFlags(int openMode) {
            if (openMode & OPEN_CREATE)
                if (openMode & OPEN_APPEND)
                    return OPEN_ALWAYS;
                else
                    return CREATE_ALWAYS;
            else
                return OPEN_EXISTING;
        }

        DWORD getExtraFlags(int openMode) {
            DWORD extra = FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS;// | FILE_FLAG_WRITE_THROUGH;
            if (openMode & OPEN_ASYNC) {
                extra |= FILE_FLAG_OVERLAPPED;
                #ifdef SEQAN_DIRECTIO
                    extra |= FILE_FLAG_NO_BUFFERING;
                #endif
            }
            if (openMode & OPEN_TEMPORARY)  extra |= FILE_FLAG_DELETE_ON_CLOSE;
            return extra;
        }

    };


    //////////////////////////////////////////////////////////////////////////////
    // (SeqAn adaption)
    //////////////////////////////////////////////////////////////////////////////

    struct aiocb_win32 {
        OVERLAPPED  overlapped;
        Event       xmitDone;
    };

	template <typename TSpec>
    struct aRequest<File<Async<TSpec> > >
    {
        typedef aiocb_win32 Type;
    };
/*
	template <typename TSpec>
    struct aEvent<File<Async<TSpec> > >
    {
        typedef Event Type;
    };


	template <typename TSpec>
    struct aQueue<File<Async<TSpec> > >
    {
        typedef IOQueue Type;
    };

	template <typename TSpec>
    struct aHint<File<Async<TSpec> > >
    {
        typedef typename aQueue<File<Async<TSpec> > >::Type::aHint Type;
    };

	template <typename TSpec>
    struct aCallback<File<Async<TSpec> > >
    {
        typedef typename aQueue<File<Async<TSpec> > >::Type::aCallback Type;
    };*/


	template <typename TSpec>
    inline typename Size<File<Async<TSpec> > >::Type size(File<Async<TSpec> > &me) {
        return me.size();
    }

	template <typename TSpec>
    inline bool setEOF(File<Async<TSpec> > &me) {
        return me.setEOF();
    }

	template <typename TSpec>
    inline unsigned sectorSize(File<Async<TSpec> > const &) {
        DWORD SpC, nofC, tnoC, aligning;
        if (GetDiskFreeSpace(NULL, &SpC, &aligning, &nofC, &tnoC) == 0)  {
            ::std::cerr << "Error " << GetLastError() << " while querying cluster size" << ::std::endl;
            return 4096;
        }
        return aligning;
    }


    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    inline bool areadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb_win32 &request)
    {
        SEQAN_PROTIMESTART(tw);
        LARGE_INTEGER ofs;
        ofs.QuadPart = fileOfs;
        ofs.QuadPart *= sizeof(TValue);
        request.overlapped.Offset = ofs.LowPart;
        request.overlapped.OffsetHigh = ofs.HighPart;
        if (!request.xmitDone) open(request.xmitDone);
        request.overlapped.hEvent = request.xmitDone.hEvent;
        if (ReadFile(
            me.hFileAsync, 
            memPtr, 
            count * sizeof(TValue),
            &ofs.LowPart,
            &request.overlapped) || (me.error() == ERROR_IO_PENDING))
        {
            SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
            return true;
        }
        if (me.error() == ERROR_NO_SYSTEM_RESOURCES) {  // read synchronoulsy instead
            #ifdef SEQAN_DEBUG_OR_TEST_
            	::std::cerr << "Warning: Falling back to sync. read. :( " << ::std::endl;
            #endif
			signal(request.xmitDone);
            return readAt(me, memPtr, count, fileOfs);
        }
        return false;
    }
    
    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    inline bool awriteAt(File<Async<TSpec> > & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aiocb_win32 &request)
    {
        SEQAN_PROTIMESTART(tw);
        LARGE_INTEGER ofs;
        ofs.QuadPart = fileOfs;
        ofs.QuadPart *= sizeof(TValue);
        request.overlapped.Offset = ofs.LowPart;
        request.overlapped.OffsetHigh = ofs.HighPart;
        if (!request.xmitDone) open(request.xmitDone);
        request.overlapped.hEvent = request.xmitDone.hEvent;
        if (WriteFile(
            me.hFileAsync, 
            memPtr, 
            count * sizeof(TValue),
            &ofs.LowPart,
            &request.overlapped) || (me.error() == ERROR_IO_PENDING))
        {
            SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
            return true;
        }
        if (me.error() == ERROR_NO_SYSTEM_RESOURCES) {  // write synchronoulsy instead
            #ifdef SEQAN_DEBUG_OR_TEST_
            	::std::cerr << "Warning: Falling back to sync. write. :( " << ::std::endl;
            #endif
			signal(request.xmitDone);
            return writeAt(me, memPtr, count, fileOfs);
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////
    // queue specific functions

    inline bool waitFor(aiocb_win32 &request) {
        SEQAN_PROTIMESTART(tw);
		if (!waitFor(request.xmitDone, 60000))
            ::std::cerr << "waitFor timeout" << ::std::endl;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return true;
	}

    template < typename TTime >
    inline bool waitFor(aiocb_win32 &request, TTime timeout_millis) {
        SEQAN_PROTIMESTART(tw);
		bool result = waitFor(request.xmitDone, timeout_millis);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb_win32 const * const contexts[], TSize count, DWORD timeout_millis = Event::Infinite) {
        Event::Handle handles[count];
        for(TSize i = 0; i < count; ++i)
            handles[i] = contexts[i]->xmitDone.hEvent;

        SEQAN_PROTIMESTART(tw);
        DWORD result = WaitForMultipleObjects(count, &handles, false, timeout_millis);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        if (/*result >= WAIT_OBJECT_0 && */result < WAIT_OBJECT_0 + count)
    		return result - WAIT_OBJECT_0;
        return count;
	}

	template <typename TSpec>
    inline bool cancel(File<Async<TSpec> > & me, aiocb_win32 const &request) {
        return CancelIo(me.hFileAsync);
    }

	template <typename TSpec>
    inline bool flush(File<Async<TSpec> > & me) {
		if (me.hFile != me.hFileAsync)	// in case of equality no direct access was done -> no flush needed
        	return FlushFileBuffers(me.hFile);
        else
            return true;
    }

    template < typename TSpec, typename aRequest >
    inline void release(File<Async<TSpec> > & me, aRequest & request) { }


/*        
    //////////////////////////////////////////////////////////////////////
    // callback based read/write

    template < typename TSpec, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TSpec> > >::Type
    aread(File<Async<TSpec> > & me, TValue *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TSpec> > >::Type request = 
            me.queue->areadAt(
                me.hFileAsync,
                me.position,
                memPtr,
                bsize,
                cb,
                hint);
        me.position += bsize;
        return request;
    }
    
    template < typename TSpec, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TSpec> > >::Type
    awrite(File<Async<TSpec> > & me, TValue const *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TSpec> > >::Type request = 
            me.queue->awriteAt(
                memPtr,
                me.hFileAsync,
                me.position,
                bsize,
                cb,
                hint);
        me.position += bsize;
        return request;
    }

    template < typename TSpec, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TSpec> > >::Type
    areadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->areadAt(
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            memPtr,
            bsize,
            cb,
            hint);
    }
    
    template < typename TSpec, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TSpec> > >::Type
    awriteAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->awriteAt(
            memPtr,
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            bsize,
            cb,
            hint);
    }


    //////////////////////////////////////////////////////////////////////
    // event based read/write

    template < typename TSpec, typename TValue, typename TSize,
               typename aEvent >
    inline typename aRequest<File<Async<TSpec> > >::Type
    aread(File<Async<TSpec> > & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TSpec> > >::Type request = 
            me.queue->areadAt(
                me.hFileAsync,
                me.position,
                memPtr,
                bsize,
                event);
        me.position += bsize;
        return request;
    }
    
    template < typename TSpec, typename TValue, typename TSize,
               typename aEvent >
    inline typename aRequest<File<Async<TSpec> > >::Type
    awrite(File<Async<TSpec> > & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TSpec> > >::Type request =  
            me.queue->awriteAt(
                memPtr,
                me.hFileAsync,
                me.position,
                bsize,
                event);
        me.position += bsize;
        return request;
    }

    template < typename TSpec, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename aRequest<File<Async<TSpec> > >::Type
    areadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->areadAt(
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            memPtr,
            bsize,
            event);
    }
    
    template < typename TSpec, TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename aRequest<File<Async<TSpec> > >::Type
    awriteAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->awriteAt(
            memPtr,
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            bsize,
            event);
    }


    //////////////////////////////////////////////////////////////////////
    // queue specific functions

	template <typename TSpec>
    inline void flush(File<Async<TSpec> > & me) {
        me.queue->flush();
    }

    template < typename TSpec, typename aRequest >
    inline void release(File<Async<TSpec> > & me, aRequest & request) {
        me.queue->release(request);
    }
*/

	//////////////////////////////////////////////////////////////////////////////
	// page aligned allocate for direct file io

    struct TagAllocateAligned_;	//< allocate page aligned memory for direct i/o access
    typedef Tag<TagAllocateAligned_> const TagAllocateAligned;

	template <typename T, typename TValue, typename TSize>
	inline void
	allocate(T const &, 
			 TValue * & data,
			 TSize count,
			 TagAllocateAligned const)
	{
		data = (TValue *) VirtualAlloc(NULL, count * sizeof(TValue), MEM_COMMIT, PAGE_READWRITE);
        if (data)
            SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
        else
			::std::cerr << "AlignAllocator: Could not allocate memory of size " << ::std::hex << count * sizeof(TValue) << ::std::dec << ". (ErrNo=" << GetLastError() << ")" << ::std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////
	// page aligned deallocate for direct file io

	template <typename T, typename TValue, typename TSize>
	inline void 
	deallocate( T const &,
				TValue * data, 
				TSize count,
				TagAllocateAligned const)
	{
		if (data) {
			VirtualFree(data, 0, MEM_RELEASE);
			if (count)	// .. to use count if SEQAN_PROFILE is not defined
				SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
		}
	}

#else

    
	template <typename TSpec>
    class File<Async<TSpec> > : public File<Sync<TSpec> >
    {
    public:

        typedef File<Sync<TSpec> >  Base;

        typedef off_t			FilePtr;
		typedef off_t           SizeType;   // type of file size
        typedef size_t          _SizeType;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handleAsync;
		using Base::handle;

		File(void * = NULL): 	// to be compatible with the FILE*(NULL) constructor
			handleAsync(-1) {}

        bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) {
            handle = ::open(fileName, Base::_getOFlag(openMode & ~OPEN_ASYNC), S_IREAD | S_IWRITE);
			if (handle == -1) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Open failed on file " << fileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
				return false;
			}

			if (Base::_getOFlag(openMode | OPEN_ASYNC) & O_DIRECT) {
				handleAsync = ::open(fileName, Base::_getOFlag(openMode | OPEN_ASYNC & ~OPEN_CREATE), S_IREAD | S_IWRITE);
				if (handleAsync == -1 || errno == EINVAL) {	// fall back to cached access
					#ifdef SEQAN_DEBUG_OR_TEST_
						if (!(openMode & OPEN_QUIET))
							::std::cerr << "Warning: Direct access openening failed. (" << ::strerror(errno) << ")" << ::std::endl;
					#endif
					handleAsync = handle;
				}
				#ifdef SEQAN_DEBUG_OR_TEST_
				    else
						if (!(openMode & OPEN_QUIET))
							::std::cerr << "Direct access successfully initiated" << ::std::endl;
				#endif
			} else
				handleAsync = handle;

			if (sizeof(FilePtr) < 8 && !(openMode & OPEN_QUIET))
				// To remove this warning, you have to options:
				// 1. include the following line before including anything in your application
				//    #define _FILE_OFFSET_BITS 64
				// 2. include <seqan/platform.h> or <seqan/sequence.h> before any other include
				::std::cerr << "WARNING: FilePtr is not 64bit wide" << ::std::endl;


			SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
            return true;
        }

        bool openTemp(int openMode = DefaultOpenTempMode<File>::VALUE) {
#ifdef SEQAN_DEFAULT_TMPDIR
			char tmpFileName[] = SEQAN_DEFAULT_TMPDIR "/GNDXXXXXXX";
#else
			char tmpFileName[] = "/var/tmp/GNDXXXXXXX";
#endif
			if ((handle = handleAsync = ::mkstemp(tmpFileName)) == -1) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Couldn't create temporary file " << tmpFileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
				return false;
			}
			if (!(close() && open(tmpFileName, openMode))) return false;
            #ifdef SEQAN_DEBUG
				if (::unlink(tmpFileName) == -1 && !(openMode & OPEN_QUIET))
					::std::cerr << "Couldn't unlink temporary file " << tmpFileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
            #else
				::unlink(tmpFileName);
			#endif
			return true;
        }

        inline bool close() {
			bool result = true;
			if (handle != handleAsync)
	            result &= (::close(handleAsync) == 0);
            result &= (::close(handle) == 0);
            handleAsync = -1;
            handle = -1;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return result;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // (SeqAn adaption)
    //////////////////////////////////////////////////////////////////////////////
/*
	template <typename TSpec>
    struct aQueue<File<Async<TSpec> > >
    {
        typedef void* Type;
    };
*/
	template <typename TSpec>
    struct aRequest<File<Async<TSpec> > >
    {
        typedef aiocb Type;
    };
/*
	template <typename TSpec>
    struct aEvent<File<Async<TSpec> > >
    {
        typedef aiocb Type;
    };
*/

    //////////////////////////////////////////////////////////////////////
    // event based read/write

//    enum { _AsyncIOSignal = SIGIO };

	inline void printRequest(aiocb &request, const char *_hint = NULL) {
		::std::cerr << ::std::hex;
		if (_hint)
			::std::cerr << _hint << ::std::endl;
		::std::cerr << "fildes:  " << request.aio_fildes << ::std::endl;
		::std::cerr << "buffer:  " << (unsigned long)request.aio_buf << ::std::endl;
		::std::cerr << "offset:  " << request.aio_offset<< ::std::endl;
		::std::cerr << "nbytes:  " << request.aio_nbytes << ::std::endl;
		::std::cerr << "event:   " << request.aio_sigevent.sigev_notify << ::std::endl;
		::std::cerr << "Raddr:   " << &request << ::std::endl;
		::std::cerr << ::std::dec;
	}

    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    bool areadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb &request)
    {
        SEQAN_PROTIMESTART(tw);
        memset(&request, 0, sizeof(aiocb));
        request.aio_fildes = me.handleAsync;
        request.aio_buf = memPtr;
        request.aio_offset = fileOfs;
        request.aio_offset *= sizeof(TValue);
        request.aio_nbytes = count * sizeof(TValue);
        request.aio_sigevent.sigev_notify = SIGEV_NONE;
/*      request.aio_sigevent.sigev_notify = SIGEV_SIGNAL;
        request.aio_sigevent.sigev_signo = _AsyncIOSignal;
        request.aio_sigevent.sigev_value.sival_ptr = &request;
		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_read():");
		#endif
*/      SEQAN_PROADD(SEQAN_PROIO, (request.aio_nbytes + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
		int result = aio_read(&request);
        SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
        #ifdef SEQAN_DEBUG
			if (result) ::std::cerr << "areadAt returned " << result << " and errno is " << ::strerror(errno) << ::std::endl;
		#endif
		return result == 0;
    }
    
    template < typename TSpec, typename TValue, typename TSize, typename TPos >
    bool awriteAt(File<Async<TSpec> > & me, const TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb &request)
    {
        SEQAN_PROTIMESTART(tw);
        memset(&request, 0, sizeof(aiocb));
        request.aio_fildes = me.handleAsync;
        request.aio_buf = const_cast<TValue*>(memPtr);
        request.aio_offset = fileOfs;
        request.aio_offset *= sizeof(TValue);
        request.aio_nbytes = count * sizeof(TValue);
        request.aio_sigevent.sigev_notify = SIGEV_NONE;
/*      request.aio_sigevent.sigev_notify = SIGEV_SIGNAL;
        request.aio_sigevent.sigev_signo = _AsyncIOSignal;
        request.aio_sigevent.sigev_value.sival_ptr = &request;
		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_write():");
		#endif
*/      SEQAN_PROADD(SEQAN_PROIO, (request.aio_nbytes + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
		int result = aio_write(&request);
        SEQAN_PROADD(SEQAN_PROIWAIT, SEQAN_PROTIMEDIFF(tw));
        #ifdef SEQAN_DEBUG
			if (result) ::std::cerr << "awriteAt returned " << result << " and errno is " << ::strerror(errno) << ::std::endl;
		#endif
        return result == 0;
    }

	template <typename TSpec>
    inline bool flush(File<Async<TSpec> > & me) {
		#if _POSIX_SYNCHRONIZED_IO > 0
			return me.handle == me.handleAsync || fdatasync(me.handle) == 0;
		#else
			return me.handle == me.handleAsync || fsync(me.handle) == 0;
		#endif
    }

    //////////////////////////////////////////////////////////////////////
    // queue specific functions

	inline bool waitFor(aiocb &request) {
/*		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_suspend():");
		#endif
*/
		aiocb * cblist = &request;
        SEQAN_PROTIMESTART(tw);
		int result = aio_suspend(&cblist, 1, NULL);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        #ifdef SEQAN_DEBUG
			if (result) {
	 			int eno = aio_error(&request);
				if (eno != EINPROGRESS)
					::std::cerr << "waitFor: aio_error returned " << ::strerror(eno) << " and errno=" << errno << " " << ::strerror(errno) << ::std::endl;
			}
		#endif
		return result == 0;
	}

	inline bool waitFor(aiocb &request, long timeout_millis) {
/*		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_suspend_timeout():");
		#endif
*/
		int result;
		if (timeout_millis == 0)
			result = aio_error(&request);
		else {
			aiocb * cblist = &request;
			timespec ts;
			ts.tv_sec = timeout_millis / 1000;
			ts.tv_nsec = (timeout_millis % 1000) * 1000;
			SEQAN_PROTIMESTART(tw);
			result = aio_suspend(&cblist, 1, &ts);
			SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
		}

        #ifdef SEQAN_DEBUG
			if (result) {
	 			int eno = aio_error(&request);
				if (eno != EINPROGRESS)
					::std::cerr << "waitFor(timeOut=" << timeout_millis << "): aio_error returned " << ::strerror(eno) << " and errno=" << errno << " " << ::strerror(errno) << ::std::endl;
			}
		#endif
        return result == 0;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb const * const contexts[], TSize count) {
        SEQAN_PROTIMESTART(tw);
		bool result = aio_suspend(contexts, count, NULL);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result == 0;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb const * const contexts[], TSize count, long timeout_millis) {
        timespec ts;
        ts.tv_sec = timeout_millis / 1000;
        ts.tv_nsec = (timeout_millis % 1000) * 1000;
        SEQAN_PROTIMESTART(tw);
		bool result = aio_suspend(contexts, count, &ts);
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result == 0;
	}

	template <typename TSpec>
    inline bool cancel(File<Async<TSpec> > & me, aiocb &request) {
/*		#ifdef SEQAN_VVERBOSE
			printRequest(request, "aio_cancel():");
		#endif
*/      return aio_cancel(me.handleAsync, &request) == 0;
    }

    inline int error(aiocb const & request) {
        return aio_error(&request);
    }

    inline int return_value(aiocb & request) {
        return aio_return(&request);
    }

	template <typename TSpec>
    inline void release(File<Async<TSpec> > & /*me*/, aiocb const & /*request*/) {}

/*
    typedef void (*sighandler_t)(int);
    static unsigned _AsyncIOHandlerRefCount = 0;
    static struct sigaction _AsyncIOOldSig;

    inline void _AsyncIOHandler(int sigNo, siginfo_t *info, void *hint) {
        SEQAN_ASSERT(sigNo == _AsyncIOSignal);
        // TODO: signal respective event
        // currently we don't need async IO handlers because
        // we only wait for single events
    }

    static sighandler_t _addAsyncIOHandler() {
        struct sigaction newSig, oldSig;
        newSig.sa_sigaction = _AsyncIOHandler;
        sigemptyset(&newSig.sa_mask);
        newSig.sa_flags = SA_RESTART + SA_SIGINFO;
        if (sigaction(_AsyncIOSignal, &newSig, &oldSig) < 0)
            return SIG_ERR;
        return oldSig.sa_handler;
    }
*/
    
	//////////////////////////////////////////////////////////////////////////////
	// page aligned allocate for direct file io

    struct TagAllocateAligned_;	//< allocate page aligned memory for direct i/o access
    typedef Tag<TagAllocateAligned_> const TagAllocateAligned;

	template <typename T, typename TValue, typename TSize>
	inline void
	allocate(T const & /*me*/, 
			 TValue * & data,
			 TSize count,
			 TagAllocateAligned const)
	{
		data = (TValue *) valloc(count * sizeof(TValue));
        if (data)
			SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
		else
			::std::cerr << "AlignAllocator: Could not allocate memory of size " << ::std::hex << 
				count * sizeof(TValue) << " with page alignment. (ErrNo=" << ::std::dec <<
				errno << ")" << ::std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////
	// page aligned deallocate for direct file io

	template <typename T, typename TValue, typename TSize>
	inline void 
	deallocate( T const & /*me*/,
				TValue * data, 
				TSize count,
				TagAllocateAligned const)
	{
        if (data && count)	// .. to use count if SEQAN_PROFILE is not defined
			SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
		free(data);
	}

#endif

    //////////////////////////////////////////////////////////////////////////////
    // global functions

	template <typename TSpec>
    struct Size< File<Async<TSpec> > >
    {
        typedef typename File<Async<TSpec> >::SizeType Type;
    };

	template <typename TSpec>
    struct Position< File<Async<TSpec> > >
    {
        typedef typename File<Async<TSpec> >::FilePtr Type;
    };

	template <typename TSpec>
    struct Difference< File<Async<TSpec> > >
    {
        typedef typename File<Async<TSpec> >::FilePtr Type;
    };



    template < typename TSpec, typename TValue, typename TSize>
	inline void
	allocate( File<Async<TSpec> > const & me,
			  TValue * & data, 
			  TSize count)
	{
		allocate(me, data, count, TagAllocateAligned());
	}

    template <typename TSpec, typename TValue, typename TSize>
	inline void
	deallocate( File<Async<TSpec> > const & me,
				TValue * data, 
				TSize count)
	{
		deallocate(me, data, count, TagAllocateAligned());
	}


}

#endif
