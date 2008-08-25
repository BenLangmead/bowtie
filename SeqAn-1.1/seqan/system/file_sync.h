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
  $Id: file_sync.h,v 1.1 2008/08/25 16:20:08 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_FILE_SIMPLE_H
#define SEQAN_HEADER_FILE_SIMPLE_H

#include <fcntl.h>          // O_CREAT ..
#include <sys/stat.h>       // 
#include <cstdio>           // tmpnam(..)

#ifdef PLATFORM_WINDOWS
# include <io.h>            // read(..) ..
#else
# include <cstdlib>
# include <cerrno>
# include <unistd.h>
#endif


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


	template <typename TSpec /* = void */>
	struct Sync;


#ifdef PLATFORM_WINDOWS

    //////////////////////////////////////////////////////////////////////////////
    // Windows rtl file access
	template <typename TSpec>
	class File<Sync<TSpec> >
    {
    public:

		typedef __int64			FilePtr;
		typedef __int64         SizeType;   // type of file size
        typedef unsigned int    _SizeType;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handle;

        File(void *dummy = NULL): // to be compatible with the FILE*(NULL) constructor
            handle(-1) {}

        inline int _getOFlag(int openMode) const {
			int result;

			switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    result = _O_RDONLY;
					break;
                case OPEN_WRONLY:
                    result = _O_WRONLY;
					break;
                case OPEN_RDWR:
                    result = _O_RDWR;
					break;
			}

			if (openMode & OPEN_CREATE)     result |= _O_CREAT;
			if (!(openMode & OPEN_APPEND))	result |= _O_TRUNC;
            if (openMode & OPEN_TEMPORARY)  result |= _O_TEMPORARY;
			return result | _O_BINARY;
        }

        bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) {
            handle = _open(fileName, _getOFlag(openMode), _S_IREAD | _S_IWRITE);
			if (handle == -1) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Open failed on file " << fileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
				return false;
			}
			SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
            return true;
        }

        bool openTemp(int openMode = DefaultOpenTempMode<File>::VALUE) {
#ifdef SEQAN_DEFAULT_TMPDIR
			char *fileName = _tempnam(SEQAN_DEFAULT_TMPDIR, "GNDX");
#else
			char *fileName = _tempnam(NULL, "GNDX");
#endif
			if (!fileName) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Cannot create a unique temporary filename" << ::std::endl;
				return false;
			}
            bool result = open(fileName, openMode | OPEN_TEMPORARY);
			free(fileName);
			return result;
        }

        inline bool close() {
            if (_close(handle) != 0)
                return false;
            handle = -1;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return true;
        }

		inline int read(void *buffer, _SizeType count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    int result = _read(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline int write(void const *buffer, _SizeType count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    int result = _write(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline FilePtr seek(FilePtr pos, int origin = SEEK_SET) const {
			return _lseeki64(handle, pos, origin);
		}

		inline FilePtr tell() const {
			return _telli64(handle);
		}

		static int error() {
			return errno;
		}

        operator bool () const {
            return handle != -1;
        }
    };

	inline bool fileExists(const char *fileName) {
		struct _stat buf;
		return _stat(fileName, &buf) == 0;
	}

	inline bool fileUnlink(const char *fileName) {
		return _unlink(fileName) == 0;
	}

#else

    //////////////////////////////////////////////////////////////////////////////
    // Unix file access
	template <typename TSpec>
	class File<Sync<TSpec> >
    {
    public:

		typedef off_t			FilePtr;
		typedef off_t           SizeType;   // type of file size
        typedef size_t          _SizeType;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handle;

        File(void */*dummy*/ = NULL): // to be compatible with the FILE*(NULL) constructor
            handle(-1) {}

        inline int _getOFlag(int openMode) const {
			int result = O_LARGEFILE;

			switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    result |= O_RDONLY;
					break;
                case OPEN_WRONLY:
                    result |= O_WRONLY;
					if (!(openMode & OPEN_APPEND))	result |= O_TRUNC;
					break;
                case OPEN_RDWR:
                    result |= O_RDWR;
					if (!(openMode & OPEN_APPEND))	result |= O_TRUNC;
					break;
			}

			if (openMode & OPEN_CREATE)     result |= O_CREAT;
//			if (openMode & OPEN_TEMPORARY)  result |= O_TEMPORARY;
        #ifdef SEQAN_DIRECTIO
    		if (openMode & OPEN_ASYNC)		result |= O_DIRECT;
        #endif
			return result;
        }

        bool open(char const *fileName, int openMode = DefaultOpenMode<File>::VALUE) {
            handle = ::open(fileName, _getOFlag(openMode), S_IREAD | S_IWRITE);
			if (handle == -1 && errno == EINVAL) {	// fall back to cached access
	            #ifdef SEQAN_DEBUG_OR_TEST_
					if (!(openMode & OPEN_QUIET))
						::std::cerr << "Warning: Direct access openening failed: " << fileName << "." << ::std::endl;
				#endif			
          	    handle = ::open(fileName, _getOFlag(openMode & ~OPEN_ASYNC), S_IREAD | S_IWRITE);
			}
			
			if (handle == -1) {
				if (!(openMode & OPEN_QUIET))
					::std::cerr << "Open failed on file " << fileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
				return false;
			}

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
			char tmpFileName[] = "/GNDtmp";
			//if ((handle = ::mkstemp(tmpFileName)) == -1) {
			//	if (!(openMode & OPEN_QUIET))
			//		::std::cerr << "Cannot create temporary file " << tmpFileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
			//	return false;
			//}
			if (!(close() && open(tmpFileName, openMode))) return false;
			#ifdef SEQAN_DEBUG
				int result = 
            #endif
			::unlink(tmpFileName);
            #ifdef SEQAN_DEBUG
				if (result == -1 && !(openMode & OPEN_QUIET))
					::std::cerr << "Cannot unlink temporary file " << tmpFileName << ". (" << ::strerror(errno) << ")" << ::std::endl;
            #endif
			return true;
        }

        inline bool close() {
            if (::close(this->handle) == -1) return false;
            handle = -1;
			SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
            return true;
        }

		inline ssize_t read(void *buffer, _SizeType count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    ssize_t result = ::read(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline ssize_t write(void const *buffer, _SizeType count) const {
            SEQAN_PROADD(SEQAN_PROIO, (count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    ssize_t result = ::write(handle, buffer, count);
            SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline FilePtr seek(FilePtr pos, int origin = SEEK_SET) const {
            FilePtr result = ::lseek(handle, pos, origin);
//			#ifdef SEQAN_DEBUG
				if (result < 0)
					::std::cerr << "lseek returned " << result << ". (" << ::strerror(errno) << ")" << ::std::endl;
//			#endif
			return result;
		}

		inline FilePtr tell() const {
            return seek(0, SEEK_CUR);
        }

		static int error() {
            return errno;
		}

        operator bool () const {
            return handle != -1;
        }
    };

	inline bool fileExists(const char *fileName) {
		struct stat buf;
		return stat(fileName, &buf) != -1;
	}

	inline bool fileUnlink(const char *fileName) {
		return unlink(fileName) == 0;
	}

#endif

    //////////////////////////////////////////////////////////////////////////////
    // global functions

	template <typename TSpec>
    struct Size< File<Sync<TSpec> > >
    {
        typedef typename File<Sync<TSpec> >::SizeType Type;
    };

	template <typename TSpec>
    struct Position< File<Sync<TSpec> > >
    {
        typedef typename File<Sync<TSpec> >::FilePtr Type;
    };

	template <typename TSpec>
    struct Difference< File<Sync<TSpec> > >
    {
        typedef typename File<Sync<TSpec> >::FilePtr Type;
    };

    template < typename TSpec, typename TValue, typename TSize >
    inline bool read(File<Sync<TSpec> > & me, TValue *memPtr, TSize const count) {
		return (int) me.read(memPtr, count * sizeof(TValue)) == (int) (count * sizeof(TValue));
    }
    
    template < typename TSpec, typename TValue, typename TSize >
    inline bool write(File<Sync<TSpec> > & me, TValue const *memPtr, TSize const count) {
		return (int) me.write(memPtr, count * sizeof(TValue)) == (int) (count * sizeof(TValue));
    }

}

#endif
