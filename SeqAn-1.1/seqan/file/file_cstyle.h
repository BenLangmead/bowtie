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
  $Id: file_cstyle.h,v 1.1 2008/08/25 16:20:04 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_CSTYLE_H
#define SEQAN_HEADER_FILE_CSTYLE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/*    template <>
    struct Value< FILE* >
    {
	    typedef unsigned char Type;
    };
*/
/*  already defined in "seqan/cstream.h" ...

    template <>
    struct Position< FILE* >
    {
	    typedef long Type;
    };
*/
/*
    template <>
    struct Size< FILE* >
    {
	    typedef size_t Type;
    };
*/
/*

    template <>
    struct Position< FILE* >
    {
	    typedef long Type;
    };
*/
    template <>
    struct Difference< FILE* >
    {
	    typedef long Type;
    };


	inline const char * _getCStyleOpenMode(int openMode) {
		switch (openMode & OPEN_MASK) {
            case OPEN_WRONLY:
                if (!(openMode & OPEN_APPEND))
                    if (openMode & OPEN_CREATE)
                        return "w";
                    else
                        return "r+";
                else
                    return "a";
            case OPEN_RDWR:
                if (!(openMode & OPEN_APPEND))
                    if (openMode & OPEN_CREATE)
                        return "w+";
                    else
                        return "r+";
                else
                    return "a+";
				break;
		}
        return "r";
    }

    inline bool open(FILE* &me, const char *fileName, int openMode) {
		SEQAN_PROADD(SEQAN_PROOPENFILES, 1);
        return (me = fopen(fileName, _getCStyleOpenMode(openMode))) != NULL;
    }

    inline bool open(FILE* &me, const char *fileName) {
		return open(me, fileName, DefaultOpenMode<FILE*>::VALUE);
	}

    inline bool openTemp(FILE* &me) {
		SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
        return (me = tmpfile()) != NULL;
    }

    inline bool close(FILE* me) {
		SEQAN_PROSUB(SEQAN_PROOPENFILES, 1);
        return fclose(me) == 0;
    }

    inline unsigned sectorSize(FILE* const &) {
        return 4096;
    }

    template < typename TPos >
    inline Size<FILE*>::Type seek(FILE* me, TPos const fileOfs, int origin) {
        fseek(me, fileOfs, origin);
		return ftell(me);
    }
    template < typename TPos >
    inline Size<FILE*>::Type seek(FILE* me, TPos const fileOfs) {
		return seek(me, fileOfs, SEEK_BEGIN);
    }

    inline Size<FILE*>::Type tell(FILE* me) {
		return ftell(me);
    }

    template < typename TValue, typename TSize >
    inline bool read(FILE* me, TValue *memPtr, TSize const count) {
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fread(memPtr, sizeof(TValue), count, me) == count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    template < typename TValue, typename TSize >
    inline bool write(FILE* me, TValue const *memPtr, TSize const count) {
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fwrite(memPtr, sizeof(TValue), count, me) == count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    template < typename TValue, typename TSize, typename TPos >
    inline bool readAt(FILE* me, TValue *memPtr, TSize const count, TPos const fileOfs) {
		typedef typename Position<FILE*>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fread(memPtr, sizeof(TValue), count, me) == count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }
    
    template < typename TValue, typename TSize, typename TPos >
    inline bool writeAt(FILE* me, TValue const *memPtr, TSize const count, TPos const fileOfs) {
		typedef typename Position<FILE*>::Type pos_t;
		seek(me, (pos_t)fileOfs * (pos_t)sizeof(TValue));
        SEQAN_PROADD(SEQAN_PROIO, (sizeof(TValue) * count + SEQAN_PROPAGESIZE - 1) / SEQAN_PROPAGESIZE);
        SEQAN_PROTIMESTART(tw);
        bool result = fwrite(memPtr, sizeof(TValue), count, me) == count;
        SEQAN_PROADD(SEQAN_PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
    }

    inline Size<FILE*>::Type size(FILE* me) {
        Size<FILE*>::Type old_pos = tell(me);
        Size<FILE*>::Type result = 0;
        if (seek(me, 0, SEEK_END) == 0)
            result = tell(me);
        seek(me, old_pos, SEEK_BEGIN);
        return result;
    }

    template < typename TSize >
    inline void resize(FILE* me, TSize new_length) {
        Size<FILE*>::Type old_pos = tell(me);
        seek(me, new_length, SEEK_BEGIN);
        seek(me, old_pos, SEEK_BEGIN);
    }

	inline bool flush(FILE*) {
		return true; 
	}

    template < typename aRequest >
	inline void release(FILE*, aRequest &) {
	}

}

#endif
