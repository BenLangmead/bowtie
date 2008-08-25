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
  $Id: file_manual_forwards.h,v 1.1 2008/08/25 16:20:08 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_MANUAL_FORWARDS_H 
#define SEQAN_HEADER_FILE_MANUAL_FORWARDS_H 

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

//____________________________________________________________________________
// TagAllocateAligned_

struct TagAllocateAligned_;       	// "projects/library/seqan/file/file_async.h"(283)
//struct aiocb;

//////////////////////////////////////////////////////////////////////////////
// TYPEDEFS
//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
// (*sighandler_t)(int)

typedef void (*sighandler_t)(int);       	// "projects/library/seqan/file/file_async.h"(258)

//____________________________________________________________________________
// TagAllocateAligned

typedef Tag<TagAllocateAligned_> const TagAllocateAligned;       	// "projects/library/seqan/file/file_async.h"(284)


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
// allocate

template <typename T, typename TValue, typename TSize> inline void allocate(T const & me, TValue * & data, TSize count, TagAllocateAligned const);       	// "projects/library/seqan/file/file_async.h"(292)
template <typename TSpec, typename TValue, typename TSize> inline void allocate( File<Async<TSpec> > const & me, TValue * & data, TSize count);       	// "projects/library/seqan/file/file_async.h"(351)

/*
//____________________________________________________________________________
// areadAt

template <typename TSpec, typename TValue, typename TSize, typename TPos > bool areadAt(File<Async<TSpec> > & me, TValue *memPtr, TSize const count, TPos const fileOfs, aiocb &request);       	// "projects/library/seqan/file/file_async.h"(136)

//____________________________________________________________________________
// awriteAt

template <typename TSpec, typename TValue, typename TSize, typename TPos > bool awriteAt(File<Async<TSpec> > & me, const TValue *memPtr, TSize const count, TPos const fileOfs, aiocb &request);       	// "projects/library/seqan/file/file_async.h"(160)

//____________________________________________________________________________
// cancel

template <typename TSpec> inline bool cancel(File<Async<TSpec> > & me, aiocb &request);       	// "projects/library/seqan/file/file_async.h"(242)

//____________________________________________________________________________
// deallocate

template <typename T, typename TValue, typename TSize> inline void deallocate( T const & me, TValue * data, TSize count, TagAllocateAligned const);       	// "projects/library/seqan/file/file_async.h"(310)
template <typename TSpec, typename TValue, typename TSize> inline void deallocate( File<Async<TSpec> > const & me, TValue * data, TSize count);       	// "projects/library/seqan/file/file_async.h"(360)

//____________________________________________________________________________
// error

inline int error(aiocb const &request);       	// "projects/library/seqan/file/file_async.h"(246)

//____________________________________________________________________________
// fileExists

inline bool fileExists(const char *fileName);       	// "projects/library/seqan/file/file_sync.h"(189)

//____________________________________________________________________________
// fileUnlink

inline bool fileUnlink(const char *fileName);       	// "projects/library/seqan/file/file_sync.h"(194)

//____________________________________________________________________________
// flush

template <typename TSpec> inline bool flush(File<Async<TSpec> > & me);       	// "projects/library/seqan/file/file_async.h"(182)

//____________________________________________________________________________
// printRequest

inline void printRequest(aiocb &request);       	// "projects/library/seqan/file/file_async.h"(125)

//____________________________________________________________________________
// read

template <typename TSpec, typename TValue, typename TSize > inline bool read(File<Sync<TSpec> > & me, TValue *memPtr, TSize const count);       	// "projects/library/seqan/file/file_sync.h"(226)

//____________________________________________________________________________
// release

template <typename TSpec> inline void release(File<Async<TSpec> > & me, aiocb const &request);       	// "projects/library/seqan/file/file_async.h"(255)

//____________________________________________________________________________
// return_value

inline int return_value(aiocb &request);       	// "projects/library/seqan/file/file_async.h"(250)

//____________________________________________________________________________
// waitFor

inline bool waitFor(aiocb &request);       	// "projects/library/seqan/file/file_async.h"(189)
inline bool waitFor(aiocb &request, long timeout_millis);       	// "projects/library/seqan/file/file_async.h"(204)

//____________________________________________________________________________
// waitForAny

template <typename TSize > inline TSize waitForAny(aiocb const * const contexts[], TSize count);       	// "projects/library/seqan/file/file_async.h"(223)
template <typename TSize > inline TSize waitForAny(aiocb const * const contexts[], TSize count, long timeout_millis);       	// "projects/library/seqan/file/file_async.h"(231)

*/

//____________________________________________________________________________
// write

template <typename TSpec, typename TValue, typename TSize > inline bool write(File<Sync<TSpec> > & me, TValue const *memPtr, TSize const count);       	// "projects/library/seqan/file/file_sync.h"(231)

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif

