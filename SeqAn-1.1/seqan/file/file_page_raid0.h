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
  $Id: file_page_raid0.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_PAGE_RAID0_H
#define SEQAN_HEADER_FILE_PAGE_RAID0_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


    //////////////////////////////////////////////////////////////////////////////
    // page based read/write for striped files

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec >
	inline bool 
	readPage(
		int pageNo, 
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		File< Striped<_FileCount, TFile> > &file)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READING;
		return areadAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size(pf), 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf), 
			pf.request);
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec >
	inline bool 
	writePage(
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		int pageNo, 
		File< Striped<_FileCount, TFile> > &file)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.status = pf.WRITING;
		return awriteAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size(pf), 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf), 
			pf.request);
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec, typename TSize >
	inline bool 
	readLastPage(
		int pageNo, 
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		File< Striped<_FileCount, TFile> > &file,
		TSize size)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
		return readAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size, 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf));
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec, typename TSize >
	inline bool 
	writeLastPage(
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		int pageNo, 
		File< Striped<_FileCount, TFile> > &file,
		TSize size)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return writeAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size, 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf));
	}


	//////////////////////////////////////////////////////////////////////////////
	// bucket based read/write methods for striped files

	template < typename TValue, unsigned _FileCount, typename TFile >
	inline unsigned 
	readBucket(
		PageBucket<TValue> &b, 
		int pageNo, 
		unsigned pageSize, 
		unsigned dataSize, 
		File< Striped<_FileCount, TFile> > &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
        unsigned readSize = _min(dataSize - b.pageOfs, (unsigned)(b.end - b.begin));
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readBucket:  " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)(pageNo / _FileCount) * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << readSize << ::std::endl;
		#endif
        if (readSize && readAt(file[pageNo % _FileCount], b.begin, readSize, (pos_t)(pageNo / _FileCount) * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += readSize;
            b.cur = b.begin;
            b.end = b.begin + readSize;
            return readSize;
        } else
            return 0;
	}

	template < typename TValue, unsigned _FileCount, typename TFile >
	inline bool 
	writeBucket(
		PageBucket<TValue> &b,
		int pageNo, 
		unsigned pageSize, 
		File< Striped<_FileCount, TFile> > &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)(pageNo / _FileCount) * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << b.cur - b.begin << ::std::endl;
		#endif
        if ((b.cur == b.begin) || writeAt(file[pageNo % _FileCount], b.begin, b.cur - b.begin, (pos_t)(pageNo / _FileCount) * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += b.cur - b.begin;
            b.cur = b.begin;
            return true;
        } else
            return false;
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec >
	inline bool 
	writeBucket(
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, Dynamic<TSpec> > &pf, 
		unsigned &pageOfs, 
		File< Striped<_FileCount, TFile> > &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << pf.begin;
			::std::cerr << " from page " << ::std::dec << pf.pageNo << " at " << (pos_t)(pf.pageNo / _FileCount) * (pos_t)pageSize(pf) + pageOfs;
			::std::cerr << " size " << size(pf) << ::std::endl;
		#endif
        if (pf.end == pf.begin) return true;
        if (awriteAt(file[pf.pageNo % _FileCount], pf.begin, size(pf), (pos_t)(pf.pageNo / _FileCount) * (pos_t)pageSize(pf) + pageOfs, pf.request)) {
            pf.status = pf.WRITING;
            pageOfs += size(pf);
            return true;
        } else
            return false;
	}

}

#endif
