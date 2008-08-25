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
  $Id: file_filereaderiterator.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_FILEREADEITERATOR_H
#define SEQAN_HEADER_FILE_FILEREADEITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TFormat, typename TFile = FILE*, typename TSpec = Default>
struct FileReader;

//////////////////////////////////////////////////////////////////////////////
// FileReader: an iterator that scans through the data of a file
// note: this is not the iterator of the FileReader string (see file_filereader.h)
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
class Iter<TFile, FileReader<TFormat, TFile2, TSpec> >
{
public:
	typedef typename Value<TFile>::Type TValue;
	typedef typename Position<TFile>::Type TPosition;

	TPosition data_file_pos;	//position of the last read char relative to data begin (only valid if data_eof == false)
	TFile * data_host;		//the host file
	TValue data_char;		//the last read char
	bool data_eof;			//true if reached end of record
//	TFilePosition data_begin_pos;

	Iter(TFile & file_, bool skip_meta = true):
		data_file_pos(0),
		data_host(& file_),
		data_eof(false)
	{
		data_char = _streamGet(file_);
		goBegin(*this, skip_meta);
	}
	Iter(Iter const & other_):
		data_file_pos(other_.data_file_pos),
		data_host(other_.data_host),
		data_char(other_.data_char),
		data_eof(other_.data_eof)
	{
	}
	~Iter() 
	{
	}

	Iter const &
	operator = (Iter const & other_)
	{
		data_file_pos = other_.data_file_pos;
		data_host = other_.data_host;
		data_char = other_.data_char;
		data_eof = other_.data_eof;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
struct Value< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >:
	Value<TFile>
{
};

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
struct GetValue< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >
{
	typedef typename Value< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type Type;
};

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
struct Reference< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >
{
	typedef typename Value< Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type & Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline TFile &
host(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
	return *(it.data_host);
}


template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline typename Reference<Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type
value(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
	return it.data_char;
}

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline typename GetValue<Iter<TFile, FileReader<TFormat, TFile2, TSpec> > >::Type
getValue(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
	return it.data_char;
}

template <typename TFile, typename TFormat, typename TFile2, typename TSpec>
inline bool
atEnd(Iter<TFile, FileReader<TFormat, TFile2, TSpec> > & it)
{
	return it.data_eof;
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
