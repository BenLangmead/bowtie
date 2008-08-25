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
  $Id: file_format_genbank.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_GENBANK_H
#define SEQAN_HEADER_FILE_GENBANK_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Genbank
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Genbank:
	Genbank format for sequences from the Genbank database.
*/
struct TagGenbank_;
typedef Tag<TagGenbank_> const Genbank;


//////////////////////////////////////////////////////////////////////////////
// FileReader Iterator
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Genbank, TFile2, TSpec> > & it, bool skip_meta = true)
{
SEQAN_CHECKPOINT
	String<char> line;

	if (_streamEOF(host(it)))
	{//end of file
		it.data_eof = true;
		return;
	}

	if (skip_meta && (it.data_char != ' '))
	{//skip metadata block
		while (true)
		{
			if (it.data_char == '/')
			{//end of record
				_stream_skipLine(host(it), it.data_char);
				it.data_eof = true;
				return;
			}
			if ((it.data_char == 'O') || (it.data_char == 'o'))
			{
				clear(line);
				_stream_appendLine(host(it), line, it.data_char);
				if ((prefix(line, 6) == "ORIGIN") || (prefix(line, 6) == "origin"))
				{//end of metadata
					break;
				}
			}
			//skip meta line
			_stream_skipLine(host(it), it.data_char);

			if (_streamEOF(host(it)))
			{//end of file
				it.data_eof = true;
				return;
			}
		}
	}

	//find first character
	while (true)
	{
		if (_streamEOF(host(it)))
		{//end of file
			it.data_eof = true;
			return;
		}
		if ((it.data_char != ' ') && ((it.data_char < '0') || (it.data_char > '9')))
		{
			if ((it.data_char != '\n') && (it.data_char != '\r'))
			{//fist char found
				break;
			}

			it.data_char = _streamGet(host(it));
			if (it.data_char == '/')
			{//end of record
				_stream_skipLine(host(it), it.data_char);
				it.data_eof = true;
				return;
			}
		}
		else
		{
			it.data_char = _streamGet(host(it));
		}
	}

	it.data_file_pos = _streamTellG(host(it));
	it.data_file_pos -=1;
	it.data_eof = _streamEOF(host(it));
}


template <typename TFile, typename TFile2, typename TSpec>
inline void
goNext(Iter<TFile, FileReader<Genbank, TFile2, TSpec> > & it)
{
SEQAN_CHECKPOINT
	do
	{
		it.data_char = _streamGet(host(it));
		if (_streamEOF(host(it)))
		{
			it.data_eof = true;
			return;
		}
		it.data_file_pos += 1;

		if ((it.data_char == '\n') || (it.data_char == '\r'))
		{//linebreak detected: find begin of next line
			do
			{
				it.data_char = _streamGet(host(it));
				if (_streamEOF(host(it)))
				{
					it.data_eof = true;
					return;
				}
				it.data_file_pos += 1;
			} while ((it.data_char == '\n') || (it.data_char == '\r'));

			if (it.data_char == '/')
			{//end of record
				_stream_skipLine(host(it), it.data_char);
				_streamUnget(host(it));
				//it.data_file_pos is invalid now
				it.data_eof = true;
				return;
			}
		}
	} while ((it.data_char == ' ') || ((it.data_char >= '0') && (it.data_char <= '9')));
}

//////////////////////////////////////////////////////////////////////////////
// File Format Access Function
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
read(TFile & file,
	 TData & data,
	 Genbank)
{
SEQAN_CHECKPOINT
	Iter<TFile, FileReader<Genbank> > it(file);

	clear(data);
	while (!atEnd(it))
	{
		appendValue(data, getValue(it));
		goNext(it);
	}
}

template <typename TFile, typename TData, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Genbank)
{
SEQAN_CHECKPOINT
	typename Size<TData>::Type siz = length(data);
	Iter<TFile, FileReader<Genbank> > it(file);

	clear(data);
	while (!atEnd(it) && (siz < limit))
	{
		appendValue(data, getValue(it));
		goNext(it);
	}
	while (!atEnd(it))
	{
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TMeta>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 Genbank)
{
SEQAN_CHECKPOINT
	typedef typename Value<TMeta>::Type TValue;
	String<char> line;

	clear(meta);

	if (_streamEOF(file))
	{
		return;
	}

	TValue c = _streamGet(file);

	while (!_streamEOF(file))
	{
		clear(line);
		_stream_appendLine(file, line, c);

		if (c == '/')
		{//end of record
			_streamUnget(file);
			break;
		}

		append(meta, line);
		appendValue(meta, '\n');

		if ((prefix(line, 6) == "ORIGIN") || (prefix(line, 6) == "origin"))
		{//end of metadata
			_streamUnget(file);
			break;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
inline void
goNext(TFile & file,
	   Genbank)
{
SEQAN_CHECKPOINT
	typedef typename Value<TFile>::Type TValue;

	if (_streamEOF(file))
	{
		return;
	}

	while (!_streamEOF(file))
	{
		TValue c = _streamGet(file);
		if (c == '/')
		{//end of record
			_stream_skipLine(file, c);
			_streamUnget(file);
			return;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
/*
template <typename TFile>
inline void
length(TFile & file,
	   Genbank)
{
SEQAN_CHECKPOINT
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
write(TFile & file,
	  TData & data,
	  Genbank)
{
SEQAN_CHECKPOINT
	enum
	{
		BLOCK_SIZE = 10,
		BLOCKS_PER_LINE = 6
	};
	char const * NUM_BLOCK_FORMAT = "%9d";

	typedef typename Size<TData>::Type TSize;
	typedef typename Iterator<TData, Standard>::Type TIterator;

	TSize count = 0;
	TIterator it = begin(data, Standard());
	TIterator it_end = end(data, Standard());

	while (it != it_end)
	{
		//write count 
		_streamPutInt(file, count+1, NUM_BLOCK_FORMAT);

		int block_count = 0;
		int char_in_block_count = BLOCK_SIZE;

		//write rest of line
		while (it != it_end)
		{
			if (char_in_block_count == BLOCK_SIZE)
			{//begin new block
				if (block_count >= BLOCKS_PER_LINE)
				{//end of line
					_streamPut(file, '\n');
					break;
				}
				_streamPut(file, ' ');
				char_in_block_count = 0;
				++block_count;
			}

			//write next character
			_streamPut(file, *it);
			++it;
			++count;
			++char_in_block_count;
		}
	}

	write(file, "\n//\n");
}

template <typename TFile, typename TData, typename TMeta>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  Genbank)
{
SEQAN_CHECKPOINT
	write(file, meta);
	write(file, data, Genbank());
}


//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
