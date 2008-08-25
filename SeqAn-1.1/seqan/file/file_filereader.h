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
  $Id: file_filereader.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_FILEREADER_H
#define SEQAN_HEADER_FILE_FILEREADER_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// FileReader String
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.File Reader String:
..cat:Strings
..cat:Files
..general:Class.String
..summary:Read sequence data from file.
..signature:String<TValue, FileReader<TFormat, TFile, TSpec> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...metafunction:Metafunction.Value
..param.TFormat:A file format.
...type:Tag.File Format
..param.TFile:A file.
..param.TSpec:A further specializing type.
...default:@Tag.Default@
*/


template <typename TValue, typename TFormat, typename TFile, typename TSpec>
class String<TValue, FileReader<TFormat, TFile, TSpec> >
{
public:
	enum
	{
		BLOCK_SIZE = 0x1000
	};

	typedef typename Position<TFile>::Type TFilePosition;

	typedef typename Size<TFile>::Type TFileSize;
	typedef String<TFileSize> TABL;

	typedef typename Position<TABL>::Type TABLPosition;

	typedef String<TValue> TBuf;

	Holder<TFile, Tristate2> data_file;
	TFilePosition data_file_begin;		//file pointer to begin of data in file
	TABL data_abl;						//accumulated block lengths
	TABLPosition data_active_block;		//number of active block
	TFileSize data_active_block_begin;	//begin position of active block
	TFileSize data_active_block_end;	//end position of active block
	TBuf data_buf;						//data of active block
	bool data_scanned;					//true if the complete string was scanned


	//TODO
	//String()...

	String(TFile & fl_)
		: data_scanned(false)
	{
		reserve(data_buf, (size_t) BLOCK_SIZE, Exact());

		setValue(data_file, fl_);
		_FileReaderString_construct(*this);
	}
	template <typename TString>
	String(TString const & str_)
		: data_scanned(false)
	{
		reserve(data_buf, (size_t) BLOCK_SIZE, Exact());

		create(data_file);
		if (!_streamOpen(value(data_file), str_))
		{
			clear(data_file);
		}
		_FileReaderString_construct(*this);
	}
	~String()
	{
		if (!dependent(data_file) && !empty(data_file))
		{
			_streamClose(value(data_file));
			clear(data_file);
		}
	}
};

//////////////////////////////////////////////////////////////////////////////

//template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TIteratorSpec>
//struct Iterator<String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec>;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
struct Value<String<TValue, FileReader<TFormat, TFile, TSpec> > >
{
	typedef TValue Type;
};

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
struct GetValue<String<TValue, FileReader<TFormat, TFile, TSpec> > >
{
	typedef TValue const Type;
};
template <typename TValue, typename TFormat, typename TFile, typename TSpec>
struct Reference<String<TValue, FileReader<TFormat, TFile, TSpec> > >
{
	typedef TValue const Type;
};

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
struct Size<String<TValue, FileReader<TFormat, TFile, TSpec> > >:
	Size<TFile>
{
};
template <typename TValue, typename TFormat, typename TFile, typename TSpec>
struct Difference<String<TValue, FileReader<TFormat, TFile, TSpec> > >:
	Difference<TFile>
{
};
template <typename TValue, typename TFormat, typename TFile, typename TSpec>
struct Position<String<TValue, FileReader<TFormat, TFile, TSpec> > >:
	Position<TFile>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline TFile &
_dataFile(String<TValue, FileReader<TFormat, TFile, TSpec> > & me)
{
	return value(me.data_file);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPosition>
inline void
_FileReaderString_loadblock(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
							TPosition blocknum)
{
	typedef String<TValue, FileReader<TFormat, TFile, TSpec> > TString;

	typedef typename Size<TFile>::Type TFileSize;
	typedef String<TFileSize> TABL;
	typedef typename Position<TABL>::Type TABLPosition;
	TABLPosition blocknum2 = blocknum;

	if (blocknum2 > length(me.data_abl))
	{
		if (me.data_scanned)
		{
			_FileReaderString_loadblock(me, length(me.data_abl) - 1);
		}
		else
		{
			for (TABLPosition bp = length(me.data_abl); !me.data_scanned && (bp <= blocknum2); ++bp)
			{
				_FileReaderString_loadblock(me, bp);
			}
		}
	}
	else
	{
		typedef Iter<TFile, FileReader<TFormat> > TFileReaderIt;
		_streamSeekG(_dataFile(me), me.data_file_begin + blocknum2 * TString::BLOCK_SIZE);
		TABLPosition end_filepos = me.data_file_begin + (blocknum2 + 1) * TString::BLOCK_SIZE;
		TFileReaderIt fit(_dataFile(me), false);

		unsigned int len;
		for (len = 0; !atEnd(fit) && (fit.data_file_pos < end_filepos); goNext(fit))
		{
			me.data_buf[len] = value(fit);
			++len;
		}
		resize(me.data_buf, len);

		if (blocknum2 == length(me.data_abl))
		{
			if (blocknum2 > 0)
			{
				appendValue(me.data_abl, me.data_abl[blocknum2-1] + length(me.data_buf));
			}
			else
			{
				appendValue(me.data_abl, length(me.data_buf));
			}

			if (atEnd(fit))
			{
				me.data_scanned = true;
			}
		}

		me.data_active_block = blocknum2;
		me.data_active_block_begin = (blocknum2) ? me.data_abl[blocknum2 - 1] : 0;
		me.data_active_block_end = me.data_abl[blocknum2];
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPosition>
inline unsigned int
_FileReaderString_findblock(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
							TPosition pos)
{
	typedef typename Size<TFile>::Type TFileSize;

	while (!me.data_scanned && (me.data_abl[length(me.data_abl) - 1] <= pos))
	{
		_FileReaderString_loadblock(me, length(me.data_abl));
	}

	if (pos >= me.data_abl[length(me.data_abl) - 1])
	{//pos greater than file length
		return length(me.data_abl);
	}

	return ::std::lower_bound(begin(me.data_abl, Standard()), end(me.data_abl, Standard()), (TFileSize) pos) - begin(me.data_abl, Standard());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline void
_FileReaderString_construct(String<TValue, FileReader<TFormat, TFile, TSpec> > & me)
{
	//find begin of data in file
	typedef Iter<TFile, FileReader<TFormat> > TFileReaderIt;
	TFileReaderIt fit(_dataFile(me));
	me.data_file_begin = fit.data_file_pos;

	_FileReaderString_loadblock(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

//tests whether block_number will exist when the file was scanned completely

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TUint>
inline bool
_FileReaderString_isValidBlock(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
							   TUint block_number)
{
	typedef typename Size<TFile>::Type TFileSize;
	TFileSize block_number2 = block_number;

	while (!me.data_scanned && (length(me.data_abl) <= block_number2))
	{
		_FileReaderString_loadblock(me, length(me.data_abl));
	}
	return (length(me.data_abl) > block_number2);

}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline void
_FileReaderString_loadComplete(String<TValue, FileReader<TFormat, TFile, TSpec> > & me)
{
	if (!me.data_scanned)
	{//scan the whole sequence
		typedef typename Position<TFile>::Type TPosition;
		_FileReaderString_loadblock(me, supremumValue<TPosition>());
	}
}

//////////////////////////////////////////////////////////////////////////////
 
template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline void const * 
id(String<TValue, FileReader<TFormat, TFile, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	return &me;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPos>
inline TValue
value(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
	  TPos pos)
{
	typedef typename Size<TFile>::Type TFileSize;
	TFileSize pos2 = pos;

	if ((me.data_active_block_begin > pos2) || (me.data_active_block_end <= pos2))
	{//change block
		_FileReaderString_loadblock(me, _FileReaderString_findblock(me, pos2));
	}
	return me.data_buf[pos2 - me.data_active_block_begin];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline typename Size< String<TValue, FileReader<TFormat, TFile, TSpec> > >::Type
length(String<TValue, FileReader<TFormat, TFile, TSpec> > & me)
{
	_FileReaderString_loadComplete(me);

	return me.data_abl[length(me.data_abl) - 1];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TIteratorSpec>
inline typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec >::Type
begin(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
	  Tag<TIteratorSpec> const)
{
	typedef typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec >::Type TIterator;
	return TIterator(me);
}
template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TIteratorSpec>
inline typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> > const, TIteratorSpec >::Type
begin(String<TValue, FileReader<TFormat, TFile, TSpec> > const & me,
	  Tag<TIteratorSpec> const)
{
	typedef String<TValue, FileReader<TFormat, TFile, TSpec> > TString;
	typedef typename Iterator< TString const, TIteratorSpec >::Type TIterator;
	return TIterator(const_cast<TString &>(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TIteratorSpec>
inline typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec >::Type
end(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
	Tag<TIteratorSpec> const)
{
	typedef typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec >::Type TIterator;
	return TIterator(me, GoEnd());
}
template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TIteratorSpec>
inline typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> > const, TIteratorSpec >::Type
end(String<TValue, FileReader<TFormat, TFile, TSpec> > const & me,
	Tag<TIteratorSpec> const)
{
	typedef String<TValue, FileReader<TFormat, TFile, TSpec> > TString;
	typedef typename Iterator< TString const, TIteratorSpec >::Type TIterator;
	return TIterator(const_cast<TString &>(me), GoEnd());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPosition, typename TIteratorSpec>
inline typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec >::Type
iterator(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
		 TPosition pos,
		 Tag<TIteratorSpec> const)
{
	typedef typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec >::Type TIterator;
	return TIterator(me, pos);
}
template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPosition, typename TIteratorSpec>
inline typename Iterator< String<TValue, FileReader<TFormat, TFile, TSpec> > const, TIteratorSpec >::Type
iterator(String<TValue, FileReader<TFormat, TFile, TSpec> > const & me,
		 TPosition pos,
		 Tag<TIteratorSpec> const)
{
	typedef String<TValue, FileReader<TFormat, TFile, TSpec> > TString;
	typedef typename Iterator< TString const, TIteratorSpec >::Type TIterator;
	return TIterator(const_cast<TString &>(me), pos);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Iterator for FileReader String
// (note: do not confuse with FileReader Iterator, see file_filereaderiterator.h)
//////////////////////////////////////////////////////////////////////////////

struct FileReaderIterator;

//helper meta function for storing types associated with file reader string.
//Due to a bug in VC++, _FileReader_Types is instantiated for arbitrary TContainer types
//when instantiating Iter<TContainer, FileReaderIterator>
//"ABL" = "active buffer lengths" table
template <typename T>
struct _FileReader_Types
{// dummy implementation to make VC++ happy
	typedef int TABLPosition;
	typedef int TBuf;
};
template <typename TValue, typename TFormat, typename TFile, typename TSpec>
struct _FileReader_Types<String<TValue, FileReader<TFormat, TFile, TSpec> > >
{
	typedef typename Size<TFile>::Type TFileSize;
	typedef String<TFileSize> TABL;
	typedef typename Position<TABL>::Type TABLPosition;
	typedef String<TValue> TBuf;
};


template <typename TContainer>
class Iter<TContainer, FileReaderIterator>
{
public:
	typedef typename _FileReader_Types<TContainer>::TABLPosition TABLPosition;
	typedef typename _FileReader_Types<TContainer>::TBuf TBuf;

	typedef typename Position<TBuf>::Type TBufPosition;
	typedef typename Size<TBuf>::Type TBufSize;

	TContainer * data_container;
	TABLPosition data_abl_pos;		//number of block
	TBufPosition data_buf_pos;		//number of char in block
	TBufSize data_buf_len;			//length of block
	bool data_atEnd;				//true if iterator is atEnd

	Iter(TContainer & cont_)
		: data_container(& cont_)
	{
		goBegin(*this);
	}

	Iter(TContainer & cont_, GoEnd)
		: data_container(& cont_)
	{
		goEnd(*this);
	}

	template <typename TPos>
	Iter(TContainer & cont_, TPos pos_)
		: data_container(& cont_)
	{
		setPosition(*this, pos_);
	}

	Iter(Iter const & other_)
		: data_container(other_.data_container)
		, data_abl_pos(other_.data_abl_pos)
		, data_buf_pos(other_.data_buf_pos)
		, data_buf_len(other_.data_buf_len)
		, data_atEnd(other_.data_atEnd)
	{
	}

	Iter &
	operator = (Iter const & other_)
	{
		data_container = other_.data_container;
		data_abl_pos = other_.data_abl_pos;
		data_buf_pos = other_.data_buf_pos;
		data_buf_len = other_.data_buf_len;
		data_atEnd = other_.data_atEnd;

		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TIteratorSpec>
struct Iterator<String<TValue, FileReader<TFormat, TFile, TSpec> >, TIteratorSpec>
{
	typedef Iter<String<TValue, FileReader<TFormat, TFile, TSpec> >, FileReaderIterator> Type;
};
template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TIteratorSpec>
struct Iterator<String<TValue, FileReader<TFormat, TFile, TSpec> > const, TIteratorSpec>
{
	typedef Iter<String<TValue, FileReader<TFormat, TFile, TSpec> >, FileReaderIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
struct Value<Iter<TContainer, FileReaderIterator> >:
	Value<TContainer>
{
};
template <typename TContainer>
struct GetValue<Iter<TContainer, FileReaderIterator> >:
	Value<TContainer>
{
};
template <typename TContainer>
struct Reference<Iter<TContainer, FileReaderIterator> >:
	Value<TContainer>
{
};

template <typename TContainer>
struct Size<Iter<TContainer, FileReaderIterator> >:
	Size<TContainer>
{
};
template <typename TContainer>
struct Difference<Iter<TContainer, FileReaderIterator> >:
	Difference<TContainer>
{
};
template <typename TContainer>
struct Position<Iter<TContainer, FileReaderIterator> >:
	Position<TContainer>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline typename GetValue<Iter<TContainer, FileReaderIterator> >::Type
getValue(Iter<TContainer, FileReaderIterator> & it)
{
	TContainer & cont = *(it.data_container);
	if (cont.data_active_block != it.data_abl_pos)
	{
		_FileReaderString_loadblock(cont, it.data_abl_pos);
	}
	return cont.data_buf[it.data_buf_pos];
}
template <typename TContainer>
inline typename GetValue<Iter<TContainer, FileReaderIterator> >::Type
getValue(Iter<TContainer, FileReaderIterator> const & it)
{
	TContainer & cont = *(it.data_container);
	if (cont.data_active_block != it.data_abl_pos)
	{
		_FileReaderString_loadblock(cont, it.data_abl_pos);
	}
	return cont.data_buf[it.data_buf_pos];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline typename Reference<Iter<TContainer, FileReaderIterator> >::Type
value(Iter<TContainer, FileReaderIterator> & it)
{
	return getValue(it);
}
template <typename TContainer>
inline typename Reference<Iter<TContainer, FileReaderIterator> >::Type
value(Iter<TContainer, FileReaderIterator> const & it)
{
	return getValue(it);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline TContainer &
container(Iter<TContainer, FileReaderIterator> & it)
{
	return *(it.data_container);
}
template <typename TContainer>
inline TContainer &
container(Iter<TContainer, FileReaderIterator> const & it)
{
	return *(it.data_container);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline typename Position<Iter<TContainer, FileReaderIterator> >::Type
position(Iter<TContainer, FileReaderIterator> const & it)
{
	TContainer & cont = *(it.data_container);
	if (it.data_atEnd)
	{
		return length(cont);
	}
	else
	{
		if (it.data_abl_pos)
		{
			return cont.data_abl[it.data_abl_pos-1] + it.data_buf_pos;
		}
		else
		{
			return it.data_buf_pos;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TPos>
inline void
setPosition(Iter<TContainer, FileReaderIterator> & it,
			TPos pos)
{
	TContainer & cont = *(it.data_container);
	it.data_abl_pos = _FileReaderString_findblock(cont, pos);
	it.data_atEnd = (cont.data_scanned && (pos >= length(cont)));
	if (it.data_atEnd)
	{
		it.data_buf_pos = 0;
		it.data_buf_len = 0;
	}
	else
	{
		if (it.data_abl_pos == 0)
		{
			it.data_buf_pos = pos;
		}
		else
		{
			it.data_buf_pos = pos - cont.data_abl[it.data_abl_pos - 1];
		}
		if (cont.data_active_block != it.data_abl_pos)
		{
			_FileReaderString_loadblock(cont, it.data_abl_pos);
		}
		it.data_buf_len = length(cont.data_buf);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goNext(Iter<TContainer, FileReaderIterator> & it)
{
	++it.data_buf_pos;
	if (it.data_buf_pos >= it.data_buf_len)
	{
		if (!it.data_atEnd)
		{
			it.data_buf_pos = 0;
			++it.data_abl_pos;

			TContainer & cont = *(it.data_container);
			it.data_atEnd = !_FileReaderString_isValidBlock(cont, it.data_abl_pos);
			if (it.data_atEnd)
			{//at end
				it.data_buf_len = 0;
			}
			else
			{//not at end
				if (cont.data_active_block != it.data_abl_pos)
				{
					_FileReaderString_loadblock(cont, it.data_abl_pos);
				}
				it.data_buf_len = length(cont.data_buf);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goPrevious(Iter<TContainer, FileReaderIterator> & it)
{
	if (it.data_buf_pos > 0)
	{
		--it.data_buf_pos;
	}
	else
	{
		TContainer & cont = *(it.data_container);
		if (it.data_atEnd)
		{
			//_FileReaderString_loadComplete(cont);// length(cont) will do it
			it.data_atEnd = (length(cont) == 0);
			it.data_abl_pos = length(cont.data_abl) - 1;
		}
		else if (it.data_abl_pos)
		{
			--it.data_abl_pos;
		}
		else return; //already in begin pos

		if (cont.data_active_block != it.data_abl_pos)
		{
			_FileReaderString_loadblock(cont, it.data_abl_pos);
		}
		it.data_buf_len = length(cont.data_buf);
		it.data_buf_pos = it.data_buf_len - 1;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goBegin(Iter<TContainer, FileReaderIterator> & it)
{
	it.data_abl_pos = 0;
	it.data_buf_pos = 0;

	TContainer & cont = *(it.data_container);
	if (_FileReaderString_isValidBlock(cont, 0))
	{
		it.data_buf_len = cont.data_abl[0];
		it.data_atEnd = false;
	}
	else
	{
		it.data_buf_len = 0;
		it.data_atEnd = true;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goEnd(Iter<TContainer, FileReaderIterator> & it)
{
	it.data_abl_pos = 0;
	it.data_buf_pos = 0;
	it.data_buf_len = 0;
	it.data_atEnd = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool
atEnd(Iter<TContainer, FileReaderIterator> & it)
{
	return it.data_atEnd;
}
template <typename TContainer>
inline bool
atEnd(Iter<TContainer, FileReaderIterator> const & it)
{
	return it.data_atEnd;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool
atBegin(Iter<TContainer, FileReaderIterator> & it)
{
	return (it.data_abl_pos == 0) && (it.data_buf_pos == 0);
}
template <typename TContainer>
inline bool
atBegin(Iter<TContainer, FileReaderIterator> const & it)
{
	return (it.data_abl_pos == 0) && (it.data_buf_pos == 0);
}

//////////////////////////////////////////////////////////////////////////////
// operator ==
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator == (Iter<TContainer, FileReaderIterator> const & left,
			 Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return (atEnd(left) == atEnd(right)) && ((atEnd(left) && atEnd(right)) || position(left) == position(right));
}

//////////////////////////////////////////////////////////////////////////////
// operator !=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator != (Iter<TContainer, FileReaderIterator> const & left,
			 Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return (atEnd(left) != atEnd(right)) || (position(left) != position(right));
}

//////////////////////////////////////////////////////////////////////////////
// operator < / >
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator < (Iter<TContainer, FileReaderIterator> const & left,
			Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return (atEnd(right) && !atEnd(left)) || (position(left) < position(right));
}

template <typename TContainer>
inline bool 
operator > (Iter<TContainer, FileReaderIterator> const & left,
			Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return (atEnd(left) && !atEnd(right)) || (position(left) > position(right));
}

//////////////////////////////////////////////////////////////////////////////
// operator <= / >=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline bool 
operator <= (Iter<TContainer, FileReaderIterator> const & left,
			 Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return atEnd(right) || (position(left) <= position(right));
}

template <typename TContainer>
inline bool 
operator >= (Iter<TContainer, FileReaderIterator> const & left,
			 Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return atEnd(left) || (position(left) >= position(right));
}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, FileReaderIterator>  
operator + (Iter<TContainer, FileReaderIterator> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, FileReaderIterator>(container(left), position(left) + right);
}
template <typename TContainer, typename TIntegral>
inline Iter<TContainer, FileReaderIterator>  
operator + (TIntegral left,
			Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, FileReaderIterator>(container(right), position(right) + left);
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, FileReaderIterator> &
operator += (Iter<TContainer, FileReaderIterator> & left,
			 TIntegral right)
{
SEQAN_CHECKPOINT
	left.data_buf_pos += right;
	if (left.data_buf_pos >= left.data_buf_len)
	{
		setPosition(left, position(left) + right);
	}
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, FileReaderIterator>  
operator - (Iter<TContainer, FileReaderIterator> const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	return Iter<TContainer, FileReaderIterator>(container(left), position(left) - right);
}

//____________________________________________________________________________

template <typename TContainer>
inline typename Difference<Iter<TContainer, FileReaderIterator> >::Type  
operator - (Iter<TContainer, FileReaderIterator> const & left,
			Iter<TContainer, FileReaderIterator> const & right)
{
SEQAN_CHECKPOINT
	return position(left) - position(right);
}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, FileReaderIterator> &
operator -= (Iter<TContainer, FileReaderIterator> & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	if (left.data_buf_pos < right)
	{
		setPosition(left, position(left) - right);
	}
	else
	{
		left.data_buf_pos -= right;
	}
	return left;
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
