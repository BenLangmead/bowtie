#ifndef PACKEDSTRINGIO_H_
#define PACKEDSTRINGIO_H_

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include "endian.h"

using namespace std;
using namespace seqan;

/**
 * Exception to throw when a Fasta file is malformed.
 */
class MalformedFastaException : public runtime_error {
public:
	MalformedFastaException(const string& msg = "") : runtime_error(msg) {}
};

/**
 * Exception to throw when a type was of an unexpected size for
 * marshaling.
 */
class UnexpectedTypeSizeException : public runtime_error {
public:
	UnexpectedTypeSizeException(const string& _msg = "", int _actually = 0) :
		runtime_error(_msg), msg(_msg), actually(_actually) { }
	
    virtual const char* what() const throw() {
    	stringstream ss(msg); ss << " was: " << actually;
    	return ss.str().c_str();
    }
private:
	const string& msg;
	const int actually;
};

/**
 * Assumes the host of a packed string is a String<unsigned int>
 */
template<typename TVal>
static void
writePacked(ostream& out, const vector<String<TVal, Packed<> > >& ss, bool verbose = false)
throw(UnexpectedTypeSizeException)
{
	for(size_t i = 0; i < ss.size(); i++) {
		writePacked(out, ss[i], verbose);
	}
}

/**
 * Writes given packed string to given output stream.  Format of the
 * packed file is simply a list of packed sequences, where each packed
 * sequence is a 8-byte unsigned (unpacked) sequence length followed by
 * the packed sequence.
 * 
 * Assumes the host of a packed string is a String<unsigned int>
 */
template<typename TVal>
static void
writePacked(ostream& out, const String<TVal, Packed<> >& s, bool verbose = false)
throw(UnexpectedTypeSizeException)
{
	uint64_t tmp;
	
	// Endianness test
	static uint8_t __etest[2] = { 1, 0 };
	bool littleEndian = (*(uint16_t *)__etest == 1);

	// Get length of string and write to stream, being careful to swap
	// if this is a little-endian machine
	uint64_t len = tmp = length(s);
	if(verbose) cout << "len: " << len << endl;
	if(littleEndian) tmp = endianSwapU64(tmp);
	out.write((char*)&tmp, 8);
	
	// Determine size of unsigned int; can only handle 2-, 4-, 8-byte
	// ints for now
	int sui = sizeof(unsigned int);
	if(sui != 2 && sui != 4 && sui != 8) {
		throw UnexpectedTypeSizeException("sizeof(unsigned int) must be 2, 4 or 8", sui);
	}
	int charsPerUi = _PackedConsts<String<TVal, Packed<> > >::VALUES_PER_WORD;
	int uisPerUi64;
	switch(sui) {
		case 2: uisPerUi64 = 4; break;
		case 4: uisPerUi64 = 2; break;
		case 8: uisPerUi64 = 1; break;
		default: throw; // can't happen
	}
	int charsPerUi64 = charsPerUi * uisPerUi64;
	
	if(verbose) cout << "sui: " << sui << endl;
	if(verbose) cout << "charsPerUi: " << charsPerUi << endl;
	if(verbose) cout << "uisPerUi64: " << uisPerUi64 << endl;
	if(verbose) cout << "charsPerUi64: " << charsPerUi64 << endl;
	
	// Write the packed sequence in 64-bit chunks, being careful to
	// swap if this is a little-endian machine.
	const String<unsigned int>& h = host(s);
	// Round up to nearest 8-byte boundary to get the number of 8-byte
	// chunks
	uint64_t blen64 = (len+(charsPerUi64-1))/charsPerUi64;
	if(verbose) cout << "blen64: " << blen64 << endl;
	for(uint64_t i = 0; i < blen64; i++) {
		uint64_t tmp;
		uint64_t off = i * uisPerUi64;
		switch(uisPerUi64) {
			case 1: {
				tmp = h[off];
				break;
			}
			case 2: {
				tmp  = ((off < len-1)? h[off+1] : 0);
				tmp <<= 32;
				tmp |= h[off];
				break;
			}
			case 4: {
				tmp  = ((off < len-3)? h[off+3] : 0);
				tmp <<= 16;
				tmp |= ((off < len-2)? h[off+2] : 0);
				tmp <<= 16;
				tmp |= ((off < len-1)? h[off+1] : 0);
				tmp <<= 16;
				tmp |= h[off];
				break;
			}
			default: throw; // can't happen
		}
		if(littleEndian) tmp = endianSwapU64(tmp);
		out.write((char*)&tmp, 8);
	}
}

/**
 * Reads the sequences of a fasta file into a vector of packed Strings.
 */
template<typename TVal>
static void readAndPackFasta(istream& in, vector<String<TVal, Packed<> > >& ss, bool verbose = false)
throw(MalformedFastaException)
{
	string line;
	if(getline(in, line).eof()) {
		// Error; needs to be at least one line
		throw MalformedFastaException("EOF reading first line");
	}
	if(line[0] != '>') {
		throw MalformedFastaException("First line did not begin with '>'");
	}
	String<TVal, Packed<> > s;
	do {
		if(line[0] == '>') {
			// Push completed string onto vector
			if(!empty(s)) {
				ss.push_back(s);
				clear(s);
				assert(empty(s));
			}
		} else if(line.find('>') != string::npos) {
			throw MalformedFastaException("'>' found in the middle of a line that didn't begin with '>'");
		} else if (!line.empty()) {
			// Append line onto String in progress
			//if(verbose) cout << line;
			append(s, line.c_str()); // Wildcards are normalized to As!
			assert(!empty(s));
		}
	} while(!getline(in, line).eof());
	if(!empty(s)) {
		ss.push_back(s);
	}
}

/**
 * Read packed strings from the given input stream until it is
 * exhausted, appending each to ss.
 */
template<typename TVal>
static void readPacked(istream& in,
                       vector<String<TVal, Packed<> > >& ss,
                       int upto = -1,
                       bool verbose = false)
throw(UnexpectedTypeSizeException)
{
	// Get length
	in.seekg (0, ios::end);
	streampos flen = in.tellg();
	in.seekg (0, ios::beg);
	// Read packed strings until file is exhausted
	while(in.tellg() < flen) {
		String<TVal, Packed<> > ps;
		readPacked(in, ps, verbose);
		ss.push_back(ps);
		if(upto != -1 && ss.size() >= (size_t)upto) {
			return;
		}
	}
}

/**
 * Read one packed string from the given input stream and store it in
 * s.
 */
template<typename TVal>
static void readPacked(istream& in,
                       String<TVal, Packed<> >& s,
                       bool verbose = false)
throw(UnexpectedTypeSizeException)
{
	if(verbose) cout << "Reading a packed string" << endl;
	
	// Endianness test
	static uint8_t __etest[2] = { 1, 0 };
	int littleEndian = (*(uint16_t *) __etest == 1);

	// Read length of string
	uint64_t len;
	in.read((char*)&len, 8);
	assert(in.gcount() == 8);
	if(littleEndian) len = endianSwapU64(len);
	
	// Determine size of unsigned it; we can only handle 4- and 8-byte
	// unsigned ints for now
	int sui = sizeof(unsigned int);
	if(sui != 2 && sui != 4 && sui != 8) {
		throw UnexpectedTypeSizeException("sizeof(unsigned int) must be 2, 4 or 8", sui);
	}
	int charsPerUi = _PackedConsts<String<TVal, Packed<> > >::VALUES_PER_WORD;
	int uisPerUi64;
	switch(sui) {
		case 2: uisPerUi64 = 4; break;
		case 4: uisPerUi64 = 2; break;
		case 8: uisPerUi64 = 1; break;
		default: throw; // can't happen
	}
	int charsPerUi64 = charsPerUi * uisPerUi64;
	if(verbose) cout << "sui: " << sui << endl;
	if(verbose) cout << "charsPerUi: " << charsPerUi << endl;
	if(verbose) cout << "uisPerUi64: " << uisPerUi64 << endl;
	if(verbose) cout << "charsPerUi64: " << charsPerUi64 << endl;

	// Number of unsigned ints in the packed String<>
	uint64_t blen64 = (len+(charsPerUi64-1))/charsPerUi64;
	uint64_t blen = (len+(charsPerUi-1))/charsPerUi;
	String<unsigned int> h;
	reserve(h, blen+3);
	
	// Read each unsigned int, being careful about whether the uint is
	// 4 or 8 bytes, and being careful to swap if this is a little-
	// endian machine
	for(uint64_t i = 0; i < blen64; i++) {
		uint64_t tmp;
		in.read((char*)&tmp, 8);
		assert(in.gcount() == 8);
		if(littleEndian) tmp = endianSwapU64(tmp);
		switch(uisPerUi64) {
			case 1: {
				append(h, tmp);
				break;
			}
			case 2: {
				uint32_t tmp0 = (uint32_t)tmp;
				uint32_t tmp1 = (uint32_t)(tmp >> 32);
				append(h, tmp0);
				append(h, tmp1);
				break;
			}
			case 4: {
				uint16_t tmp0 = (uint16_t)tmp;
				uint16_t tmp1 = (uint16_t)(tmp >> 16);
				uint16_t tmp2 = (uint16_t)(tmp >> 32);
				uint16_t tmp3 = (uint16_t)(tmp >> 48);
				append(h, tmp0);
				append(h, tmp1);
				append(h, tmp2);
				append(h, tmp3);
				break;
			}
			default: throw; // can't happen
		}
	}
	assign(host(s), h);
	_setLength(s, len);
}

#endif /*PACKEDSTRINGIO_H_*/
