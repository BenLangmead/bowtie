#ifndef EBWT_H_
#define EBWT_H_

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "bitpack.h"
#include "blockwise_sa.h"
#include "endian.h"
#include "word_io.h"
#include "ccnt_lut.h"
#include "random_source.h"
#include "hit.h"
#include "ref_read.h"

using namespace std;
using namespace seqan;

#ifndef VMSG_NL
#define VMSG_NL(args...) { \
	stringstream tmp; \
	tmp << args << endl; \
	this->verbose(tmp.str()); \
}
#endif

#ifndef VMSG
#define VMSG(args...) { \
	stringstream tmp; \
	tmp << args; \
	this->verbose(tmp.str()); \
}
#endif

/**
 * Extended Burrows-Wheeler transform header.  This together with the
 * actual data arrays and other text-specific parameters defined in
 * class Ebwt constitute the entire Ebwt.
 */
class EbwtParams {

public:
	EbwtParams(uint32_t __len,
	           int32_t __lineRate,
	           int32_t __linesPerSide,
	           int32_t __offRate,
	           int32_t __ftabChars,
	           int32_t __chunkRate) :
	           _len(__len),        // # characters in the original string
	           _bwtLen(_len + 1),  // # characters in the BWT string ; extra 1 for $ suffix 
	           _sz((_len+3)/4),    // # bytes needed to hold original string
	           _bwtSz(_len/4 + 1), // # bytes needed to hold BWT string
	           _lineRate(__lineRate), 
	           _linesPerSide(__linesPerSide),
	           _origOffRate(__offRate),
	           _offRate(__offRate),
	           _offMask(0xffffffff << _offRate),
	           _ftabChars(__ftabChars),
	           _eftabLen(_ftabChars*2),
	           _eftabSz(_eftabLen*4),
	           _ftabLen((1 << (_ftabChars*2))+1),
	           _ftabSz(_ftabLen*4),
	           _offsLen((_bwtLen + (1 << _offRate) - 1) >> _offRate),
	           _offsSz(_offsLen*4),
	           _lineSz(1 << _lineRate),
	           _sideSz(_lineSz * _linesPerSide),
	           _sideBwtSz(_sideSz - 8),
	           _sideBwtLen(_sideBwtSz*4),
	           _numSidePairs((_bwtSz+(2*_sideBwtSz)-1)/(2*_sideBwtSz)), // string size rounded up to nearest line pair
	           _numSides(_numSidePairs*2),
	           _numLines(_numSides * _linesPerSide),
	           _ebwtTotLen(_numSidePairs * (2*_sideSz)),
	           _ebwtTotSz(_ebwtTotLen),
	           _chunkRate(__chunkRate),
	           _chunkLen(1 << _chunkRate),
	           _chunkMask(0xffffffff << _chunkRate),
	           _numChunks((_len + _chunkLen - 1) / _chunkLen)
	{
		assert(repOk());
	}
	
	EbwtParams(const EbwtParams& eh) :
	           _len(eh._len),
	           _bwtLen(eh._bwtLen), 
	           _sz(eh._sz),
	           _bwtSz(eh._bwtSz),
	           _lineRate(eh._lineRate), 
	           _linesPerSide(eh._linesPerSide),
	           _origOffRate(eh._origOffRate),
	           _offRate(eh._offRate),
	           _offMask(eh._offMask),
	           _ftabChars(eh._ftabChars),
	           _eftabLen(eh._eftabLen),
	           _eftabSz(eh._eftabSz),
	           _ftabLen(eh._ftabLen),
	           _ftabSz(eh._ftabSz),
	           _offsLen(eh._offsLen),
	           _offsSz(eh._offsSz),
	           _lineSz(eh._lineSz),
	           _sideSz(eh._sideSz),
	           _sideBwtSz(eh._sideBwtSz),
	           _sideBwtLen(eh._sideBwtLen),
	           _numSidePairs(eh._numSidePairs),
	           _numSides(eh._numSides),
	           _numLines(eh._numLines),
	           _ebwtTotLen(eh._ebwtTotLen),
	           _ebwtTotSz(eh._ebwtTotSz),
	           _chunkRate(eh._chunkRate),
	           _chunkLen(eh._chunkLen),
	           _chunkMask(eh._chunkMask),
	           _numChunks(eh._numChunks)
	{
		assert(repOk());
	}

	uint32_t len() const           { return _len; }
	uint32_t bwtLen() const        { return _bwtLen; }
	uint32_t sz() const            { return _sz; }
	uint32_t bwtSz() const         { return _bwtSz; }
	int32_t  lineRate() const      { return _lineRate; }
	int32_t  linesPerSide() const  { return _linesPerSide; }
	int32_t  origOffRate() const   { return _origOffRate; }
	int32_t  offRate() const       { return _offRate; }
	uint32_t offMask() const       { return _offMask; }
	int32_t  ftabChars() const     { return _ftabChars; }
	uint32_t eftabLen() const      { return _eftabLen; }
	uint32_t eftabSz() const       { return _eftabSz; }
	uint32_t ftabLen() const       { return _ftabLen; }
	uint32_t ftabSz() const        { return _ftabSz; }
	uint32_t offsLen() const       { return _offsLen; }
	uint32_t offsSz() const        { return _offsSz; }
	uint32_t lineSz() const        { return _lineSz; }
	uint32_t sideSz() const        { return _sideSz; }
	uint32_t sideBwtSz() const     { return _sideBwtSz; }
	uint32_t sideBwtLen() const    { return _sideBwtLen; }
	uint32_t numSidePairs() const  { return _numSidePairs; }
	uint32_t numSides() const      { return _numSides; }
	uint32_t numLines() const      { return _numLines; }
	uint32_t ebwtTotLen() const    { return _ebwtTotLen; }
	uint32_t ebwtTotSz() const     { return _ebwtTotSz; }
	int32_t  chunkRate() const     { return _chunkRate; }
	uint32_t chunkLen() const      { return _chunkLen; }
	uint32_t chunkMask() const     { return _chunkMask; }
	uint32_t numChunks() const     { return _numChunks; }

	void setOffRate(int __offRate) {
		_offRate = __offRate;
        _offMask = 0xffffffff << _offRate;
        _offsLen = (_bwtLen + (1 << _offRate) - 1) >> _offRate;
        _offsSz = _offsLen*4;
	}

	/// Check that this EbwtParams is internally consistent
	bool repOk() const {
		assert_gt(_len, 0);
		assert_gt(_lineRate, 3);
		assert_geq(_offRate, 0);
		assert_leq(_ftabChars, 16);
		assert_geq(_ftabChars, 1);
		assert_geq(_chunkRate, 1);
		assert_lt(_chunkRate, 32);
		assert_lt(_lineRate, 32);
		assert_lt(_linesPerSide, 32);
		assert_lt(_ftabChars, 32);
		assert_eq(0, _ebwtTotSz % (2*_lineSz));
		return true;
	}

	/**
	 * Pretty-print the header contents to the given output stream.
	 */
	void print(ostream& out) const {
		out << "Headers:" << endl
		    << "    len: "          << _len << endl
		    << "    bwtLen: "       << _bwtLen << endl
		    << "    sz: "           << _sz << endl
		    << "    bwtSz: "        << _bwtSz << endl
		    << "    lineRate: "     << _lineRate << endl
		    << "    linesPerSide: " << _linesPerSide << endl
		    << "    offRate: "      << _offRate << endl
		    << "    offMask: 0x"    << hex << _offMask << dec << endl
		    << "    ftabChars: "    << _ftabChars << endl
		    << "    eftabLen: "     << _eftabLen << endl
		    << "    eftabSz: "      << _eftabSz << endl
		    << "    ftabLen: "      << _ftabLen << endl
		    << "    ftabSz: "       << _ftabSz << endl
		    << "    offsLen: "      << _offsLen << endl
		    << "    offsSz: "       << _offsSz << endl
		    << "    lineSz: "       << _lineSz << endl
		    << "    sideSz: "       << _sideSz << endl
		    << "    sideBwtSz: "    << _sideBwtSz << endl
		    << "    sideBwtLen: "   << _sideBwtLen << endl
		    << "    numSidePairs: " << _numSidePairs << endl
		    << "    numSides: "     << _numSides << endl
		    << "    numLines: "     << _numLines << endl
		    << "    ebwtTotLen: "   << _ebwtTotLen << endl
		    << "    ebwtTotSz: "    << _ebwtTotSz << endl
		    << "    chunkRate: "    << _chunkRate << endl
		    << "    chunkLen: "     << _chunkLen << endl
	        << "    chunkMask: "    << hex << _chunkMask << dec << endl
            << "    numChunks: "    << _numChunks << endl;
	}

private:
    uint32_t _len;
    uint32_t _bwtLen;
    uint32_t _sz;
    uint32_t _bwtSz;
    int32_t  _lineRate;
    int32_t  _linesPerSide;
    int32_t  _origOffRate;
    int32_t  _offRate;
    uint32_t _offMask;
    int32_t  _ftabChars;
    uint32_t _eftabLen;
    uint32_t _eftabSz;
    uint32_t _ftabLen;
    uint32_t _ftabSz;
	uint32_t _offsLen;
	uint32_t _offsSz;
	uint32_t _lineSz;
	uint32_t _sideSz;
	uint32_t _sideBwtSz;
	uint32_t _sideBwtLen;
	uint32_t _numSidePairs;
	uint32_t _numSides;
	uint32_t _numLines;
	uint32_t _ebwtTotLen;
	uint32_t _ebwtTotSz;
	int32_t  _chunkRate;
	uint32_t _chunkLen;
	uint32_t _chunkMask;
	uint32_t _numChunks;
};

/**
 * Exception to throw when a file-realted error occurs.
 */
class EbwtFileOpenException : public std::runtime_error {
public:
	EbwtFileOpenException(const std::string& msg = "") :
		std::runtime_error(msg) { }
};

// Forward declarations for Ebwt class
class SideLocus;
template<typename TStr> class EbwtSearchParams;
template<typename TStr> class EbwtSearchState;

/**
 * Extended Burrows-Wheeler transform data.
 * 
 * An Ebwt may be transferred to and from RAM with calls to
 * evictFromMemory() and loadIntoMemory().  By default, a newly-created
 * Ebwt is not loaded into memory; if the user would like to use a
 * newly-created Ebwt to answer queries, they must first call
 * loadIntoMemory().
 */
template <typename TStr>
class Ebwt {
public:
	typedef typename Value<TStr>::Type TAlphabet;

	#define Ebwt_INITS \
	    _toBigEndian(currentlyBigEndian()), \
	    _overrideOffRate(__overrideOffRate), \
	    _verbose(__verbose), \
	    _sanity(__sanityCheck), \
	    _in1(), \
	    _in2(), \
	    _zOff(0xffffffff), \
	    _zEbwtByteOff(0xffffffff), \
	    _zEbwtBpOff(-1), \
	    _nPat(0), \
	    _plen(NULL), \
	    _pmap(NULL), \
	    _fchr(NULL), \
	    _ftab(NULL), \
	    _eftab(NULL), \
	    _offs(NULL), \
	    _ebwt(NULL)
	
	/// Construct an Ebwt from the given input file
	Ebwt(const string& in,
	     int32_t __overrideOffRate = -1,
	     bool __verbose = false,
	     bool __sanityCheck = false) :
	     Ebwt_INITS,
	     _eh(readIntoMemory(true, in + ".1.ebwt", in + ".2.ebwt"))
	{
		// Read offs from secondary stream; if the offRate has been
		// overridden, make sure to sample the offs rather than take them
		// all
		if(_overrideOffRate > _eh.offRate()) {
			_eh.setOffRate(_overrideOffRate);
			assert_eq(_overrideOffRate, _eh.offRate());
		}
		assert(repOk());
	}
	
	/// Construct an Ebwt from the given header parameters and string
	/// vector, optionally using a blockwise suffix sorter with the
	/// given 'bmax' and 'dcv' parameters.  The string vector is
	/// ultimately joined and the joined string is passed to buildToDisk().
	Ebwt(int32_t lineRate,
	     int32_t linesPerSide,
	     int32_t offRate,
	     int32_t ftabChars,
	     int32_t chunkRate,
	     const string& file,   // base filename for EBWT files
	     bool useBlockwise,
	     uint32_t bmax,
	     uint32_t bmaxSqrtMult,
	     uint32_t bmaxDivN,
	     int dcv,
	     vector<istream*>& is,
	     vector<uint32_t>& szs,
	     uint32_t sztot,
	     const RefReadInParams& refparams,
	     uint32_t seed,
	     int32_t __overrideOffRate = -1,
	     bool __verbose = false,
	     bool __sanityCheck = false) :
	     Ebwt_INITS,
	     _eh(joinedLen(szs, chunkRate), lineRate, linesPerSide, offRate, ftabChars, chunkRate)
	{
		string file1 = file + ".1.ebwt";
		string file2 = file + ".2.ebwt";
		string file3 = file + ".3.ebwt";
		// Open output files
		ofstream fout1(file1.c_str(), ios::binary);
		ofstream fout2(file2.c_str(), ios::binary);
		ofstream fout3(file3.c_str(), ios::binary);
//		{
//			Timer timer(cout, "  Time for call to writePacked: ", __verbose);
//			VMSG_NL("Writing packed representation to " << file3);
//			writePacked<Dna>(fout3, ss, __verbose);
//		}
		// Build
		initFromVector(is,
		               szs,
		               sztot,
		               refparams,
		               fout1,
		               fout2,
		               fout3,
		               useBlockwise,
		               bmax,
		               bmaxSqrtMult,
		               bmaxDivN,
		               dcv,
		               seed);
		// Close output files
		fout1.flush();
		VMSG_NL("Wrote " << fout1.tellp() << " bytes to primary EBWT file: " << file1);
		fout1.close();
		fout2.flush();
		VMSG_NL("Wrote " << fout2.tellp() << " bytes to secondary EBWT file: " << file2);
		fout2.close();
		// Reopen as input streams
		VMSG_NL("Re-opening _in1 and _in2 as input streams");
		_in1.open(file1.c_str(), ios_base::in | ios::binary);
		_in2.open(file2.c_str(), ios_base::in | ios::binary);
		assert_eq((streamoff)_in1.tellg(), ios::beg);
		assert_eq((streamoff)_in2.tellg(), ios::beg);
		assert(_in1.good());
		assert(_in2.good());
		if(_sanity) {
			VMSG_NL("Sanity-checking Ebwt");
			assert(!isInMemory());
			readIntoMemory(false);
			sanityCheckAll();
			evictFromMemory();
			assert(!isInMemory());
			assert(_in1.is_open()); assert(_in1.good());
			assert(_in2.is_open()); assert(_in2.good());
			assert_eq((streamoff)_in1.tellg(), ios::beg);
			assert_eq((streamoff)_in2.tellg(), ios::beg);
		}
		VMSG_NL("Returning from Ebwt constructor");
	}
	
	/**
	 * Helper for the two complete constructors above.  Takes a vector
	 * of text strings and renders them into a single string with a
	 * call to join, then appends a 'N' ('$'-equivalent) and constructs
	 * the suffix array over that string, then calls the main Ebwt
	 * building routing (buildToDisk()) passing the string and the
	 * suffix-array builder as input.
	 * 
	 * Note that the suffix array might be extremely large - in fact,
	 * it's exactlty 4 bytes per text character.  If the input text is
	 * the human genome (3 Gbases), the suffix array will be about 12
	 * GB.
	 * 
	 * Both the suffix array and the joined string go out of scope at
	 * the end of this function, causing them both to be freed.
	 */
	void initFromVector(vector<istream*>& is,
	                    vector<uint32_t>& szs,
	                    uint32_t sztot,
	                    const RefReadInParams& refparams,
	                    ofstream& out1,
	                    ofstream& out2,
	                    ofstream& out3,
	                    bool useBlockwise,
	                    uint32_t bmax,
	                    uint32_t bmaxSqrtMult,
	                    uint32_t bmaxDivN,
	                    int dcv,
	                    uint32_t seed) 
	{
		// Compose text strings into single string; doing so
		// initializes _plen, _pmap, _nPat
		VMSG_NL("Calculating joined length");
		TStr s; // holds the entire joined reference after call to joinToDisk
		uint32_t jlen = joinedLen(szs, _eh.chunkRate());
		assert_geq(jlen, sztot);
		VMSG_NL("  = " << jlen << " (" << (jlen-sztot) << " bytes of padding)");
		VMSG_NL("Writing header");
		writeFromMemory(true, out1, out2);
		try {
			VMSG_NL("Reserving space for joined string");
			seqan::reserve(s, jlen, Exact());
			// FIXME: 7/8/08: Why does seqan::capacity(s) return 0
			// here, causing the assert to fire?  When I run in a
			// debugger, seqan::capacity() returns the correct value.
			// Happens even with optimizations turned off.  Very, very
			// strange.
			//assert_gt(cap, seqan::capacity(s));
			VMSG_NL("Joining reference sequences");
			{
				Timer timer(cout, "  Time to join reference sequences: ", _verbose);
				joinToDisk(is, szs, sztot, refparams, s, out1, out2, out3, seed);
			}
			out3.flush();
			VMSG_NL("Wrote " << out3.tellp() << " bytes to tertiary EBWT file");
			out3.close();
		} catch(bad_alloc& e) {
			cerr << "Out of memory creating joined string in "
			     << "Ebwt::initFromVector() at " << __FILE__ << ":"
			     << __LINE__ << endl;
			throw e;
		}
		assert_geq(length(s), jlen);
		if(useBlockwise) {
			if(bmax != 0xffffffff) {
				VMSG_NL("bmax according to bmax setting: " << bmax);
			}
			else if(bmaxSqrtMult != 0xffffffff) {
				bmax *= bmaxSqrtMult;
				VMSG_NL("bmax according to bmaxSqrtMult setting: " << bmax);
			}
			else if(bmaxDivN != 0xffffffff) {
				bmax = max<uint32_t>(jlen / bmaxDivN, 1);
				VMSG_NL("bmax according to bmaxDivN setting: " << bmax);
			}
			else {
				bmax = (uint32_t)sqrt(length(s));
				VMSG_NL("bmax defaulted to: " << bmax);
			}
			VMSG("Using blockwise SA w/ bmax=" << bmax);
			if(dcv == 0) {
				VMSG_NL(" and *no difference cover*");
			} else {
				VMSG_NL(", dcv=" << dcv);
			}
			KarkkainenBlockwiseSA<TStr> bsa(s, bmax, dcv, seed, _sanity, _verbose);
			assert(bsa.suffixItrIsReset());
			assert_eq(bsa.size(), length(s)+1);
			// Build Ebwt; doing so writes everything else (p
			buildToDisk(bsa, s, out1, out2);
		} else {
			VMSG_NL("Using entire SA");
			SillyBlockwiseDnaSA<TStr> bsa(s, 32, _sanity, _verbose);
			assert(bsa.suffixItrIsReset());
			assert_eq(bsa.size(), length(s)+1);
			// Build Ebwt; doing so initializes everything else
			buildToDisk(bsa, s, out1, out2);
		}
		assert(repOk());
		VMSG_NL("Returning from initFromVector");
	}
	
	/// Return the length that the joined string of the given string list will have
	uint32_t joinedLen(vector<uint32_t>& szs, uint32_t chunkRate) {
		uint32_t ret = 0;
		uint32_t chunkLen = 1 << chunkRate;
		for(unsigned int i = 0; i < szs.size(); i++) {
			//if(i < szs.size() - 1) {
			ret += ((szs[i] + chunkLen - 1) / chunkLen) * chunkLen;
			//} else {
			//	ret += szs[i];
			//}
		}
		return ret;
	}

	/// Destruct an Ebwt
	~Ebwt() {
		// Delete everything that was allocated in read(false, ...)
		if(_fchr  != NULL) delete[] _fchr;  _fchr  = NULL;
		if(_ftab  != NULL) delete[] _ftab;  _ftab  = NULL;
		if(_eftab != NULL) delete[] _eftab; _eftab = NULL;
		if(_offs  != NULL) delete[] _offs;  _offs  = NULL;
		if(_plen  != NULL) delete[] _plen;  _plen  = NULL;
		if(_pmap  != NULL) delete[] _pmap;  _pmap  = NULL;
		if(_ebwt  != NULL) delete[] _ebwt;  _ebwt  = NULL;
		try {
			if(_in1.is_open()) _in1.close();
			if(_in2.is_open()) _in2.close();
		} catch(...) {
			VMSG_NL("~Ebwt(): Caught an exception while closing streams!");
		}
	}
	
	/// Accessors
	const EbwtParams& eh() const     { return _eh; }
	uint32_t    zOff() const         { return _zOff; }
	uint32_t    zEbwtByteOff() const { return _zEbwtByteOff; }
	int         zEbwtBpOff() const   { return _zEbwtBpOff; }
	uint32_t    nPat() const         { return _nPat; }
	uint32_t*   fchr() const         { return _fchr; }
	uint32_t*   ftab() const         { return _ftab; }
	uint32_t*   eftab() const        { return _eftab; }
	uint32_t*   offs() const         { return _offs; }
	uint32_t*   plen() const         { return _plen; }
	uint32_t*   pmap() const         { return _pmap; }
	uint8_t*    ebwt() const         { return _ebwt; }
	bool        toBe() const         { return _toBigEndian; }
	bool        verbose() const      { return _verbose; }
	bool        sanityCheck() const  { return _sanity; }
	
	/// Return true iff the Ebwt is currently in memory
	bool isInMemory() const {
		assert(_eh.repOk());
		if(_ebwt != NULL) {
			assert(_ftab != NULL);
			assert(_eftab != NULL);
			assert(_fchr != NULL);
			assert(_offs != NULL);
			assert(_pmap != NULL);
			assert_neq(_zEbwtByteOff, 0xffffffff);
			assert_neq(_zEbwtBpOff, -1);
			return true;
		} else {
			assert(_ftab == NULL);
			assert(_eftab == NULL);
			assert(_fchr == NULL);
			assert(_offs == NULL);
			assert(_pmap == NULL);
			assert_eq(_zEbwtByteOff, 0xffffffff);
			assert_eq(_zEbwtBpOff, -1);
			return false;
		}
	}

	/// Return true iff the Ebwt is currently stored on disk
	bool isEvicted() const {
		return !isInMemory();
	}

	/**
	 * Load this Ebwt into memory by reading it in from the _in1 and
	 * _in2 streams.
	 */
	void loadIntoMemory() {
		readIntoMemory(false);
	}

	/**
	 * Frees memory associated with the Ebwt.
	 */
	void evictFromMemory() {
		assert(isInMemory());
		delete[] _fchr;  _fchr  = NULL;
		delete[] _ftab;  _ftab  = NULL;
		delete[] _eftab; _eftab = NULL;
		delete[] _offs;  _offs  = NULL;
		// Keep plen; it's small and the client may want to query it
		// even when the others are evicted.
		//delete[] _plen;  _plen  = NULL;
		delete[] _pmap;  _pmap  = NULL;
		delete[] _ebwt;  _ebwt  = NULL;
		_zEbwtByteOff = 0xffffffff;
		_zEbwtBpOff = -1;
	}

	/**
	 * Non-static facade for static function ftabHi.
	 */
	uint32_t ftabHi(uint32_t i) const {
		return Ebwt::ftabHi(_ftab, _eftab, _eh.len(), _eh.ftabLen(), _eh.eftabLen(), i);
	}
	
	/**
	 * Get "high interpretation" of ftab entry at index i.  The high
	 * interpretation of a regular ftab entry is just the entry
	 * itself.  The high interpretation of an extended entry is the
	 * second correpsonding ui32 in the eftab.
	 * 
	 * It's a static member because it's convenient to ask this
	 * question before the Ebwt is fully initialized.
	 */
	static uint32_t ftabHi(uint32_t *ftab,
	                       uint32_t *eftab,
	                       uint32_t len,
	                       uint32_t ftabLen,
	                       uint32_t eftabLen,
	                       uint32_t i)
	{
		assert_lt(i, ftabLen);
		if(ftab[i] <= len) {
			return ftab[i];
		} else {
			uint32_t efIdx = ftab[i] ^ 0xffffffff;
			assert_lt(efIdx*2+1, eftabLen);
			return eftab[efIdx*2+1];
		}
	}

	/**
	 * Non-static facade for static function ftabLo.
	 */
	uint32_t ftabLo(uint32_t i) const {
		return Ebwt::ftabLo(_ftab, _eftab, _eh.len(), _eh.ftabLen(), _eh.eftabLen(), i);
	}

	/**
	 * Get "low interpretation" of ftab entry at index i.  The low
	 * interpretation of a regular ftab entry is just the entry
	 * itself.  The low interpretation of an extended entry is the
	 * first correpsonding ui32 in the eftab.
	 * 
	 * It's a static member because it's convenient to ask this
	 * question before the Ebwt is fully initialized.
	 */
	static uint32_t ftabLo(uint32_t *ftab,
	                       uint32_t *eftab,
	                       uint32_t len,
	                       uint32_t ftabLen,
	                       uint32_t eftabLen,
	                       uint32_t i)
	{
		assert_lt(i, ftabLen);
		if(ftab[i] <= len) {
			return ftab[i];
		} else {
			uint32_t efIdx = ftab[i] ^ 0xffffffff;
			assert_lt(efIdx*2+1, eftabLen);
			return eftab[efIdx*2];
		}
	}

	/**
	 * When using read() to create an Ebwt, we have to set a couple of
	 * additional fields in the Ebwt object that aren't part of the
	 * parameter list and are not stored explicitly in the file.  Right
	 * now, this just involves initializing _zEbwtByteOff and
	 * _zEbwtBpOff from _zOff.
	 */
	void postReadInit(EbwtParams& eh) {
		uint32_t sideNum     = _zOff / eh.sideBwtLen();
		uint32_t sideCharOff = _zOff % eh.sideBwtLen();
		uint32_t sideByteOff = sideNum * eh.sideSz();
		_zEbwtByteOff = sideCharOff >> 2;
		assert_lt(_zEbwtByteOff, eh.sideBwtSz());
		_zEbwtBpOff = sideCharOff & 3;
		assert_lt(_zEbwtBpOff, 4);
		if((sideNum & 1) == 0) {
			// This is an even (backward) side
			_zEbwtByteOff = eh.sideBwtSz() - _zEbwtByteOff - 1;
			_zEbwtBpOff = 3 - _zEbwtBpOff;
			assert_lt(_zEbwtBpOff, 4);
		}
		_zEbwtByteOff += sideByteOff;
		assert(repOk(eh)); // Ebwt should be fully initialized now
	}
	
	/**
	 * Pretty-print the Ebwt to the given output stream.
	 */
	void print(ostream& out) const {
		print(out, _eh);
	}

	/**
	 * Pretty-print the Ebwt and given EbwtParams to the given output
	 * stream.
	 */
	void print(ostream& out, const EbwtParams& eh) const {
		eh.print(out); // print params
		out << "Ebwt (" << (isInMemory()? "memory" : "disk") << "):" << endl
		    << "    zOff: "         << _zOff << endl
		    << "    zEbwtByteOff: " << _zEbwtByteOff << endl
		    << "    zEbwtBpOff: "   << _zEbwtBpOff << endl
		    << "    nPat: "  << _nPat << endl
		    << "    plen: ";
		if(_plen == NULL) {
			out << "NULL" << endl;
		} else {
			out << "non-NULL, [0] = " << _plen[0] << endl;
		}
		out << "    pmap: ";
		if(_pmap == NULL) {
			out << "NULL" << endl;
		} else {
			out << "non-NULL, [0] = " << _pmap[0] << endl;
		}
		out << "    ebwt: ";
		if(_ebwt == NULL) {
			out << "NULL" << endl;
		} else {
			out << "non-NULL, [0] = " << _ebwt[0] << endl;
		}
		out << "    fchr: ";
		if(_fchr == NULL) {
			out << "NULL" << endl;
		} else {
			out << "non-NULL, [0] = " << _fchr[0] << endl;
		}
		out << "    ftab: ";
		if(_ftab == NULL) {
			out << "NULL" << endl;
		} else {
			out << "non-NULL, [0] = " << _ftab[0] << endl;
		}
		out << "    eftab: ";
		if(_eftab == NULL) {
			out << "NULL" << endl;
		} else {
			out << "non-NULL, [0] = " << _eftab[0] << endl;
		}
		out << "    offs: ";
		if(_offs == NULL) {
			out << "NULL" << endl;
		} else {
			out << "non-NULL, [0] = " << _offs[0] << endl;
		}
	}
	
	// Building
	static TStr join(vector<TStr>& l, uint32_t chunkRate, uint32_t seed);
	static TStr join(vector<istream*>& l, vector<uint32_t>& szs, uint32_t sztot, const RefReadInParams& refparams, uint32_t chunkRate, uint32_t seed);
	void joinToDisk(vector<istream*>& l, vector<uint32_t>& szs, uint32_t sztot, const RefReadInParams& refparams, TStr& ret, ostream& out1, ostream& out2, ostream& out3, uint32_t seed);
	void buildToDisk(InorderBlockwiseSA<TStr>& sa, const TStr& s, ostream& out1, ostream& out2);

	// I/O
	EbwtParams readIntoMemory(bool justHeader, const string& in1, const string& in2);
	EbwtParams readIntoMemory(bool justHeader, const string& in1, const string& in2, bool& bigEndian);
	EbwtParams readIntoMemory(bool justHeader);
	EbwtParams readIntoMemory(bool justHeader, bool& be);
	void writeFromMemory(bool justHeader, ostream& out1, ostream& out2) const;
	void writeFromMemory(bool justHeader, const string& out1, const string& out2) const;

	// Sanity checking
	void printRangeFw(uint32_t begin, uint32_t end) const;
	void printRangeBw(uint32_t begin, uint32_t end) const;
	void sanityCheckUpToSide(int upToSide) const;
	void sanityCheckAll() const;
	void restore(TStr& s) const;
	
	// Searching and reporting
	inline bool report(uint32_t off, uint32_t top, uint32_t bot, EbwtSearchState<TStr>& s) const;
	inline bool reportChaseOne(uint32_t i, EbwtSearchState<TStr>& s, SideLocus *l /* = NULL */) const;
	inline bool reportChaseRange(EbwtSearchState<TStr>& s) const;
	inline bool reportChaseSample(EbwtSearchState<TStr>& s) const;
	inline bool reportMultiple(EbwtSearchState<TStr>& s) const;
	inline bool reportOne(uint32_t i, EbwtSearchState<TStr>& s, SideLocus *l /* = NULL */) const;
	inline int rowL(const SideLocus& l) const;
	inline uint32_t countUpTo(const SideLocus& l, int c) const;
	inline uint32_t countFwSide(const SideLocus& l, int c) const;
	inline uint32_t countBwSide(const SideLocus& l, int c) const;
	inline uint32_t mapLF(const SideLocus& l) const;
	inline uint32_t mapLF(const SideLocus& l, int c) const;
	inline uint32_t mapLF1(const SideLocus& l, int c) const;
	inline bool searchFinish(EbwtSearchState<TStr>& s) const;
	inline bool searchFinish1(EbwtSearchState<TStr>& s) const;
	inline void searchWithFchr(EbwtSearchState<TStr>& s) const;
	inline void searchWithFtab(EbwtSearchState<TStr>& s) const;
	void search(const TStr& qry, EbwtSearchParams<TStr>& params, uint32_t seed /* = 0 */, bool inexact /* = false */) const;
	void search(EbwtSearchState<TStr>& s, EbwtSearchParams<TStr>& params, bool inexact /* = false */) const;
	bool search1MismatchOrBetter(const TStr& qry, EbwtSearchParams<TStr>& params, bool allowExactHits /* = true*/, uint32_t seed /* = 0*/) const;
	bool search1MismatchOrBetter(EbwtSearchState<TStr>& s, EbwtSearchParams<TStr>& params, bool allowExactHits /* = true*/) const;
	
	/// Check that in-memory Ebwt is internally consistent with respect
	/// to given EbwtParams; assert if not
	bool inMemoryRepOk(const EbwtParams& eh) const {
		assert_leq(ValueSize<TAlphabet>::VALUE, 4);
		assert_geq(_zEbwtBpOff, 0);
		assert_lt(_zEbwtBpOff, 4);
		assert_lt(_zEbwtByteOff, eh.ebwtTotSz());
		assert_lt(_zOff, eh.bwtLen());
		return true;
	}

	/// Check that in-memory Ebwt is internally consistent; assert if
	/// not
	bool inMemoryRepOk() const {
		return repOk(_eh);
	}
	
	/// Check that Ebwt is internally consistent with respect to given
	/// EbwtParams; assert if not
	bool repOk(const EbwtParams& eh) const {
		assert(_eh.repOk());
		if(isInMemory()) {
			return inMemoryRepOk(eh);
		}
		return true;
	}
	
	/// Check that Ebwt is internally consistent; assert if not
	bool repOk() const {
		return repOk(_eh);
	}

private:
	
	ostream& log() const {
		return cout; // TODO: turn this into a parameter
	}
	
	/// Print a verbose message and flush (flushing is helpful for
	/// debugging)
	void verbose(const string& s) const {
		if(this->verbose()) {
			this->log() << s;
			this->log().flush();
		}
	}

	bool       _toBigEndian;
	int32_t    _overrideOffRate;
	bool       _verbose;
	bool       _sanity;
	ifstream   _in1;
	ifstream   _in2;
	uint32_t   _zOff;
	uint32_t   _zEbwtByteOff;
	int        _zEbwtBpOff;
	uint32_t   _nPat;
	uint32_t*  _plen;
	// _plen ans e
	uint32_t*  _pmap;
	// _fchr, _ftab and _eftab are expected to be relatively small
	// (usually < 1MB, perhaps a few MB if _fchr is particularly large
	// - like, say, 11).  For this reason, we don't bother with writing
	// them to disk through separate output streams; we
	uint32_t*  _fchr;
	uint32_t*  _ftab;
	uint32_t*  _eftab; // "extended" entries for _ftab
	// _offs may be extremely large.  E.g. for DNA w/ offRate=4 (one
	// offset every 16 rows), the total size of _offs is the same as
	// the total size of the input sequence
	uint32_t*  _offs;
	// _ebwt is the Extended Burrows-Wheeler Transform itself, and thus
	// is at least as large as the input sequence.
	uint8_t*   _ebwt;
	EbwtParams _eh;
};

#define SAMPLE_THRESH 20

typedef enum {
	MHP_CHASE_ALL = 0,
	MHP_MULTICHASE_ALL,
	MHP_PICK_1_RANDOM,
	MHP_CHASE_SAMPLE
} MultiHitPolicy;

/**
 * Keeps some statistics about how many patterns were tried and placed
 * (or not).
 */
template<typename TStr>
struct EbwtSearchStats {
	uint32_t read;                 // patterns read
	uint64_t readLenTot;           // patterns read
	uint32_t fwRead;               // forward patterns read
	uint32_t rcRead;               // reverse-complement patterns read
	uint32_t tries;                // total "tries" of any kind
	uint32_t exactTries;           // total exact-match "tries"
	uint32_t fwTries;              // forward "tries" of any kind
	uint32_t fwExactTries;         // forward-pattern exact-match "tries"
	uint32_t rcTries;              // reverse-complement "tries" of any kind
	uint32_t rcExactTries;         // reverse-complement exact-match "tries"
	uint64_t exactHitTotalCnt;     // # total hits from all tries
	uint32_t exactHitOneOrMore;    // # tried patterns w/ 1 or more hits
	uint32_t exactHitOneOrMoreFw;  // # tried fw pats w/ 1 or more hits
	uint32_t exactHitOneOrMoreRc;  // # tried rw pats w/ 1 or more hits
	// # pushthroughs (a double-pushthrough counts as 1) while exact matching
	uint64_t pushthroughs;
	// # double pushthroughs (of top and bottom) while exact matching
	uint64_t pushthroughsDouble;
	// # exact double pushthroughs where both original arrows pointed into the
	// same line
	uint64_t pushthroughsDoubleSingleLine;
	// # exact double pushthroughs where both original arrows pointed into
	// distinct lines
	uint64_t pushthroughsDoubleDoubleLine;
	// # exact single pushthroughs (where just one candidate remains, and
	// therefore bot == top+1)
	uint64_t pushthroughsSingle;
	uint64_t pushthroughsSingleChase;
	uint64_t pushthroughsSingleMatch;
	uint64_t stops;
	uint64_t narrowingBackHalfAdvance[64];
	uint64_t backHalfAdvance[64];
	uint64_t stopsAtDepth[64];
	uint64_t stopsAtDepthFw[64];
	uint64_t stopsAtDepthRc[64];
	EbwtSearchStats() :
		read(0), readLenTot(0), fwRead(0), rcRead(0), tries(0), exactTries(0),
		fwTries(0), fwExactTries(0), rcTries(0), rcExactTries(0),
		exactHitTotalCnt(0), exactHitOneOrMore(0), exactHitOneOrMoreFw(0),
		exactHitOneOrMoreRc(0), pushthroughs(0), pushthroughsDouble(0),
		pushthroughsDoubleSingleLine(0), pushthroughsDoubleDoubleLine(0),
		pushthroughsSingle(0), pushthroughsSingleChase(0),
		pushthroughsSingleMatch(0), stops(0)
	{
		bzero(narrowingBackHalfAdvance, 64 * sizeof(uint64_t));
		bzero(backHalfAdvance,          64 * sizeof(uint64_t));
		bzero(stopsAtDepth,             64 * sizeof(uint64_t));
		bzero(stopsAtDepthFw,           64 * sizeof(uint64_t));
		bzero(stopsAtDepthRc,           64 * sizeof(uint64_t));
	}
	void incRead(const EbwtSearchState<TStr> state, const TStr& pat) {
		read++;
		readLenTot += length(pat);
		if(state.params().fw()) fwRead++;
		else rcRead++;
	}
	void incTry(const EbwtSearchState<TStr> state) {
		tries++;
		exactTries++;
		if(state.params().fw()) { fwTries++; fwExactTries++; }
		else { rcTries++; rcExactTries++; }
	}
	void addExactHits(const EbwtSearchState<TStr> state) {
		exactHitTotalCnt += state.spread();
		exactHitOneOrMore++;
		if(state.params().fw()) exactHitOneOrMoreFw++;
		else exactHitOneOrMoreRc++;
	}
	void incPushthrough(const EbwtSearchState<TStr> state,
	                    bool chase = false,
	                    bool singleOverride = false)
	{
		pushthroughs++;
		if(state.spread() > 1 && !singleOverride) {
			pushthroughsDouble++;
		} else {
			pushthroughsSingle++;
			if(chase) {
				pushthroughsSingleChase++;
			} else {
				pushthroughsSingleMatch++;
			}
		}
	}
	void incStopAt(const EbwtSearchState<TStr> state) {
		int i = (int)(state.qlen() - state.qidx() - 1);
		assert_lt(i, 64);
		assert_geq(i, 0);
		stops++;
		stopsAtDepth[i]++;
		if(state.params().fw()) stopsAtDepthFw[i]++;
		else stopsAtDepthRc[i]++;
	}
	void incNarrowHalfAdvance(const EbwtSearchState<TStr> state, bool narrowing = false) {
		int i = (int)(state.qlen() - state.qidx() - 1);
		assert_lt(i, 64);
		assert_geq(i, 0);
		backHalfAdvance[i]++;
		if(narrowing) narrowingBackHalfAdvance[i]++;
	}
	/// Write statistics to given output stream.
	void write(ostream& out) const {
		#define FLOAT_OUT(f) setprecision(1) << fixed << (f)
		out << "Statistics:" << endl;
		out << "  Patterns read: " << read << endl;
		out << "    Forward patterns read: " << fwRead << endl;
		out << "    Reverse-complement patterns read: " << rcRead << endl;
		out << "  Average pattern length: " << FLOAT_OUT((double)readLenTot/read) << endl;
		out << "  Tries: " << tries << endl;
		out << "    Exact tries: " << exactTries << endl;
		out << "  Tries: " << tries << endl;
		out << "    Forward tries: " << fwTries
			<< " (" << FLOAT_OUT(fwTries*100.0/tries) << "% of tries)" << endl;
		out << "      Forward exact tries: " << fwExactTries << endl;
		out << "    Reverse-complement tries: " << rcTries
			<< " (" << FLOAT_OUT(rcTries*100.0/tries) << "% of tries)" << endl;
		out << "      Reverse-complement exact tries: " << rcExactTries << endl;
		// Stops
		out << "  Stops: " << stops << endl;
		out << "    At depth: " << endl;
		out << "      Both          Just fw       Just rc" << endl;
		for(int i = 0; i < 64; i++) {
			if(stopsAtDepth[i] > 0) {
				stringstream ss;
				out << "      ";
				ss << i << ":" << stopsAtDepth[i];
				while(ss.str().length() < 14) {
					ss << " ";
				}
				ss << i << ":" << stopsAtDepthFw[i];
				while(ss.str().length() < 28) {
					ss << " ";
				}
				ss << i << ":" << stopsAtDepthRc[i];
				while(ss.str().length() < 42) {
					ss << " ";
				}
				out << ss.str() << endl;
			}
		}
		// Back-half advances
		out << "  Back-half advances:" << endl;
		out << "    At depth: " << endl;
		out << "      All           Just narrowing" << endl;
		for(int i = 0; i < 64; i++) {
			if(backHalfAdvance[i] > 0) {
				stringstream ss;
				out << "      ";
				ss << i << ":" << backHalfAdvance[i];
				while(ss.str().length() < 14) {
					ss << " ";
				}
				ss << i << ":" << narrowingBackHalfAdvance[i] << " (";
				ss << FLOAT_OUT(narrowingBackHalfAdvance[i]*100.0/backHalfAdvance[i]) << "%)";
				out << ss.str() << endl;
			}
		}
		out << "  Total # of exact hits: " << exactHitTotalCnt
		    << " (" << FLOAT_OUT((float)exactHitTotalCnt/exactTries) << " per exact try)" << endl;
		out << "  Exact tries w/ >=1 hits: " << exactHitOneOrMore
		    << " (" << FLOAT_OUT(exactHitOneOrMore*100.0/exactTries) << "% of exact tries)"
		    << " (" << FLOAT_OUT(exactHitOneOrMore*100.0/read) << "% of patterns)" << endl;
		out << "    Exact fw tries w/ >=1 hits: " << exactHitOneOrMoreFw
		    << " (" << FLOAT_OUT(exactHitOneOrMoreFw*100.0/fwExactTries) << "% of exact fw tries) " << endl;
		out << "    Exact rc tries w/ >=1 hits: " << exactHitOneOrMoreRc
		    << " (" << FLOAT_OUT(exactHitOneOrMoreRc*100.0/rcExactTries) << "% of exact rc tries) " << endl;
		out << "  Exact pushthroughs: " << pushthroughs << endl;
		out << "    Exact single-pushthroughs: " << pushthroughsSingle
		    << " (" << FLOAT_OUT(pushthroughsSingle*100.0/pushthroughs) << "%)" << endl;
		out << "      Exact single-pushthroughs for matching: " << pushthroughsSingleMatch
		    << " (" << FLOAT_OUT(pushthroughsSingleMatch*100.0/pushthroughsSingle) << "% of single pushthroughs)" << endl;
		out << "      Exact single-pushthroughs for chasing results: " << pushthroughsSingleChase
		    << " (" << FLOAT_OUT(pushthroughsSingleChase*100.0/pushthroughsSingle) << "% of single pushthroughs)" << endl;
		out << "    Exact double-pushthroughs: " << pushthroughsDouble
		    << " (" << FLOAT_OUT(pushthroughsDouble*100.0/pushthroughs) << "%)" << endl;
		out << "      Exact single-line double-pushthroughs: " << pushthroughsDoubleSingleLine << endl;
		out << "      Exact double-line double-pushthroughs: " << pushthroughsDoubleDoubleLine << endl;
	}
};

/**
 * Structure encapsulating search parameters, such as whether and how
 * to backtrack and how to deal with multiple equally-good hits.
 */
template<typename TStr>
class EbwtSearchParams {
public:
	EbwtSearchParams(HitSink& __sink,
	                 EbwtSearchStats<TStr>& __stats,
	                 MultiHitPolicy __mhp,
	                 const vector<TStr>& __texts,
	                 bool __revcomp = true,
	                 bool __fw = true,
	                 bool __ebwtFw = true,
	                 bool __arrowMode = false,
	                 bool __comprehensiveBacktrack = false) :
		_sink(__sink),
		_stats(__stats),
		_mhp(__mhp),
		_texts(__texts),
		_patid(0xffffffff),
		_revcomp(__revcomp),
		_fw(__fw),
		_ebwtFw(__ebwtFw),
		_backtracking(false),
        _arrowMode(__arrowMode),
		_comprehensiveBacktrack(__comprehensiveBacktrack),
		_suppress(false) { }
	MultiHitPolicy multiHitPolicy() const { return _mhp; }
	HitSink& sink() const            { return _sink; }
	void setPatId(uint32_t __patid)  { _patid = __patid; }
	uint32_t patId() const           { return _patid; }
	void setFw(bool __fw)            { _fw = __fw; }
	bool fw() const                  { return _fw; }
	void setEbwtFw(bool __ebwtFw)    { _ebwtFw = __ebwtFw; }
	bool ebwtFw() const              { return _ebwtFw; }
	EbwtSearchStats<TStr>& stats() const { return _stats; }
	/**
	 * Report a hit 
	 */
	void reportHit(EbwtSearchState<TStr>& s,
	               U32Pair h,          // hit locus
	               U32Pair a,          // arrow pair
	               uint32_t tlen,      // length of text
	               uint32_t len,       // length of query
	               uint32_t oms) const 
	{
		// The search functions should not have allowed us to get here
		assert(!_suppress);
		//uint32_t mm = _backtracking ? 1 : 0;
		bitset<max_read_bp> mm = 0;
		hit_pat_t pat;
		string patQuals;
		string patName = s.query_name();
		if (_ebwtFw)
		{
			pat = s.query();
			patQuals = s.query_quals();
		}
		else
		{
			for(size_t i = 0; i < len; i++) 
			{
				appendValue(pat, s.query()[len-i-1]);
				patQuals.push_back(s.query_quals()[len-i-1]);
			}
		}
		
		if (s.mismatch() != 0xffffffff)
		{	
			
			if ((_ebwtFw && !_fw) || (_fw && !_ebwtFw))
				mm.set(len - s.mismatch() - 1);
			else
				mm.set(s.mismatch());
			
			//else
//			{
//				if ((_ebwtFw && !_fw) || (!_ebwtFw && _fw))
//					mm.set(len - s.mismatch() - 1);
//				else
//					mm.set(s.mismatch());
//			}
		}
		
		bool provisional = (_backtracking && _mhp == MHP_PICK_1_RANDOM && _fw && _revcomp);
		if(!_ebwtFw && !_arrowMode) {
			h.second = tlen - h.second - 1;
			h.second -= (s.qlen()-1);
		}
		// Check the hit against the original text, if it's available
		if(_texts.size() > 0 && !_arrowMode) {
			assert_lt(h.first, _texts.size());
			assert_eq(tlen, length(_texts[h.first]));
			bitset<max_read_bp> diffs = 0;
			// This type of check assumes that only mismatches are
			// possible.  If indels are possible, then we either need
			// the caller to provide information about indel locations,
			// or we need to extend this to a more complicated check.
			for(size_t i = 0; i < len; i++) {
				assert_lt(h.second + i, length(_texts[h.first]));
				if(_ebwtFw) {
					// Forward pattern appears at h
					if(s.query()[i] != _texts[h.first][h.second + i]) {
						diffs.set(i);
					}
				} else {
					// Reverse of pattern appears at h
					if(s.query()[len-i-1] != _texts[h.first][h.second + i]) {
						// Text and query are already reversed, so we
						// don't need to reverse our entries in diffs
						diffs.set(i);
					}
				}
			}
			if(diffs != mm) {
				cerr << "Expected " << mm << " mismatches, got " << diffs << endl;
				cerr << "  Pat:  ";
				for(size_t i = 0; i < len; i++) {
					if(_ebwtFw) cerr << s.query()[i];
					else cerr << s.query()[len-i-1];
				}
				cerr << endl;
				cerr << "  Tseg: ";
				for(size_t i = 0; i < len; i++) {
					cerr << _texts[h.first][h.second + i];
				}
				cerr << endl;
				if(length(_texts[h.first]) < 80) {
					cerr << "  Text: " << _texts[h.first] << endl;
				}
				cerr << "  FW: " << _fw << endl;
				cerr << "  Ebwt FW: " << _ebwtFw << endl;
				cerr << "  Provisional: " << provisional << endl;
			}
			assert_eq(diffs, mm);
		}
		if(provisional) {
			// Provisional hits may or may not be 'accepted' later on;
			// this might happen if we find a hit with one mismatch 
			// but haven't yet tried to exact-match the pattern's
			// reverse complement.  If the reverse complement does
			// eventually match, then we'll reject this provisional
			// 1-mismatch hit.  Otherwise we'll accept it.
			sink().reportProvisionalHit(_arrowMode? a : h, _patid, patName, pat, patQuals, _fw, mm, oms);
		} else {
			sink().reportHit(_arrowMode? a : h, _patid, patName, pat, patQuals, _fw, mm, oms);
		}
	}
	void write(ostream& out) const {
		const char *mhpToStr[] = {
			"Chase all",
			"Multichase all",
			"Random pick 1"
		};
		out << "Multi-hit policy: " << mhpToStr[_mhp] << endl;
	}
	bool arrowMode() const {
		return _arrowMode;
	}
	bool comprehensiveBacktrack() const {
		return _comprehensiveBacktrack;
	}
	void setComprehensiveBacktrack(bool __comprehensiveBacktrack) {
		_comprehensiveBacktrack = __comprehensiveBacktrack;
	}
	void setSuppressHits(bool __suppress) {
		_suppress = __suppress;
	}
	bool suppressHits() const { return _suppress; }
	/// Return true iff we're in backtracking mode
	bool backtracking() const { return _backtracking; }
	/// Set whether we're in backtracking mode; when in backtracking
	/// mode, we don't change the backtracking state at all
	void setBacktracking(bool __backtracking) {
		_backtracking = __backtracking;
	}
	const vector<TStr>& texts() const { return _texts; }
private:
	HitSink& _sink;
	EbwtSearchStats<TStr>& _stats;
	MultiHitPolicy _mhp;    // policy for when read hits multiple spots
    const vector<TStr>& _texts; // original texts, if available (if not
                                // available, _texts.size() == 0)
	uint32_t _patid;      // id of current read
	bool _revcomp;        // whether reverse complements are enabled
	bool _fw;             // current read is forward-oriented
	bool _ebwtFw;         // current Ebwt is forward-oriented
	bool _backtracking;   // we're currently backtracking
	bool _arrowMode;      // report arrows
	bool _comprehensiveBacktrack; // Allow backtracking to points after
	                           // the read's halfway mark (for
	                           // 1-mismatch)
	bool _suppress;         // suppress hits?
};

struct SideLocus {
	SideLocus() :
	_sideByteOff(0),
	_charOff(0),
	_fw(-1),
	_by(-1),
	_bp(-1),
	_side(NULL) { }
	
	/**
	 * Construct from row and other relevant information about the Ebwt.
	 */
	SideLocus(uint32_t row, const EbwtParams& ep, uint8_t* ebwt) {
		initFromRow(row, ep, ebwt);
	}
	
	/**
	 * Calculate SideLocus based on a row and other relevant
	 * information about the Ebwt.
	 */
	void initFromRow(uint32_t row, const EbwtParams& ep, uint8_t* ebwt) {
		const uint32_t sideBwtLen = ep.sideBwtLen();
		const uint32_t sideBwtSz  = ep.sideBwtSz();
		const uint32_t sideSz     = ep.sideSz();
		const uint32_t sideNum    = row / sideBwtLen;
		_charOff                  = row % sideBwtLen;
		_sideByteOff              = sideNum * sideSz;
		assert_leq(_sideByteOff + sideSz, ep.ebwtTotSz());
		_side = ebwt + _sideByteOff;
		__builtin_prefetch((const void *)_side,
		                   0 /* prepare for read */,
		                   0 /* no locality */);
		_fw = sideNum & 1;   // odd-numbered sides are forward
		_by = _charOff >> 2; // byte within side
		assert_lt(_by, (int)sideBwtSz);
		if(!_fw) _by = sideBwtSz - _by - 1;
		_bp = _charOff & 3;  // bit-pair within byte
		if(!_fw) _bp ^= 3;
	}
	
    uint32_t _sideByteOff; // offset of top side within ebwt[]
    uint32_t _charOff;     // character offset within side
    int _fw;               // side is forward or backward?
    int _by;               // byte within side (not adjusted for bw sides)
    int _bp;               // bitpair within byte (not adjusted for bw sides)
    uint8_t *_side;        // ptr to beginning of top side
};

/**
 * Data structure that maintains all state associated with a simple
 * single-query search (without any sort of backtracking).
 * 
 * This state gets passed around and mutated by the various
 * Ebwt<>::search*() functions.
 */
template<typename TStr>
class EbwtSearchState {
public:
	typedef typename Value<TStr>::Type TVal;
	EbwtSearchState(const Ebwt<TStr>& __ebwt,
	                const TStr& __query,
					const string& __query_name,
					const string& __query_quals,
	                const EbwtSearchParams<TStr>& __params,
	                uint32_t seed = 0) :
	    _ebwt(__ebwt),
		_params(__params),
		_rand(seed),
		_top(0xffffffff),
		_bot(0xffffffff),
		_query(__query),
		_query_name(__query_name),
		_query_quals(__query_quals),
		_remainders(),
		_tried(),
		_tops(),
		_bots(),
		_firstNzRemainder(-1),
		_qlen(length(__query)),
		_qidx(length(_query)-1),
		_topSideLocus(),
		_botSideLocus(),
		_mism(0xffffffff)
	{
		_narrowHalfLen = _qlen;
		// Let the forward Ebwt take the middle character
		if(_params.ebwtFw()) _narrowHalfLen++;
		_narrowHalfLen >>= 1;
		fill(_remainders, _narrowHalfLen, 0);
		fill(_tried,      _narrowHalfLen, 0);
		fill(_tops,       _narrowHalfLen, 0);
		fill(_bots,       _narrowHalfLen, 0);
	}
	const EbwtSearchParams<TStr>& params() const { return _params; }
	RandomSource& rand()                         { return _rand;   }
	const TStr& query() const                    { return _query;  }
	const string& query_quals() const			 { return _query_quals; }
	const string& query_name() const			 { return _query_name; }
	uint32_t top() const                         { return _top;    }
	uint32_t bot() const                         { return _bot;    }
	void setTopBot(uint32_t __top, uint32_t __bot) {
		assert_geq(__bot, __top);
		_top = __top; _bot = __bot;
		if(_bot > _top) {
			uint32_t diff = _bot - _top;
			// Calculate top from scratch and (possibly) prefetch
			_topSideLocus.initFromRow(_top, _ebwt.eh(), _ebwt.ebwt());
			// Is bot within the same side?
			if(_topSideLocus._charOff + diff < _ebwt.eh().sideBwtLen()) {
				// Yes; copy most of top's info; don't prefetch a 2nd time
				SideLocus& tsl = _topSideLocus;
				SideLocus& bsl = _botSideLocus;
				bsl._charOff     = tsl._charOff + diff;
				bsl._sideByteOff = tsl._sideByteOff;
				assert_leq(bsl._sideByteOff + _ebwt.eh().sideSz(), _ebwt.eh().ebwtTotSz());
				bsl._side        = tsl._side;
				bsl._fw          = tsl._fw;
				bsl._by          = bsl._charOff >> 2; // byte within side
				if(!bsl._fw) bsl._by = _ebwt.eh().sideBwtSz() - bsl._by - 1;
				bsl._bp          = bsl._charOff & 3;  // bit-pair within byte
				if(!bsl._fw) bsl._bp ^= 3;
			} else {
				// No; calculate bot from scratch and (possibly) prefetch
				_botSideLocus.initFromRow(_bot, _ebwt.eh(), _ebwt.ebwt());
			}
		}
	}
	uint32_t spread() const     { return _bot - _top; }
	bool repOk() const {
		assert_geq(_bot, _top);
		return true;
	}
	uint32_t qlen() const       { return _qlen; }
	uint32_t qidx() const       { return _qidx; }
	uint32_t mismatch() const	{ return _mism; }
	/// Return true iff we're currently matching in the "back half" of the read
	bool inNarrowHalf() const {
		return _qidx < _narrowHalfLen;
	}
	bool closed() const         { return _bot == _top; }
	bool qAtBeginning() const   { return _qidx == _qlen-1; }
	bool qExhausted() const     { return _qidx == 0xffffffff; }
	void decQidx()              { assert(_qidx != 0xffffffff); _qidx--; }
	void subQidx(uint32_t s)    { assert_leq(s-1, _qidx); _qidx -= s; }
	/// Get char at the current query position
	int chr() const {
		return chr(_qidx);
	}
	/// Get char at a query position
	int chr(uint32_t qidx) const {
		assert(!qExhausted());
		int c = (TVal)_query[qidx];
		assert_lt(c, 4);  // for sanity
		assert_geq(c, 0); // for sanity
		return c;
	}
	/// Assert that state is initialized
	bool initialized() const {
		assert_eq(0xffffffff, _top);
		assert_eq(0xffffffff, _bot);
		assert_eq(_qlen-1, _qidx);
		assert_eq(_mism, 0xffffffff);
		return true;
	}
	/**
	 * 
	 */
	uint32_t mapLF1(const Ebwt<TStr>& ebwt) {
		assert_eq(1, spread());
		int c = chr();
		uint32_t r = ebwt.mapLF1(_topSideLocus, c);
		// no need to update _tried
		if(r == 0xffffffff) {
			_params.stats().incStopAt(*this);
			if(!_params.backtracking() && inNarrowHalf()) {
				assert(_firstNzRemainder == -1 || _firstNzRemainder > (int)_qidx);
				assert_eq(0, _remainders[_qidx]);
				assert_eq(0, _tried[_qidx]);
				assert_eq(0, _tops[_qidx]);
				assert_eq(0, _bots[_qidx]);
				_firstNzRemainder = _qidx;
				_remainders[_qidx] = 1;
				_tried[_qidx] = c;
				_tops[_qidx] = _top;
				_bots[_qidx] = _bot;
			}
		} else {
			_params.stats().incPushthrough(*this);
			if(inNarrowHalf()) _params.stats().incNarrowHalfAdvance(*this);
			setTopBot(r, r+1);
		}
		assert(_qidx != 0xffffffff); _qidx--;
		return r;
	}
	/**
	 * Advance top and bottom arrows according to the next character in
	 * the pattern.  Return true iff there are match candidates
	 * remaining (i.e., iff bot > top).
	 */
	bool mapLF(const Ebwt<TStr>& ebwt, int c = -1) {
		if(c == -1) c = chr();
		assert_lt(c, 4);
		assert_geq(c, 0);
		uint32_t oldTop = _top;
		uint32_t oldBot = _bot;
		uint32_t diff = _bot - _top;
		uint32_t top = ebwt.mapLF(_topSideLocus, c);
		uint32_t bot = ebwt.mapLF(_botSideLocus, c);
		assert_geq(bot, top);
		if(bot == top) {
			if(!_params.backtracking() && inNarrowHalf()) {
				// Arrows moved closer together; might want to
				// backtrack to here
				assert(_firstNzRemainder == -1 || _firstNzRemainder > (int)_qidx);
				assert_eq(0, _remainders[_qidx]);
				assert_eq(0, _tried[_qidx]);
				assert_eq(0, _tops[_qidx]);
				assert_eq(0, _bots[_qidx]);
				_firstNzRemainder = _qidx;
				_remainders[_qidx] = diff;
				_tried[_qidx] = c;
				_tops[_qidx] = oldTop;
				_bots[_qidx] = oldBot;
			}
			// Set _top and _bot but don't calculate their side loci
			// (that's a waste of time - no one cares).  Just leave the
			// loci as-is.
			_bot = bot; _top = top;
			assert(_qidx != 0xffffffff); _qidx--;
			return false;
		}
		setTopBot(top, bot);
		// Calculate old difference minus new difference
		uint32_t diffDiff = diff - (_bot - _top);
		assert_leq(diffDiff, diff);
		if(!_params.backtracking() && inNarrowHalf() && diffDiff > 0) {
			// Arrows moved closer together; might want to backtrack
			// to here
			assert(_firstNzRemainder == -1 || _firstNzRemainder > (int)_qidx);
			assert_eq(0, _remainders[_qidx]);
			assert_eq(0, _tried[_qidx]);
			assert_eq(0, _tops[_qidx]);
			assert_eq(0, _bots[_qidx]);
			_firstNzRemainder = _qidx;
			_remainders[_qidx] = diffDiff;
			_tried[_qidx] = c;
			_tops[_qidx] = oldTop;
			_bots[_qidx] = oldBot;
		} else if(!_params.backtracking() && inNarrowHalf()) {
			assert_eq(0, _remainders[_qidx]);
		}
		if(inNarrowHalf()) {
			_params.stats().incNarrowHalfAdvance(*this, diffDiff > 0);
		}
		assert(_qidx != 0xffffffff); _qidx--;
		return true;
	}
	/**
	 * Try all as-yet-untaken branches at query position qidx; return
	 * true iff one or more 1-mismatch hits were found.
	 */
	bool backtrack1From(const Ebwt<TStr>& ebwt, uint32_t qidx) {
		typedef typename Value<TStr>::Type TVal;
		assert(_params.backtracking());
		assert_geq(qidx, (uint32_t)_firstNzRemainder);
		assert_lt(qidx, _narrowHalfLen);
		assert_lt(_tried[qidx], 4);
		assert_gt(_bots[qidx], _tops[qidx]);
		bool oneHit = (_params.multiHitPolicy() == MHP_PICK_1_RANDOM);
		// Restore top and bot from this branch point; also set the top
		// and bot loci
		setTopBot(_tops[qidx], _bots[qidx]);
		bool gotHits = false;
		assert_lt(qidx, _qlen);
		_qidx = qidx; // Restore qidx
		// Randomly select a character
		unsigned int st = _rand.nextU32();
		ASSERT_ONLY(int skipped = 0);
		// For each character starting with the randomly-selected one
		for(uint32_t i = 0; i < 4; i++) {
			unsigned int imod = (i + st) & 3;
			// Have I tried this one yet?
			if(_tried[qidx] == imod) {
				// Yes...
				ASSERT_ONLY(skipped++);
				continue;
			}
			// Haven't tried this character from this position yet;
			// proceed one step with this character
			ASSERT_ONLY(int savedNz = _firstNzRemainder);
			// Save loci, since the calls to mapLF and searchFinish
			// may mutate them
		    SideLocus savedTopSideLocus = _topSideLocus;
		    SideLocus savedBotSideLocus = _botSideLocus;
			if(!mapLF(ebwt, imod)) {
				// Top/bot loci should not have changed, so there is no
				// need to restore them
				_top = _tops[qidx];
				_bot = _bots[qidx];
				_qidx = qidx;
				assert_eq(savedNz, _firstNzRemainder);
				continue;
			}
			assert_eq(savedNz, _firstNzRemainder);
			assert_eq(0, params().sink().numProvisionalHits());
			// Now finish off the search
			bool ret = ebwt.searchFinish(*this);
			// Did the call to searchFinish generate one or more hits?
			if(ret || params().sink().numProvisionalHits() > 0) {
				// Yes...
				if(oneHit) {
					// Got one 1-mismatch hit; stop now
					return true;
				}
				gotHits = true;
			}
			// Return state to the backtrack point
			_top = _tops[qidx];
			_bot = _bots[qidx];
			_qidx = qidx;
			_topSideLocus = savedTopSideLocus;
			_botSideLocus = savedBotSideLocus;
		}
		if(oneHit) {
			assert_eq(0, params().sink().numProvisionalHits());
		}
		assert_eq(1, skipped);
		return gotHits;
	}
	/**
	 * Try all as-yet-untaken branches at all "narrow-half" query
	 * positions after which the arrows narrowed.
	 */
	void backtrack1(const Ebwt<TStr>& ebwt) {
		bool oneHit = (_params.multiHitPolicy() == MHP_PICK_1_RANDOM);
		assert(_params.backtracking());
		if(_firstNzRemainder == -1) {
			// No backtracking opportunities
			return;
		}
		// It must be the case that the arrows narrowed between this
		// query position and the one "after" (to the left of) it
		assert_gt(_remainders[_firstNzRemainder], 0);
		assert_gt(_bots[_firstNzRemainder], _tops[_firstNzRemainder]);
		#ifndef NDEBUG
		// It can't be the case that arrows narrowed after any query
		// position "after" (to the left of) _firstNzRemainder
		for(int i = 0; i < _firstNzRemainder; i++) {
			assert_eq(0, _remainders[i]);
		}
		// If we have the texts...
		if(_firstNzRemainder > 0 && _params.texts().size() > 0) {
			// ...then we can double-check that each position "before" 
			// (to the right of) and including _firstNzRemainder match
			// their text counterparts
			for(size_t i = _firstNzRemainder; i < _qlen; i++) {
				// ...
			}
		}
		#endif
		for(uint32_t i = _firstNzRemainder; i < _narrowHalfLen; i++) {
			if(_remainders[i] == 0) {
				#ifndef NDEBUG
				#endif
				continue;
			}
			
			_mism = i;
			
			if(backtrack1From(ebwt, i) && oneHit) {
				// We got one or more 1-mismatch hits; we can stop now
				return;
			}
		}
		// no hits
		assert(!oneHit || params().sink().numProvisionalHits() == 0);
		return;
	}
	const SideLocus& topSideLocus() const { return _topSideLocus; }
	const SideLocus& botSideLocus() const { return _botSideLocus; }
	/// Next two accessors break encapsulation for _top/_botSideLocus;
	/// might be worth fixing
	SideLocus* topSideLocusPtr() { return &_topSideLocus; }
	SideLocus* botSideLocusPtr() { return &_botSideLocus; }

private:
	const Ebwt<TStr>& _ebwt; // EBWT
	const EbwtSearchParams<TStr>& _params; // search params
	RandomSource _rand; // random source for picking a random hit
	uint32_t _top;      // top index
	uint32_t _bot;      // bot index
	const TStr& _query; // query string
	const string& _query_name; // query name
	const string& _query_quals; // query qualities string
	String<uint32_t> _remainders; // space left to explore
	String<uint8_t> _tried; // denotes characters originally tried at each pos
	String<uint32_t> _tops; // tops we might want to backtrack to
	String<uint32_t> _bots; // bots we might want to backtrack to
	int _firstNzRemainder; // index of first non-zero element in _remainders
	uint32_t _qlen;     // length of query string
    uint32_t _qidx;     // offset into query (moves from high to low offsets)
    SideLocus _topSideLocus; // top EBWT side coordinates
    SideLocus _botSideLocus; // bot EBWT side coordinates
    uint32_t _narrowHalfLen;
	uint32_t _mism;
};

///////////////////////////////////////////////////////////////////////
//
// Functions for printing and sanity-checking Ebwts
//
///////////////////////////////////////////////////////////////////////

/**
 * Given a range of positions in the EBWT array within the BWT portion
 * of a forward side, print the characters at those positions along
 * with a summary occ[] array. 
 */
template<typename TStr>
void Ebwt<TStr>::printRangeFw(uint32_t begin, uint32_t end) const {
	assert(isInMemory());
	uint32_t occ[] = {0, 0, 0, 0};
	assert_gt(end, begin);
	for(uint32_t i = begin; i < end; i++) {
		uint8_t by = this->ebwt()[i];
		for(int j = 0; j < 4; j++) {
			// Unpack from lowest to highest bit pair
			int twoBit = unpack_2b_from_8b(by, j);
			occ[twoBit]++;
			cout << toDna(twoBit);
		}
		assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) & 3);
	}
	cout << ":{" << occ[0] << "," << occ[1] << "," << occ[2] << "," << occ[3] << "}" << endl;
}

/**
 * Given a range of positions in the EBWT array within the BWT portion
 * of a backward side, print the characters at those positions along
 * with a summary occ[] array. 
 */
template<typename TStr>
void Ebwt<TStr>::printRangeBw(uint32_t begin, uint32_t end) const {
	assert(isInMemory());
	uint32_t occ[] = {0, 0, 0, 0};
	assert_gt(end, begin);
	for(uint32_t i = end-1; i >= begin; i--) {
		uint8_t by = this->ebwt()[i];
		for(int j = 3; j >= 0; j--) {
			// Unpack from lowest to highest bit pair
			int twoBit = unpack_2b_from_8b(by, j);
			occ[twoBit]++;
			cout << toDna(twoBit);
		}
		assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) & 3);
		if(i == 0) break;
	}
	cout << ":{" << occ[0] << "," << occ[1] << "," << occ[2] << "," << occ[3] << "}" << endl;
}

/**
 * Check that the ebwt array is internally consistent up to (and not
 * including) the given side index by re-counting the chars and
 * comparing against the embedded occ[] arrays.
 */
template<typename TStr>
void Ebwt<TStr>::sanityCheckUpToSide(int upToSide) const {
	assert(isInMemory());
	uint32_t occ[] = {0, 0, 0, 0};
	uint32_t occ_save[] = {0, 0};
	uint32_t cur = 0; // byte pointer
	const EbwtParams& eh = this->eh();
	bool fw = false;
	while(cur < (upToSide * eh.sideSz())) {
		assert_leq(cur + eh.sideSz(), eh.ebwtTotLen());
		for(uint32_t i = 0; i < eh.sideBwtSz(); i++) {
			uint8_t by = this->ebwt()[cur + (fw ? i : eh.sideBwtSz()-i-1)];
			for(int j = 0; j < 4; j++) {
				// Unpack from lowest to highest bit pair
				int twoBit = unpack_2b_from_8b(by, fw ? j : 3-j);
				occ[twoBit]++;
				//if(_verbose) cout << toDna(twoBit);
			}
			assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) % 4);
		}
		assert_eq(0, (occ[0] + occ[1] + occ[2] + occ[3]) % eh.sideBwtLen());
		if(fw) {
			// Finished forward bucket; check saved [G] and [T]
			// against the two uint32_ts encoded here
			ASSERT_ONLY(uint32_t *u32ebwt = reinterpret_cast<uint32_t*>(&this->ebwt()[cur + eh.sideBwtSz()]));
			ASSERT_ONLY(uint32_t gs = u32ebwt[0]);
			ASSERT_ONLY(uint32_t ts = u32ebwt[1]);
			assert_eq(gs, occ_save[0]);
 			assert_eq(ts, occ_save[1]);
			fw = false;
		} else {
			// Finished backward bucket; check current [A] and [C]
			// against the two uint32_ts encoded here
			ASSERT_ONLY(uint32_t *u32ebwt = reinterpret_cast<uint32_t*>(&this->ebwt()[cur + eh.sideBwtSz()]));
			ASSERT_ONLY(uint32_t as = u32ebwt[0]);
			ASSERT_ONLY(uint32_t cs = u32ebwt[1]);
			assert(as == occ[0] || as == occ[0]-1); // one 'a' is a skipped '$' and doesn't count toward occ[]
			assert_eq(cs, occ[1]);
 			occ_save[0] = occ[2]; // save gs
 			occ_save[1] = occ[3]; // save ts
			fw = true;
		}
		cur += eh.sideSz();
	}
}

/**
 * Sanity-check various pieces of the Ebwt
 */
template<typename TStr>
void Ebwt<TStr>::sanityCheckAll() const {
	const EbwtParams& eh = this->eh();
	assert(isInMemory());
	// Check ftab
	for(uint32_t i = 1; i < eh.ftabLen(); i++) {
		assert_geq(this->ftabHi(i), this->ftabLo(i-1));
		assert_geq(this->ftabLo(i), this->ftabHi(i-1));
		assert_leq(this->ftabHi(i), eh.bwtLen()+1);
	}
	assert_eq(this->ftabHi(eh.ftabLen()-1), eh.bwtLen());

	// Check offs
	int seenLen = (eh.bwtLen() + 31) >> 5; 
	uint32_t *seen;
	try {
		seen = new uint32_t[seenLen]; // bitvector marking seen offsets
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating seen[] at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}
	bzero(seen, 4 * seenLen);
	uint32_t offsLen = eh.offsLen();
	for(uint32_t i = 0; i < offsLen; i++) {
		assert_lt(this->offs()[i], eh.bwtLen());
		int w = this->offs()[i] >> 5;
		int r = this->offs()[i] & 31;
		assert_eq(0, (seen[w] >> r) & 1); // shouldn't have been seen before
		seen[w] |= (1 << r);
	}
	delete[] seen;
	
	// Check nPat
	assert_gt(this->nPat(), 0);
	
	// Check plen
	for(uint32_t i = 0; i < this->nPat(); i++) {
		assert_gt(this->plen()[i], 0);
	}

	// Check pmap/plen
	for(uint32_t i = 0; i < eh.numChunks()*2; i += 2) {
		assert_lt(this->pmap()[i], this->nPat());             // valid pattern id
		if(i > 0) { assert_geq(this->pmap()[i], this->pmap()[i-2]); } // pattern id in order
		assert_lt(this->pmap()[i+1], this->plen()[this->pmap()[i]]); // valid offset into that pattern
		assert_eq(this->pmap()[i+1] & eh.chunkMask(), this->pmap()[i+1]);
		if(i >= 2 && this->pmap()[i] == this->pmap()[i-2]) {  // same pattern as last entry?
			assert_eq(this->pmap()[i+1], this->pmap()[i-2+1] + eh.chunkLen())
		}
	}

	// Check ebwt
	sanityCheckUpToSide(eh.numSides());
	VMSG_NL("Ebwt::sanityCheck passed");
}

///////////////////////////////////////////////////////////////////////
//
// Functions for searching Ebwts
//
///////////////////////////////////////////////////////////////////////

/**
 * Return the final character in row i (i.e. the i'th character in the
 * BWT transform).  Note that the 'L' in the name of the function
 * stands for 'last', as in the literature.
 */
template<typename TStr>
inline int Ebwt<TStr>::rowL(const SideLocus& l) const {
	// Extract and return appropriate bit-pair
	return unpack_2b_from_8b(l._side[l._by], l._bp);
}

/**
 * Tricky-bit-bashing population count function for 64-bit argument.
 */
inline static int pop(uint64_t x) {
   x = x - ((x >> 1) & 0x5555555555555555llu);
   x = (x & 0x3333333333333333llu) + ((x >> 2) & 0x3333333333333333llu);
   x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Fllu;
   x = x + (x >> 8);
   x = x + (x >> 16);
   x = x + (x >> 32);
   return x & 0x3F;
} 

/**
 * Tricky-bit-bashing bitpair counting for given two-bit value (0-3)
 * within a 64-bit argument.
 */
inline static int countInU64(int c, uint64_t dw) {
	uint64_t dwA  = dw &  0xAAAAAAAAAAAAAAAAllu;
	uint64_t dwNA = dw & ~0xAAAAAAAAAAAAAAAAllu;
	int ret;
	switch(c) {
	case 0:
		ret = 32 - pop((dwA >> 1) | dwNA);
		break;
	case 1:
		ret = pop(~(dwA >> 1) & dwNA);
		break;
	case 2:
		ret = pop((dwA >> 1) & ~dwNA);
		break;
	case 3:
		ret = pop((dwA >> 1) & dwNA);
		break;
	default:
		throw;
	}
	assert_leq(ret, 32);
	assert_geq(ret, 0);
	return ret;
}

/**
 * Counts the number of occurrences of character 'c' in the given Ebwt
 * side up to (but not including) the given byte/bitpair (by/bp).
 * 
 * This is a performance-critical function.  This is the top search-
 * related hit in the time profile.
 */
template<typename TStr>
inline uint32_t Ebwt<TStr>::countUpTo(const SideLocus& l, int c) const {
	uint32_t cCnt = 0;
	int i = 0;
	// Count occurrences of c in each 64-bit (using bit trickery);
	// note: this seems does not seem to lend a significant boost to
	// performance.  If you comment out this whole loop (which won't
	// affect correctness - it will just cause the following loop to
	// take up the slack) then runtime does not change noticeably.
	// Someday the countInU64() and pop() functions should be
	// vectorized/SSE-ized in case that helps.
	for(; i+7 < l._by; i += 8) {
		cCnt += countInU64(c, *(uint64_t*)&l._side[i]);
	}
	// Count occurences of c in the rest of the side (using LUT)
	for(; i < l._by; i++) {
		cCnt += cCntLUT[c][l._side[i]];
	}
	// Count occurences of c in the rest of the byte
	for(i = 0; i < l._bp; i++) {
		if(unpack_2b_from_8b(l._side[l._by], i) == c) cCnt++;
	}
	return cCnt;
}

/**
 * Count all occurrences of character c from the beginning of the
 * forward side to <by,bp> and add in the occ[] count up to the side
 * break just prior to the side.
 */
template<typename TStr>
inline uint32_t Ebwt<TStr>::countFwSide(const SideLocus& l, int c) const {
	const EbwtParams& eh = this->eh();
	int by = l._by;
	int bp = l._bp;
	uint8_t *ebwtSide = l._side;
	uint32_t sideByteOff = l._sideByteOff;
	assert_lt(c, 4);
	assert_geq(c, 0);
	assert_lt(by, (int)eh.sideBwtSz());
	assert_geq(by, 0);
	assert_lt(bp, 4);
	assert_geq(bp, 0);
	uint32_t cCnt = countUpTo(l, c);
	assert_leq(cCnt, eh.sideBwtLen());
	// Now factor in the occ[] count at the side break
	uint32_t *ac = reinterpret_cast<uint32_t*>(ebwtSide - 8);                // prev
	uint32_t *gt = reinterpret_cast<uint32_t*>(ebwtSide + eh.sideSz() - 8); // cur
	if(c == 0 && sideByteOff <= _zEbwtByteOff && sideByteOff + by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if(sideByteOff + by > _zEbwtByteOff ||
		   sideByteOff + by == _zEbwtByteOff && bp > _zEbwtBpOff)
		{
			cCnt--; // Adjust for '$' looking like an 'A'
		}
	}
	uint32_t ret;
	switch(c) {
		case 0: ret = cCnt + ac[0] + this->_fchr[c]; break;
		case 1: ret = cCnt + ac[1] + this->_fchr[c]; break;
		case 2: ret = cCnt + gt[0] + this->_fchr[c]; break;
		case 3: ret = cCnt + gt[1] + this->_fchr[c]; break;
		default: throw;
	}
#ifndef NDEBUG
	assert_leq(ret, this->_fchr[c+1]); // can't have jumpded into next char's section
	if(c == 0) {
		assert_leq(cCnt, eh.sideBwtLen());
	} else {
		assert_lt(ret, eh.bwtLen());
	}
#endif
	return ret;
}

/**
 * Count all instances of character c from <by,bp> to the logical end
 * (actual beginning) of the backward side, and subtract that from the
 * occ[] count up to the side break.
 */
template<typename TStr>
inline uint32_t Ebwt<TStr>::countBwSide(const SideLocus& l, int c) const {
	const EbwtParams& eh = this->eh();
	int by = l._by;
	int bp = l._bp;
	uint8_t *ebwtSide = l._side;
	uint32_t sideByteOff = l._sideByteOff;
	assert_lt(c, 4);
	assert_geq(c, 0);
	assert_lt(by, (int)eh.sideBwtSz());
	assert_geq(by, 0);
	assert_lt(bp, 4);
	assert_geq(bp, 0);
	uint32_t cCnt = countUpTo(l, c);
	if(unpack_2b_from_8b(ebwtSide[by], bp) == c) cCnt++;
	assert_leq(cCnt, eh.sideBwtLen());
	// Now factor in the occ[] count at the side break
	uint32_t *ac = reinterpret_cast<uint32_t*>(ebwtSide + eh.sideSz() - 8);     // cur
	uint32_t *gt = reinterpret_cast<uint32_t*>(ebwtSide + (2*eh.sideSz()) - 8); // next
	assert_leq(ac[0], eh.numSides() * eh.sideBwtLen()); // b/c it's used as padding
	assert_lt(ac[1], eh.len()); assert_lt(gt[0], eh.len()); assert_lt(gt[1], eh.len());
	if(c == 0 && sideByteOff <= _zEbwtByteOff && sideByteOff + by >= _zEbwtByteOff) {
		// Adjust for the fact that we represented $ with an 'A', but
		// shouldn't count it as an 'A' here
		if(sideByteOff + by > _zEbwtByteOff ||
		   sideByteOff + by == _zEbwtByteOff && bp >= _zEbwtBpOff)
		{
			if(_verbose) cout << "Adjusting for $" << endl;
			cCnt--;
		}
	}
	uint32_t ret;
	switch(c) {
		case 0: ret = ac[0] - cCnt + this->_fchr[c]; break;
		case 1: ret = ac[1] - cCnt + this->_fchr[c]; break;
		case 2: ret = gt[0] - cCnt + this->_fchr[c]; break;
		case 3: ret = gt[1] - cCnt + this->_fchr[c]; break;
		default: throw;
	}
#ifndef NDEBUG
	assert_leq(ret, this->_fchr[c+1]); // can't have jumpded into next char's section
	if(c == 0) {
		assert_leq(cCnt, eh.sideBwtLen());
	} else {
		assert_lt(ret, eh.bwtLen());
	}
#endif
	return ret;
}

/**
 * Given row i, return the row that the LF mapping maps i to.
 */
template<typename TStr>
inline uint32_t Ebwt<TStr>::mapLF(const SideLocus& l) const {
	uint32_t ret;
	int c = unpack_2b_from_8b(l._side[l._by], l._bp);
	assert_lt(c, 4);
	assert_geq(c, 0);
	if(l._fw) ret = countFwSide(l, c); // Forward side
	else      ret = countBwSide(l, c); // Backward side
	assert_lt(ret, this->eh().bwtLen());
	return ret;
}

/**
 * Given row i and character c, return the row that the LF mapping maps
 * i to on character c.
 */
template<typename TStr>
inline uint32_t Ebwt<TStr>::mapLF(const SideLocus& l, int c) const {
	uint32_t ret;
	assert_lt(c, 4);
	assert_geq(c, 0);
	if(l._fw) ret = countFwSide(l, c); // Forward side
	else      ret = countBwSide(l, c); // Backward side
	assert_lt(ret, this->eh().bwtLen());
	return ret;
}

/**
 * Like mapLF1, except it returns 0xffffffff if the character int the
 * final column of row i is not c.  This is an optimization for the
 * (hopefully common) case where we're matching and the top and bottom
 * arrows are separated by a single row. 
 * 
 * This was at the top of the profile for a recent run:
 * 
 *   %    self              
 *  time  seconds   calls   name
 * 
 *  24.13  38.45 151119757  Ebwt<>::mapLF1(uint, int) const
 *  21.42  34.14 267596260  Ebwt<>::countUpTo(uchar*, int, int, int) const
 *   8.04  12.81 133825619  Ebwt<>::countBwSide(uchar*, uint, int, int, int) const
 *   5.35   8.53 133770641  Ebwt<>::countFwSide(uchar*, uint, int, int, int) const
 * 
 * Is it because this is the first place to read from the side, thus
 * taking the cache miss?
 */
template<typename TStr>
inline uint32_t Ebwt<TStr>::mapLF1(const SideLocus& l, int c) const {
	uint32_t ret;
	if(unpack_2b_from_8b(l._side[l._by], l._bp) != c) { // L2 miss?
		return 0xffffffff;
	}
	assert_lt(c, 4);
	assert_geq(c, 0);
	if(l._fw) ret = countFwSide(l, c); // Forward side
	else      ret = countBwSide(l, c); // Backward side
	assert_lt(ret, this->eh().bwtLen());
	return ret;
}

/**
 * Report a potential match at offset 'off' with pattern length
 * 'matchlen'.  We must be careful to filter out spurious matches that
 * fall partially within the padding that separates texts.
 */
template<typename TStr>
inline bool Ebwt<TStr>::report(uint32_t off,
                               uint32_t top,
                               uint32_t bot,
                               EbwtSearchState<TStr>& s) const
{
	VMSG_NL("In report");
	assert_lt(off, this->eh().len());
	if(s.params().arrowMode()) {
		// Call reportHit with a bogus genome position; in this mode,
		// all we care about are the top and bottom arrows
		s.params().reportHit(
				s,                   // state
				make_pair(0, 0),     // (bogus) position
				make_pair(top, bot), // arrows
				0,                   // (bogus) tlen
				s.qlen(),            // qlen
				bot-top-1);          // # other hits
		return true;
	}
	// Check whether our match overlaps with the padding between two
	// texts, in which case the match is spurious
	uint32_t ptabOff = (off >> this->eh().chunkRate())*2;
	uint32_t coff = off & ~(this->eh().chunkMask()); // offset into chunk
	uint32_t tidx = this->pmap()[ptabOff];           // id of text matched
	uint32_t toff = this->pmap()[ptabOff+1];         // hit's offset in text
	uint32_t tlen = this->plen()[tidx];              // length of text
	assert_lt(toff, tlen);
	assert_lt(tidx, this->nPat());
	if(toff + coff + s.qlen() > tlen) {
		// Spurious result
		return false;
	} else {
		// Genuine result
		if(_verbose) {
			cout << "report tidx=" << tidx << ", off=" << (toff+coff) << ", absoff=" << off << ", toff=" << toff << endl;
		}
		s.params().reportHit(
				s,                            // state
				make_pair(tidx, toff + coff), // position
				make_pair(top, bot),          // arrows
				tlen,                         // tlen
				s.qlen(),                     // qlen
				bot-top-1);                   // # other hits
		return true;
	}
}

/**
 * Report a result.  Involves walking backwards along the original
 * string by way of the LF-mapping until we reach a marked SA row or
 * the row corresponding to the 0th suffix.  A marked row's offset
 * into the original string can be read directly from the this->offs()[]
 * array.
 */
template<typename TStr>
inline bool Ebwt<TStr>::reportChaseOne(uint32_t i,
                                       EbwtSearchState<TStr>& s,
                                       SideLocus *l = NULL) const
{
	VMSG_NL("In reportChaseOne");
	assert(!s.params().arrowMode());
	uint32_t off;
	uint32_t jumps = 0;
	bool own_locus = false;
	if(l == NULL) {
		own_locus = true;
		l = new SideLocus(i, this->_eh, this->_ebwt);
	}
	assert(l != NULL);
	// Walk along until we reach the next marked row to the left
	while(((i & this->_eh.offMask()) != i) && i != _zOff) {
		// Not a marked row; walk left one more char
		s.params().stats().incPushthrough(s, true, true);
		uint32_t newi = mapLF(*l); // calc next row
		assert_neq(newi, i);
		i = newi;                                  // update row
		l->initFromRow(i, this->_eh, this->_ebwt); // update locus
		jumps++;
	}
	// This is a marked row
	if(i == _zOff) {
		off = jumps;
		VMSG_NL("reportChaseOne found zoff off=" << off << " (jumps=" << jumps << ")");
	} else {
		off = this->offs()[i >> this->_eh.offRate()] + jumps;
		VMSG_NL("reportChaseOne found off=" << off << " (jumps=" << jumps << ")");
	}
	if (own_locus) delete l;
	return report(off, s.top(), s.bot(), s);
}

#define NUM_SAMPLES 10

/**
 * Report a randomly-selected sample from a range of hits.  The goal is to 
 * reduce the time spent enumerating a large range of hits and is predicated
 * on the assumption that the vast majority of those hits have similar
 * extensions.
 */
template<typename TStr>
inline bool Ebwt<TStr>::reportChaseSample(EbwtSearchState<TStr>& s) const
{	
	assert(!s.params().arrowMode());
	assert_gt(s.spread(), NUM_SAMPLES);
	set<uint32_t> sampled_hits;
	bool reported = false;
	while (sampled_hits.size() < NUM_SAMPLES)
	{
		uint32_t r = s.top() + (s.rand().nextU32() % s.spread());
		if (sampled_hits.find(r) == sampled_hits.end())
		{
			sampled_hits.insert(r);
			if(_verbose) cout << "reportChaseRange r: " << r << endl;
			// TODO: Could calculate i's locus here, as an offset from
			// s.top()'s locus
			if (reportChaseOne(r, s)) reported = true;
		}
	}
	return reported;
}

/**
 * Report a range of results by chasing down each one individually.
 */
template<typename TStr>
inline bool Ebwt<TStr>::reportChaseRange(EbwtSearchState<TStr>& s) const
{
	assert(!s.params().arrowMode());
	assert_gt(s.spread(), 1);
	bool reported = false;
	for(uint32_t i = s.top(); i < s.bot(); i++)  {
		if(_verbose) cout << "reportChaseRange i: " << i << endl;
		// TODO: Could calculate i's locus here, as an offset from
		// s.top()'s locus
		if(reportChaseOne(i, s)) reported = true;
	}
	return reported;
}

/**
 * Report one result.  
 * 
 * Returns true if one or more results were reported.  If false is
 * returned, that means that one or more results were spurious becuase
 * they overlapped padding between two text sequences.
 */
template<typename TStr>
inline bool Ebwt<TStr>::reportOne(uint32_t i,
                                  EbwtSearchState<TStr>& s,
                                  SideLocus *l = NULL) const
{
	VMSG_NL("In reportOne");
	assert_eq(1, s.spread());
	if(s.params().suppressHits()) return true; // return without reporting
	s.params().stats().addExactHits(s);
	if(s.params().arrowMode()) {
		// Don't chase anything; just report arrows and bail
		report(i, s.top(), s.bot(), s);
		return true;
	}
	return reportChaseOne(i, s, l);
}

/**
 * Given that our search is finished and there are multiple results
 * between top and bot, decide how to report them.  One option is to
 * report all the results by chasing each one down 1-by-1, using
 * reportChaseRange().  Another option is to report all using a
 * "multichase" (iterative 4-way push-throughs).  Another option is
 * to report a single hit chosen at random.
 */
template<typename TStr>
inline bool Ebwt<TStr>::reportMultiple(EbwtSearchState<TStr>& s) const
{	
	VMSG_NL("In reportMultiple");
	assert_gt(s.spread(), 1);
	if(s.params().suppressHits()) return true; // return without reporting
	s.params().stats().addExactHits(s);
	if(s.params().arrowMode()) {
		// Don't chase anything; just report arrows and bail
		report(0, s.top(), s.bot(), s);
		return true;
	}
	switch(s.params().multiHitPolicy()) {
		case MHP_CHASE_ALL: {
			return reportChaseRange(s);
			break;
		}
		case MHP_CHASE_SAMPLE: {
			if (s.spread() > SAMPLE_THRESH)
				return reportChaseSample(s);
			else
				return reportChaseRange(s);
			break;
		}
		case MHP_MULTICHASE_ALL: {
			throw runtime_error("MHP_MULTICHASE_ALL not yet implemented");
			return false;
		}
		case MHP_PICK_1_RANDOM: {
			uint32_t r = s.top() + (s.rand().nextU32() % s.spread());
			for(uint32_t i = 0; i < s.spread(); i++) {
				uint32_t ri = r + i;
				if(ri >= s.bot()) ri -= s.spread();
				// TODO: Could calculate r's locus here, as an offset
				// from s.top()'s locus
				if(reportChaseOne(ri, s)) return true;
			}
			return false;
		}
		default: {
			throw runtime_error("Bad multi-hit policy");
			return false;
		}
	}
}

/**
 * Search the Ebwt starting with the given row.
 */
template<typename TStr>
inline bool Ebwt<TStr>::searchFinish1(EbwtSearchState<TStr>& s) const
{
	VMSG_NL("In searchFinish1");
	assert_eq(1, s.spread());
	int jumps = 0;
	uint32_t off = 0xffffffff;
	assert_lt(s.top(), this->eh().len());
	while(!s.closed() && !s.qExhausted()) {
		VMSG_NL("searchFinish1: qidx: " << s.qidx() << ", top: " << s.top());
		uint32_t r = s.mapLF1(*this);
		if(r == 0xffffffff) return false; // no match
		// Did we stumble into a row with a known offset?  If so, make
		// note so that don't have to chase it later.
		if(!s.params().arrowMode()) {
			if(off == 0xffffffff && ((r & this->eh().offMask()) == r || r == _zOff)) {
				off = ((r == _zOff) ? 0 : this->offs()[r >> this->eh().offRate()]);
				assert_lt(off, this->eh().len());
				jumps = 0; // reset
			} else {
				jumps++;
			}
		}
	}
	assert(!s.closed());
	assert_neq(s.top(), 0xffffffff);
	assert_geq(jumps, 0);
	if(s.params().suppressHits()) return true; // return without reporting
	if(s.params().arrowMode()) {
		// Just return arrows; we haven't tried to set off and we're
		// not going to try to chase the result
		return report(0, s.top(), s.bot(), s); 
	}
	if(off != 0xffffffff) {
		// Good luck - we previously made note of a marked row and
		// need not chase it any further 
		off -= jumps;
		assert_lt(off, this->eh().len());
		assert_lt((unsigned int)jumps, s.qlen());
		assert_geq(jumps, 0);
		VMSG_NL("searchFinish1 pre-found off=" << off);
		s.params().stats().addExactHits(s);
		return report(off, s.top(), s.bot(), s); 
		assert_leq(s.params().sink().numProvisionalHits(), 1);
	} else {
		// Haven't yet encountered a marked row for this result;
		// chase it
		return reportOne(s.top(), s, s.topSideLocusPtr());
	}
}

/**
 * Search the Ebwt starting with the given bounds of top/bot.  If and
 * when the difference between top and bot shrinks to 1, hand control
 * to searchFinish1, which is optimized for that (hopefully common)
 * case.  If the entire pattern is consumed and there is bot is still
 * greater top, then that range represents a set of potential matches
 * and we hand it off to reportChaseRange.
 */
template<typename TStr>
inline bool Ebwt<TStr>::searchFinish(EbwtSearchState<TStr>& s) const
{
	VMSG_NL("In searchFinish");
	assert(s.repOk());
	while(!s.closed() && !s.qExhausted()) {
		s.params().stats().incPushthrough(s);
		s.mapLF(*this);
		assert(s.repOk());
		if(s.spread() == 1) {
			// Switch to searching loop specialized for the case where
			// one potential result remains
			return searchFinish1(s);
		}
	}
	if(!s.closed()) {
		// Multiple potential results to report
		if(s.spread() == 1) {
			s.params().stats().addExactHits(s);
			return reportOne(s.top(), s, s.topSideLocusPtr());
		} else {
			return reportMultiple(s);
		}
	} else {
		s.params().stats().incStopAt(s);
		return false;
	}
}

/**
 * Do an ftab lookup to find initial top and bottom then call
 * searchFinish().
 * 
 * The Ebwt must be in memory.
 */
template<typename TStr>
inline void Ebwt<TStr>::searchWithFtab(EbwtSearchState<TStr>& s) const
{
	typedef typename Value<TStr>::Type TVal;
	assert(isInMemory());
	assert(s.initialized());
	assert(s.qAtBeginning());
	s.params().stats().incTry(s);
	int ftabChars = this->eh().ftabChars();
	assert_geq(s.qlen(), (unsigned int)ftabChars);
	// Rightmost char gets least significant bit-pair
	const TStr& qry = s.query();
	uint32_t ftabOff = qry[s.qlen() - ftabChars];
	for(int i = ftabChars - 1; i > 0; i--) {
		ftabOff <<= 2;
		ftabOff |= (int)(TVal)qry[s.qlen()-i];
	}
	assert_lt(ftabOff, this->eh().ftabLen()-1);
	if(_verbose) cout << "Looking up in ftab with: " << u32ToDna(ftabOff, ftabChars) << endl;
	s.setTopBot(ftabHi(ftabOff), ftabLo(ftabOff+1));
	s.subQidx(ftabChars);
	searchFinish(s);
}

/**
 * Do an fchr lookup to find (coarse) initial top and bottom then call
 * searchFinish().
 * 
 * The Ebwt must be in memory.
 */
template<typename TStr>
inline void Ebwt<TStr>::searchWithFchr(EbwtSearchState<TStr>& s) const
{
	assert(isInMemory());
	assert(s.initialized());
	assert(s.qAtBeginning());
	s.params().stats().incTry(s);
	s.setTopBot(this->_fchr[s.chr()], this->_fchr[s.chr()+1]);
	s.decQidx();
	searchFinish(s);
}

template<typename TStr>
void Ebwt<TStr>::search(EbwtSearchState<TStr>& s,
                        EbwtSearchParams<TStr>& params,
                        bool inexact = false) const
{
	assert(isInMemory());
	params.setBacktracking(false);
	if(s.qlen() >= (unsigned int)this->eh().ftabChars() &&
	   // If we're backtracking, don't let ftab skip over any of the
	   // revisitable region
	   (!inexact ||
	    (unsigned int)this->eh().ftabChars() <= (s.qlen()>>1)))
	{
		searchWithFtab(s);
	} else {
		searchWithFchr(s);
	}
}

/**
 * Search the Ebwt; if the query is long enough, start by looking up
 * the initial top and bottom in the ftab.  If the input string is too
 * short for the ftab, get the initial top and bottom from the fchr.
 * 
 * The Ebwt must be in memory.
 */
template<typename TStr>
void Ebwt<TStr>::search(const TStr& qry,
                        EbwtSearchParams<TStr>& params,
                        uint32_t seed = 0,
                        bool inexact = false) const
{
	assert(isInMemory());
	assert(!empty(qry));
	EbwtSearchState<TStr> s(*this, qry, params, seed);
	search(s, params, inexact);
}

/**
 * Search the Ebwt; if the query is long enough, start by looking up
 * the initial top and bottom in the ftab.  If the input string is too
 * short for the ftab, get the initial top and bottom from the fchr.
 * 
 * Returns true iff one or more exact hits were found
 * 
 * The Ebwt must be in memory.
 */
template<typename TStr>
bool Ebwt<TStr>::search1MismatchOrBetter(EbwtSearchState<TStr>& s,
                                         EbwtSearchParams<TStr>& params,
                                         bool allowExactHits = true) const
{
	assert(isInMemory());
	params.setBacktracking(false);
	uint64_t hits = s.params().sink().numHits();
	bool oneHit = (params.multiHitPolicy() == MHP_PICK_1_RANDOM);
	// If exact hits are not allowed, then suppress hits here; the call
	// to search() will still fill the EbwtSearchState with the
	// appropriate backtrack points
	params.setSuppressHits(!allowExactHits);
	search(s, params, true);
	// Un-suppress hits
	params.setSuppressHits(false);
	assert(s.params().sink().numHits() == hits || allowExactHits);
	if(s.params().sink().numHits() > hits && oneHit) {
		// One or more exact hits found
		return true;
	}
	// No exact hits found; try backtracking for 1-mismatch hits
	params.setBacktracking(true);
	s.backtrack1(*this);
	return false;
}

/**
 * Transform this Ebwt into the original string in linear time by using
 * the LF mapping to walk backwards starting at the row correpsonding
 * to the end of the string.  The result is written to s.  The Ebwt
 * must be in memory.
 */
template<typename TStr>
void Ebwt<TStr>::restore(TStr& s) const {
	assert(isInMemory());
	resize(s, this->eh().len(), Exact());
	uint32_t jumps = 0;
	uint32_t i = this->_eh.len(); // should point to final SA elt (starting with '$')
	SideLocus l(i, this->_eh, this->_ebwt);
	while(i != _zOff) {
		assert_lt(jumps, this->_eh.len());
		if(_verbose) cout << "restore: i: " << i << endl;
		// Not a marked row; go back a char in the original string
		uint32_t newi = mapLF(l);
		assert_neq(newi, i);
		s[this->eh().len() - jumps - 1] = rowL(l);
		i = newi;
		l.initFromRow(i, this->_eh, this->_ebwt);
		jumps++;
	}
	assert_eq(jumps, this->_eh.len());
}

///////////////////////////////////////////////////////////////////////
//
// Functions for reading and writing Ebwts
//
///////////////////////////////////////////////////////////////////////


/**
 * Read an Ebwt from file with given filename.
 */
template<typename TStr>
EbwtParams Ebwt<TStr>::readIntoMemory(bool justHeader,
                                      const string& in1,
                                      const string& in2)
{
	bool bigEndian; // dummy; caller doesnn't care
	return readIntoMemory(justHeader, in1, in2, bigEndian);
}

/**
 * Read an Ebwt from a file with given filename.  The endianness of the
 * data read is installed in the bigEndian out parameter (true=big).
 */
template<typename TStr>
EbwtParams Ebwt<TStr>::readIntoMemory(bool justHeader,
                                      const string& in1,
                                      const string& in2,
                                      bool& bigEndian)
{
	// Initialize our primary and secondary input-stream fields
	if(!_in1.is_open()) {
		if(this->verbose()) cout << "Opening \"" << in1 << "\"" << endl;
		_in1.open(in1.c_str(), ios_base::in | ios::binary);
		if(!_in1.is_open()) {
			throw EbwtFileOpenException("Cannot open file " + in1);
		}
	}
	assert(_in1.is_open());
	assert(_in1.good());
	assert_eq((streamoff)_in1.tellg(), ios::beg);
	if(!_in2.is_open()) {
		if(this->verbose()) cout << "Opening \"" << in2 << "\"" << endl;
		_in2.open(in2.c_str(), ios_base::in | ios::binary);
		if(!_in2.is_open()) {
			throw EbwtFileOpenException("Cannot open file " + in2);
		}
	}
	assert(_in2.is_open());
	assert(_in2.good());
	assert_eq((streamoff)_in2.tellg(), ios::beg);
	EbwtParams eh(readIntoMemory(justHeader, bigEndian));
	return eh;
}

/**
 * Read an Ebwt into memory (specifically, into the fields of this Ebwt
 * object) from the object's input-stream pair fields.
 */
template<typename TStr>
EbwtParams Ebwt<TStr>::readIntoMemory(bool justHeader) {
	bool bigEndian; // dummy; caller doesn't care
	return readIntoMemory(justHeader, bigEndian);
}

/**
 * Read an Ebwt from an input-stream pair.  The endianness of the data
 * read is installed in the bigEndian out parameter (true=big).  The
 * _in1 and _in2 istreams must already be initialized and open. 
 * 
 * The caller is responsible for ensuring that both istreams are set up
 * to point just before the Ebwt records.  In most cases, this means
 * that they will both be completely rewound.
 * 
 * This function rewinds the streams before returning.
 * 
 */
template<typename TStr>
EbwtParams Ebwt<TStr>::readIntoMemory(bool justHeader, bool& be) {
	// _in1 and _in2 must already be open with the get cursor at the
	// beginning and no error flags set.
	assert(_in1.is_open()); assert(_in1.good());
	assert(_in2.is_open()); assert(_in2.good());
	assert_eq((streamoff)_in1.tellg(), ios::beg);
	assert_eq((streamoff)_in2.tellg(), ios::beg);

	if(_verbose) cout << "  Reading header" << endl;

	// Read endianness hints from both streams
	be = false;
	uint32_t one = readU32(_in1, be); // 1st word of primary stream
	#ifndef NDEBUG
	assert_eq(one, readU32(_in2, be)); // should match!
	#else
	readU32(_in2, be);
	#endif
	if(one != 1) {
		assert_eq((1u<<24), one);
		assert_eq(1, endianSwapU32(one));
		be = true;
	}

	// Reads header entries one by one from primary stream
	uint32_t len          = readU32(_in1, be);
	int32_t  lineRate     = readI32(_in1, be);
	int32_t  linesPerSide = readI32(_in1, be);
	int32_t  offRate      = readI32(_in1, be);
	int32_t  ftabChars    = readI32(_in1, be);
	int32_t  chunkRate    = readI32(_in1, be);

	// Create a new EbwtParams from the entries read from primary stream
	EbwtParams eh(len, lineRate, linesPerSide, offRate, ftabChars, chunkRate);
	if(_verbose) eh.print(cout);
	uint32_t offsLen = eh.offsLen();
	uint32_t offRateDiff = 0;
	uint32_t offsLenSampled = offsLen;
	if(_overrideOffRate > offRate) {
		offRateDiff = _overrideOffRate - offRate;
	}
	offsLenSampled >>= offRateDiff;

	// Read nPat from primary stream
	this->_nPat = readI32(_in1, be);
	try {
		// Read plen from primary stream
		if(_verbose) cout << "Reading plen (" << this->nPat() << ")" << endl;
		this->_plen = new uint32_t[this->nPat()];
		if(be) {
			for(uint32_t i = 0; i < this->nPat(); i++) {
				this->_plen[i] = readU32(_in1, be);
			}
		} else {
			_in1.read((char *)this->_plen, this->nPat()*4);
			assert_eq(this->_nPat*4, (uint32_t)_in1.gcount());
		}
		for(uint32_t i = 0; i < this->nPat(); i++) {
			assert_leq(this->_plen[i], len);
			assert_gt(this->_plen[i], 0);
		}
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating plen[] in Ebwt::read()"
		     << " at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}

	// TODO: I'm not consistent on what "header" means.  Here I'm using
	// "header" to mean everything that would exist in memory if we
	// started to build the Ebwt but stopped short of the build*() step
	// (i.e. everything up to and including join()).
	if(justHeader) goto done;

	// Read pmap from primary stream
	try {
		uint32_t pmapEnts = eh.numChunks()*2;
		if(_verbose) cout << "Reading pmap (" << pmapEnts << ")" << endl;
		this->_pmap = new uint32_t[pmapEnts];
		if(be) {
			for(uint32_t i = 0; i < pmapEnts; i += 2) {
				this->_pmap[i]   = readU32(_in1, be); // pat #
				this->_pmap[i+1] = readU32(_in1, be); // pat offset
			}
		} else {
			_in1.read((char *)this->_pmap, pmapEnts*4);
			assert_eq(pmapEnts*4, (uint32_t)_in1.gcount());
		}
		for(uint32_t i = 0; i < pmapEnts; i += 2) {
			assert_lt(this->_pmap[i], this->_nPat);
			assert_lt(this->_pmap[i+1], this->_plen[this->_pmap[i]]);
		}
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating pmap[] in Ebwt::read()"
		     << " at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}

	// Allocate ebwt (big allocation) 
	if(_verbose) cout << "Reading ebwt (" << eh.ebwtTotLen() << ")" << endl;
	try {
		this->_ebwt = new uint8_t[eh.ebwtTotLen()];
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating ebwt[] in Ebwt::read()"
		     << " at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}
	// Read ebwt from primary stream
	_in1.read((char *)this->_ebwt, eh.ebwtTotLen());
	assert_eq(eh.ebwtTotLen(), (uint32_t)_in1.gcount());

	// Read zOff from primary stream
	_zOff = readU32(_in1, be);
	assert_lt(_zOff, len);
	
	try {
		// Read fchr from primary stream
		if(_verbose) cout << "Reading fchr (5)" << endl;
		this->_fchr = new uint32_t[5];
		for(int i = 0; i < 5; i++) {
			this->_fchr[i] = readU32(_in1, be);
			assert_leq(this->_fchr[i], len);
			if(i > 0) assert_geq(this->_fchr[i], this->_fchr[i-1]);
		}
		// Read ftab from primary stream
		if(_verbose) cout << "Reading ftab (" << eh.ftabLen() << ")" << endl;
		this->_ftab = new uint32_t[eh.ftabLen()];
		if(be) {
			for(uint32_t i = 0; i < eh.ftabLen(); i++)
				this->_ftab[i] = readU32(_in1, be);
		} else {
			_in1.read((char *)this->_ftab, eh.ftabLen()*4);
			assert_eq(eh.ftabLen()*4, (uint32_t)_in1.gcount());
		}
		// Read etab from primary stream
		if(_verbose) cout << "Reading eftab (" << eh.eftabLen() << ")" << endl;
		this->_eftab = new uint32_t[eh.eftabLen()];
		if(be) {
			for(uint32_t i = 0; i < eh.eftabLen(); i++)
				this->_eftab[i] = readU32(_in1, be);
		} else {
			_in1.read((char *)this->_eftab, eh.eftabLen()*4);
			assert_eq(eh.eftabLen()*4, (uint32_t)_in1.gcount());
		}
		for(uint32_t i = 0; i < eh.eftabLen(); i++) {
			if(i > 0 && this->_eftab[i] > 0) {
				assert_geq(this->_eftab[i], this->_eftab[i-1]);
			} else if(i > 0 && this->_eftab[i-1] == 0) {
				assert_eq(0, this->_eftab[i]);
			}
		}
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating fchr[], ftab[] or eftab[] in "
		     << "Ebwt::read()  at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}

	// Allocate offs (big allocation)
	try {
		if(_verbose) cout << "Reading offs (" << offsLenSampled << ")" << endl;
		this->_offs = new uint32_t[offsLenSampled];
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating offs[] in Ebwt::read()"
		     << " at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}
	if(be || offRateDiff > 0) {
		for(uint32_t i = 0; i < offsLen; i++) {
			if((i & ~(0xffffffff << offRateDiff)) != 0) {
				char tmp[4];
				_in2.read(tmp, 4);
			} else {
				this->_offs[i >> offRateDiff] = readU32(_in2, be);
			}
		}
	} else {
		_in2.read((char *)this->_offs, offsLen*4);
		assert_eq(offsLen*4, (uint32_t)_in2.gcount());
	}
	for(uint32_t i = 0; i < offsLenSampled; i++) {
		assert_leq(this->_offs[i], len);
	}

	postReadInit(eh); // Initialize fields of Ebwt not read from file
	if(_verbose) print(cout, eh);

	// The fact that _ebwt and friends actually point to something
	// (other than NULL) now signals to other member functions that the
	// Ebwt is loaded into memory.
	
  done: // Exit hatch for both justHeader and !justHeader
  
	// Be kind
	_in1.clear(); _in1.seekg(0, ios::beg); 
	_in2.clear(); _in2.seekg(0, ios::beg); 
	assert(_in1.is_open()); assert(_in1.good());
	assert(_in2.is_open()); assert(_in2.good());
	return eh;
}

/**
 * Write an extended Burrows-Wheeler transform to a pair of output
 * streams.
 * 
 * @param out1 output stream to primary file
 * @param out2 output stream to secondary file
 * @param be   write in big endian?
 */
template<typename TStr>
void Ebwt<TStr>::writeFromMemory(bool justHeader,
                                 ostream& out1,
                                 ostream& out2) const
{
	const EbwtParams& eh = this->eh();
	assert(eh.repOk());
	uint32_t be = this->toBe();
	assert(out1.good());
	assert(out2.good());
	
	// When building an Ebwt, these header parameters are known
	// "up-front", i.e., they can be written to disk immediately,
	// before we join() or buildToDisk() 
	writeI32(out1, 1, be); // endian hint for priamry stream
	writeI32(out2, 1, be); // endian hint for secondary stream
	writeU32(out1, eh.len(),          be); // length of string (and bwt and suffix array)
	writeI32(out1, eh.lineRate(),     be); // 2^lineRate = size in bytes of 1 line
	writeI32(out1, eh.linesPerSide(), be); // not used
	writeI32(out1, eh.offRate(),      be); // every 2^offRate chars is "marked"
	writeI32(out1, eh.ftabChars(),    be); // number of 2-bit chars used to address ftab
	writeI32(out1, eh.chunkRate(),    be);

	if(!justHeader) {
		assert(isInMemory());
		// These Ebwt parameters are known after the inputs strings have
		// been joined() but before they have been built().  These can
		// written to the disk next and then discarded from memory.
		writeU32(out1, this->nPat(),      be);
		for(uint32_t i = 0; i < this->nPat(); i++)
		writeU32(out1, this->plen()[i], be);
		for(uint32_t i = 0; i < eh.numChunks()*2; i++)
			writeU32(out1, this->pmap()[i], be);
	
		// These Ebwt parameters are discovered only as the Ebwt is being
		// built (in buildToDisk()).  Of these, only 'offs' and 'ebwt' are
		// terribly large.  'ebwt' is written to the primary file and then
		// discarded from memory as it is built; 'offs' is similarly
		// written to the secondary file and discarded.
		out1.write((const char *)this->ebwt(), eh.ebwtTotLen());
		writeU32(out1, this->zOff(), be);
		uint32_t offsLen = eh.offsLen();
		for(uint32_t i = 0; i < offsLen; i++)
			writeU32(out2, this->offs()[i], be);
	
		// 'fchr', 'ftab' and 'eftab' are not fully determined until the
		// loop is finished, so they are written to the primary file after
		// all of 'ebwt' has already been written and only then discarded
		// from memory.
		for(int i = 0; i < 5; i++)
			writeU32(out1, this->_fchr[i], be);
		for(uint32_t i = 0; i < eh.ftabLen(); i++)
			writeU32(out1, this->ftab()[i], be);
		for(uint32_t i = 0; i < eh.eftabLen(); i++)
			writeU32(out1, this->eftab()[i], be);
	}
}

/**
 * Given a pair of strings representing output filenames, and assuming
 * this Ebwt object is currently in memory, write out this Ebwt to the
 * specified files.
 * 
 * If sanity-checking is enabled, then once the streams have been
 * fully written and closed, we reopen them and read them into a 
 * (hopefully) exact copy of this Ebwt.  We then assert that the
 * current Ebwt and the copy match in all of their fields.
 */
template<typename TStr>
void Ebwt<TStr>::writeFromMemory(bool justHeader,
                                 const string& out1,
                                 const string& out2) const
{
	const EbwtParams& eh = this->eh();
	assert(isInMemory());
	assert(eh.repOk());
    
	ofstream fout1(out1.c_str(), ios::binary);
	ofstream fout2(out2.c_str(), ios::binary);
	writeFromMemory(justHeader, fout1, fout2);
	fout1.close();
	fout2.close();

	// Read the file back in and assert that all components match
	if(_sanity) {
		if(_verbose)
			cout << "Re-reading \"" << out1 << "\"/\"" << out2 << "\" for sanity check" << endl;
		Ebwt copy(out1, out2, _verbose, _sanity);
		assert(!isInMemory());
		copy.loadIntoMemory();
		assert(isInMemory());
	    assert_eq(eh.lineRate(),     copy.eh().lineRate());
	    assert_eq(eh.linesPerSide(), copy.eh().linesPerSide());
	    assert_eq(eh.offRate(),      copy.eh().offRate());
	    assert_eq(eh.ftabChars(),    copy.eh().ftabChars());
	    assert_eq(eh.len(),          copy.eh().len());
	    assert_eq(eh.chunkRate(),    copy.eh().chunkRate());
	    assert_eq(_zOff,             copy.zOff());
	    assert_eq(_zEbwtBpOff,       copy.zEbwtBpOff());
	    assert_eq(_zEbwtByteOff,     copy.zEbwtByteOff());
		assert_eq(_nPat,             copy.nPat());
		for(uint32_t i = 0; i < _nPat; i++)
			assert_eq(this->plen()[i], copy.plen()[i]);
		for(uint32_t i = 0; i < eh.numChunks()*2; i++)
			assert_eq(this->pmap()[i], copy.pcap()[i]);
		for(uint32_t i = 0; i < 5; i++)
			assert_eq(this->_fchr[i], copy.fchr()[i]);
		for(uint32_t i = 0; i < eh.ftabLen(); i++)
			assert_eq(this->ftab()[i], copy.ftab()[i]);
		for(uint32_t i = 0; i < eh.eftabLen(); i++)
			assert_eq(this->eftab()[i], copy.eftab()[i]);
		for(uint32_t i = 0; i < eh.offsLen(); i++)
			assert_eq(this->offs()[i], copy.offs()[i]);
		for(uint32_t i = 0; i < eh.ebwtTotLen(); i++)
			assert_eq(this->ebwt()[i], copy.ebwt()[i]);
		copy.sanityCheckAll();
		if(_verbose)
			cout << "Read-in check passed for \"" << out1 << "\"/\"" << out2 << "\"" << endl;
	}
}

///////////////////////////////////////////////////////////////////////
//
// Functions for building Ebwts
//
///////////////////////////////////////////////////////////////////////

/**
 * Join several text strings together in a way that's compatible with
 * the text-chunking scheme dictated by chunkRate parameter.
 * 
 * The non-static member Ebwt::join additionally builds auxilliary
 * arrays that maintain a mapping between chunks in the joined string
 * and the original text strings.
 */
template<typename TStr>
TStr Ebwt<TStr>::join(vector<TStr>& l, uint32_t chunkRate, uint32_t seed) {
	RandomSource rand(seed); // reproducible given same seed
	uint32_t chunkLen = 1 << chunkRate;
	uint32_t chunkMask = 0xffffffff << chunkRate;
	TStr ret;
	size_t guessLen = 0;
	for(size_t i = 0; i < l.size(); i++) {
		guessLen += length(l[i]) + chunkLen;
	}
	reserve(ret, guessLen, Exact());
	for(size_t i = 0; i < l.size(); i++) {
		TStr& s = l[i];
		assert_gt(length(s), 0);
		append(ret, s);
		// s isn't the last pattern; padding between s and the next
		// pattern may be necessary
		uint32_t diff = 0;
		uint32_t leftover = length(ret) & ~chunkMask;
		if(leftover > 0) {
			// The joined string currently ends in the middle of a
			// chunk, so we have to pad it by 'diff'
			diff = chunkLen - leftover;
			assert_gt(diff, 0);
		}
		for(size_t j = 0; j < diff; j++) {
			// Append random characters to fill the gap; note that
			// the randomness is reproducible as long as the 'seed'
			// argument is the same, which helps us to sanity-check
			// the result.
			appendValue(ret, (Dna)(rand.nextU32() & 3)); // append random junk
			assert_lt((uint8_t)ret[length(ret)-1], 4);
		}
		// Padded pattern ends on a chunk boundary
		assert_eq(length(ret), length(ret) & chunkMask);
	}
	return ret;
}

/**
 * Join several text strings together in a way that's compatible with
 * the text-chunking scheme dictated by chunkRate parameter.
 * 
 * The non-static member Ebwt::join additionally builds auxilliary
 * arrays that maintain a mapping between chunks in the joined string
 * and the original text strings.
 */
template<typename TStr>
TStr Ebwt<TStr>::join(vector<istream*>& l,
                      vector<uint32_t>& szs,
                      uint32_t sztot,
                      const RefReadInParams& refparams,
                      uint32_t chunkRate,
                      uint32_t seed)
{
	RandomSource rand(seed); // reproducible given same seed
	RefReadInParams rpcp = refparams;
	uint32_t chunkLen = 1 << chunkRate;
	uint32_t chunkMask = 0xffffffff << chunkRate;
	TStr ret;
	size_t guessLen = sztot + (szs.size() * chunkLen);
	reserve(ret, guessLen, Exact());
	for(size_t i = 0; i < l.size(); i++) {
		// For each sequence we can pull out of istream l[i]...
		assert(l[i]->good());
		assert_geq(rpcp.numSeqCutoff, -1);
		assert_geq(rpcp.baseCutoff, -1);
		bool first = true;
		while(l[i]->good() && rpcp.numSeqCutoff != 0 && rpcp.baseCutoff != 0) {
			size_t bases = fastaRefReadAppend(*l[i], ret, rpcp, first);
			if(bases == 0) continue;
			first = false;
			if(rpcp.numSeqCutoff != -1) rpcp.numSeqCutoff--;
			if(rpcp.baseCutoff != -1)   rpcp.baseCutoff -= bases;
			assert_geq(rpcp.numSeqCutoff, -1);
			assert_geq(rpcp.baseCutoff, -1);
		    // insert padding
			uint32_t diff = 0;
			uint32_t rlen = length(ret);
			uint32_t leftover = rlen & ~chunkMask;
			if(leftover > 0) {
				// The joined string currently ends in the middle of a
				// chunk, so we have to pad it by 'diff'
				diff = chunkLen - leftover;
				assert_gt(diff, 0);
			}
			for(uint32_t i = 0; i < diff; i++) {
				// Append random characters to fill the gap; note that
				// the randomness is reproducible as long as the 'seed'
				// argument is the same, which helps us to sanity-check
				// the result.
				appendValue(ret, (Dna)(rand.nextU32() & 3));
				assert_lt((uint8_t)(Dna)ret[length(ret)-1], 4);
				//appendNE(ret, rand.nextU32() & 3);
			}
			// Pattern now ends on a chunk boundary
			assert_eq(length(ret), length(ret) & chunkMask);
		}
	}
	return ret;
}

/**
 * Join several text strings together according to the text-chunking
 * scheme specified in the EbwtParams.  Ebwt fields calculated in this
 * function (_nPat, _plen, _pmap) are written directly to disk.  The
 * _nPat and _plen fields are also retained in the Ebwt, but the _pmap
 * field is not.  _pmap is relatively big, so we avoid keeping it in
 * memory unless the user specifically loads the Ebwt.
 * 
 * It is assumed, but not required, that the header values have already
 * been written to 'out1' before this function is called.
 * 
 * The static member Ebwt::join just returns a joined version of a
 * list of strings without building any of the auxilliary arrays.
 * Because the pseudo-random number generator is the same, we expect
 * this function and the static function to give the same result given
 * the same seed.
 */
template<typename TStr>
void Ebwt<TStr>::joinToDisk(vector<istream*>& l,
                            vector<uint32_t>& szs,
                            uint32_t sztot,
                            const RefReadInParams& refparams,
                            TStr& ret,
                            ostream& out1,
                            ostream& out2,
                            ostream& out3,
                            uint32_t seed = 0)
{
	RandomSource rand(seed); // reproducible given same seed
	const EbwtParams& eh = this->eh();
	RefReadInParams rpcp = refparams;
	assert_gt(szs.size(), 0);
	assert_gt(sztot, 0);
	this->_nPat = szs.size(); // store this in memory
	this->_pmap = NULL;
	writeU32(out1, this->_nPat, this->toBe());
	ASSERT_ONLY(uint32_t pmapEnts = eh.numChunks()*2);
	uint32_t pmapOff = 0;
	// Allocate plen[]
	try {
		this->_plen = new uint32_t[this->_nPat];
	} catch(bad_alloc& e) {
		cerr << "Out of memory allocating plen[] in Ebwt::join()"
		     << " at " << __FILE__ << ":" << __LINE__ << endl;
		throw e;
	}
	// For each pattern, set plen
	for(size_t i = 0; i < this->_nPat; i++) {
		this->_plen[i] = szs[i];
		writeU32(out1, this->_plen[i], this->toBe());
	}
	size_t seqsRead = 0;
	for(unsigned int i = 0; i < l.size(); i++) {
		assert(l[i]->good());
		streampos pos = l[i]->tellg();
		assert_geq(rpcp.numSeqCutoff, -1);
		assert_geq(rpcp.baseCutoff, -1);
		bool first = true;
		// For each sequence we can pull out of istream l[i]...
		while(l[i]->good() && rpcp.numSeqCutoff != 0 && rpcp.baseCutoff != 0) {
			size_t ilen = length(ret);
			size_t bases = fastaRefReadAppend(*l[i], ret, rpcp, first);
			size_t nlen = length(ret);
			if(bases == 0) continue;
			first = false;
			assert_eq(bases, this->_plen[seqsRead]);
			seqsRead++;
			if(rpcp.numSeqCutoff != -1) rpcp.numSeqCutoff--;
			if(rpcp.baseCutoff != -1)   rpcp.baseCutoff -= bases;
			assert_geq(rpcp.numSeqCutoff, -1);
			assert_geq(rpcp.baseCutoff, -1);
			// Write reference bases (but not padding) to .3.ebwt index
			{
				String<Dna, Packed<> > pchunk;
				reserve(pchunk, 4096 + 32);
				writePackedLen(out3, bases, _verbose);
				for(size_t j = 0; j < bases; j += 4096) {
					clear(pchunk);
					size_t upper = min<size_t>(ilen + j + 4096, nlen);
					append(pchunk, infix(ret, ilen + j, upper));
					size_t len = length(pchunk);
					assert_gt(len, 0);
					assert_leq(len, 4096);
					if(upper == nlen) {
						// Add a bunch of 0s on the end, clobbering the
						// existing random junk - this ensures the
						// content of the .3.ebwt file is deterministic 
						for(size_t j = 0; j < 31; j++) {
							appendValue(pchunk, (Dna)0);
						}
					}
					writePacked(out3, pchunk, len, false, _verbose);
				}
			}
		    // insert padding
			uint32_t diff = 0;
			uint32_t rlen = length(ret);
			uint32_t leftover = rlen & ~(eh.chunkMask());
			if(leftover > 0) {
				// The joined string currently ends in the middle of a
				// chunk, so we have to pad it by 'diff'
				diff = eh.chunkLen() - leftover;
				assert_gt(diff, 0);
			}
			for(uint32_t i = 0; i < diff; i++) {
				// Append random characters to fill the gap; note that
				// the randomness is reproducible as long as the 'seed'
				// argument is the same, which helps us to sanity-check
				// the result.
				appendValue(ret, (Dna)(rand.nextU32() & 3)); // append random junk
				assert_lt((uint8_t)(Dna)ret[length(ret)-1], 4);
			}
			// Pattern now ends on a chunk boundary
			assert_eq(length(ret), length(ret) & eh.chunkMask());
			// Initialize elements of the pmap that cover this pattern
			for(unsigned int j = 0; j < bases; j += eh.chunkLen()) {
				assert_lt(pmapOff+1, pmapEnts);
				pmapOff += 2;
				writeU32(out1, seqsRead-1, this->toBe()); // pattern id
				writeU32(out1, j, this->toBe()); // offset into pattern
			}
		}
		l[i]->clear();
		l[i]->seekg(pos);
		assert(!l[i]->bad());
		assert(!l[i]->fail());
		l[i]->clear();
		assert(l[i]->good());
		assert(!l[i]->eof());
		#ifndef NDEBUG
		int c = l[i]->get();
		assert_eq('>', c);
		assert(l[i]->good());
		assert(!l[i]->eof());
		l[i]->seekg(pos);
		l[i]->clear();
		assert(l[i]->good());
		assert(!l[i]->eof());
		#endif
	}
	assert_eq(pmapOff, pmapEnts); // initialized every pmap element
}


/**
 * Build an Ebwt from a string 's' and its suffix array 'sa' (which
 * might actually be a suffix array *builder* that builds blocks of the
 * array on demand).  The bulk of the Ebwt, i.e. the ebwt and offs
 * arrays, is written directly to disk.  This is by design: keeping
 * those arrays in memory needlessly increases the footprint of the
 * building process.  Instead, we prefer to build the Ebwt directly 
 * "to disk" and then read it back into memory later as necessary.
 * 
 * It is assumed that the header values and join-related values (nPat,
 * plen, pmap) have already been written to 'out1' before this function
 * is called.  When this function is finished, it will have
 * additionally written ebwt, zOff, fchr, ftab and eftab to the primary
 * file and offs to the secondary file. 
 * 
 * Assume DNA/RNA/any alphabet with 4 or fewer elements.
 * Assume occ array entries are 32 bits each.
 * 
 * @param sa            the suffix array to convert to a Ebwt 
 * @param s             the original string
 * @param out           
 */
template<typename TStr>
void Ebwt<TStr>::buildToDisk(InorderBlockwiseSA<TStr>& sa,
                             const TStr& s,
                             ostream& out1,
                             ostream& out2)
{
	const EbwtParams& eh = this->eh();
	
	assert(eh.repOk());
	assert_eq(length(s)+1, sa.size());
	assert_eq(length(s), eh.len());
	assert_gt(eh.lineRate(), 3);
	assert(sa.suffixItrIsReset());
	assert_leq((int)ValueSize<Dna>::VALUE, 4);
	
	uint32_t  len = eh.len();
	uint32_t  ftabLen = eh.ftabLen();
	uint32_t  sideSz = eh.sideSz();
	uint32_t  ebwtTotSz = eh.ebwtTotSz();
	uint32_t  fchr[] = {0, 0, 0, 0, 0};
	uint32_t* ftab = NULL;
	uint32_t  zOff = 0xffffffff;

	// Save # of occurrences of each character as we walk along the bwt
	uint32_t occ[4] = {0, 0, 0, 0};
	// Save 'G' and 'T' occurrences between backward and forward buckets
	uint32_t occSave[2] = {0, 0};

	// Record rows that should "absorb" adjacent rows in the ftab.
	// The absorbed rows represent suffixes shorter than the ftabChars
	// cutoff.
	uint8_t absorbCnt = 0;
	uint8_t *absorbFtab;
	try {
		VMSG_NL("Allocating ftab, absorbFtab");
		ftab = new uint32_t[ftabLen];
		bzero(ftab, 4 * ftabLen);
		absorbFtab = new uint8_t[ftabLen];
		bzero(absorbFtab, ftabLen);
	} catch(bad_alloc &e) {
		cerr << "Out of memory allocating ftab[] or absorbFtab[] "
		     << "in Ebwt::buildToDisk() at " << __FILE__ << ":"
		     << __LINE__ << endl;
		throw e;
	}
	assert(ftab != NULL);
	assert(absorbFtab != NULL);

	// Allocate the side buffer; holds a single side as its being
	// constructed and then written to disk.  Reused across all sides.
	uint8_t *ebwtSide = NULL;
	try {
		ebwtSide = new uint8_t[sideSz];
	} catch(bad_alloc &e) {
		cerr << "Out of memory allocating ebwtSide[] in "
		     << "Ebwt::buildToDisk() at " << __FILE__ << ":"
		     << __LINE__ << endl;
		throw e;
	}
	assert(ebwtSide != NULL);
	
	// Points to the base offset within ebwt for the side currently
	// being written
	uint32_t side = 0;
	// Points to a byte offset from 'side' within ebwt[] where next
	// char should be written
	int sideCur = eh.sideBwtSz() - 1;

	// Whether we're assembling a forward or a reverse bucket
	bool fw = false;
	
	// Did we just finish writing a forward bucket?  (Must be true when
	// we exit the loop.)
	bool wroteFwBucket = false;
	
	// Have we skipped the '$' in the last column yet?
	ASSERT_ONLY(bool dollarSkipped = false);
	
	uint32_t si = 0;   // string offset (chars)
	ASSERT_ONLY(uint32_t lastSufInt = 0);
	ASSERT_ONLY(bool inSA = true); // true iff saI still points inside suffix
	                               // array (as opposed to the padding at the
	                               // end)
	// Iterate over packed bwt bytes 
	VMSG_NL("Entering Ebwt loop");
	ASSERT_ONLY(uint32_t beforeEbwtOff = (uint32_t)out1.tellp());
	while(side < ebwtTotSz) {
		wroteFwBucket = false;
		// Sanity-check our cursor into the side buffer
		assert_geq(sideCur, 0);
		assert_lt(sideCur, (int)eh.sideBwtSz());
		assert_eq(0, side % sideSz); // 'side' must be on side boundary
		ebwtSide[sideCur] = 0; // clear 
		assert_lt(side + sideCur, ebwtTotSz);
		// Iterate over bit-pairs in the si'th character of the BWT
		for(int bpi = 0; bpi < 4; bpi++, si++) {
			int bwtChar;
			bool count = true;
			if(si <= len) {
				// Still in the SA; extract the bwtChar
				uint32_t saElt = sa.nextSuffix();
				// (that might have triggered sa to calc next suf block)
				if(saElt == 0) {
					// Don't add the '$' in the last column to the BWT
					// transform; we can't encode a $ (only A C T or G)
					// and counting it as, say, an A, will mess up the
					// LR mapping
					bwtChar = 0; count = false;
					ASSERT_ONLY(dollarSkipped = true);
					zOff = si; // remember the SA row that
					           // corresponds to the 0th suffix
				} else {
					bwtChar = (int)(Dna)(s[saElt-1]);
					assert_lt(bwtChar, 4);
					// Update the fchr
					fchr[bwtChar]++;
				}
				// Update ftab
				if((len-saElt) >= (uint32_t)eh.ftabChars()) {
					// Turn the first ftabChars characters of the
					// suffix into an integer index into ftab
					uint32_t sufInt = 0;
					for(int i = 0; i < eh.ftabChars(); i++) {
						sufInt <<= 2;
						assert_lt(i, (int)(len-saElt));
						sufInt |= (unsigned char)(Dna)(s[saElt+i]);
					}
					// Assert that this prefix-of-suffix is greater
					// than or equal to the last one (true b/c the
					// suffix array is sorted) 
					#ifndef NDEBUG
					if(lastSufInt > 0) assert_geq(sufInt, lastSufInt);
					lastSufInt = sufInt;
					#endif
					// Update ftab
					assert_lt(sufInt+1, ftabLen);
					ftab[sufInt+1]++;
					if(absorbCnt > 0) {
						// Absorb all short suffixes since the last
						// transition into this transition
						absorbFtab[sufInt] = absorbCnt;
						absorbCnt = 0;
					}
				} else {
					// Otherwise if suffix is fewer than ftabChars
					// characters long, then add it to the 'absorbCnt';
					// it will be absorbed into the next transition
					assert_lt(absorbCnt, 255);
					absorbCnt++;
				}
				// Offset boundary? - update offset array
				if((si & eh.offMask()) == si) {
					assert_lt((si >> eh.offRate()), eh.offsLen());
					// Write offsets directly to the secondary output
					// stream, thereby avoiding keeping them in memory
					writeU32(out2, saElt, this->toBe());
				}
			} else {
				// Strayed off the end of the SA, now we're just
				// padding out a bucket
				#ifndef NDEBUG
				if(inSA) {
					// Assert that we wrote all the characters in the
					// string before now
					assert_eq(si, len+1);
					inSA = false;
				}
				#endif
				// 'A' used for padding; important that padding be
				// counted in the occ[] array
				bwtChar = 0; 
			}
			if(count) occ[bwtChar]++;
			// Append BWT char to bwt section of current side
			if(fw) {
				// Forward bucket: fill from least to most
				pack_2b_in_8b(bwtChar, ebwtSide[sideCur], bpi);
				assert_eq((ebwtSide[sideCur] >> (bpi*2)) & 3, bwtChar);
			} else {
				// Backward bucket: fill from most to least
				pack_2b_in_8b(bwtChar, ebwtSide[sideCur], 3-bpi);
				assert_eq((ebwtSide[sideCur] >> ((3-bpi)*2)) & 3, bwtChar);
			}
		} // end loop over bit-pairs
		assert_eq(dollarSkipped ? 3 : 0, (occ[0] + occ[1] + occ[2] + occ[3]) & 3);
		assert_eq(0, si & 3);
		if(fw) sideCur++;
		else   sideCur--;
		if(sideCur == (int)eh.sideBwtSz()) {
			// Forward side boundary
			assert_eq(0, si % eh.sideBwtLen());
			sideCur = eh.sideBwtSz() - 1;
			assert(fw); fw = false; wroteFwBucket = true;
			// Write 'G' and 'T'
			assert_leq(occSave[0], occ[2]);
			assert_leq(occSave[1], occ[3]);
			uint32_t *u32side = reinterpret_cast<uint32_t*>(ebwtSide);
			side += sideSz;
			assert_leq(side, eh.ebwtTotSz());
			u32side[(sideSz >> 2)-2] = endianizeU32(occSave[0], this->toBe());
			u32side[(sideSz >> 2)-1] = endianizeU32(occSave[1], this->toBe());
			// Write forward side to primary file
			out1.write((const char *)ebwtSide, sideSz);
		} else if (sideCur == -1) {
			// Backward side boundary
			assert_eq(0, si % eh.sideBwtLen());
			sideCur = 0;
			assert(!fw); fw = true;
			// Write 'A' and 'C'
			uint32_t *u32side = reinterpret_cast<uint32_t*>(ebwtSide);
			side += sideSz;
			assert_leq(side, eh.ebwtTotSz());
			u32side[(sideSz >> 2)-2] = endianizeU32(occ[0], this->toBe());
			u32side[(sideSz >> 2)-1] = endianizeU32(occ[1], this->toBe());
			occSave[0] = occ[2]; // save 'G' count
			occSave[1] = occ[3]; // save 'T' count
			// Write backward side to primary file
			out1.write((const char *)ebwtSide, sideSz);
		}
	}
	VMSG_NL("Exited Ebwt loop");
	assert(ftab != NULL);
	assert_neq(zOff, 0xffffffff);
	if(absorbCnt > 0) {
		// Absorb any trailing, as-yet-unabsorbed short suffixes into
		// the last element of ftab
		absorbFtab[ftabLen-1] = absorbCnt;
	}
	// Assert that our loop counter got incremented right to the end
	assert_eq(side, eh.ebwtTotSz());
	// Assert that we wrote the expected amount to out1
	assert_eq(((uint32_t)out1.tellp() - beforeEbwtOff), eh.ebwtTotSz());
	// assert that the last thing we did was write a forward bucket
	assert(wroteFwBucket);
	
	//
	// Write zOff to primary stream
	//
	writeU32(out1, zOff, this->toBe());
	
	//
	// Finish building fchr
	//
	// Exclusive prefix sum on fchr
	for(int i = 1; i < 4; i++) {
		fchr[i] += fchr[i-1];
	}
	assert_eq(fchr[3], len);
	// Shift everybody up by one
	for(int i = 4; i >= 1; i--) {
		fchr[i] = fchr[i-1]; 
	}
	fchr[0] = 0;
	if(_verbose) {
		for(int i = 0; i < 5; i++)
			cout << "fchr[" << toDna5(i) << "]: " << fchr[i] << endl;
	}
	// Write fchr to primary file
	for(int i = 0; i < 5; i++) {
		writeU32(out1, fchr[i], this->toBe());
	}
	
	//
	// Finish building ftab and build eftab
	//
	// Prefix sum on ftable
	uint32_t eftabLen = 0;
	assert_eq(0, absorbFtab[0]);
	for(uint32_t i = 1; i < ftabLen; i++) {
		if(absorbFtab[i] > 0) eftabLen += 2;
	}
	assert_leq(eftabLen, (uint32_t)eh.ftabChars()*2);
	eftabLen = eh.ftabChars()*2;
	uint32_t *eftab = NULL;
	try {
		eftab = new uint32_t[eftabLen];
		bzero(eftab, 4 * eftabLen);
	} catch(bad_alloc &e) {
		cerr << "Out of memory allocating eftab[] "
		     << "in Ebwt::buildToDisk() at " << __FILE__ << ":"
		     << __LINE__ << endl;
		throw e;
	}
	assert(eftab != NULL);
	uint32_t eftabCur = 0;
	for(uint32_t i = 1; i < ftabLen; i++) {
		uint32_t lo = ftab[i] + Ebwt::ftabHi(ftab, eftab, len, ftabLen, eftabLen, i-1);
		if(absorbFtab[i] > 0) {
			// Skip a number of short pattern indicated by absorbFtab[i]
			uint32_t hi = lo + absorbFtab[i];
			assert_lt(eftabCur*2+1, eftabLen);
			eftab[eftabCur*2] = lo;
			eftab[eftabCur*2+1] = hi;
			ftab[i] = (eftabCur++) ^ 0xffffffff; // insert pointer into eftab
			assert_eq(lo, Ebwt::ftabLo(ftab, eftab, len, ftabLen, eftabLen, i));
			assert_eq(hi, Ebwt::ftabHi(ftab, eftab, len, ftabLen, eftabLen, i));
		} else {
			ftab[i] = lo;
		}
	}
	assert_eq(Ebwt::ftabHi(ftab, eftab, len, ftabLen, eftabLen, ftabLen-1), len+1);
	// Write ftab to primary file
	for(uint32_t i = 0; i < ftabLen; i++) {
		writeU32(out1, ftab[i], this->toBe());
	}
	// Write eftab to primary file
	for(uint32_t i = 0; i < eftabLen; i++) {
		writeU32(out1, eftab[i], this->toBe());
	}
	delete[] ftab;
	delete[] eftab;
	delete[] absorbFtab;

	// Note: if you'd like to sanity-check the Ebwt, you'll have to
	// read it back into memory first!
	assert(!isInMemory());
	VMSG_NL("Exiting Ebwt::buildToDisk()");
}

#endif /*EBWT_H_*/
