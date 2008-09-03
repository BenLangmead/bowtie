#ifndef HIT_H_
#define HIT_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <bitset>
#include <seqan/sequence.h>
#include "assert_helpers.h"

/**
 * Classes for dealing with reporting alignments.
 */

using namespace std;
using namespace seqan;

typedef pair<uint32_t,uint32_t> U32Pair;
// For now, we support reads up to 63 bp long, which is the same as Maq, as 
// of 0.6.7
static const int max_read_bp = 63;

/**
 * Encapsulates a hit, including a text-id/text-offset pair, a pattern
 * id, and a boolean indicating whether it matched as its forward or
 * reverse-complement version.
 */
struct Hit {	
	Hit(U32Pair _h, 
		uint32_t _patId,
		const String<char>& _patName,
		const String<Dna5>& _patSeq,
		const String<char>& _patQualities,
		bool _fw, 
		const bitset<max_read_bp>& _mms,
		uint32_t _oms = 0) : h(_h), 
							 patId(_patId),
							 patName(_patName),
							 patSeq(_patSeq),
							 patQualities(_patQualities),
							 mms(_mms),
							 oms(_oms),
							 fw(_fw) {}
	
	U32Pair  h;
	uint32_t patId;
	String<char> patName;
	String<Dna5> patSeq;
	String<char> patQualities;
	bitset<max_read_bp> mms;
	uint32_t oms;   // # of other possible mappings; 0 -> this is unique
	bool fw;
	Hit& operator = (const Hit &other) {
	    this->h = other.h;
	    this->patId = other.patId;
		this->patName = other.patName;
		this->patSeq = other.patSeq;
		this->patQualities = other.patQualities;
	    this->mms = other.mms;
	    this->oms = other.oms;
		this->fw  = other.fw;
	    return *this;
	}
};

/// Sort by text-id then by text-offset
bool operator< (const Hit& a, const Hit& b);

/**
 * Encapsulates an object that accepts hits, optionally retains them in
 * a vector, and does something else with them according to
 * descendent's implementation of pure virtual member reportHitImpl().
 */
class HitSink {
public:
	HitSink(ostream& __out = cout,
	        bool __keep = false,
	        vector<string>* __refnames = NULL) :
		_out(__out),
		_hits(),
		_provisionalHits(),
		_keep(__keep),
		_numHits(0llu),
		_refnames(__refnames) { }
	virtual ~HitSink() { }
	virtual void reportHit(const U32Pair& h,
						   uint32_t patId,
						   const String<char>& patName,
						   const String<Dna5>& patSeq,
						   const String<char>& patQualities,
						   bool fw,
						   const bitset<max_read_bp>& mms,
						   uint32_t oms) = 0;
	/**
	 * A provisional hit is a hit that we might want to report, but we
	 * aren't sure yet.
	 */
	virtual void reportProvisionalHit(
			const U32Pair& h,
			uint32_t patId,
		    const String<char>& patName,
			const String<Dna5>& patSeq,
			const String<char>& patQualities,
            bool fw,
            const bitset<max_read_bp>& mms,
            uint32_t oms) = 0;
	/// Report and purge all accumulated provisional hits
	virtual void acceptProvisionalHits() {
		// Save the old _keep and set it to false while we report the
		// provisional hits; this is because we already retained them
		// when they were reported provisionally.
		bool keep = _keep;
		_keep = false;
		for(size_t i = 0; i < _provisionalHits.size(); i++) {
			const Hit& h = _provisionalHits[i];
			reportHit(h.h, h.patId, h.patName, h.patSeq, h.patQualities, h.fw, h.mms, h.oms);
		}
		_keep = keep; // restore _keep
		_provisionalHits.clear();
	}
	/// Purge all accumulated provisional hits without reporting them
	virtual void rejectProvisionalHits() {
		_provisionalHits.clear();
	}
	virtual void finishImpl() { }
	virtual void finish()               { finishImpl(); }
	virtual void flush()                { _out.flush(); }
	virtual void setRetainHits(bool r)  { _keep = r; }
	virtual bool retainHits()           { return _keep; }
	virtual void clearRetainedHits()    { _hits.clear(); }
	virtual vector<Hit>& retainedHits() { return _hits; }
	virtual ostream& out()              { return _out; }
	virtual uint64_t numHits()          { return _numHits; } 
	virtual size_t numProvisionalHits() { return _provisionalHits.size(); } 
protected:
	ostream& _out;
	vector<Hit> _hits;
	vector<Hit> _provisionalHits;
	bool _keep;
	uint64_t _numHits;
	vector<string>* _refnames;
};

class HitBucket : public HitSink
{
public:
	HitBucket() : HitSink(cout, true) { }
	virtual void reportHit(const U32Pair& h,
						   uint32_t patId,
						   const String<char>& patName,
						   const String<Dna5>& patSeq,
						   const String<char>& patQualities,
						   bool fw,
						   const bitset<max_read_bp>& mms,
						   uint32_t oms) 
	{
		_hits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
		_numHits++;
	}
	
	virtual void reportProvisionalHit(const U32Pair& h,
									  uint32_t patId,
									  const String<char>& patName,
									  const String<Dna5>& patSeq,
									  const String<char>& patQualities,
									  bool fw,
									  const bitset<max_read_bp>& mms,
									  uint32_t oms)
	{
		_hits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
		_provisionalHits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
	}
};

/**
 * Sink that prints lines like this:
 * (pat-id)[-|+]:<hit1-text-id,hit2-text-offset>,<hit2-text-id...
 */
class PrettyHitSink : public HitSink {
public:
	PrettyHitSink(
			ostream& __out,
			bool __revcomp = false,
			bool __reportOpps = false,
			bool __keep = false,
			vector<string>* __refnames = NULL) :
		HitSink(__out, __keep, __refnames),
		_revcomp(__revcomp),
		_reportOpps(__reportOpps),
		_lastPat(0xffffffff),
		_lastFw(false),
		_first(true) { }

	virtual void reportHit(
			const U32Pair& h,
			uint32_t patId,
			const String<char>& patName,
			const String<Dna5>& patSeq,
		    const String<char>& patQualities,
			bool fw,
			const bitset<max_read_bp>& mms,
			uint32_t oms)
	{
		assert(!_revcomp || (patId & 1) == 0 || !fw);
		assert(!_revcomp || (patId & 1) == 1 || fw);
		if(_revcomp) patId >>= 1;
		if(patId != _lastPat || fw != _lastFw) {
			// First hit on a new line
			_lastPat = patId;
			_lastFw = fw;
			if(_first) _first = false;
			else       out() << endl;
			out() << patId << (fw? "+":"-") << ":";
		} else {
			// Not the first hit on the line
			out() << ",";
		}
		assert(!_first);
    	// .first is text id, .second is offset
		out() << "<" << h.first << "," << h.second << "," << mms.count();
		if(_reportOpps) out() << "," << oms;
		out() << ">";
		if(_keep) {
			_hits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
		}
		_numHits++;
	}
	/**
	 * A provisional hit is a hit that we might want to report, but we
	 * aren't sure yet.
	 */
	virtual void reportProvisionalHit(
			const U32Pair& h,
			uint32_t patId,
			const String<char>& patName,
			const String<Dna5>& patSeq,
			const String<char>& patQualities,
			bool fw,
			const bitset<max_read_bp>& mms,
			uint32_t oms)
	{
		if(_keep) {
			_hits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
		}
		_provisionalHits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
	}
	virtual void finishImpl() {
		if(_first) {
			out() << "No results" << endl;
		} else {
			out() << endl;
		}
	}

private:
	bool _revcomp;
	bool _reportOpps;
	uint32_t _lastPat;
	bool _lastFw;
	bool _first;
};

/**
 * Sink that prints lines like this:
 * (pat-id)[-|+]:<hit1-text-id,hit2-text-offset>,<hit2-text-id...
 */
class VerboseHitSink : public HitSink {
public:
	VerboseHitSink(ostream& __out,
				   bool __revcomp = false,
				   bool __keep = false,
				   vector<string>* __refnames = NULL) :
	HitSink(__out, __keep, __refnames),
	_revcomp(__revcomp),
	_lastPat(0xffffffff),
	_lastFw(false),
	_first(true) { }
	
	virtual void reportHit(const U32Pair& h,
						   uint32_t patId,
						   const String<char>& patName,
						   const String<Dna5>& patSeq,
						   const String<char>& patQualities,
						   bool fw,
						   const bitset<max_read_bp>& mms,
						   uint32_t oms)
	{
		assert(!_revcomp || (patId & 1) == 0 || !fw);
		assert(!_revcomp || (patId & 1) == 1 || fw);
		if(_revcomp) patId >>= 1;
		
		_first = false;
		
		out() << patName <<" \t" << (fw? "+":"-") << "\t";
		
    	// .first is text id, .second is offset
		
		if(this->_refnames != NULL && h.first < this->_refnames->size()) {
			out() << (*this->_refnames)[h.first];
		} else {
			out() << h.first;
		}
		out() << "\t" << h.second;
		out() << "\t" << patSeq;
		
		out() << "\t" << patQualities;
		
		out() << "\t" << oms;
		out() << "\t";
		
		bool firstmiss = true;
		for (unsigned int i = 0; i < mms.size(); ++ i)
		{
			if (mms.test(i))
			{
				if (!firstmiss)
					out() <<",";
				out() << i;
				firstmiss = false;
			}
		}
		
		out () << endl;
		if(_keep) {
			_hits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
		}
		_numHits++;
	}
	/**
	 * A provisional hit is a hit that we might want to report, but we
	 * aren't sure yet.
	 */
	virtual void reportProvisionalHit(const U32Pair& h,
									  uint32_t patId,
									  const String<char>& patName,
									  const String<Dna5>& patSeq,
									  const String<char>& patQualities,
									  bool fw,
									  const bitset<max_read_bp>& mms,
									  uint32_t oms)
	{
		if(_keep) {
			_hits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
		}
		_provisionalHits.push_back(Hit(h, patId, patName, patSeq, patQualities, fw, mms, oms));
	}
	virtual void finishImpl() {
		if(_first) {
			out() << "No results" << endl;
		} else {
			//out() << endl;
		}
	}
	
private:
	bool _revcomp;
	bool _reportOpps;
	uint32_t _lastPat;
	bool _lastFw;
	bool _first;
};


/**
 * Sink that writes 16-byte binary records for each hit.  It also
 * buffers those records in 4K _buf.
 */
class BufferedBinaryHitSink : public HitSink {
public:
	BufferedBinaryHitSink(
			ostream& __out,
			bool __revcomp = false,
			bool __reportOpps = false,
			bool __keep = false,
			vector<string>* __refnames = NULL) :
		HitSink(__out, __keep, __refnames),
		_revcomp(__revcomp),
		_reportOpps(__reportOpps),
		_cur(0) { }
	
	virtual void reportHit(
			const U32Pair& h,
			uint32_t patId,
			const String<char>& patName,
			const String<Dna5>& patSeq,
		    const String<char>& patQualities,
			bool fw,
			const bitset<max_read_bp>& mms,
			uint32_t oms)
	{
		//FIXME:
/*
		assert((pat & 1) == 0 || !fw);
		assert((pat & 1) == 1 || fw);
		if(_revcomp) patId >>= 1;
		assert_eq(8, sizeof(U32Pair));
		assert_eq(4, sizeof(uint32_t));
		*((U32Pair*)&_buf[_cur]) = h;
		_cur += 8;
		*((uint32_t*)&_buf[_cur]) = patId; // pattern 
		_cur += 4;
		*((uint32_t*)&_buf[_cur]) = fw? 1 : 0; // orientation
		_cur += 4;
		*((uint32_t*)&_buf[_cur]) = mms; // # mismatches
		_cur += 4;
		if(_reportOpps) {
			*((uint32_t*)&_buf[_cur]) = oms; // # other mappings
			_cur += 4;
		}
		assert_leq(_cur, BUFSZ);
		if(_cur + 20 >= BUFSZ) {
			out().write((const char *)_buf, _cur);
			assert(!out().bad());
			_cur = 0;
		}
		if(_keep) {
			_hits.push_back(Hit(h, patId, fw, mms, oms));
		}
		_numHits++;
 */
	}
	/**
	 * A provisional hit is a hit that we might want to report, but we
	 * aren't sure yet.
	 */
	virtual void reportProvisionalHit(
			const U32Pair& h,
			uint32_t patId,
			const String<char>& patName,
			const String<Dna5>& patSeq,
			const String<char>& patQualities,
			bool fw,
			const bitset<max_read_bp>& mms,
			uint32_t oms)
	{
		/*
		if(_keep) {
			_hits.push_back(Hit(h, patId, fw, mms, oms));
		}
		_provisionalHits.push_back(Hit(h, patId, fw, mms, oms));
		 */
	}
	virtual void finish() {
		finishImpl();
		char zeros[] = {0, 0, 0, 0, 0, 0, 0, 0};
		_out.write((const char *)zeros, 8);
		// Write a bunch of newlines to make it easier to 'tail' the
		// hits file
		for(int i = 0; i < 256; i++) _out << endl;
		_out.flush();
	}
	virtual void finishImpl() {
		if(_cur > 0) {
			out().write((const char *)_buf, _cur);
			assert(!out().bad());
		}
	}

private:
	static const int BUFSZ = 4096;
	bool _revcomp;
	bool _reportOpps;
	int _cur;
	char _buf[BUFSZ];
};


#endif /*HIT_H_*/

