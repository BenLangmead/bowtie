#ifndef HIT_H_
#define HIT_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <bitset>
#include <sstream>
#include <seqan/sequence.h>
#include "assert_helpers.h"
#include "spinlock.h"
#include "threading.h"

/**
 * Classes for dealing with reporting alignments.
 */

using namespace std;
using namespace seqan;

typedef pair<uint32_t,uint32_t> U32Pair;
// For now, we support reads up to 63 bp long, which is the same as Maq, as
// of 0.6.7
static const int max_read_bp = 1023;

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
		const String<char>& _quals,
		bool _fw,
		const bitset<max_read_bp>& _mms,
		uint32_t _oms = 0) :
		h(_h),
		patId(_patId),
		patName(_patName),
		patSeq(_patSeq),
		quals(_quals),
		mms(_mms),
		oms(_oms),
		fw(_fw) {}

	U32Pair             h;       /// reference index & offset
	uint32_t            patId;   /// read index
	String<char>        patName; /// read name
	String<Dna5>        patSeq;  /// read sequence
	String<char>        quals;   /// read qualities
	bitset<max_read_bp> mms;     /// mismatch mask
	uint32_t            oms;     /// # of other possible mappings; 0 -> this is unique
	bool                fw;      /// orientation of read in alignment

	Hit& operator = (const Hit &other) {
	    this->h       = other.h;
	    this->patId   = other.patId;
		this->patName = other.patName;
		this->patSeq  = other.patSeq;
		this->quals   = other.quals;
	    this->mms     = other.mms;
	    this->oms     = other.oms;
		this->fw      = other.fw;
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
	HitSink(ostream&        __out = cout,
	        vector<string>* __refnames = NULL) :
		_out(__out),
		_refnames(__refnames),
		_lock()
	{
#ifdef USE_SPINLOCK
		// No initialization
#else
   		MUTEX_INIT(_lock);
#endif
	}
	virtual ~HitSink() { }

	/// Implementation of hit-report
	virtual void reportHit(const Hit& h) = 0;

	/// Called when all alignments are complete
	virtual void finish()               { }
	/// Flushes the alignment output stream
	virtual void flush()                { _out.flush(); }
	/// Returns the alignment output stream
	virtual ostream& out()              { return _out; }
protected:
	void mylock() {
#ifdef USE_SPINLOCK
		_lock.Enter();
#else
		MUTEX_LOCK(_lock);
#endif
	}
	void myunlock() {
#ifdef USE_SPINLOCK
		_lock.Leave();
#else
		MUTEX_UNLOCK(_lock);
#endif
	}
	ostream&        _out;      /// the alignment output stream
	vector<string>* _refnames; /// map from reference indexes to names
#ifdef USE_SPINLOCK
	SpinLock _lock;
#else
	MUTEX_T _lock;     /// mutex for locking critical regions
#endif
};

/**
 * A per-thread wrapper for a HitSink.  Incorporates state that a
 * single search thread cares about, i.e., this thread's provisional
 * hits and reported-hit count.
 */
class HitSinkPerThread {
public:
	HitSinkPerThread(HitSink& sink, bool __keep = false) :
		_sink(sink),
		_keep(__keep),
		_numHits(0llu),
		_hits(),
		_provHits() { }

	/// Purge all accumulated provisional hits without reporting them
	void rejectProvisionalHits() {
		_provHits.clear();
	}

	/// Report and purge all accumulated provisional hits
	void acceptProvisionalHits() {
		// Save the old _keep and set it to false while we report the
		// provisional hits; this is because we already retained them
		// when they were reported provisionally.
		bool keep = _keep;
		_keep = false;
		for(size_t i = 0; i < _provHits.size(); i++) {
			const Hit& h = _provHits[i];
			reportHit(h);
		}
		_keep = keep; // restore _keep
		_provHits.clear();
	}

	/// Set whether to retain hits in a vector or not
	void setRetainHits(bool r)  { _keep = r; }
	/// Return whether we're retaining hits or not
	bool retainHits()           { return _keep; }
	/// Clear all hits in the retained-hits vector
	void clearRetainedHits()    { _hits.clear(); }
	/// Return the vector of retained hits
	vector<Hit>& retainedHits() { return _hits; }

	/**
	 * Make note of a hit but don't pass it on to the HitSink yet
	 * because it may be canceled in the future.
	 */
	void reportProvisionalHit(const Hit& h)	{
		if(_keep) _hits.push_back(h);
		_provHits.push_back(h);
	}

	/**
	 * Implementation for hit reporting; update per-thread _hits and
	 * _numHits variables and call the master HitSink to do the actual
	 * reporting
	 */
	void reportHit(const Hit& h) {
		_sink.reportHit(h);
		if(_keep) _hits.push_back(h);
		_numHits++;
	}

	/// Return the number of hits reported so far
	uint64_t numHits()          { return _numHits; }
	/// Return the number of provisional hits stored as of now
	size_t numProvisionalHits() { return _provHits.size(); }

private:
	HitSink&    _sink;
	bool        _keep;
	uint64_t    _numHits;
	vector<Hit> _hits;
	vector<Hit> _provHits;
};

/**
 * Sink that prints lines like this:
 * (pat-id)[-|+]:<hit1-text-id,hit2-text-offset>,<hit2-text-id...
 */
class ConciseHitSink : public HitSink {
public:
	ConciseHitSink(
			ostream&        __out,
			bool            __reportOpps = false,
			vector<string>* __refnames = NULL) :
		HitSink(__out, __refnames),
		_reportOpps(__reportOpps),
		_lastPat(0xffffffff),
		_lastFw(false),
		_first(true) { }

	/**
	 * TODO: This output method is inherently non-threadsafe if we're
	 * ever reporting more than one result per line.  Need more
	 * buffering logic in that case.
	 */
	virtual void reportHit(const Hit& h) {
		ostringstream ss;
		if(h.patId != _lastPat || h.fw != _lastFw) {
			// First hit on a new line
			_lastPat = h.patId;
			_lastFw  = h.fw;
			if(_first) _first = false;
			else ss << endl;
			ss << h.patId << (h.fw? "+" : "-") << ":";
		} else {
			// Not the first hit on the line
			ss << ",";
		}
		assert(!_first);
    	// .first is text id, .second is offset
		ss << "<" << h.h.first << "," << h.h.second << "," << h.mms.count();
		if(_reportOpps) ss << "," << h.oms;
		ss << ">";
		mylock();
		out() << ss.str();
		myunlock();
	}

	virtual void finish() {
		if(_first) out() << "No results";
		out() << endl;
	}

private:
	bool     _reportOpps;
	uint32_t _lastPat;
	bool     _lastFw;
	bool     _first; /// true -> first hit hasn't yet been reported
};

/**
 * Sink that prints lines like this:
 * pat-name \t [-|+] \t ref-name \t ref-off \t pat \t qual \t #-alt-hits \t mm-list
 */
class VerboseHitSink : public HitSink {
public:
	VerboseHitSink(ostream&        __out,
				   vector<string>* __refnames = NULL) :
	HitSink(__out, __refnames),
	_first(true) { }

	virtual void reportHit(const Hit& h) {
		_first = false;
		ostringstream ss;
		ss << h.patName << "\t" << (h.fw? "+":"-") << "\t";
    	// .first is text id, .second is offset
		if(this->_refnames != NULL && h.h.first < this->_refnames->size()) {
			ss << (*this->_refnames)[h.h.first];
		} else {
			ss << h.h.first;
		}
		ss << "\t" << h.h.second;
		ss << "\t" << h.patSeq;
		ss << "\t" << h.quals;
		ss << "\t" << h.oms;
		ss << "\t";
		bool firstmiss = true;
		for (unsigned int i = 0; i < h.mms.size(); ++ i) {
			if (h.mms.test(i)) {
				if (!firstmiss) ss << ",";
				ss << i;
				firstmiss = false;
			}
		}
		ss << endl;
		mylock();
		out() << ss.str();
		myunlock();
	}

	virtual void finish() {
		if(_first) out() << "No results" << endl;
	}

private:
	bool _first; /// true -> first hit hasn't yet been reported
};

#endif /*HIT_H_*/

