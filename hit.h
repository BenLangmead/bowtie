#ifndef HIT_H_
#define HIT_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "spinlock.h"
#include "threading.h"
#include "bitset.h"

/**
 * Classes for dealing with reporting alignments.
 */

using namespace std;
using namespace seqan;

typedef enum output_types {
	FULL = 1,
	CONCISE,
	BINARY,
	NONE
};

static const std::string output_type_names[] = {
	"Invalid!",
	"Full",
	"Concise",
	"Binary",
	"None"
};

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
		const FixedBitset<max_read_bp>& _mms,
		const vector<char>& _refcs,
		uint32_t _oms = 0) :
		h(_h),
		patId(_patId),
		patName(_patName),
		patSeq(_patSeq),
		quals(_quals),
		mms(_mms),
		refcs(_refcs),
		oms(_oms),
		fw(_fw)
	{
		// Enforce the constraints imposed by the binary output format
		if(length(patName) > 0xffff) {
			cerr << "Error: One or more read names are 2^16 characters or longer.  Please truncate" << endl
			     << "read names and re-run bowtie." << endl;
			exit(1);
		}
		if(mms.count() > 0xff) {
			cerr << "Error: The alignment contains 256 or more mismatches.  bowtie cannot handle" << endl
			     << "alignments with this many alignments.  Please provide smaller reads or consider" << endl
			     << "using a different tool." << endl;
			exit(1);
		}
		if(length(quals) > 0xffff) {
			cerr << "Error: One or more quality strings are 2^16 characters or longer.  Please" << endl
			     << "truncate reads and re-run bowtie." << endl;
			exit(1);
		}
		if(length(patSeq) > 0xffff) {
			cerr << "Error: One or more read sequences are 2^16 characters or longer.  Please" << endl
			     << "truncate reads and re-run bowtie." << endl;
			exit(1);
		}
	}

	U32Pair             h;       /// reference index & offset
	uint32_t            patId;   /// read index
	String<char>        patName; /// read name
	String<Dna5>        patSeq;  /// read sequence
	String<char>        quals;   /// read qualities
	FixedBitset<max_read_bp> mms;     /// mismatch mask
	vector<char>        refcs;   /// reference characters for mms
	uint32_t            oms;     /// # of other possible mappings; 0 -> this is unique
	bool                fw;      /// orientation of read in alignment

	Hit& operator = (const Hit &other) {
	    this->h       = other.h;
	    this->patId   = other.patId;
		this->patName = other.patName;
		this->patSeq  = other.patSeq;
		this->quals   = other.quals;
	    this->mms     = other.mms;
	    this->refcs   = other.refcs;
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
	void lock() {
#ifdef USE_SPINLOCK
		_lock.Enter();
#else
		MUTEX_LOCK(_lock);
#endif
	}
	void unlock() {
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
 * single search thread cares about.
 */
class HitSinkPerThread {
public:
	HitSinkPerThread(HitSink& sink, bool __keep = false) :
		_sink(sink),
		_bestRemainingStratum(0),
		_numHits(0llu),
		_keep(__keep),
		_hits() { }
	
	virtual ~HitSinkPerThread() { }

	/// Set whether to retain hits in a vector or not
	void setRetainHits(bool r)  { _keep = r; }
	/// Return whether we're retaining hits or not
	bool retainHits()           { return _keep; }
	/// Clear all hits in the retained-hits vector
	void clearRetainedHits()    { _hits.clear(); }
	/// Return the vector of retained hits
	vector<Hit>& retainedHits() { return _hits; }
	
	/// Finalize current read
	virtual void finishRead() = 0;

	/**
	 * Implementation for hit reporting; update per-thread _hits and
	 * _numHits variables and call the master HitSink to do the actual
	 * reporting
	 */
	bool reportHit(const Hit& h) {
		if(_keep) _hits.push_back(h);
		_numHits++;
		return reportHitImpl(h);
	}
	
	/**
	 * Concrete subclasses override this to (possibly) report a hit and
	 * return true iff the caller should continue to report more hits. 
	 */
	virtual bool reportHitImpl(const Hit& h) = 0;

	/// Return the number of hits reported so far
	uint64_t numHits()          { return _numHits; }
	
	/**
	 * The search routine is informing us that it will not be reporting
	 * any more hits at the given stratum.
	 */
	void finishedWithStratum(int stratum) {
		assert_eq(_bestRemainingStratum, stratum);
		_bestRemainingStratum++;
		finishedWithStratumImpl(stratum);
	}
	
	/**
	 * Concrete subclasses override this to determine whether the
	 * search routine should keep searching after having finished
	 * reporting all alignments at the given stratum.
	 */
	virtual bool finishedWithStratumImpl(int stratum) = 0;
	
protected:
	HitSink&    _sink; /// Ultimate destination of reported hits
	/// Least # mismatches in alignments that will be reported in the
	/// future.  Updated by the search routine.
	int         _bestRemainingStratum;
	/// # hits reported to this HitSink so far (not all of which were
	/// necesssary reported to _sink)
	uint64_t    _numHits;
private:
	bool        _keep; /// Whether to retain all reported hits in _hits
	vector<Hit> _hits; /// Repository for retained hits
};

/**
 * Report first N good alignments encountered; trust search routine
 * to try alignments in something approximating a best-first order.
 * Best used in combination with a stringent alignment policy.
 */
class FirstNGoodHitSinkPerThread : public HitSinkPerThread {

public:
	FirstNGoodHitSinkPerThread(
			HitSink& sink,
	        int __n,
	        bool __keep = false) : 
	        HitSinkPerThread(sink, __keep),
	        _hitsForThisRead(0),
	        _n(__n)
	{
		assert_gt(_n, 0);
	}

	/// Finalize current read
	virtual void finishRead() {
		_hitsForThisRead = 0;
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHitImpl(const Hit& h) {
		_sink.reportHit(h);
		_hitsForThisRead++;
		assert_leq(_hitsForThisRead, _n);
		if(_hitsForThisRead == _n) {
			return false; // already reported N good hits; stop!
		}
		return true; // not at N yet; keep going
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits. 
	 */
	virtual bool finishedWithStratumImpl(int stratum) { return true; }

private:
	int _hitsForThisRead; /// # hits for this read so far
	int _n;               /// max # hits to report per read
};

/**
 * Report the first N best alignments encountered.  Aggregator and
 * search routine collaborate to ensure that there exists no
 * alignment that's better than the ones reported.
 */
class FirstNBestHitSinkPerThread : public HitSinkPerThread {

public:
	FirstNBestHitSinkPerThread(
			HitSink& sink,
	        int __n,
	        bool __keep = false) : 
	        HitSinkPerThread(sink, __keep),
	        _n(__n)
	{
		assert_gt(_n, 0);
	}

private:
	int _n; // max # hits to report
};

/**
 * Report the first N best alignments encountered in a single
 * alignment stratum (i.e., the best stratum in which there were
 * any hits).  Aggregator and search routine collaborate to ensure
 * that there exists no alignment that's better than the ones
 * reported.
 */
class FirstNBestStratifiedHitSinkPerThread : public HitSinkPerThread {

public:
	FirstNBestStratifiedHitSinkPerThread(
			HitSink& sink,
	        int __n,
	        bool __keep = false) : 
	        HitSinkPerThread(sink, __keep),
	        _n(__n)
	{
		assert_gt(_n, 0);
	}

private:
	int _n; // max # hits to report
};

/**
 * Report all valid alignments.
 */
class AllHitSinkPerThread : public HitSinkPerThread {
public:
	AllHitSinkPerThread(
			HitSink& sink,
	        bool __keep = false) : 
	        HitSinkPerThread(sink, __keep) { }
};

/**
 * Report all alignments encountered in the best alignment stratum
 * for which there were any alignments.
 */
class AllStratifiedHitSinkPerThread : public HitSinkPerThread {
public:
	AllStratifiedHitSinkPerThread(
			HitSink& sink,
	        bool __keep = false) : 
	        HitSinkPerThread(sink, __keep) { }
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
			ss << endl;
			ss << h.patId << (h.fw? "+" : "-") << ":";
		} else {
			// Not the first hit on the line
			ss << ",";
		}
    	// .first is text id, .second is offset
		ss << "<" << h.h.first << "," << h.h.second << "," << h.mms.count();
		if(_reportOpps) ss << "," << h.oms;
		ss << ">";
		lock();
		_first = false;
		out() << ss.str();
		unlock();
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
		// Output mismatch column
		bool firstmiss = true;
		size_t c = 0;
		for (unsigned int i = 0; i < h.mms.size(); ++ i) {
			if (h.mms.test(i)) {
				if (!firstmiss) ss << ",";
				ss << i;
				if(h.refcs.size() > 0) {
					assert_gt(h.refcs.size(), i);
					ASSERT_ONLY(char cc = toupper(h.refcs[i]));
					assert(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T');
					char refChar = toupper(h.refcs[i]);
					char qryChar = (h.fw ? h.patSeq[i] : h.patSeq[length(h.patSeq)-i-1]);
					assert_neq(refChar, qryChar);
					ss << ":" << refChar << ">" << qryChar;
				}
				firstmiss = false;
				c++;
			}
		}
		ss << endl;
		// Make sure to grab lock before writing to output stream
		lock();
		_first = false;
		out() << ss.str();
		unlock();
	}

	virtual void finish() {
		if(_first) out() << "No results" << endl;
	}

private:
	bool _first; /// true -> first hit hasn't yet been reported
};

/**
 * Sink for binary output:
 */
class BinaryHitSink : public HitSink {
public:
	BinaryHitSink(ostream&        __out,
				  vector<string>* __refnames = NULL) :
	HitSink(__out, __refnames) { }

	virtual void reportHit(const Hit& h) {
		lock();
		// true iff we're going to print the reference sequence name
		bool refName = this->_refnames != NULL &&
		                h.h.first < this->_refnames->size();
		uint16_t pnamelen = (uint16_t)length(h.patName);
		// Write read name
		_out.write((const char *)&pnamelen, 2);
		_out.write(begin(h.patName), pnamelen);
		// Write fw/refname flags
		uint8_t flags = (h.fw ? 1 : 0) | (refName? 2 : 0);
		_out.write((const char *)&flags, 1);
		if(refName) {
			// Write reference name as string
			uint16_t rnamelen = (uint16_t)(*this->_refnames)[h.h.first].length();
			_out.write((const char *)&rnamelen, 2);
			_out.write((*this->_refnames)[h.h.first].c_str(), rnamelen);
		} else {
			// Write reference name as index into global reference name
			// list
			_out.write((const char *)&h.h.first, 4);
		}
		// Write reference offset
		_out.write((const char *)&h.h.second, 4);
		// Write pattern sequence
		uint16_t plen = (uint16_t)length(h.patSeq);
		for(size_t i = 0; i < plen; i += 2) {
			uint8_t twoChrs = (uint8_t)h.patSeq[i];
			if(i+1 < plen) {
				twoChrs |= ((uint8_t)h.patSeq[i+1] << 4);
			}
			_out.write((const char *)&twoChrs, 1);
		}
		// Write quals sequence
		uint16_t qlen = (uint16_t)length(h.quals);
		_out.write(begin(h.quals), qlen);
		// Write oms
		_out.write((const char *)&h.oms, 4);
		// Write # mismatches
		uint8_t numMms = h.mms.count();
		_out.write((const char *)&numMms, 1);
		// Output mismatches
		size_t c = 0;
		for (uint8_t i = 0; i < h.mms.size(); ++ i) {
			if (h.mms.test(i)) {
				_out.write((const char *)&i, 1);
				assert_gt(h.refcs.size(), i);
				assert_eq(1, dna4Cat[(int)h.refcs[i]]);
				uint8_t refChar = charToDna5[(int)h.refcs[i]];
				assert_leq(refChar, 4);
				uint8_t qryChar = (h.fw ? (int)h.patSeq[i] :
				                          (int)h.patSeq[length(h.patSeq)-i-1]);
				assert_leq(refChar, 4);
				assert_neq(refChar, qryChar);
				uint8_t both = refChar | (qryChar << 4);
				_out.write((const char *)&both, 1);
				c++;
			}
		}
		unlock();
	}

	virtual void finish() { }
};

/**
 * Sink that does nothing:
 */
class StubHitSink : public HitSink {
public:
	StubHitSink() : HitSink(cout, NULL) { }
	virtual void reportHit(const Hit& h) { }
	virtual void finish() {	}
};

#endif /*HIT_H_*/

