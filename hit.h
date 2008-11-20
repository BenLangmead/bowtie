#ifndef HIT_H_
#define HIT_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <seqan/sequence.h>
#include "alphabet.h"
#include "assert_helpers.h"
#include "spinlock.h"
#include "threading.h"
#include "bitset.h"
#include "tokenize.h"
#include "pat.h"

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
	Hit() :
		h(make_pair(0, 0)),
		patId(0),
		patName(String<char>("")),
		patSeq(String<char>("")),
		quals(String<char>("")),
		mms(),
		refcs(),
		oms(0),
		fw(true) { }

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
		if(seqan::length(patName) > 0xffff) {
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
		if(seqan::length(quals) > 0xffff) {
			cerr << "Error: One or more quality strings are 2^16 characters or longer.  Please" << endl
			     << "truncate reads and re-run bowtie." << endl;
			exit(1);
		}
		if(seqan::length(patSeq) > 0xffff) {
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
	FixedBitset<max_read_bp> mms;/// mismatch mask
	vector<char>        refcs;   /// reference characters for mms
	uint32_t            oms;     /// # of other possible mappings; 0 -> this is unique
	bool                fw;      /// orientation of read in alignment

	size_t length() const { return seqan::length(patSeq); }

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
		_outs(),
		_deleteOuts(false),
		_refnames(__refnames),
		_locks()
	{
	    _outs.push_back(&__out);
		_locks.resize(1);
#ifdef USE_SPINLOCK
		// No initialization
#else
   		MUTEX_INIT(_locks[0]);
#endif
	}

	/**
	 * Open a number of output streams; usually one per reference
	 * sequence.  For now, we give then names refXXXXX.map where XXXXX
	 * is the 0-padded reference index.  Someday we may want to include
	 * the name of the reference sequence in the filename somehow.
	 */
	HitSink(size_t numOuts, vector<string>* __refnames = NULL) :
		_outs(),
		_deleteOuts(true),
		_refnames(__refnames),
		_locks()
	{
		// Open all files for writing and initialize all locks
		for(size_t i = 0; i < numOuts; i++) {
			_outs.push_back(NULL); // we open output streams lazily
			_locks.resize(i+1);
#ifdef USE_SPINLOCK
			// No initialization
#else
			MUTEX_INIT(_locks[i]);
#endif
		}
	}

	virtual ~HitSink() {
		if(_deleteOuts) {
			for(size_t i = 0; i < _outs.size(); i++) {
				if(_outs[i] != NULL) delete _outs[i];
			}
		}
	}

	/**
	 * Called by concrete subclasses to figure out which elements of
	 * the _outs/_locks array to use when outputting the alignment.
	 */
	size_t refIdxToStreamIdx(size_t refIdx) {
		if(refIdx >= _outs.size()) return 0;
		return refIdx;
	}

	/// Implementation of hit-report
	virtual void reportHit(const Hit& h) = 0;

	virtual void reportHits(vector<Hit>& hs) {
		for(size_t i = 0; i < hs.size(); i++) {
#ifndef NDEBUG
			for(size_t j = i+1; j < hs.size(); j++) {
				assert_eq(hs[i].h.first, hs[j].h.first);
			}
#endif
			reportHit(hs[i]);
		}
	}

	/// Called when all alignments are complete
	virtual void finish()               { }
	/// Flushes the alignment output stream
	virtual void flush() {
		for(size_t i = 0; i < _outs.size(); i++) {
			if(_outs[i] != NULL) _outs[i]->flush();
		}
	}
	/// Returns the alignment output stream; if the stream needs to be
	/// created, create it
	virtual ostream& out(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
		if(_outs[strIdx] == NULL) {
			assert(_deleteOuts);
			ostringstream oss;
			oss << "ref";
			if     (strIdx < 10)    oss << "0000";
			else if(strIdx < 100)   oss << "000";
			else if(strIdx < 1000)  oss << "00";
			else if(strIdx < 10000) oss << "0";
			oss << strIdx << ".map";
			_outs[strIdx] = new ofstream(oss.str().c_str());
		}
		assert(_outs[strIdx] != NULL);
		return *(_outs[strIdx]);
	}

protected:
	void lock(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
#ifdef USE_SPINLOCK
		_locks[strIdx].Enter();
#else
		MUTEX_LOCK(_locks[strIdx]);
#endif
	}
	void mainlock() {
#ifdef USE_SPINLOCK
		_mainlock.Enter();
#else
		MUTEX_LOCK(_mainlock);
#endif
	}
	void unlock(size_t refIdx) {
		size_t strIdx = refIdxToStreamIdx(refIdx);
#ifdef USE_SPINLOCK
		_locks[strIdx].Leave();
#else
		MUTEX_UNLOCK(_locks[strIdx]);
#endif
	}
	void mainunlock() {
#ifdef USE_SPINLOCK
		_mainlock.Leave();
#else
		MUTEX_UNLOCK(_mainlock);
#endif
	}
	vector<ostream*> _outs;     /// the alignment output stream(s)
	bool             _deleteOuts; /// Whether to delete elements of _outs upon exit
	vector<string>*  _refnames; /// map from reference indexes to names
#ifdef USE_SPINLOCK
	vector<SpinLock> _locks;    /// spinlocks for per-file critical sections
	SpinLock         _mainlock; /// spinlocks for fields of this object
#else
	vector<MUTEX_T>  _locks;    /// pthreads mutexes for per-file critical sections
	MUTEX_T          _mainlock; /// pthreads mutexes for fields of this object
#endif
};

/**
 * A per-thread wrapper for a HitSink.  Incorporates state that a
 * single search thread cares about.
 */
class HitSinkPerThread {
public:
	HitSinkPerThread(HitSink& sink, uint32_t __max = 0xffffffff, bool __keep = false) :
		_sink(sink),
		_bestRemainingStratum(0),
		_bufferHits(true),
		_numReportableHits(0llu),
		_numValidHits(0llu),
		_keep(__keep),
		_hits(),
		_bufferedHits(),
		_strata(),
		_hitsForThisRead(),
		_max(__max) { }

	virtual ~HitSinkPerThread() { }

	/// Set whether to retain hits in a vector or not
	void setRetainHits(bool r)  { _keep = r; }
	/// Return whether we're retaining hits or not
	bool retainHits()           { return _keep; }

	/// Return the vector of retained hits
	vector<Hit>& retainedHits()   { return _hits; }
	vector<int>& retainedStrata() { return _strata; }

	/// Finalize current read
	virtual void finishRead(PatternSourcePerThread& p,
	                        ostream *dumpUnalignFa = NULL,
	                        ostream *dumpUnalignFq = NULL)
	{
		_bestRemainingStratum = 0;
		if(_hitsForThisRead == 0) {
			// Dump the unaligned read to one or more of the unaligned-
			// read output streams
			if(dumpUnalignFa != NULL) {
				(*dumpUnalignFa) << ">" << p.name() << endl << p.patFw() << endl;
			}
			if(dumpUnalignFq != NULL) {
				(*dumpUnalignFq) << "@" << p.name() << endl << p.patFw() << endl << "+" << endl << p.qualFw() << endl;
			}
		}
		finishReadImpl();
		if(_bufferHits && _bufferedHits.size() > 0) {
			// Flush buffered hits
			_sink.reportHits(_bufferedHits);
			_bufferedHits.clear();
		}
	}

	virtual void finishReadImpl() = 0;

	/**
	 * Implementation for hit reporting; update per-thread _hits and
	 * _numReportableHits variables and call the master HitSink to do the actual
	 * reporting
	 */
	virtual void reportHitImpl(const Hit& h, int stratum) {
		_numReportableHits++;
		if(_keep) {
#ifndef NDEBUG
			// Check whether any of the previous 256 hits are identical
			// to this one.  Really, we should check them all but that
			// would be too slow for long runs.
			size_t len = _hits.size();
			size_t lookBack = 256;
			if(lookBack > len) lookBack = len;
			for(size_t i = 0; i < lookBack; i++) {
				if(h.patId    == _hits[len-i-1].patId &&
				   h.h.first  == _hits[len-i-1].h.first &&
				   h.h.second == _hits[len-i-1].h.second &&
				   h.fw       == _hits[len-i-1].fw)
				{
					cerr << "Repeat hit: " << h.patId << (h.fw? "+":"-")
					     << ":<" << h.h.first << "," << h.h.second << ","
					     << h.mms.count() << ">" << endl;
					assert(false); // same!
				}
			}
#endif
			_hits.push_back(h);
			_strata.push_back(stratum);
		}
		if(_bufferHits) {
#ifndef NDEBUG
			// Ensure all buffered hits have the same patid
			for(size_t i = 1; i < _bufferedHits.size(); i++) {
				assert_eq(_bufferedHits[0].patId, _bufferedHits[i].patId);
			}
#endif
			_bufferedHits.push_back(h);
		} else {
			_sink.reportHit(h);
		}
	}

	/**
	 * Concrete subclasses override this to (possibly) report a hit and
	 * return true iff the caller should continue to report more hits.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		_numValidHits++;
		return true;
	}

	/// Return the number of valid hits so far (not necessarily
	/// reported).  It's up to the concrete subclasses
	uint64_t numValidHits()    { return _numValidHits; }

	/// Return the number of hits reported so far
	uint64_t numReportedHits() { return _numReportableHits; }

	/**
	 * The search routine is informing us that it will not be reporting
	 * any more hits at the given stratum.
	 */
	bool finishedWithStratum(int stratum) {
		bool ret = finishedWithStratumImpl(stratum);
		_bestRemainingStratum = stratum+1;
		return ret;
	}

	/**
	 * Concrete subclasses override this to determine whether the
	 * search routine should keep searching after having finished
	 * reporting all alignments at the given stratum.
	 */
	virtual bool finishedWithStratumImpl(int stratum) = 0;

	/// Return the maximum number of hits allowed per read
	virtual uint32_t maxHits() = 0;

	/// The mhits maximum
	uint32_t overThresh() { return _max; }

	/// Whether this thread, for this read, knows that we have already
	/// exceeded the mhits maximum
	bool exceededOverThresh() { return _hitsForThisRead > _max; }

	/// Return whether we span strata
	virtual bool spanStrata() = 0;

	/// Return whether we report only the best possible hits
	virtual bool best() = 0;

protected:
	HitSink&    _sink; /// Ultimate destination of reported hits
	/// Least # mismatches in alignments that will be reported in the
	/// future.  Updated by the search routine.
	int         _bestRemainingStratum;
	/// Whether to buffer-and-flush such that all hits for a given read
	/// appear consecutively (per phase)
	bool        _bufferHits;
	/// # hits reported to this HitSink so far (not all of which were
	/// necesssary reported to _sink)
	uint64_t    _numReportableHits;
	uint64_t    _numValidHits;
	bool        _keep; /// Whether to retain all reported hits in _hits
	vector<Hit> _hits; /// Repository for retained hits
	/// Buffered hits, to be reported and flushed at end of read-phase
	vector<Hit> _bufferedHits;
	vector<int> _strata; /// Repository for retained strata

	// Following variables are declared in the parent but maintained in
	// the concrete subcalsses
	uint32_t _hitsForThisRead; /// # hits for this read so far
	uint32_t _max;             /// don't report any hits if there were > _max
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
	        uint32_t __n,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
	        HitSinkPerThread(sink, __max, __keep),
	        _n(__n)
	{
		assert_gt(_n, 0);
	}

	virtual uint32_t maxHits() { return _n; }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return false; // we settle for "good" hits
	}

	/// Finalize current read
	virtual void finishReadImpl() {
		_hitsForThisRead = 0;
	}

	/**
	 * Report and then return true if we've already reported N good
	 * hits.  Ignore the stratum - it's not relevant for finding "good"
	 * hits.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		_hitsForThisRead++;
		if(_hitsForThisRead <= _n) {
			// Only report hit if we haven't
			reportHitImpl(h, stratum);
		}
		if(_hitsForThisRead > _max) {
			_bufferedHits.clear();
			return true; // done - report nothing
		}
		if(_hitsForThisRead == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits and max isn't set; stop!
		}
		return false; // not at N or max yet; keep going
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) { return false; }

private:
	uint32_t _n;               /// max # hits to report per read
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
	        uint32_t __n,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
	        HitSinkPerThread(sink, __max, __keep),
	        _n(__n)
	{
		assert_gt(_n, 0);
	}

	virtual uint32_t maxHits() { return _n; }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		assert_geq(stratum, _bestRemainingStratum);
		if(stratum == _bestRemainingStratum) {
			// This hit is within th best possible remaining stratum,
			// so it should definitely count
			_hitsForThisRead++;
			if(_hitsForThisRead > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			if(_hitsForThisRead <= _n) {
				reportHitImpl(h, stratum);
			}
			if(_hitsForThisRead == _n &&
			   (_max == 0xffffffff || _max < _n))
			{
				return true; // already reported N good hits; stop!
			}
		} else {
			// The alignment is not guaranteed to be in the best
			// stratum, but it might be, so save it for possible
			// future reporting
			_hitStrata[stratum].push_back(h);
		}
		return false; // not at N yet; keep going
	}

	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual void finishReadImpl() {
		if(_hitsForThisRead == _n || _hitsForThisRead > _max) {
			reset();
			return; // already reported N good hits; stop!
		}
		for(int j = 0; j < 4; j++) {
			for(size_t i = 0; i < _hitStrata[j].size(); i++) {
				// This hit is within th best possible remaining stratum,
				// so it should definitely count
				_hitsForThisRead++;
				if(_hitsForThisRead > _max) {
					_bufferedHits.clear();
					return; // done - report nothing
				}
				if(_hitsForThisRead <= _n) {
					reportHitImpl(_hitStrata[j][i], j);
				}
				if(_hitsForThisRead == _n &&
				   (_max == 0xffffffff || _max < _n))
				{
					reset();
					return; // already reported N good hits; stop!
				}
			}
			_hitStrata[j].clear();
		}
		reset();
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
#ifndef NDEBUG
		for(int j = 0; j < stratum; j++) {
			assert_eq(0, _hitStrata[j].size());
		}
#endif
		if(_hitsForThisRead == _n) {
			return true; // already reported N good hits; stop!
		}
		for(int j = stratum; j <= stratum+1; j++) {
			for(size_t i = 0; i < _hitStrata[j].size(); i++) {
				// This hit is within the best possible remaining stratum,
				// so it should definitely count
				_hitsForThisRead++;
				if(_hitsForThisRead > _max) {
					_bufferedHits.clear();
					return true; // done - report nothing
				}
				if(_hitsForThisRead <= _n) {
					reportHitImpl(_hitStrata[j][i], j);
				}
				if(_hitsForThisRead == _n &&
				   (_max == 0xffffffff || _max < _n))
				{
					_hitStrata[j].clear();
					return true; // already reported N good hits; stop!
				}
			}
			_hitStrata[j].clear();
		}
		return false; // keep going
	}

private:
	void reset() {
		clearAll();
		_hitsForThisRead = 0;
	}

	void clearAll() {
		for(int j = 0; j < 4; j++) {
			_hitStrata[j].clear();
		}
	}

	uint32_t _n;               /// max # hits to report
	vector<Hit> _hitStrata[4]; /// lower numbered strata are better
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
	        uint32_t __n,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
	        HitSinkPerThread(sink, __max, __keep),
	        _n(__n),
	        _bestStratumReported(999)
	{
		assert_gt(_n, 0);
	}

	virtual uint32_t maxHits() { return _n; }

	virtual bool spanStrata() {
		return false; // we do not span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		assert_geq(stratum, _bestRemainingStratum);
		if(stratum == _bestRemainingStratum &&
		   stratum <= _bestStratumReported)
		{
			// This hit is within th best possible remaining stratum,
			// so it should definitely count
			_hitsForThisRead++;
			if(_hitsForThisRead > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			if(_hitsForThisRead <= _n) {
				reportHitImpl(h, stratum);
			}
			_bestStratumReported = stratum;
			if(_hitsForThisRead == _n &&
			   (_max == 0xffffffff || _max < _n))
			{
				return true; // already reported N good hits; stop!
			}
		} else if(stratum <= _bestStratumReported) {
			_hitStrata[stratum].push_back(h);
			_bestStratumReported = stratum;
		} else {
			// No need to buffer the hit because we're guaranteed not
			// to report it
		}
		return false; // not at N yet; keep going
	}

	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual void finishReadImpl() {
		if(_bestStratumReported < 999 &&
		   _hitsForThisRead < _n &&
		   _hitsForThisRead <= _max)
		{
			assert_leq(_bestStratumReported, 3);
#ifndef NDEBUG
			// All better strata should be empty
			for(int j = 0; j < _bestStratumReported; j++) {
				assert_eq(0, _hitStrata[j].size());
			}
#endif
			for(size_t i = 0; i < _hitStrata[_bestStratumReported].size(); i++) {
				// This hit is within th best possible remaining stratum,
				// so it should definitely count
				_hitsForThisRead++;
				if(_hitsForThisRead > _max) {
					_bufferedHits.clear();
					return; // done - report nothing
				}
				if(_hitsForThisRead <= _n) {
					reportHitImpl(_hitStrata[_bestStratumReported][i], _bestStratumReported);
				}
				if(_hitsForThisRead == _n &&
				   (_max == 0xffffffff || _max < _n))
				{
					_hitStrata[_bestStratumReported].clear();
					break; // already reported N good hits; stop!
				}
			}
			_hitStrata[_bestStratumReported].clear();
		} else if(_bestStratumReported == 999) {
#ifndef NDEBUG
			// All strata should be empty
			for(int j = 0; j < 4; j++) {
				assert_eq(0, _hitStrata[j].size());
			}
#endif
		}
		reset();
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
		if(_hitsForThisRead == _n &&
		   (_max == 0xffffffff || _max < _n))
		{
			return true; // already reported N good hits; stop!
		}
		for(size_t i = 0; i < _hitStrata[stratum].size(); i++) {
			// This hit is within the best possible remaining stratum,
			// so it should definitely count
			_hitsForThisRead++;
			if(_hitsForThisRead > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			if(_hitsForThisRead <= _n) {
				reportHitImpl(_hitStrata[stratum][i], stratum);
			}
			if(_hitsForThisRead == _n &&
			   (_max == 0xffffffff || _max < _n))
			{
				_hitStrata[stratum].clear();
				return true; // already reported N good hits; stop!
			}
		}
		_hitStrata[stratum].clear();
		// reported at least once in this stratum; do not move to the
		// next stratum
		if(_hitsForThisRead > 0) return true;
		// didn't report; move to next stratum
		return false;
	}

private:

	void reset() {
		clearAll();
		_hitsForThisRead = 0;
		_bestStratumReported = 999;
	}

	void clearAll() {
		for(int j = 0; j < 4; j++) {
			_hitStrata[j].clear();
		}
	}

	uint32_t _n;               /// max # hits to report
	int _bestStratumReported;  /// stratum of best reported hit thus far
	vector<Hit> _hitStrata[4];
};

/**
 * Report all valid alignments.
 */
class AllHitSinkPerThread : public HitSinkPerThread {

public:
	AllHitSinkPerThread(
			HitSink& sink,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
		    HitSinkPerThread(sink, __max, __keep) { }

	virtual uint32_t maxHits() { return 0xffffffff; }

	virtual bool spanStrata() {
		return true; // we span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and always return true; we're finiding all hits so that
	 * search routine should always continue.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		_hitsForThisRead++;
		if(_hitsForThisRead > _max) {
			_bufferedHits.clear();
			return true; // done - report nothing
		}
		reportHitImpl(h, stratum);
		return false; // reporting all; always keep going
	}

	/**
	 * Finalize; do nothing because we haven't buffered anything
	 */
	virtual void finishReadImpl() {
		_hitsForThisRead = 0;
	}

	/**
	 * Always return false; search routine should not stop.
	 */
	virtual bool finishedWithStratumImpl(int stratum) { return false; }
};

/**
 * Report all alignments encountered in the best alignment stratum
 * for which there were any alignments.
 */
class AllStratifiedHitSinkPerThread : public HitSinkPerThread {

public:
	AllStratifiedHitSinkPerThread(
			HitSink& sink,
	        uint32_t __max = 0xffffffff,
	        bool __keep = false) :
	        HitSinkPerThread(sink, __max, __keep),
	        _bestStratumReported(999),
	        _reported(false) { }

	virtual uint32_t maxHits() { return 0xffffffff; }

	virtual bool spanStrata() {
		return false; // we do not span strata
	}

	virtual bool best() {
		return true; // we report "best" hits
	}

	/**
	 * Report and then return false if we've already reported N.
	 */
	virtual bool reportHit(const Hit& h, int stratum) {
		HitSinkPerThread::reportHit(h, stratum);
		assert_geq(stratum, _bestRemainingStratum);
		if(stratum == _bestRemainingStratum &&
		   stratum <= _bestStratumReported)
		{
			// This hit is within the best possible remaining stratum,
			// so it should definitely be reported
			_hitsForThisRead++;
			if(_hitsForThisRead > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			reportHitImpl(h, stratum);
			_reported = true;
			_bestStratumReported = stratum;
		} else if(stratum <= _bestStratumReported) {
			_hitStrata[stratum].push_back(h);
			_bestStratumReported = stratum;
		} else {
			// No need to buffer the hit because we're guaranteed not
			// to report it
		}
		return false; // keep going
	}

	/**
	 * Finalize current read by reporting any buffered hits from best
	 * to worst until they're all reported or until we've reported all
	 * N
	 */
	virtual void finishReadImpl() {
		if(_bestStratumReported < 999) {
#ifndef NDEBUG
			// All better strata should be empty
			for(int j = 0; j < _bestStratumReported; j++) {
				assert_eq(0, _hitStrata[j].size());
			}
#endif
			for(size_t i = 0; i < _hitStrata[_bestStratumReported].size(); i++) {
				// This hit is within th best possible remaining stratum,
				// so it should definitely count
				_hitsForThisRead++;
				if(_hitsForThisRead > _max) {
					_bufferedHits.clear();
					return; // done - report nothing
				}
				reportHitImpl(_hitStrata[_bestStratumReported][i], _bestStratumReported);
			}
		} else {
#ifndef NDEBUG
			// All strata should be empty
			for(int j = 0; j < 4; j++) {
				assert_eq(0, _hitStrata[j].size());
			}
#endif
		}
		reset();
	}

	/**
	 * Always return true; search routine should only stop if it's
	 * already reported N hits.
	 */
	virtual bool finishedWithStratumImpl(int stratum) {
#ifndef NDEBUG
		// All better strata should be empty
		for(int j = 0; j < stratum; j++) {
			assert_eq(0, _hitStrata[j].size());
		}
#endif
		for(size_t i = 0; i < _hitStrata[stratum].size(); i++) {
			// This hit is within the best possible remaining stratum,
			// so it should definitely count
			_hitsForThisRead++;
			if(_hitsForThisRead > _max) {
				_bufferedHits.clear();
				return true; // done - report nothing
			}
			reportHitImpl(_hitStrata[stratum][i], stratum);
			_reported = true;
		}
		_hitStrata[stratum].clear();
		// reported at least once in this stratum; do not move to the
		// next stratum
		if(_reported) return true; // stop
		// didn't report; move to next stratum
		return false;
	}

private:

	void reset() {
		clearAll();
		_reported = false;
		_bestStratumReported = 999;
	}

	void clearAll() {
		for(int j = 0; j < 4; j++) {
			_hitStrata[j].clear();
		}
	}

	int _bestStratumReported;          /// stratum of best reported hit thus far
	bool _reported;
	vector<Hit> _hitStrata[4];
};

/**
 * Sink that prints lines like this:
 * (pat-id)[-|+]:<hit1-text-id,hit2-text-offset>,<hit2-text-id...
 *
 * Activated with --concise
 */
class ConciseHitSink : public HitSink {
public:
	/**
	 * Construct a single-stream ConciseHitSink (default)
	 */
	ConciseHitSink(
			ostream&        __out,
			bool            __reportOpps = false,
			vector<string>* __refnames = NULL) :
		HitSink(__out, __refnames),
		_reportOpps(__reportOpps),
		_first(true),
		_numReported(0llu) { }

	/**
	 * Construct a multi-stream ConciseHitSink with one stream per
	 * reference string (see --refout)
	 */
	ConciseHitSink(
	        size_t          __numOuts,
			bool            __reportOpps = false,
			vector<string>* __refnames = NULL) :
		HitSink(__numOuts, __refnames),
		_reportOpps(__reportOpps),
		_first(true),
		_numReported(0llu) { }

	/**
	 * Append a verbose, readable hit to the given output stream.
	 */
	static void append(ostream& ss, const Hit& h, bool reportOpps) {
		ss << h.patId << (h.fw? "+" : "-") << ":";
    	// .first is text id, .second is offset
		ss << "<" << h.h.first << "," << h.h.second << "," << h.mms.count();
		if(reportOpps) ss << "," << h.oms;
		ss << ">" << endl;
	}

	/**
	 * Append a verbose, readable hit to the given output stream.
	 */
	void append(ostream& ss, const Hit& h) {
		ConciseHitSink::append(ss, h, this->_reportOpps);
	}

	/**
	 * Report a concise alignment to the appropriate output stream.
	 */
	virtual void reportHit(const Hit& h) {
		ostringstream ss;
		append(ss, h);
		lock(h.h.first);
		out(h.h.first) << ss.str();
		unlock(h.h.first);
		mainlock();
		_first = false;
		_numReported++;
		mainunlock();
	}

	/**
	 * Report a set of hits in one big critical section, so that there
	 * is no interlacing of hits for different reads.
	 */
	virtual void reportHits(vector<Hit>& hs) {
		if(hs.size() == 0) return;
		if(_outs.size() > 1 && hs.size() > 2) {
			sort(hs.begin(), hs.end());
		}
		for(size_t i = 0; i < hs.size(); i++) {
			const Hit& h = hs[i];
			if(i == 0) {
				lock(h.h.first);
			} else if(refIdxToStreamIdx(h.h.first) != refIdxToStreamIdx(hs[i-1].h.first)) {
				unlock(hs[i-1].h.first);
				lock(h.h.first);
			}
			append(out(h.h.first), h);
		}
		unlock(hs[hs.size()-1].h.first);
		mainlock();
		_first = false;
		_numReported += hs.size();
		mainunlock();
	}

	/**
	 * Wrap up and report a short summary.
	 */
	virtual void finish() {
		if(_first) {
			assert_eq(0llu, _numReported);
			cout << "No results" << endl;
		}
		else cout << "Reported " << _numReported << " alignments to " << _outs.size() << " output stream(s)" << endl;
	}

private:
	bool     _reportOpps;
	bool     _first;       /// true -> first hit hasn't yet been reported
	uint64_t _numReported;
};

/**
 * Sink that prints lines like this:
 * pat-name \t [-|+] \t ref-name \t ref-off \t pat \t qual \t #-alt-hits \t mm-list
 */
class VerboseHitSink : public HitSink {
public:
	/**
	 * Construct a single-stream VerboseHitSink (default)
	 */
	VerboseHitSink(ostream&        __out,
				   vector<string>* __refnames = NULL,
				   int             __partition = 0) :
	HitSink(__out, __refnames),
	_first(true),
	_numReported(0llu),
	_partition(__partition) { }

	/**
	 * Construct a multi-stream VerboseHitSink with one stream per
	 * reference string (see --refout)
	 */
	VerboseHitSink(size_t          __numOuts,
				   vector<string>* __refnames = NULL,
				   int             __partition = 0) :
	HitSink(__numOuts, __refnames),
	_first(true),
	_numReported(0llu),
	_partition(__partition) { }

	/**
	 * Given a line of output from the VerboseHitSink, parse it into a
	 * Hit object and return the object (by value).  If the reference
	 * and/or read are identified in the output by name, then the patid
	 * and/or h.first fields may be invalid.  TODO: load the id/name
	 * mapping from the .ebwt file so that we can handle either case.
	 */
	static bool readHit(Hit& h,
	                    istream& in,
	                    vector<string>* refnames,
	                    bool verbose = false)
	{
		char buf[4096];
		in.getline(buf, 4096);
		if(in.eof()) {
			return false;
		}
		if(in.bad()) {
			cerr << "Alignment file set \"bad\" bit" << endl;
			exit(1);
		}
		if(in.fail()) {
			cerr << "A line from the alignment file was longer than 4K" << endl;
			exit(1);
		}
		string line(buf);
		vector<string> toks;
		tokenize(line, "\t", toks);
		// Skip over a partition key, if one exists
		if(toks[0].find_first_of(" ", 0) != string::npos) {
			if(verbose) cout << "Erased partition string" << endl;
			toks.erase(toks.begin());
		}
		if(verbose) {
			for(size_t i = 0; i < toks.size(); i++) {
				cout << toks[i] << ", ";
			}
			cout << endl;
		}
		// Parse read name
		bool readNameIsIdx = true;
		const string& readName = toks[0];
		for(size_t i = 0; i < readName.length(); i++) {
			if(readName[i] < '0' || readName[i] > '9') readNameIsIdx = false;
		}
		uint32_t patid = 0;
		if(readNameIsIdx) {
			istringstream readIdxSs(readName);
			readIdxSs >> patid;
		}
		// Parse orientation
		bool orientation = (toks[1] == "+");
		// Parse reference sequence id
		bool refIsIdx = true;
		for(size_t i = 0; i < toks[2].length(); i++) {
			if(toks[2][i] < '0' || toks[2][i] > '9') {
				refIsIdx = false;
				break;
			}
		}
		istringstream refIdxSs(toks[2]);
		uint32_t refIdx;
		if(refIsIdx) refIdxSs >> refIdx;
		else if(refnames != NULL) {
			bool found = false;
			for(size_t i = 0; i < refnames->size(); i++) {
				if((*refnames)[i] == toks[2]) {
					found = true;
					refIdx = i;
					break;
				}
			}
			if(!found) {
				refIdx = refnames->size();
				refnames->push_back(toks[2]);
			}
		} else {
			// reference was named and we have no way of mapping it to
			// an index, so we ignore it
			cerr << "Could not find an id to map reference name \"" << toks[2] << "\" to." << endl;
			exit(1);
		}
		// Parse reference sequence offset
		istringstream refOffSs(toks[3]);
		uint32_t refOff; refOffSs >> refOff;
		// Parse read sequence (5' is on RHS iff !orientation)
		const string& readSeq = toks[4];
		// Parse read qualities (5' is on RHS iff !orientation)
		const string& readQual = toks[5];
		assert_eq(readSeq.length(), readQual.length());
		// Parse the # other hits at this stratum estimate
		istringstream omsSs(toks[6]);
		uint32_t oms; omsSs >> oms;
		vector<char> refcs;
		refcs.resize(readSeq.length(), 0);
		FixedBitset<max_read_bp> mms;
		// toks.size() == 8 iff there are one or more mismatches
		if(toks.size() == 8) {
			vector<string> mmStrs;
			// Separate comma-delimited mismatch tokens
			tokenize(toks[7], ",", mmStrs);
			for(size_t i = 0; i < mmStrs.size(); i++) {
				vector<string> cs;
				// Split the offset from the characters
				tokenize(mmStrs[i], ":", cs);
				assert_eq(2, cs.size());
				// first character of the "A>C" string is the reference
				// character
				istringstream offSs(cs[0]);
				uint32_t off; offSs >> off;
				assert_lt(off, readSeq.length());
				// 'off' is an offset from the 5' end of the read.  We
				// want to translate it into an offset from the left
				// end of the read as it aligns to the forward strand
				// of the reference.  This means that we need to invert
				// it if the 5' end is on the right (i.e., if the
				// orientation is "-")
				//
				// BTL: never mind; I think this is too confusing.  We
				// shouldn't overload the Hit structure such that there
				// are multiple interpretations of its fields.
				//
				size_t noff = off;
				//if(!orientation) {
				//	// read is reversed; invert offset
				//	noff = readSeq.length() - off - 1;
				//}
				mms.set(noff);
				refcs[noff] = cs[1][0]; // reference char is before the >
				//assert_eq(readSeq[noff], cs[1][2]);
				if(orientation) {
					assert_eq(readSeq[noff], cs[1][2]);
				} else {
					assert_eq(readSeq[readSeq.length() - noff - 1], cs[1][2]);
				}
				if(verbose) {
					cout << "  Set mm at offset " << noff << " to "
					     << refcs[noff] << endl;
				}
			}
			if(verbose) cout << "Parsing mismatches" << endl;
		}
		h.h = make_pair<uint32_t>(refIdx, refOff);
		h.patId = (uint32_t)patid; // patid
		h.patName = String<char>(readName);
		h.patSeq = String<Dna5>(readSeq);
		h.quals = String<char>(readQual);
		h.fw = orientation; // fw
		h.mms = mms;    // mms
		h.refcs = refcs;  // refcs
		h.oms = oms;   // oms
		return true;
	}

	/**
	 * Append a verbose, readable hit to the given output stream.
	 */
	static void append(ostream& ss,
	                   const Hit& h,
	                   const vector<string>* refnames,
	                   int partition)
	{
		if(partition > 0) {
			// Output a partitioning key
			ss << h.h.first;
			ostringstream ss2; ss2 << (h.h.second / partition);
			string s2 = ss2.str();
			while(s2.length() < 10) {
				s2 = "0" + s2;
			}
			assert_eq(10, s2.length());
			ss << " " << s2.c_str() << "\t";
		}
		ss << h.patName << "\t" << (h.fw? "+":"-") << "\t";
    	// .first is text id, .second is offset
		if(refnames != NULL && h.h.first < refnames->size()) {
			ss << (*refnames)[h.h.first];
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
	}

	/**
	 * Append a verbose, readable hit to the output stream
	 * corresponding to the hit.
	 */
	void append(ostream& ss, const Hit& h) {
		VerboseHitSink::append(ss, h, this->_refnames, this->_partition);
	}


	/**
	 * Report a verbose, human-readable alignment to the appropriate
	 * output stream.
	 */
	virtual void reportHit(const Hit& h) {
		ostringstream ss;
		append(ss, h);
		// Make sure to grab lock before writing to output stream
		lock(h.h.first);
		out(h.h.first) << ss.str();
		unlock(h.h.first);
		mainlock();
		_first = false;
		_numReported++;
		mainunlock();
	}

	/**
	 * Report a list of verbose, human-readable alignments to the
	 * appropriate output stream.
	 */
	virtual void reportHits(vector<Hit>& hs) {
		if(hs.size() == 0) return;
		if(_outs.size() > 1 && hs.size() > 2) {
			sort(hs.begin(), hs.end());
		}
		for(size_t i = 0; i < hs.size(); i++) {
			const Hit& h = hs[i];
			if(i == 0) {
				lock(h.h.first);
			} else if(refIdxToStreamIdx(h.h.first) != refIdxToStreamIdx(hs[i-1].h.first)) {
				unlock(hs[i-1].h.first);
				lock(h.h.first);
			}
			append(out(h.h.first), h);
		}
		unlock(hs[hs.size()-1].h.first);
		mainlock();
		_first = false;
		_numReported += hs.size();
		mainunlock();
	}

	/**
	 * Finalize the alignment output by printing a summary message to
	 * stdout.
	 */
	virtual void finish() {
		if(_first) {
			assert_eq(0llu, _numReported);
			cout << "No results" << endl;
		} else {
			cout << "Reported " << _numReported << " alignments to "
			     << _outs.size() << " output stream(s)" << endl;
		}
	}

private:
	bool     _first;       /// true iff this object hasn't yet reported a hit
	uint64_t _numReported; /// number of hits reported
	int      _partition;   /// partition size, or 0 if partitioning is disabled
};

/**
 * Sink for outputting alignments in a binary format.
 */
class BinaryHitSink : public HitSink {
public:

	/**
	 * Construct a single-stream BinaryHitSink (default)
	 */
	BinaryHitSink(ostream&        __out,
				  vector<string>* __refnames = NULL) :
	HitSink(__out, __refnames),
	_first(true), _numReported(0llu) { }

	/**
	 * Construct a multi-stream BinaryHitSink with one stream per
	 * reference string (see --refout)
	 */
	BinaryHitSink(size_t          __numOuts,
				  vector<string>* __refnames = NULL) :
	HitSink(__numOuts, __refnames),
	_first(true), _numReported(0llu) { }

	/**
	 * Append a binary alignment to the output stream corresponding to
	 * the reference sequence involved.
	 */
	static void append(ostream& o,
					   const Hit& h,
					   const vector<string>* refnames)
	{
		// true iff we're going to print the reference sequence name
		bool refName = refnames != NULL &&
		                h.h.first < refnames->size();
		uint16_t pnamelen = (uint16_t)length(h.patName);
		// Write read name
		o.write((const char *)&pnamelen, 2);
		o.write(begin(h.patName), pnamelen);
		// Write fw/refname flags
		uint8_t flags = (h.fw ? 1 : 0) | (refName? 2 : 0);
		o.write((const char *)&flags, 1);
		if(refName) {
			// Write reference name as string
			uint16_t rnamelen = (uint16_t)(*refnames)[h.h.first].length();
			o.write((const char *)&rnamelen, 2);
			o.write((*refnames)[h.h.first].c_str(), rnamelen);
		} else {
			// Write reference name as index into global reference name
			// list
			o.write((const char *)&h.h.first, 4);
		}
		// Write reference offset
		o.write((const char *)&h.h.second, 4);
		// Write pattern sequence
		uint16_t plen = (uint16_t)length(h.patSeq);
		o.write((const char *)&plen, 2);
		for(size_t i = 0; i < plen; i += 2) {
			uint8_t twoChrs = (uint8_t)h.patSeq[i];
			if(i+1 < plen) {
				twoChrs |= ((uint8_t)h.patSeq[i+1] << 4);
			}
			o.write((const char *)&twoChrs, 1);
		}
		// Write quals sequence
		uint16_t qlen = (uint16_t)length(h.quals);
		assert_eq(plen, qlen);
		o.write(begin(h.quals), qlen);
		// Write oms
		o.write((const char *)&h.oms, 4);
		// Write # mismatches
		uint8_t numMms = h.mms.count();
		o.write((const char *)&numMms, 1);
		// Output mismatches
		size_t c = 0;
		for (uint8_t i = 0; i < h.mms.size(); ++ i) {
			if (h.mms.test(i)) {
				o.write((const char *)&i, 1);
				assert_gt(h.refcs.size(), i);
				assert_eq(1, dna4Cat[(int)h.refcs[i]]);
				uint8_t refChar = charToDna5[(int)h.refcs[i]];
				assert_leq(refChar, 4);
				uint8_t qryChar = (h.fw ? (int)h.patSeq[i] :
				                          (int)h.patSeq[length(h.patSeq)-i-1]);
				assert_leq(refChar, 4);
				assert_neq(refChar, qryChar);
				uint8_t both = refChar | (qryChar << 4);
				o.write((const char *)&both, 1);
				c++;
			}
		}
	}

	/**
	 * Append a binary alignment to the output stream corresponding to
	 * the reference sequence involved.
	 */
	void append(ostream& o, const Hit& h) {
		BinaryHitSink::append(o, h, this->_refnames);
	}

	/**
	 * Read a binary-encoded hit (written by append() above) from an
	 * input stream.
	 */
	static bool readHit(Hit& h,
	                    istream& in,
	                    vector<string>* refnames,
	                    bool verbose = false)
	{
		if(!in.good()) return false;
		uint16_t pnamelen;
		in.read((char *)&pnamelen, 2);
		if(in.eof()) return false;
		if(verbose) cout << "name(" << pnamelen << ")=\"";
		seqan::resize(h.patName, pnamelen);
		in.read((char*)begin(h.patName), pnamelen);
		if(verbose) cout << h.patName << "\", ";
		// Parse read name
		h.patId = 0;
		for(int i = 0; i < pnamelen; i++) {
			if(h.patName[i] < '0' || h.patName[i] > '9') {
				h.patId = 0;
				break;
			}
			h.patId *= 10;
			h.patId += (h.patName[i] - '0');
		}
		if(verbose) cout << "patid(" << h.patId << ")=\"";
		// Write fw/refname flags
		uint8_t flags;
		in.read((char *)&flags, 1);
		h.fw         = ((flags & 1) != 0);
		bool refName = ((flags & 2) != 0);
		if(verbose) cout << "fw=" << (h.fw ? "true":"false") << ", ";
		if(refName) {
			// Read and ignore reference name
			uint16_t rnamelen;
			char buf[2048];
			in.read((char *)&rnamelen, 2);
			in.read(buf, rnamelen);
			buf[rnamelen] = '\0';
			if(verbose) cout << "refname=\"" << buf << "\", ";
			string name(buf);
			// Add to the refnames vector; this isn't efficient - if
			// this ends up being a problem, we should use a map
			// instead of a vector
			if(refnames != NULL) {
				bool found = true;
				h.h.first = 0;
				for(int i = 0; i < rnamelen; i++) {
					if(buf[i] < '0' || buf[i] > '9') {
						h.h.first = 0;
						found = false;
						break;
					}
					h.h.first *= 10;
					h.h.first += (buf[i] - '0');
				}
				if(!found) {
					for(size_t i = 0; i < refnames->size(); i++) {
						if((*refnames)[i] == name) {
							h.h.first = i;
							found = true;
							break;
						}
					}
					if(!found) {
						h.h.first = refnames->size();
						refnames->push_back(name);
					}
				}
			} else {
				h.h.first = 0;
				for(int i = 0; i < rnamelen; i++) {
					if(buf[i] < '0' || buf[i] > '9') {
						h.h.first = 0;
						break;
					}
					h.h.first *= 10;
					h.h.first += (buf[i] - '0');
				}
			}
			if(verbose) cout << "refidx=\"" << h.h.first << "\", ";
		} else {
			// Read reference id as index into global reference name list
			in.read((char*)&h.h.first, 4);
			if(verbose) cout << "refidx=" << h.h.first << ", ";
		}
		// Read reference offset
		in.read((char*)&h.h.second, 4);
		if(verbose) cout << "refoff=" << h.h.second << ", ";
		// Read pattern length
		uint16_t plen;
		in.read((char*)&plen, 2);
		if(verbose) cout << "plen=" << plen << ", ";
		h.refcs.resize(plen, 0);
		assert_gt(plen, 1);
		assert_lt(plen, 1024);
		// Read pattern sequence
		seqan::resize(h.patSeq, plen);
		for(size_t i = 0; i < plen; i += 2) {
			uint8_t twoChars;
			in.read((char *)&twoChars, 1);
			assert_lt((twoChars & 15), 5);
			assert_lt((twoChars >> 4), 5);
			h.patSeq[i] = (twoChars & 15);
			if(i+1 < plen) {
				h.patSeq[i+1] = (twoChars >> 4);
			}
		}
		if(verbose) cout << "pat=\"" << h.patSeq << "\", ";
		// Read quals sequence
		uint16_t qlen = plen;
		seqan::resize(h.quals, qlen);
		in.read((char*)begin(h.quals), qlen);
		if(verbose) cout << "quals=\"" << h.quals << "\", ";
#ifndef NDEBUG
		for(size_t i = 0; i < length(h.quals); i++) {
			assert_geq(h.quals[i], 33);
			assert_leq(h.quals[i], 126);
		}
#endif
		// Read oms
		in.read((char *)&h.oms, 4);
		if(verbose) cout << "oms=" << h.oms << ", ";
		// Read # mismatches
		uint8_t numMms;
		in.read((char *)&numMms, 1);
		if(verbose) cout << "nummms=" << ((int)numMms) << ", ";
		// Read mismatches
		for(uint8_t i = 0; i < numMms; i++) {
			uint8_t ii;
			// Read the read offset (from 5' end)
			in.read((char*)&ii, 1);
			h.mms.set(ii);
			uint8_t both;
			// Read the reference and query characters involved
			in.read((char*)&both, 1);
			uint8_t refChar = both & 15;
			uint8_t qryChar = both >> 4;
			assert_neq(refChar, qryChar);
			h.refcs[ii] = "ACGTN"[refChar];
			if(verbose) cout << ((int)ii) << ":" << "ACGTN"[refChar] << ">" << "ACGTN"[qryChar] << ", ";
		}
		if(verbose) cout << endl;
		// Done
		return true;
	}

	/**
	 * Report a single hit to the appropriate output stream.
	 */
	virtual void reportHit(const Hit& h) {
		lock(h.h.first);
		append(out(h.h.first), h);
		unlock(h.h.first);
		mainlock();
		_first = false;
		_numReported++;
		mainunlock();
	}

	/**
	 * Report a list of hits to the appropriate output stream.
	 */
	virtual void reportHits(vector<Hit>& hs) {
		if(hs.size() == 0) return;
		// Sort the hits in order of
		if(_outs.size() > 1 && hs.size() > 2) {
			sort(hs.begin(), hs.end());
		}
		for(size_t i = 0; i < hs.size(); i++) {
			const Hit& h = hs[i];
			if(i == 0) {
				// Lock the first stream
				lock(h.h.first);
			} else if(refIdxToStreamIdx(h.h.first) != refIdxToStreamIdx(hs[i-1].h.first)) {
				// Move to the next stream; be sure to unlock before
				// locking the next one
				unlock(hs[i-1].h.first);
				lock(h.h.first);
			}
			// Actually write the hit
			append(out(h.h.first), h);
		}
		// Unlock the last stream
		unlock(hs[hs.size()-1].h.first);
		mainlock();
		_first = false;
		_numReported += hs.size();
		mainunlock();
	}

	/**
	 * Finalize the alignment output by printing a summary message to
	 * stdout.
	 */
	virtual void finish() {
		if(_first) {
			assert_eq(0llu, _numReported);
			cout << "No results" << endl;
		} else {
			cout << "Reported " << _numReported << " alignments to "
			     << _outs.size() << " output stream(s)" << endl;
		}
	}
private:
	bool _first;           /// true iff this object hasn't yet reported a hit
	uint64_t _numReported; /// number of hits reported
};

/**
 * Sink that does nothing.
 */
class StubHitSink : public HitSink {
public:
	StubHitSink() : HitSink(cout, NULL) { }
	virtual void reportHit(const Hit& h) { }
	virtual void finish() {	}
};

#endif /*HIT_H_*/

