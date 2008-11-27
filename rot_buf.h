#ifndef ROT_BUF_H_
#define ROT_BUF_H_

#include "alphabet.h"

/**
 * Abstract parent for classes that analyze information aggregated in
 * multiple-alignment columns.  Individual column elements have type T.
 */
template<typename T>
class ColumnAnalyzer {
public:
	virtual ~ColumnAnalyzer() { }
	virtual void analyze(uint32_t refidx,
	                     uint32_t refoff,
	                     const vector<T>& col) = 0;
private:
};

/**
 * A concrete column-analysis class where each column is a vector of
 * char, char pairs where the first char contains the reference
 * character at that position (hi 4 bits) along with the read character
 * for that column and row (lo 4 bits).  The second char contains the
 * quality value.
 */
class SNPColumnCharPairAnalyzer : public ColumnAnalyzer<pair<char, char> > {
public:
	SNPColumnCharPairAnalyzer(ostream& out) :
		ColumnAnalyzer<pair<char, char> >(), out_(out) { }

	virtual ~SNPColumnCharPairAnalyzer() { }

	/**
	 * Analyze a column of ref char, qry char, qual values and output a
	 * homozygous SNP call where evidence is sufficient.
	 */
	virtual void analyze(uint32_t refidx,
	                     uint32_t refoff,
	                     const vector<pair<char, char> >& col)
	{
		// Quality "mass" (total) for all of A/C/G/T
		int totalMass = 0;
		// Quality "mass" (total) for each of A/C/G/T
		int snpMass[4] = {0, 0, 0, 0};
		int rc = -1;
		// Iterate through rows in column
		for(size_t i = 0; i < col.size(); i++) {
			// Update total mass
			totalMass += (col[i].second - 33);
			// If this row is different from the reference...
			if(col[i].first != -1) {
				// Pick out the reference character
				assert_lt((col[i].first >> 4), 4);
				if(rc == -1) {
					rc = col[i].first >> 4; // ref
					assert_lt(rc, 4);
				} else {
					// Assert that reference character is the same
					assert_eq(rc, (col[i].first >> 4));
				}
				// Pick out the read character
				int c = col[i].first & 15;
				if(c == 4) continue; // If it's an N, ignore it
				assert_lt(c, 4);
				assert_neq(c, rc); // Read char must not match ref char
				// Update the "mass" for the read character
				snpMass[c] += ((int)col[i].second - 33);
			}
		}
		// Bail if we didn't encounter any rows
		if(rc == -1) return;
		assert_eq(0, snpMass[rc]);
		int maxSnpMass = 0;
		// If the total mass exceeds a sufficient threshold
		if(totalMass > 300) {
			// Sufficient evidence to try to make a call
			for(int i = 0; i < 4; i++) {
				// Find the largest per-A/C/G/T mass
				if(snpMass[i] > maxSnpMass) {
					maxSnpMass = snpMass[i];
				}
				// Take ratio of highest-mass character to the total
				// mass
				double ratio = snpMass[i] / (double)totalMass;
				if(ratio > 0.85) {
					// 85% is the threshold
					assert_neq(-1, rc);
					out_ << refidx << " " << refoff << " "
					     << "ACGT"[rc] << " " << "ACGT"[i] << endl;
				}
			}
		}
	}

private:
	ostream& out_;
};

/**
 * User may add more elements to the conceptual right-hand side,
 * "ash" elements off the conceptual left-hand side, and access
 * elements between the sides.
 *
 * T is the type of element stored in the column vectors.
 */
template<typename T, int S>
class RotatingBuf {
public:

	/**
	 * Initialize new rotating buf.
	 */
	RotatingBuf() : buf_(), refidx_(0), lpos_(0), rpos_(0) {
		buf_.resize(S);
	}

	/**
	 * Potentially slide the window of interest to the right so that it
	 * encompasses reference positions newLhs through newRhs.  This may
	 * involve "ashing" and analyzing finished columns from the left-
	 * hand-side using the provided ColumnAnalyzer.
	 */
	void ashAndExtendRange(uint32_t newLhs,
	                       uint32_t newRhs,
	                       ColumnAnalyzer<T> *analyzer = NULL)
	{
		assert_leq(rpos_ - lpos_, S); // size of window of interest can't exceed S
		assert_leq(newLhs, newRhs);   // left-hand-side can't be to the right of right-hand
		// New left-hand side must be greater than or equal to previous
		// left-hand side.
		if(newLhs < lpos_) {
			cerr << "Error: unsorted references supplied to RotatingBuf" << endl;
			exit(1);
		}
		// Check if there are columns we're finished adding rows to
		if(newLhs > lpos_ && analyzer != NULL) {
			for(size_t i = lpos_; i < newLhs; i++) {
				// Column i is done, send it to the column analyzer
				if(buf_[i % S].size() > 0) {
					// Analyze this column
					analyzer->analyze(refidx_, i, buf_[i % S]);
				}
				// Free it up
				buf_[i % S].clear();
			}
		}
		// Update left-hand side
		lpos_ = newLhs;
		if(newRhs > rpos_) {
#ifndef NDEBUG
			for(size_t i = rpos_; i < newRhs; i++) {
				assert_eq(0, buf_[i % S].size());
			}
#endif
			rpos_ = newRhs;
		}
		assert_leq(rpos_ - lpos_, S); // size of window of interest can't exceed S
	}

	/**
	 * Add an element to the back of the column vector at position pos.
	 * pos is an absolute reference coordinate and must lie in the
	 * current window of interest.
	 */
	void add(uint32_t pos, const T& elt) {
		assert_leq(rpos_ - lpos_, S);
		assert_geq(pos, lpos_);
		assert_lt(pos, rpos_);
		buf_[pos % S].push_back(elt);
	}

	/**
	 * Return a reverence to the the column vector at absolute
	 * reference position pos.  pos must lie in the current window of
	 * interest.
	 */
	const vector<T>& get(uint32_t pos) {
		assert_geq(pos, lpos_);
		assert_lt(pos, rpos_);
		return buf_[pos % S];
	}

	/**
	 * Finish analyzing the rotating buffer.  Any as-yet-unanalyzed
	 * non-empty columns are submitted to the given ColumnAnalyzer then
	 * cleared.
	 */
	void finalize(ColumnAnalyzer<T> *analyzer = NULL) {
#ifndef NDEBUG
		size_t mn = min(lpos_ % S, rpos_ % S);
		size_t mx = max(lpos_ % S, rpos_ % S);
		for(size_t i = 0; i < mn; i++) {
			assert_eq(0, buf_[i % S].size());
		}
		for(size_t i = mx; i < S; i++) {
			assert_eq(0, buf_[i % S].size());
		}
#endif
		for(size_t i = lpos_; i < rpos_; i++) {
			if(buf_[i % S].size() > 0) {
				// Analyze this column
				if(analyzer != NULL) {
					analyzer->analyze(refidx_, i, buf_[i % S]);
				}
			}
			buf_[i % S].clear();
		}
	}

	/**
	 * Reset the rotating buffer in preparation for analyzing a new
	 * reference.
	 */
	void reset(uint32_t refidx, ColumnAnalyzer<T> *analyzer = NULL) {
		finalize(analyzer);
		refidx_ = refidx;
		lpos_ = 0;
		rpos_ = 0;
	}

private:
	vector<vector<T> > buf_; // vector of columns, where each column is a vector of T
	uint32_t refidx_; // id of the reference being scanned using the buffer
	uint32_t lpos_; // in absolute coordinates
	uint32_t rpos_; // in absolute coordinates
};

/**
 * An abstract parent class for classes that consume alignments.
 */
template<typename T, int S>
class AlignmentSink {

public:

	AlignmentSink(bool verbose = false) :
	buf_(),
	lastRefIdx_(0),
	lastRefOff_(0),
	lastLeft_(0),
	furthestRight_(0),
	verbose_(verbose) { }

	virtual ~AlignmentSink() { }

	/**
	 * Add a new alignment to the rotating buffer.
	 */
	virtual void addAlignment(const Hit& h,
	                          ColumnAnalyzer<T> *analyzer = NULL)
	{
		if(this->verbose_) {
			cout << "    Entered AlignmentSink::addAlignmentImpl" << endl;
			cout << "      Analyzer is " << (analyzer ? "non-null" : "null") << endl;
		}
		// Check if we encountered a new reference
		if(h.h.first != lastRefIdx_) {
			// Reset state and buffer for new reference sequence
			assert_gt(h.h.first, lastRefIdx_); // must have greater idx
			lastRefIdx_ = h.h.first;
			lastRefOff_ = 0;
			lastLeft_ = 0;
			furthestRight_ = 0;
			// Call reset() to finalize the old buffer and reset for
			// the new buffer.
			buf_.reset(h.h.first, analyzer);
		}
		uint32_t len = h.length();
		uint32_t rhs = h.h.second + len;
		// Update the underlying rotating buffer so that it's ready to
		// accept the column-wise evidence contributed by this
		// alignment
		buf_.ashAndExtendRange(h.h.second, rhs, analyzer);
		addAlignmentImpl(h, analyzer);
	}

	/**
	 * Analyze any remaining columns in the current buffer.
	 */
	virtual void finalize(ColumnAnalyzer<T> *analyzer = NULL) {
		this->buf_.finalize(analyzer);
		finalizeImpl(analyzer);
	}

	/**
	 * Reset the rotating buffer to deal with a new reference sequence.
	 * May involve analyzing any remaining columns in the current
	 * buffer.
	 */
	virtual void reset(uint32_t refidx,
	                   ColumnAnalyzer<T> *analyzer = NULL)
	{
		this->buf_.reset(refidx, analyzer);
		resetImpl(refidx, analyzer);
	}

protected:

	virtual void addAlignmentImpl(
			const Hit& h,
	        ColumnAnalyzer<T> *analyzer = NULL) { }

	virtual void finalizeImpl(
			ColumnAnalyzer<T> *analyzer = NULL) { }

	virtual void resetImpl(
			uint32_t refidx,
			ColumnAnalyzer<T> *analyzer = NULL) { }

	RotatingBuf<T, S> buf_;
	uint32_t lastRefIdx_;
	uint32_t lastRefOff_;
	uint32_t lastLeft_;
	uint32_t furthestRight_;
	bool verbose_;
};

/**
 * A concrete consumer of alignments that
 *
 * Assumes that alignments are coming in sorted by left-hand position
 * of alignments as they occur in the forward strand of the reference.
 *
 * We instantiate the RotatingBuf template with pair<char> on the
 * assumption that analyses consuming the alignments will only care
 * about keeping each character along with its quality value.  This
 * might not be true for more sophisticated analyses.
 */
template<int S>
class RotatingCharPairAlignmentBuf : public AlignmentSink<pair<char, char>, S> {

public:
	RotatingCharPairAlignmentBuf(bool verbose = false) :
		AlignmentSink<pair<char, char>, S>(verbose) { }

protected:

	/**
	 * Add a new alignment to the rotating buffer.
	 */
	virtual void addAlignmentImpl(const Hit& h,
	                              ColumnAnalyzer<pair<char, char> > *analyzer = NULL)
	{
		if(this->verbose_) {
			cout << "    Entered RotatingCharPairAlignmentBuf::addAlignmentImpl" << endl;
			cout << "      Analyzer is " << (analyzer ? "non-null" : "null") << endl;
		}
		uint32_t len = h.length();
		// Add char-pair evidence to each column
		for(size_t i = 0; i < len; i++) {
			if(h.mms.test(i)) {
				// The i'th character from the 5' end has a mismatch.
				// Let ii = the offset from the left-hand side of the
				// alignment (not necessarily the 5' end of the read).
				size_t ii = i;
				if(!h.fw) ii = len - i - 1;
				char q = h.quals[ii];
				int readc = (int)h.patSeq[ii];
				if(readc == 4) continue; // no evidence inherent in Ns
				char c = h.refcs[i]; // refcs also indexed from 5' end of read
				assert_eq(1, dna4Cat[(int)c]);
				c = charToDna5[(int)c];
				assert_lt(c, 4);
				assert_neq((int)c, readc);
#ifndef NDEBUG
				// Ensure that the reference character for the evidence
				// we're adding matches the reference character for all
				// evidence we've already added
				const vector<pair<char, char> >& col = this->buf_.get(h.h.second + ii);
				for(size_t j = 0; j < col.size(); j++) {
					int cc = (col[j].first >> 4);
					assert_eq((int)c, cc);
				}
#endif
				c = (c << 4) | readc;
				this->buf_.add(h.h.second + ii, make_pair(c, q));
			}
		}
	}
};

#endif /*ROT_BUF_H_*/
