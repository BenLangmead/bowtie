#ifndef ROT_BUF_H_
#define ROT_BUF_H_

#include "alphabet.h"

/**
 * Abstract parent for classes that analyze information aggregated in
 * multiple-alignment columns.
 */
template<typename T>
class ColumnAnalyzer {
public:
	virtual void analyze(uint32_t refidx,
	                     uint32_t refoff,
	                     const vector<T>& col) = 0;
private:
};

/**
 *
 */
class SNPColumnCharPairAnalyzer : public ColumnAnalyzer<pair<char, char> > {
public:
	SNPColumnCharPairAnalyzer(ostream& out) :
		ColumnAnalyzer<pair<char, char> >(), out_(out) { }

	virtual void analyze(uint32_t refidx,
	                     uint32_t refoff,
	                     const vector<pair<char, char> >& col)
	{
		int totalMass = 0;
		int snpMass[4] = {0, 0, 0, 0};
		int rc = -1;
		for(size_t i = 0; i < col.size(); i++) {
			totalMass += (col[i].second - 33);
			if(col[i].first != -1) {
				if(rc == -1) {
					rc = col[i].first >> 4; // ref
				} else {
					assert_eq(rc, (col[i].first >> 4));
				}
				int c = col[i].first & 15; // read
				if(c == 4) continue;
				assert_lt(c, 4);
				assert_lt(rc, 4);
				assert_neq(c, rc);
				snpMass[c] += ((int)col[i].second - 33);
			}
		}
		assert(snpMass[0] == 0 || snpMass[1] == 0 ||
		       snpMass[2] == 0 || snpMass[3] == 0);
		int maxSnpMass = 0;
		if(totalMass > 300) {
			// Sufficient evidence to try to make a call
			for(int i = 0; i < 4; i++) {
				if(snpMass[i] > maxSnpMass) {
					maxSnpMass = snpMass[i];
				}
				double ratio = snpMass[i] / (double)totalMass;
				if(ratio > 0.85) {
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
 */
template<typename T, int S>
class RotatingBuf {
public:
	RotatingBuf() : buf_(), refidx_(0), lpos_(0), rpos_(0) {
		buf_.resize(S);
	}

	/**
	 *
	 */
	void extendRhs(uint32_t newRhs) {
		for(size_t i = rpos_; i < newRhs; i++) {
			buf_[i % S].clear();
		}
		rpos_ = newRhs;
		assert_leq(rpos_ - lpos_, S);
	}

	/**
	 *
	 */
	void ashUpTo(uint32_t newLhs,
	             ColumnAnalyzer<T> *analyzer = NULL)
	{
		assert_leq(rpos_ - lpos_, S);
		if(newLhs > lpos_ && analyzer != NULL) {
			// We've gathered all the data that we're going to for
			// these columns, so analyze them before
			for(size_t i = lpos_; i < newLhs; i++) {
				if(buf_[i % S].size() > 0) {
					// Analyze this column
					analyzer->analyze(refidx_, i, buf_[i % S]);
				}
				buf_[i % S].clear();
			}
		}
		lpos_ = newLhs;
		assert_leq(lpos_, rpos_);
	}

	/**
	 *
	 */
	void ashAndExtendRange(uint32_t newLhs,
	                       uint32_t newRhs,
	                       ColumnAnalyzer<T> *analyzer = NULL)
	{
		assert_leq(rpos_ - lpos_, S);
		assert_leq(newLhs, newRhs);
		assert_geq(newLhs, lpos_);
		if(newLhs < lpos_) {
			cerr << "Error: alignments seem not to be sorted in reference order.  Please sort before" << endl
			     << "providing alignments to bowtie-asm." << endl;
			exit(1);
		}
		if(newLhs > lpos_ && analyzer != NULL) {
			// We've gathered all the data that we're going to for
			// these columns, so analyze them before
			for(size_t i = lpos_; i < newLhs; i++) {
				if(buf_[i % S].size() > 0) {
					// Analyze this column
					analyzer->analyze(refidx_, i, buf_[i % S]);
				}
				buf_[i % S].clear();
			}
		}
		lpos_ = newLhs;
		if(newRhs > rpos_) {
			for(size_t i = max(lpos_, rpos_); i < newRhs; i++) {
				buf_[i % S].clear();
			}
			rpos_ = newRhs;
		}
		assert_leq(rpos_ - lpos_, S);
		assert_leq(lpos_, rpos_);
	}

	/**
	 *
	 */
	void add(uint32_t pos, const T& elt) {
		assert_leq(rpos_ - lpos_, S);
		assert_geq(pos, lpos_);
		assert_lt(pos, rpos_);
		buf_[pos % S].push_back(elt);
	}

	/**
	 * Return the vector-of-Ts at absolute position 'pos'.
	 */
	const vector<T>& get(uint32_t pos) {
		assert_geq(pos, lpos_);
		assert_lt(pos, rpos_);
		return buf_[pos % S];
	}

	/**
	 *
	 */
	void finalize(ColumnAnalyzer<T> *analyzer = NULL) {
		if(analyzer == NULL) {
			// nothing to do
			return;
		}
		for(size_t i = lpos_; i < rpos_; i++) {
			if(buf_[i % S].size() > 0) {
				// Analyze this column
				analyzer->analyze(refidx_, i, buf_[i % S]);
			}
			buf_[i % S].clear();
		}
	}

	/**
	 *
	 */
	void reset(uint32_t refidx, ColumnAnalyzer<T> *analyzer = NULL) {
		finalize(analyzer);
		refidx_ = refidx;
		lpos_ = 0;
		rpos_ = 0;
	}

private:
	vector<vector<T> > buf_;
	uint32_t refidx_;
	uint32_t lpos_; // in absolute coordinates
	uint32_t rpos_; // in absolute coordinates
};

/**
 * An abstract parent class for consumers of alignments.
 */
template<typename T, int S>
class AlignmentSink {
public:
	AlignmentSink() { }
	virtual ~AlignmentSink() { }
	virtual void addAlignment(const Hit& h, ColumnAnalyzer<T> *analyzer = NULL) = 0;
	virtual void finalize(ColumnAnalyzer<T> *analyzer = NULL) = 0;
	virtual void reset(uint32_t refidx, ColumnAnalyzer<T> *analyzer = NULL) = 0;
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
	RotatingCharPairAlignmentBuf() :
		AlignmentSink<pair<char, char>, S>(),
		buf_(),
		lastRefIdx_(0),
		lastRefOff_(0),
		lastLeft_(0),
		furthestRight_(0) { }

	/**
	 * Add a new alignment to the rotating buffer.
	 */
	virtual void addAlignment(const Hit& h,
	                          ColumnAnalyzer<pair<char, char> > *analyzer = NULL)
	{
		//cout << "Adding an alignment" << endl;
		if(h.h.first != lastRefIdx_) {
			assert_gt(h.h.first, lastRefIdx_);
			// Reset state and buffer for new reference sequence
			lastRefIdx_ = h.h.first;
			lastRefOff_ = 0;
			lastLeft_ = 0;
			furthestRight_ = 0;
			buf_.reset(h.h.first, analyzer);
		}
		uint32_t len = h.length();
		uint32_t rhs = h.h.second + len;
		buf_.ashAndExtendRange(h.h.second, rhs, analyzer);
		for(size_t i = 0; i < len; i++) {
			char c = -1; // signal same-as-reference
			char q = (h.fw ? h.quals[i]  : h.quals[len - i - 1]);
			if(h.mms.test(i)) {
				c = h.refcs[i];
				assert_eq(1, dna4Cat[(int)c]);
				c = charToDna5[(int)c];
				assert_lt(c, 4);
				assert_neq((int)c, (int)h.patSeq[i]);
				c = (c << 4) | (int)h.patSeq[i];
			}
			buf_.add(h.h.second + i, make_pair(c, q));
		}
	}

	/**
	 * Analyze any remaining columns in the current buffer.
	 */
	virtual void finalize(ColumnAnalyzer<pair<char, char> > *analyzer = NULL) {
		buf_.finalize(analyzer);
	}

	/**
	 * Reset the rotating buffer to deal with a new reference sequence.
	 * May involve analyzing any remaining columns in the current
	 * buffer.
	 */
	virtual void reset(uint32_t refidx,
	                   ColumnAnalyzer<pair<char, char> > *analyzer = NULL)
	{
		buf_.reset(refidx, analyzer);
	}

private:
	RotatingBuf<pair<char, char>, S> buf_;
	uint32_t lastRefIdx_;
	uint32_t lastRefOff_;
	uint32_t lastLeft_;
	uint32_t furthestRight_;
};

#endif /*ROT_BUF_H_*/
