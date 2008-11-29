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
	RotatingBuf() : buf_(), refidx_(0), lpos_(0xffffffff), rpos_(0xffffffff) {
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
		if(lpos_ != 0xffffffff) {
			if(newLhs < lpos_) {
				cerr << "Error: alignments were not presented to RotatingBuf in sorted order" << endl;
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
		if(rpos_ == 0xffffffff) {
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
		assert_neq(0xffffffff, lpos_);
		assert_neq(0xffffffff, rpos_);
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
		assert_neq(0xffffffff, lpos_);
		assert_neq(0xffffffff, rpos_);
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
		if(lpos_ != 0xffffffff) {
#ifndef NDEBUG
			assert_neq(0xffffffff, rpos_);
			size_t mn = min(lpos_ % S, rpos_ % S);
			size_t mx = max(lpos_ % S, rpos_ % S);
			if(mn == (rpos_ % S)) {
				// The range lpos_ to rpos_ wraps
				for(size_t i = mn; i < mx; i++) {
					assert_eq(0, buf_[i % S].size());
				}
			} else {
				// The range lpos_ to rpos_ doesn't wrap
				for(size_t i = 0; i < mn; i++) {
					assert_eq(0, buf_[i % S].size());
				}
				for(size_t i = mx; i < S; i++) {
					assert_eq(0, buf_[i % S].size());
				}
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
		} else {
			assert_eq(0xffffffff, rpos_);
		}
	}

	/**
	 * Reset the rotating buffer in preparation for analyzing a new
	 * reference.
	 */
	void reset(uint32_t refidx, ColumnAnalyzer<T> *analyzer = NULL) {
		finalize(analyzer);
		refidx_ = refidx;
		lpos_ = 0xffffffff;
		rpos_ = 0xffffffff;
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

	AlignmentSink(uint32_t partitionLen = 0xffffffff, bool verbose = false) :
	buf_(),
	partition_(0xffffffff),
	partitionLen_(partitionLen),
	verbose_(verbose) { }

	virtual ~AlignmentSink() { }

	/**
	 * Add a new alignment to the rotating buffer.
	 */
	virtual void addAlignment(const Hit& h,
	                          ColumnAnalyzer<T> *analyzer = NULL)
	{
		assert_gt(h.length(), 0);
		if(this->verbose_) {
			cout << "    Entered AlignmentSink::addAlignmentImpl" << endl;
			cout << "      Analyzer is " << (analyzer ? "non-null" : "null") << endl;
		}
#ifndef NDEBUG
		// Assert that the hit overlaps with this partition
		if(partitionLen_ != 0xffffffff) {
			assert_neq(0xffffffff, partition_);
			// Overlaps LHS?
			bool left  = (h.h.second <= partition_) &&
			             (h.h.second + h.length() > partition_);
			// Overlaps RHS?
			bool right = (h.h.second < (partition_ + partitionLen_)) &&
			             (h.h.second + h.length() >= partition_ + partitionLen_);
			// Occurs inside?
			bool inside = (h.h.second >= partition_) &&
			              (h.h.second + h.length() <= partition_ + partitionLen_);
			assert(left || right || inside);
		}
#endif
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
	                   uint32_t partition,
	                   ColumnAnalyzer<T> *analyzer = NULL)
	{
		partition_ = partition;
		if(partitionLen_ != 0xffffffff) {
			assert_eq(0, partition_ % partitionLen_);
		}
		this->buf_.reset(refidx, analyzer);
		resetImpl(refidx, partition, analyzer);
	}

protected:

	virtual void addAlignmentImpl(
			const Hit& h,
	        ColumnAnalyzer<T> *analyzer = NULL) { }

	virtual void finalizeImpl(
			ColumnAnalyzer<T> *analyzer = NULL) { }

	virtual void resetImpl(
			uint32_t refidx,
			uint32_t partition,
			ColumnAnalyzer<T> *analyzer = NULL) { }

	RotatingBuf<T, S> buf_;
	uint32_t partition_;
	uint32_t partitionLen_;
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
	RotatingCharPairAlignmentBuf(uint32_t partitionLen = 0xffffffff, bool verbose = false) :
		AlignmentSink<pair<char, char>, S>(partitionLen, verbose) { }

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
		uint32_t lo = 0;
		uint32_t hi = len;
		if(this->partitionLen_ != 0xffffffff) {
			uint32_t parti = this->partition_;
			uint32_t partf = this->partition_ + this->partitionLen_;
			// If LHS of alignment is before the beginning of the
			// partition, be sure to only consider those columns
			// beginning at the partition boundary.
			if(h.h.second < parti) {
				lo = parti - h.h.second;
				assert_lt(lo, h.length());
			}
			// If RHS of alignment is after the end of the partition,
			// be sure to only consider those columns up to the
			// partition boundary.
			if((h.h.second + h.length()) > partf) {
				hi -= ((h.h.second + h.length()) - partf);
			}
		}
		assert_lt(lo, hi);
		// Add char-pair evidence to each column; it's important that
		// we proceed left-to-right along the reference
		for(size_t ii = lo; ii < hi; ii++) {
			// i = offset from 5' end, ii = offset from "left" end
			// w/r/t reference
			size_t i = ii;
			if(!h.fw) i = len - ii - 1;
			if(h.mms.test(i)) {
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
