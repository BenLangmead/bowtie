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
 * Encapsulates a summarized "SNP mass" view of the evidence in a
 * column of the multiple alignment.
 */
struct SNPMass {
	SNPMass() : totalMass(0), totalNonRefMass(0) {
		perCharMass[0] = perCharMass[1] = perCharMass[2] = perCharMass[3] = 0;
		perCharCoverage[0] = perCharCoverage[1] = perCharCoverage[2] = perCharCoverage[3] = 0;
		perCharAmbCoverage[0] = perCharAmbCoverage[1] = perCharAmbCoverage[2] = perCharAmbCoverage[3] = 0;
	}
	uint8_t refChar;      // reference character at this position
	uint32_t totalMass;   // total Phred mass covering this position
	uint32_t totalNonRefMass; // total Phtred mass from characters other than refChar covering this position
	uint32_t perCharMass[4];     // Phred mass supporting A/C/G/T
	uint32_t perCharCoverage[4]; // # reads supporting A/C/G/T
	uint32_t perCharAmbCoverage[4]; // # ambiguous reads supporting A/C/G/T
};

/**
 * A concrete column-analysis class where each column is a vector
 * containing a single pre-reduced summary of "mass" information in the
 * column.  Mass information is generally based on the bases in the
 * column as a function of their quality values and the mapping quality
 * of the alignment they came from.
 */
class SNPColumnMassAnalyzer : public ColumnAnalyzer<SNPMass> {
public:
	SNPColumnMassAnalyzer(ostream& out) :
		ColumnAnalyzer<SNPMass>(), out_(out) { }

	virtual ~SNPColumnMassAnalyzer() { }

	/**
	 * Analyze a column with a single SNPMass and output a homozygous
	 * SNP call where evidence is sufficient.
	 */
	virtual void analyze(uint32_t refidx,
	                     uint32_t refoff,
	                     const vector<SNPMass>& col)
	{
		assert_eq(1, col.size());
		const SNPMass& m = col.front();
		const uint32_t totalMassThresh = m.totalMass >> 2;
		// If the total mass exceeds a sufficient threshold
		if(m.totalMass > 120 && m.totalNonRefMass > totalMassThresh) {
			// Sufficient evidence to try to make a call
//			uint8_t maxi = 0xff;
//			uint8_t maxi2 = 0xff;
//			for(uint8_t i = 0; i < 4; i++) {
//				// Find the largest per-A/C/G/T mass
//				if(m.perCharMass[i] > maxSnpMass) {
//					maxSnpMass = m.perCharMass[i];
//					maxi = i;
//				} else if(m.perCharMass[i] > maxSnpMass) {
//
//				}
//			}
//			if(maxi != 0xff && maxi != m.refChar) {
//				// Take ratio of highest-mass character to the total
//				// mass
//				double ratio = m.perCharMass[maxi] / (double)m.totalMass;
//				if(ratio > 0.85) {
//					// 85% is the threshold
//					out_ << refidx << " " << refoff << " "
//					     << "ACGT"[m.refChar] << " " << "ACGT"[maxi] << endl;
//				}
//			}

			// Identify reference position
			out_ << refidx << " " << refoff << " " << "ACGT"[m.refChar] << " ";
			for(int i = 0; i < 4; i++) {
				if(m.perCharCoverage[i] == 0) {
					assert_eq(0, m.perCharAmbCoverage[i]);
					assert_eq(0, m.perCharMass[i]);
					out_ << "0";
				} else {
					out_ << m.perCharCoverage[i] << ","
						 << m.perCharAmbCoverage[i] << ","
						 << m.perCharMass[i];
				}
				if(i < 3) {
					out_ << " ";
				} else {
					out_ << endl;
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
	RotatingBuf() :
		buf_(), refidx_(0), lpos_(0xffffffff), rpos_(0xffffffff)
	{
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
	vector<T>& get(uint32_t pos) {
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
 * Abstract parent class of all AlignmentSinks that respect partition
 * boundaries.
 */
template<typename T, int S>
class RotatingPartitionedAlignmentBuf : public AlignmentSink<T, S> {

public:
	RotatingPartitionedAlignmentBuf(uint32_t partitionLen = 0xffffffff, bool verbose = false) :
		AlignmentSink<T, S>(partitionLen, verbose) { }

protected:

	/**
	 * Add a new alignment to the rotating buffer.
	 */
	virtual void addAlignmentImpl(const Hit& h,
	                              ColumnAnalyzer<T> *analyzer = NULL)
	{
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
		this->addAlignmentCols(h, lo, hi, analyzer);
	}

	/**
	 *
	 */
	virtual void addAlignmentCols(const Hit& h, uint32_t lo, uint32_t hi,
	                              ColumnAnalyzer<T> *analyzer) = 0;
};

/**
 * An alignment sink that considers all columns of an alignment that
 * lie within a partition and factors the evidence into a "SNPMass"
 * object.
 */
template<int S>
class RotatingSNPMassAlignmentBuf : public RotatingPartitionedAlignmentBuf<SNPMass, S> {

public:
	RotatingSNPMassAlignmentBuf(uint32_t partitionLen = 0xffffffff, bool verbose = false) :
		RotatingPartitionedAlignmentBuf<SNPMass, S>(partitionLen, verbose) { }

protected:

	/**
	 * Add a new alignment to the rotating buffer.
	 */
	virtual void addAlignmentCols(const Hit& h, uint32_t lo, uint32_t hi,
	                              ColumnAnalyzer<SNPMass> *analyzer)
	{
		uint32_t len = h.length();
		for(size_t ii = lo; ii < hi; ii++) {
			// i = offset from 5' end, ii = offset from "left" end
			// w/r/t reference
			int readc = (int)h.patSeq[ii];
			if(readc == 4) continue; // no evidence inherent in Ns
			size_t i = ii;
			if(!h.fw) {
				i = len - ii - 1;
			}
			char q = h.quals[ii];
			if(h.oms > 0) {
				// Divide Phred mass by the number of other elements in
				// the Burrows-Wheeler range
				q /= (h.oms+1);
			}
			vector<SNPMass>& col = this->buf_.get(h.h.second + ii);
			bool snp = h.mms.test(i);
			char c;
			if(snp) {
				c = h.refcs[i]; // refcs also indexed from 5' end of read
				assert_eq(1, dna4Cat[(int)c]);
				c = charToDna5[(int)c];
				assert_neq((int)c, readc);
			} else {
				c = readc;
			}
			// 'c' is the character we observe
			assert_lt(c, 4);
			if(col.empty()) {
				col.resize(1);
				// Add the Phred qualiy of this position to the total
				// mass
				col.front().totalMass = q;
				if(snp) {
					col.front().totalNonRefMass = q;
				} else {
					col.front().totalNonRefMass = 0;
				}
				col.front().perCharMass[readc] = q;
				col.front().perCharCoverage[readc] = 1;
				col.front().perCharAmbCoverage[readc] = (h.oms > 0) ? 1 : 0;
				col.front().refChar = c;
			} else {
				assert_eq(1, col.size());
				col.front().totalMass += q;
				if(snp) {
					col.front().totalNonRefMass += q;
				}
				col.front().perCharMass[readc] += q;
				col.front().perCharCoverage[readc]++;
				if(h.oms > 0) {
					col.front().perCharAmbCoverage[readc]++;
				}
				assert_eq(col.front().refChar, c);
			}
			// Done
		}
	}
};

#endif /*ROT_BUF_H_*/
