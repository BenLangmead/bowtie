/*
 * aligner_metrics.h
 */

#ifndef ALIGNER_METRICS_H_
#define ALIGNER_METRICS_H_

#include <math.h>
#include <iostream>
#include <seqan/sequence.h>
#include "alphabet.h"

using namespace std;

/**
 * Borrowed from http://www.johndcook.com/standard_deviation.html,
 * which in turn is borrowed from Knuth.
 */
class RunningStat {
public:
	RunningStat() : m_n(0) { }

	void clear() {
		m_n = 0;
	}

	void push(float x) {
		m_n++;
		// See Knuth TAOCP vol 2, 3rd edition, page 232
		if (m_n == 1) {
			m_oldM = m_newM = x;
			m_oldS = 0.0;
		} else {
			m_newM = m_oldM + (x - m_oldM)/m_n;
			m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
			// set up for next iteration
			m_oldM = m_newM;
			m_oldS = m_newS;
		}
	}

	int num() const {
		return m_n;
	}

	double mean() const {
		return (m_n > 0) ? m_newM : 0.0;
	}

	double variance() const {
		return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
	}

	double stddev() const {
		return sqrt(variance());
	}

private:
	int m_n;
	double m_oldM, m_newM, m_oldS, m_newS;
};

/**
 * Encapsulates a set of metrics that we would like an aligner to keep
 * track of, so that we can possibly use it to diagnose performance
 * issues.
 */
class AlignerMetrics {

public:

	AlignerMetrics() :
		curBacktracks_(0),
		curBwtOps_(0),
		first_(true),
		curIsLowEntropy_(false),
		curHadRanges_(false),
		reads_(0),
		lowEntReads_(0),
		hiEntReads_(0),
		alignedReads_(0),
		unalignedReads_(0),
		bwtOpsPerRead_(),
		backtracksPerRead_(),
		bwtOpsPerLoEntRead_(),
		backtracksPerLoEntRead_(),
		bwtOpsPerHiEntRead_(),
		backtracksPerHiEntRead_(),
		bwtOpsPerAlignedRead_(),
		backtracksPerAlignedRead_(),
		bwtOpsPerUnalignedRead_(),
		backtracksPerUnalignedRead_()
		{ }

	void printSummary() {
		cout << "AlignerMetrics:" << endl;
		cout << "  # Reads:         " << reads_ << endl;
		float lopct = (reads_ > 0) ? ((float)lowEntReads_/((float)lowEntReads_+hiEntReads_)) : (0.0);
		lopct *= 100.0;
		cout << "  % low-entropy:   " << (lopct) << endl;
		float unpct = (reads_ > 0) ? ((float)unalignedReads_/((float)unalignedReads_+alignedReads_)) : (0.0);
		unpct *= 100.0;
		cout << "  % unaligned:     " << (unpct) << endl;
		cout << "  BWT ops:    avg: " << bwtOpsPerRead_.mean() << ", stddev: " << bwtOpsPerRead_.stddev() << endl;
		cout << "  Backtracks: avg: " << backtracksPerRead_.mean() << ", stddev: " << backtracksPerRead_.stddev() << endl;
		cout << "  Low-entropy:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerLoEntRead_.mean() << ", stddev: " << bwtOpsPerLoEntRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerLoEntRead_.mean() << ", stddev: " << backtracksPerLoEntRead_.stddev() << endl;
		cout << "  High-entropy:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerHiEntRead_.mean() << ", stddev: " << bwtOpsPerHiEntRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerHiEntRead_.mean() << ", stddev: " << backtracksPerHiEntRead_.stddev() << endl;
		cout << "  Unaligned:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerUnalignedRead_.mean() << ", stddev: " << bwtOpsPerUnalignedRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerUnalignedRead_.mean() << ", stddev: " << backtracksPerUnalignedRead_.stddev() << endl;
		cout << "  Aligned:" << endl;
		cout << "    BWT ops:    avg: " << bwtOpsPerAlignedRead_.mean() << ", stddev: " << bwtOpsPerAlignedRead_.stddev() << endl;
		cout << "    Backtracks: avg: " << backtracksPerAlignedRead_.mean() << ", stddev: " << backtracksPerAlignedRead_.stddev() << endl;
	}

	/**
	 *
	 */
	void nextRead(const seqan::String<seqan::Dna5>& read) {
		if(!first_) {
			finishRead();
		}
		first_ = false;
		curIsLowEntropy_ = (entropyDna5(read) < 0.9f);
		curHadRanges_ = false;
		curBwtOps_ = 0;
		curBacktracks_ = 0;
	}

	/**
	 *
	 */
	void setReadHasRange() {
		curHadRanges_ = true;
	}

	/**
	 * Commit the running statistics for this read to
	 */
	void finishRead() {
		reads_++;
		if(curIsLowEntropy_) lowEntReads_++;
		else hiEntReads_++;
		if(curHadRanges_) alignedReads_++;
		else unalignedReads_++;
		bwtOpsPerRead_.push((float)curBwtOps_);
		backtracksPerRead_.push((float)curBacktracks_);
		if(curIsLowEntropy_) {
			bwtOpsPerLoEntRead_.push((float)curBwtOps_);
			backtracksPerLoEntRead_.push((float)curBacktracks_);
		} else {
			bwtOpsPerHiEntRead_.push((float)curBwtOps_);
			backtracksPerHiEntRead_.push((float)curBacktracks_);
		}
		if(curHadRanges_) {
			bwtOpsPerAlignedRead_.push((float)curBwtOps_);
			backtracksPerAlignedRead_.push((float)curBacktracks_);
		} else {
			bwtOpsPerUnalignedRead_.push((float)curBwtOps_);
			backtracksPerUnalignedRead_.push((float)curBacktracks_);
		}
	}

	// Running-total of the number of backtracks and BWT ops for the
	// current read
	uint32_t curBacktracks_;
	uint32_t curBwtOps_;

protected:

	bool first_;

	// true iff the current read is low entropy
	bool curIsLowEntropy_;
	// true iff the current read has had one or more ranges reported
	bool curHadRanges_;

	// # reads
	uint32_t reads_;
	// # low-entropy reads
	uint32_t lowEntReads_;
	// # high-entropy reads
	uint32_t hiEntReads_;
	// # reads with alignments
	uint32_t alignedReads_;
	// # reads without alignments
	uint32_t unalignedReads_;

	// Distribution of BWT operations per read
	RunningStat bwtOpsPerRead_;
	RunningStat backtracksPerRead_;

	// Distribution of BWT operations per low-entropy read
	RunningStat bwtOpsPerLoEntRead_;
	RunningStat backtracksPerLoEntRead_;

	// Distribution of BWT operations per high-entropy read
	RunningStat bwtOpsPerHiEntRead_;
	RunningStat backtracksPerHiEntRead_;

	// Distribution of BWT operations per read that "aligned" (for
	// which a range was arrived at - range may not have necessarily
	// lead to an alignment)
	RunningStat bwtOpsPerAlignedRead_;
	RunningStat backtracksPerAlignedRead_;

	// Distribution of BWT operations per read that didn't align
	RunningStat bwtOpsPerUnalignedRead_;
	RunningStat backtracksPerUnalignedRead_;
};

#endif /* ALIGNER_METRICS_H_ */
