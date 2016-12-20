#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "pat.h"
#include "filebuf.h"

using namespace std;
using namespace seqan;

/**
 * Calculate a per-read random seed based on a combination of
 * the read data (incl. sequence, name, quals) and the global
 * seed in '_randSeed'.
 */
static uint32_t genRandSeed(
	const String<Dna5>& qry,
    const String<char>& qual,
    const String<char>& name,
    uint32_t seed)
{
	// Calculate a per-read random seed based on a combination of
	// the read data (incl. sequence, name, quals) and the global
	// seed
	uint32_t rseed = (seed + 101) * 59 * 61 * 67 * 71 * 73 * 79 * 83;
	size_t qlen = seqan::length(qry);
	// Throw all the characters of the read into the random seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qry[i];
		assert_leq(p, 4);
		size_t off = ((i & 15) << 1);
		rseed ^= (p << off);
	}
	// Throw all the quality values for the read into the random
	// seed
	for(size_t i = 0; i < qlen; i++) {
		int p = (int)qual[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed ^= (p << off);
	}
	// Throw all the characters in the read name into the random
	// seed
	size_t namelen = seqan::length(name);
	for(size_t i = 0; i < namelen; i++) {
		int p = (int)name[i];
		assert_leq(p, 255);
		size_t off = ((i & 3) << 3);
		rseed ^= (p << off);
	}
	return rseed;
}

/**
 * Once name/sequence/qualities have been parsed for an
 * unpaired read, set all the other key fields of the Read
 * struct.
 */
void PatternSourcePerThread::finalize(Read& ra) {
	ra.mate = 0;
	ra.constructRevComps();
	ra.constructReverses();
	ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
}

/**
 * Once name/sequence/qualities have been parsed for a
 * paired-end read, set all the other key fields of the Read
 * structs.
 */
void PatternSourcePerThread::finalizePair(Read& ra, Read& rb) {
	ra.mate = 1;
	ra.constructRevComps();
	ra.constructReverses();
	ra.fixMateName(1);
	ra.seed = genRandSeed(ra.patFw, ra.qual, ra.name, seed_);
	
	rb.mate = 2;
	rb.constructRevComps();
	rb.constructReverses();
	rb.fixMateName(2);
	rb.seed = genRandSeed(rb.patFw, rb.qual, rb.name, seed_);
}

/**
 * Get the next paired or unpaired read from the wrapped
 * PatternComposer.  Returns a pair of bools; first indicates
 * whether we were successful, second indicates whether we're
 * done.
 */
pair<bool, bool> PatternSourcePerThread::nextReadPair() {
	// Prepare batch
	if(buf_.exhausted()) {
		pair<bool, int> res = nextBatch();
		if(res.first && res.second == 0) {
			return make_pair(false, true);
		}
		last_batch_ = res.first;
		last_batch_size_ = res.second;
		assert_eq(0, buf_.cur_buf_);
	} else {
		buf_.next(); // advance cursor
		assert_gt(buf_.cur_buf_, 0);
	}
	bool this_is_last = buf_.cur_buf_ == last_batch_size_-1;
	if(buf_.rdid() < skip_) {
		return make_pair(false, this_is_last ? last_batch_ : false);
	}
	// Parse read/pair
	assert(strlen(buf_.read_a().readOrigBuf) != 0);
	assert(buf_.read_a().empty());
	if(!parse(buf_.read_a(), buf_.read_b())) {
		return make_pair(false, false);
	}
	// Finalize read/pair
	if(paired()) {
		finalizePair(buf_.read_a(), buf_.read_b());
	} else {
		finalize(buf_.read_a());
	}
	return make_pair(true, this_is_last ? last_batch_ : false);
}

/**
 * The main member function for dispensing pairs of reads or
 * singleton reads.  Returns true iff ra and rb contain a new
 * pair; returns false if ra contains a new unpaired read.
 */
pair<bool, int> SoloPatternComposer::nextBatch(PerThreadReadBuf& pt) {
	uint32_t cur = cur_;
	while(cur < src_.size()) {
		// Patterns from srca_[cur_] are unpaired
		pair<bool, int> res;
		do {
			res = src_[cur]->nextBatch(
				pt,
				true,  // batch A (or pairs)
				true); // grab lock below
		} while(!res.first && res.second == 0);
		if(res.second == 0) {
			ThreadSafe ts(&mutex_m);
			if(cur + 1 > cur_) {
				cur_++;
			}
			cur = cur_;
			continue; // on to next pair of PatternSources
		}
		return res;
	}
	assert_leq(cur, src_.size());
	return make_pair(true, 0);
}

/**
 * The main member function for dispensing pairs of reads or
 * singleton reads.  Returns true iff ra and rb contain a new
 * pair; returns false if ra contains a new unpaired read.
 */
pair<bool, int> DualPatternComposer::nextBatch(PerThreadReadBuf& pt) {
	// 'cur' indexes the current pair of PatternSources
	uint32_t cur = cur_;
	while(cur < srca_.size()) {
		if(srcb_[cur] == NULL) {
			pair<bool, int> res = srca_[cur]->nextBatch(
				pt,
				true,  // batch A (or pairs)
				true); // grab lock below
			bool done = res.first;
			if(!done && res.second == 0) {
				ThreadSafe ts(&mutex_m);
				if(cur + 1 > cur_) cur_++;
				cur = cur_; // Move on to next PatternSource
				continue; // on to next pair of PatternSources
			}
			return make_pair(done, res.second);
		} else {
			pair<bool, int> resa, resb;
			// Lock to ensure that this thread gets parallel reads
			// in the two mate files
			{
				ThreadSafe ts(&mutex_m);
				resa = srca_[cur]->nextBatch(
					pt,
					true,   // batch A
					false); // don't grab lock below
				resb = srcb_[cur]->nextBatch(
					pt,
					false,  // batch B
					false); // don't grab lock below
				assert_eq(srca_[cur]->readCount(),
					  srcb_[cur]->readCount());
			}
			if(resa.second < resb.second) {
				cerr << "Error, fewer reads in file specified with -1 "
					 << "than in file specified with -2" << endl;
				throw 1;
			} else if(resa.second == 0 && resb.second == 0) {
				ThreadSafe ts(&mutex_m);
				if(cur + 1 > cur_) {
					cur_++;
				}
				cur = cur_; // Move on to next PatternSource
				continue; // on to next pair of PatternSources
			} else if(resb.second < resa.second) {
				cerr << "Error, fewer reads in file specified with -2 "
					 << "than in file specified with -1" << endl;
				throw 1;
			}
			assert_eq(resa.first, resb.first);
			assert_eq(resa.second, resb.second);
			return make_pair(resa.first, resa.second);
		}
	}
	assert_leq(cur, srca_.size());
	return make_pair(true, 0);
}

/**
 * Fill Read with the sequence, quality and name for the next
 * read in the list of read files.  This function gets called by
 * all the search threads, so we must handle synchronization.
 *
 * Returns pair<bool, int> where bool indicates whether we're
 * completely done, and int indicates how many reads were read.
 */
pair<bool, int> CFilePatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	bool done = false;
	int nread = 0;
	
	// synchronization at this level because both reading and manipulation of
	// current file pointer have to be protected
	ThreadSafe ts(&mutex, lock);
	pt.setReadId(readCnt_);
	while(true) { // loop that moves on to next file when needed
		do {
			pair<bool, int> ret = nextBatchFromFile(pt, batch_a);
			done = ret.first;
			nread = ret.second;
		} while(!done && nread == 0); // not sure why this would happen
		if(done && filecur_ < infiles_.size()) { // finished with this file
			open();
			resetForNextFile(); // reset state to handle a fresh file
			filecur_++;
			if(nread == 0) {
				continue;
			}
		}
		break;
	}
	assert_geq(nread, 0);
	readCnt_ += nread;
	return make_pair(done, nread);
}

/**
 * Open the next file in the list of input files.
 */
void CFilePatternSource::open() {
	if(is_open_) {
		is_open_ = false;
		fclose(fp_);
		fp_ = NULL;
		if(qfp_ != NULL) {
			fclose(qfp_);
			qfp_ = NULL;
		}
	}
	while(filecur_ < infiles_.size()) {
		// Open read
		if(infiles_[filecur_] == "-") {
			fp_ = stdin;
		} else if((fp_ = fopen(infiles_[filecur_].c_str(), "rb")) == NULL) {
			if(!errs_[filecur_]) {
				cerr << "Warning: Could not open read file \""
				     << infiles_[filecur_] << "\" for reading; skipping..."
					 << endl;
				errs_[filecur_] = true;
			}
			filecur_++;
			continue;
		}
		is_open_ = true;
		setvbuf(fp_, buf_, _IOFBF, 64*1024);
		if(!qinfiles_.empty()) {
			if(qinfiles_[filecur_] == "-") {
				qfp_ = stdin;
			} else if((qfp_ = fopen(qinfiles_[filecur_].c_str(), "rb")) == NULL) {
				if(!errs_[filecur_]) {
					cerr << "Warning: Could not open quality file \""
					     << qinfiles_[filecur_] << "\" for reading; skipping..."
						 << endl;
					errs_[filecur_] = true;
				}
				filecur_++;
				continue;
			}
			assert(qfp_ != NULL);
			setvbuf(qfp_, qbuf_, _IOFBF, 64*1024);
		}
		return;
	}
	throw 1;
}

/**
 * Constructor for vector pattern source, used when the user has
 * specified the input strings on the command line using the -c
 * option.
 */
VectorPatternSource::VectorPatternSource(
	const vector<string>& seqs,
	bool color,
	const char *dumpfile,
	int trim3,
	int trim5) :
	TrimmingPatternSource(dumpfile, trim3, trim5),
	color_(color),
	cur_(0),
	paired_(false),
	tokbuf_(),
	bufs_()
{
	// Install sequences in buffers, ready for immediate copying in
	// nextBatch().  Formatting of the buffer is just like
	// TabbedPatternSource.
	const size_t seqslen = seqs.size();
	for(size_t i = 0; i < seqslen; i++) {
		tokbuf_.clear();
		tokenize(seqs[i], ":", tokbuf_, 2);
		assert_gt(tokbuf_.size(), 0);
		assert_leq(tokbuf_.size(), 2);
		// Get another buffer ready
		bufs_.resize(bufs_.size()+1);
		bufs_.back().clear();
		// Install name
		itoa10<TReadId>(static_cast<TReadId>(i), nametmp_);
		bufs_.back() = nametmp_;
		bufs_.back().push_back('\t');
		// Install sequence
		bufs_.back().append(tokbuf_[0].c_str());
		bufs_.back().push_back('\t');
		// Install qualities
		if(tokbuf_.size() > 1) {
			bufs_.back().append(tokbuf_[1].c_str());
		} else {
			const size_t len = tokbuf_[0].length();
			for(size_t i = 0; i < len; i++) {
				bufs_.back().push_back('I');
			}
		}
	}
}

/**
 * Read next batch.  However, batch concept is not very applicable for this
 * PatternSource where all the info has already been parsed into the fields
 * in the contsructor.  This essentially modifies the pt as though we read
 * in some number of patterns.
 */
pair<bool, int> VectorPatternSource::nextBatch(
	PerThreadReadBuf& pt,
	bool batch_a,
	bool lock)
{
	ThreadSafe ts(&mutex, lock);
	pt.setReadId(cur_);
	vector<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	size_t readi = 0;
	for(; readi < pt.max_buf_ && cur_ < bufs_.size(); readi++, cur_++) {
		const size_t len = bufs_[cur_].length();
		memcpy(readbuf[readi].readOrigBuf, bufs_[cur_].c_str(), len);
		readbuf[readi].readOrigBufLen = len;
	}
	readCnt_ += readi;
	return make_pair(cur_ == bufs_.size(), readi);
}

/**
 * Finishes parsing outside the critical section.
 */
bool VectorPatternSource::parse(Read& ra, Read& rb, TReadId rdid) const {
	// Very similar to TabbedPatternSource

	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(ra.empty());
	assert_gt(ra.readOrigBufLen, 0); // raw data for read/pair is here
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = ra.readOrigBufLen;
	
	// Loop over the two ends
	for(int endi = 0; endi < 2 && c == '\t'; endi++) {
		Read& r = ((endi == 0) ? ra : rb);
		assert_eq(0, seqan::length(r.name));
		// Parse name if (a) this is the first end, or
		// (b) this is tab6
		size_t nameoff = 0;
		if(endi < 1 || paired_) {
			// Parse read name
			c = ra.readOrigBuf[cur++];
			while(c != '\t' && cur < buflen) {
				r.nameBuf[nameoff++] = c;
				c = ra.readOrigBuf[cur++];
			}
			assert_eq('\t', c);
			if(cur >= buflen) {
				return false; // record ended prematurely
			}
		} else if(endi > 0) {
			// if this is the second end and we're parsing
			// tab5, copy name from first end
			rb.name = ra.name;
		}
		r.nameBuf[nameoff] = '\0';
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, nameoff);

		// Parse sequence
		assert_eq(0, seqan::length(r.patFw));
		c = ra.readOrigBuf[cur++];
		int nchar = 0, seqoff = 0;
		if(color_ && asc2dnacat[c] > 0) {
			// First char is a DNA char (primer)
			if(asc2colcat[toupper(r.readOrigBuf[cur++])] <= 0) {
				// 2nd char isn't a color, so don't assume 'c' is primer
				cur -= 2;
			} else {
				// 'c' is primer
				r.primer = c;
			}
			c = r.readOrigBuf[cur++];
		}
		if(color_) {
			while(c != '\t' && cur < buflen) {
				if(c >= '0' && c < '4') {
					c = "ACGTN"[(int)c - '0'];
				}
				if(c == '.') {
					c = 'N';
				}
				if(isalpha(c)) {
					assert_in(toupper(c), "ACGTN");
					if(nchar++ >= this->trim5_) {
						assert_neq(0, asc2dnacat[c]);
						r.patBufFw[seqoff++] = charToDna5[c]; // ascii to int
					}
				}
				c = ra.readOrigBuf[cur++];
			}
			ra.color = true;
		} else {
			while(c != '\t' && cur < buflen) {
				if(isalpha(c)) {
					assert_in(toupper(c), "ACGTN");
					if(nchar++ >= this->trim5_) {
						assert_neq(0, asc2dnacat[c]);
						r.patBufFw[seqoff++] = charToDna5[c]; // ascii to int
					}
				}
				c = ra.readOrigBuf[cur++];
			}
		}
		assert_eq('\t', c);
		if(cur >= buflen) {
			return false; // record ended prematurely
		}
		r.patBufFw[seqoff] = '\0';
		_setBegin(r.patFw, (Dna5*)r.patBufFw);
		// record amt trimmed from 5' end due to --trim5
		r.trimmed5 = (int)(nchar - seqoff);
		// record amt trimmed from 3' end due to --trim3
		int trim3 = (seqoff < this->trim3_) ? seqoff : this->trim3_;
		_setLength(r.patFw, seqoff - trim3);
		r.trimmed3 = trim3;
		
		// Parse qualities
		assert_eq(0, seqan::length(r.qual));
		c = ra.readOrigBuf[cur++];
		int nqual = 0;
		size_t qualoff = 0;
		while(c != '\t' && c != '\n' && c != '\r') {
			if(c == ' ') {
				wrongQualityFormat(r.name);
				return false;
			}
			char cadd = charToPhred33(c, false, false);
			if(++nqual > this->trim5_) {
				r.qualBuf[qualoff++] = cadd;
			}
			if(cur >= buflen) break;
			c = ra.readOrigBuf[cur++];
		}
		if(nchar > nqual) {
			tooFewQualities(r.name);
			return false;
		} else if(nqual > nchar) {
			tooManyQualities(r.name);
			return false;
		}
		r.qualBuf[seqoff] = '\0';
		_setBegin(r.qual, r.qualBuf);
		_setLength(r.qual, seqan::length(r.patFw));
		assert(c == '\t' || c == '\n' || c == '\r' || cur >= buflen);
	}
	ra.parsed = true;
	if(!rb.parsed && rb.readOrigBufLen > 0) {
		return parse(rb, ra, rdid);
	}
	return true;
}

/**
 * Light-parse a FASTA batch into the given buffer.
 */
pair<bool, int> FastaPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	int c;
	vector<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	if(first_) {
		c = getc_unlocked(fp_);
		while(c == '\r' || c == '\n') {
			c = getc_unlocked(fp_);
		}
		if(c != '>') {
			cerr << "Error: reads file does not look like a FASTA file" << endl;
			throw 1;
		}
		first_ = false;
	}
	bool done = false;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && !done; readi++) {
		readbuf[readi].readOrigBuf[0] = '>';
		size_t bufoff = 1;
		while(true) {
			c = getc_unlocked(fp_);
			if(c < 0 || c == '>') {
				done = c < 0;
				break;
			}
			readbuf[readi].readOrigBuf[bufoff++] = c;
		}
		readbuf[readi].readOrigBufLen = bufoff;
	}
	return make_pair(done, readi);
}

/**
 * Finalize FASTA parsing outside critical section.
 */
bool FastaPatternSource::parse(Read& r, Read& rb, TReadId rdid) const {
	// We assume the light parser has put the raw data for the separate ends
	// into separate Read objects.  That doesn't have to be the case, but
	// that's how we've chosen to do it for FastqPatternSource
	assert_gt(r.readOrigBufLen, 0);
	assert(r.empty());
	int c = -1;
	size_t cur = 1;
	const size_t buflen = r.readOrigBufLen;
	
	// Parse read name
	assert_eq(0, seqan::length(r.name));
	int nameoff = 0;
	while(cur < buflen) {
		c = r.readOrigBuf[cur++];
		if(c == '\n' || c == '\r') {
			do {
				c = r.readOrigBuf[cur++];
			} while((c == '\n' || c == '\r') && cur < buflen);
			break;
		}
		r.nameBuf[nameoff++] = c;
	}
	if(cur >= buflen) {
		return false; // FASTA ended prematurely
	}
	if(nameoff > 0) {
		r.nameBuf[nameoff] = '\0';
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, nameoff);
	}
	
	// Parse sequence
	int nchar = 0, seqoff = 0;
	assert_eq(0, seqan::length(r.patFw));
	assert(c != '\n' && c != '\r');
	assert_lt(cur, buflen);

	if(color_ && asc2dnacat[c] > 0) {
		// First char is a DNA char (primer)
		if(asc2colcat[toupper(r.readOrigBuf[cur++])] <= 0) {
			// 2nd char isn't a color, so don't assume 'c' is primer
			cur -= 2;
		} else {
			// 'c' is primer
			r.primer = c;
		}
		c = r.readOrigBuf[cur++];
	}
	if(color_) {
		while(c != '\n' && cur < buflen) {
			if(c >= '0' && c < '4') {
				c = "ACGTN"[(int)c - '0'];
			}
			if(c == '.') {
				c = 'N';
			}
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(nchar++ >= this->trim5_) {
					assert_neq(0, asc2dnacat[c]);
					r.patBufFw[seqoff++] = charToDna5[c]; // ascii to int
				}
			}
			c = r.readOrigBuf[cur++];
		}
		r.color = true;
	} else {
		while(c != '\n' && cur < buflen) {
			if(c == '.') {
				c = 'N';
			}
			if(isalpha(c)) {
				// If it's past the 5'-end trim point
				if(nchar++ >= this->trim5_) {
					r.patBufFw[seqoff++] = charToDna5[c];
				}
			}
			assert_lt(cur, buflen);
			c = r.readOrigBuf[cur++];
		}
	}
	r.patBufFw[seqoff] = '\0';
	_setBegin(r.patFw, (Dna5*)r.patBufFw);
	// record amt trimmed from 5' end due to --trim5
	r.trimmed5 = (int)(nchar - seqoff);
	// record amt trimmed from 3' end due to --trim3
	int trim3 = (seqoff < this->trim3_) ? seqoff : this->trim3_;
	_setLength(r.patFw, seqoff - trim3);
	r.trimmed3 = trim3;
	
	for(size_t i = 0; i < seqoff - trim3; i++) {
		r.qualBuf[i] = 'I';
	}
	r.qualBuf[seqoff - trim3] = '\0';
	_setBegin(r.qual, r.qualBuf);
	_setLength(r.qual, seqoff - trim3);

	// Set up a default name if one hasn't been set
	if(nameoff == 0) {
		itoa10<TReadId>(static_cast<TReadId>(rdid), r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, strlen(r.nameBuf));
	}
	r.parsed = true;
	if(!rb.parsed && rb.readOrigBufLen > 0) {
		return parse(rb, r, rdid);
	}
	return true;
}

/**
 * Light-parse a FASTA-continuous batch into the given buffer.
 * This is trickier for FASTA-continuous than for other formats,
 * for several reasons:
 *
 * 1. Reads are substrings of a longer FASTA input string
 * 2. Reads may overlap w/r/t the longer FASTA string
 * 3. Read names depend on the most recently observed FASTA
 *    record name
 */
pair<bool, int> FastaContinuousPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	int c = -1;
	vector<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	size_t readi = 0;
	int nameoff = 0;
	while(readi < pt.max_buf_) {
		c = getc_unlocked(fp_);
		if(c < 0) {
			break;
		}
		if(c == '>') {
			resetForNextFile();
			nameoff = 0;
			c = getc_unlocked(fp_);
			bool sawSpace = false;
			while(c != '\n' && c != '\r') {
				if(!sawSpace) {
					sawSpace = isspace(c);
				}
				if(!sawSpace) {
					// Put it in the name prefix buffer so we
					// can re-use this prefix for all the reads
					// that are substrings of this FASTA sequence
					name_prefix_buf_[nameoff++] = c;
				}
				c = getc_unlocked(fp_);
			}
			while(c == '\n' || c == '\r') {
				c = getc_unlocked(fp_);
			}
			if(c < 0) {
				break;
			}
			name_prefix_buf_[nameoff++] = '_';
		}
		int cat = asc2dnacat[c];
		if(cat >= 2) c = 'N';
		if(cat == 0) {
			// Non-DNA, non-IUPAC char; skip
			continue;
		} else {
			// DNA char
			buf_[bufCur_++] = c;
			if(bufCur_ == 1024) {
				bufCur_ = 0; // wrap around circular buf
			}
			if(eat_ > 0) {
				eat_--;
				// Try to keep readCnt_ aligned with the offset
				// into the reference; that lets us see where
				// the sampling gaps are by looking at the read
				// name
				if(!beginning_) {
					readCnt_++;
				}
				continue;
			}
			// install name
			memcpy(readbuf[readi].readOrigBuf, name_prefix_buf_, nameoff);
			itoa10<TReadId>(readCnt_ - subReadCnt_, name_int_buf_);
			size_t bufoff = nameoff;
			memcpy(readbuf[readi].readOrigBuf + bufoff, name_int_buf_, strlen(name_int_buf_));
			bufoff += strlen(name_int_buf_);
			readbuf[readi].readOrigBuf[bufoff++] = '\t';
			// install sequence
			for(size_t i = 0; i < length_; i++) {
				if(length_ - i <= bufCur_) {
					c = buf_[bufCur_ - (length_ - i)];
				} else {
					// Rotate
					c = buf_[bufCur_ - (length_ - i) + 1024];
				}
				readbuf[readi].readOrigBuf[bufoff++] = c;
			}
			readbuf[readi].readOrigBufLen = bufoff;
			eat_ = freq_-1;
			readCnt_++;
			beginning_ = false;
			readi++;
		}
	}
	return make_pair(c < 0, readi);
}

/**
 * Finalize FASTA-continuous parsing outside critical section.
 */
bool FastaContinuousPatternSource::parse(
	Read& ra,
	Read& rb,
	TReadId rdid) const
{
	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(ra.empty());
	assert(rb.empty());
	assert_gt(ra.readOrigBufLen, 0); // raw data for read/pair is here
	assert_eq(0, rb.readOrigBufLen);
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = ra.readOrigBufLen;
	
	// Parse read name
	assert_eq(0, seqan::length(ra.name));
	int nameoff = 0;
	c = ra.readOrigBuf[cur++];
	while(c != '\t' && cur < buflen) {
		ra.nameBuf[nameoff++] = c;
		c = ra.readOrigBuf[cur++];
	}
	assert_eq('\t', c);
	if(cur >= buflen) {
		return false; // record ended prematurely
	}
	ra.nameBuf[nameoff] = '\0';
	_setBegin(ra.name, ra.nameBuf);
	_setLength(ra.name, nameoff);

	// Parse sequence
	assert_eq(0, seqan::length(ra.patFw));
	c = ra.readOrigBuf[cur++];
	int nchar = 0, seqoff = 0;
	while(cur < buflen) {
		if(isalpha(c)) {
			assert_in(toupper(c), "ACGTN");
			if(nchar++ >= this->trim5_) {
				assert_neq(0, asc2dnacat[c]);
				ra.patBufFw[seqoff++] = charToDna5[c]; // ascii to int
			}
		}
		c = ra.readOrigBuf[cur++];
	}
	ra.patBufFw[seqoff] = '\0';
	_setBegin(ra.patFw, (Dna5*)ra.patBufFw);
	// record amt trimmed from 5' end due to --trim5
	ra.trimmed5 = (int)(nchar - seqoff);
	// record amt trimmed from 3' end due to --trim3
	int trim3 = (seqoff < this->trim3_) ? seqoff : this->trim3_;
	_setLength(ra.patFw, seqoff - trim3);
	ra.trimmed3 = trim3;
	
	// Make fake qualities
	assert_eq(0, seqan::length(ra.qual));
	int qualoff = 0;
	for(size_t i = 0; i < seqoff; i++) {
		ra.qualBuf[qualoff++] = 'I';
	}
	ra.qualBuf[qualoff] = '\0';
	_setBegin(ra.qual, ra.qualBuf);
	_setLength(ra.qual, qualoff);
	ra.parsed = true;
	return true;
}

/**
 * "Light" parser.  This is inside the critical section, so the key is to do
 * just enough parsing so that another function downstream (finalize()) can do
 * the rest of the parsing.  Really this function's only job is to stick every
 * for lines worth of the input file into a buffer (r.readOrigBuf).  finalize()
 * then parses the contents of r.readOrigBuf later.
 */
pair<bool, int> FastqPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	int c = 0;
	vector<Read>& readBuf = batch_a ? pt.bufa_ : pt.bufb_;
	if(first_) {
		c = getc_unlocked(fp_);
		while(c == '\r' || c == '\n') {
			c = getc_unlocked(fp_);
		}
		if(c != '@') {
			cerr << "Error: reads file does not look like a FASTQ file" << endl;
			throw 1;
		}
		first_ = false;
		readBuf[0].readOrigBuf[0] = c;
		readBuf[0].readOrigBufLen = 1;
	}
	bool done = false, aborted = false;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && !done; readi++) {
		char* buf = readBuf[readi].readOrigBuf;
		assert(readi == 0 || readBuf[readi].readOrigBufLen == 0);
		int newlines = 4;
		while(newlines) {
			c = getc_unlocked(fp_);
			done = c < 0;
			if(c == '\n' || (done && newlines == 1)) {
				// Saw newline, or EOF that we're
				// interpreting as final newline
				newlines--;
				c = '\n';
			} else if(done) {
				aborted = true; // Unexpected EOF
				break;
			}
			buf[readBuf[readi].readOrigBufLen++] = c;
		}
	}
	if(aborted) {
		readi--;
	}
	return make_pair(done, readi);
}

/**
 * Finalize FASTQ parsing outside critical section.
 */
bool FastqPatternSource::parse(Read &r, Read& rb, TReadId rdid) const {
	// We assume the light parser has put the raw data for the separate ends
	// into separate Read objects.  That doesn't have to be the case, but
	// that's how we've chosen to do it for FastqPatternSource
	assert_gt(r.readOrigBufLen, 0);
	assert(r.empty());
	int c;
	size_t cur = 1;
	const size_t buflen = r.readOrigBufLen;

	// Parse read name
	assert_eq(0, seqan::length(r.name));
	int nameoff = 0;
	while(true) {
		assert_lt(cur, buflen);
		c = r.readOrigBuf[cur++];
		if(c == '\n' || c == '\r') {
			do {
				c = r.readOrigBuf[cur++];
			} while(c == '\n' || c == '\r');
			break;
		}
		r.nameBuf[nameoff++] = c;
	}
	r.nameBuf[nameoff] = '\0';
	_setBegin(r.name, r.nameBuf);
	_setLength(r.name, nameoff);
	
	// Parse sequence
	int nchar = 0, seqoff = 0;
	assert_eq(0, seqan::length(r.patFw));
	if(color_ && asc2dnacat[c] > 0) {
		// First char is a DNA char (primer)
		if(asc2colcat[toupper(r.readOrigBuf[cur++])] <= 0) {
			// 2nd char isn't a color, so don't assume 'c' is primer
			cur -= 2;
		} else {
			// 'c' is primer
			r.primer = c;
		}
		c = r.readOrigBuf[cur++];
	}
	if(color_) {
		while(c != '+' && cur < buflen) {
			if(c >= '0' && c < '4') {
				c = "ACGTN"[(int)c - '0'];
			}
			if(c == '.') {
				c = 'N';
			}
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(nchar++ >= this->trim5_) {
					assert_neq(0, asc2dnacat[c]);
					r.patBufFw[seqoff++] = charToDna5[c]; // ascii to int
				}
			}
			c = r.readOrigBuf[cur++];
		}
		r.color = true;
	} else {
		while(c != '+' && cur < buflen) {
			if(c == '.') {
				c = 'N';
			}
			if(isalpha(c)) {
				// If it's past the 5'-end trim point
				if(nchar++ >= this->trim5_) {
					r.patBufFw[seqoff++] = charToDna5[c];
				}
			}
			assert_lt(cur, buflen);
			c = r.readOrigBuf[cur++];
		}
	}
	_setBegin(r.patFw, (Dna5*)r.patBufFw);
	// record amt trimmed from 5' end due to --trim5
	r.trimmed5 = (int)(nchar - seqoff);
	// record amt trimmed from 3' end due to --trim3
	int trim3 = (seqoff < this->trim3_) ? seqoff : this->trim3_;
	_setLength(r.patFw, seqoff - trim3);
	r.patBufFw[seqan::length(r.patFw)] = '\0';
	r.trimmed3 = trim3;

	assert_eq('+', c);
	do {
		assert_lt(cur, buflen);
		c = r.readOrigBuf[cur++];
	} while(c != '\n' && c != '\r');
	while(cur < buflen && (c == '\n' || c == '\r')) {
		c = r.readOrigBuf[cur++];
	}
	
	assert_eq(0, seqan::length(r.qual));
	int nqual = 0, qualoff = 0;
	if (intQuals_) {
		int cur_int = 0;
		while(c != '\t' && c != '\n' && c != '\r') {
			cur_int *= 10;
			cur_int += (int)(c - '0');
			c = r.readOrigBuf[cur++];
			if(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
				char cadd = intToPhred33(cur_int, solQuals_);
				cur_int = 0;
				if (c == ' ')
					c = r.readOrigBuf[cur++];
				assert_geq(cadd, 33);
				if(++nqual > this->trim5_) {
					r.qualBuf[qualoff++] = cadd;
				}
			}
		}
	} else {
		c = charToPhred33(c, solQuals_, phred64Quals_);
		if(nqual++ >= r.trimmed5) {
			r.qualBuf[qualoff++] = c;
		}
		while(cur < buflen) {
			c = r.readOrigBuf[cur++];
			if (c == ' ') {
				wrongQualityFormat(r.name);
				return false;
			}
			if(c == '\r' || c == '\n') {
				break;
			}
			c = charToPhred33(c, solQuals_, phred64Quals_);
			if(nqual++ >= r.trimmed5) {
				r.qualBuf[qualoff++] = c;
			}
		}
		if(qualoff < seqoff) {
			tooFewQualities(r.name);
			return false;
		} else if(qualoff > seqoff) {
			tooManyQualities(r.name);
			return false;
		}
	}
	r.qualBuf[seqan::length(r.patFw)] = '\0';
	_setBegin(r.qual, r.qualBuf);
	_setLength(r.qual, seqan::length(r.patFw));

	// Set up a default name if one hasn't been set
	if(seqan::length(r.name) == 0) {
		itoa10<TReadId>(static_cast<TReadId>(readCnt_), r.nameBuf);
		_setBegin(r.name, r.nameBuf);
		_setLength(r.name, nameoff);
	}
	r.parsed = true;
	if(!rb.parsed && rb.readOrigBufLen > 0) {
		return parse(rb, r, rdid);
	}
	return true;
}

/**
 * Light-parse a batch of tabbed-format reads into given buffer.
 */
pair<bool, int> TabbedPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	int c = getc_unlocked(fp_);
	while(c >= 0 && (c == '\n' || c == '\r')) {
		c = getc_unlocked(fp_);
	}
	vector<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && c >= 0; readi++) {
		readbuf[readi].readOrigBufLen = 0;
		while(c >= 0 && c != '\n' && c != '\r') {
			readbuf[readi].readOrigBuf[readbuf[readi].readOrigBufLen++] = c;
			c = getc_unlocked(fp_);
		}
		while(c >= 0 && (c == '\n' || c == '\r')) {
			c = getc_unlocked(fp_);
		}
	}
	return make_pair(c < 0, readi);
}

/**
 * Finalize tabbed parsing outside critical section.
 */
bool TabbedPatternSource::parse(Read& ra, Read& rb, TReadId rdid) const {
	// Light parser (nextBatchFromFile) puts unparsed data
	// into Read& r, even when the read is paired.
	assert(ra.empty());
	assert(rb.empty());
	assert_gt(ra.readOrigBufLen, 0); // raw data for read/pair is here
	assert_eq(0, rb.readOrigBufLen);
	int c = '\t';
	size_t cur = 0;
	const size_t buflen = ra.readOrigBufLen;
	bool paired = false;
	
	// Loop over the two ends
	for(int endi = 0; endi < 2 && c == '\t'; endi++) {
		Read& r = ((endi == 0) ? ra : rb);
		assert_eq(0, seqan::length(r.name));
		// Parse name if (a) this is the first end, or
		// (b) this is tab6
		int nameoff = 0;
		if(endi < 1 || secondName_) {
			// Parse read name
			c = ra.readOrigBuf[cur++];
			while(c != '\t' && cur < buflen) {
				r.name[nameoff++] = c;
				c = ra.readOrigBuf[cur++];
			}
			assert_eq('\t', c);
			if(cur >= buflen) {
				return false; // record ended prematurely
			}
			r.nameBuf[nameoff] = '\0';
			_setBegin(r.name, r.nameBuf);
			_setLength(r.name, nameoff);
		} else if(endi > 0) {
			// if this is the second end and we're parsing
			// tab5, copy name from first end
			rb.name = ra.name; // not a deep copy
		}
		
		paired = endi > 0;

		// Parse sequence
		assert_eq(0, seqan::length(r.patFw));
		c = ra.readOrigBuf[cur++];
		int nchar = 0, seqoff = 0;
		if(color_ && asc2dnacat[c] > 0) {
			// First char is a DNA char (primer)
			if(asc2colcat[toupper(r.readOrigBuf[cur++])] <= 0) {
				// 2nd char isn't a color, so don't assume 'c' is primer
				cur -= 2;
			} else {
				// 'c' is primer
				r.primer = c;
			}
			c = r.readOrigBuf[cur++];
		}
		if(color_) {
			while(c != '\t' && cur < buflen) {
				if(c >= '0' && c < '4') {
					c = "ACGTN"[(int)c - '0'];
				}
				if(c == '.') {
					c = 'N';
				}
				if(isalpha(c)) {
					assert_in(toupper(c), "ACGTN");
					if(nchar++ >= this->trim5_) {
						assert_neq(0, asc2dnacat[c]);
						r.patBufFw[seqoff++] = charToDna5[c]; // ascii to int
					}
				}
				c = r.readOrigBuf[cur++];
			}
			r.color = true;
		} else {
			while(c != '\t' && cur < buflen) {
				if(isalpha(c)) {
					assert_in(toupper(c), "ACGTN");
					if(nchar++ >= this->trim5_) {
						assert_neq(0, asc2dnacat[c]);
						r.patBufFw[seqoff++] = charToDna5[c];
					}
				}
				c = ra.readOrigBuf[cur++];
			}
		}
		assert_eq('\t', c);
		if(cur >= buflen) {
			return false; // record ended prematurely
		}
		_setBegin(r.patFw, (Dna5*)r.patBufFw);
		// record amt trimmed from 5' end due to --trim5
		r.trimmed5 = (int)(nchar - seqoff);
		// record amt trimmed from 3' end due to --trim3
		int trim3 = (seqoff < this->trim3_) ? seqoff : this->trim3_;
		_setLength(r.patFw, seqoff - trim3);
		r.patBufFw[seqan::length(r.patFw)] = '\0';
		r.trimmed3 = trim3;
		
		// Parse qualities
		assert_eq(0, seqan::length(r.qual));
		c = ra.readOrigBuf[cur++];
		int nqual = 0, qualoff = 0;
		if (intQuals_) {
			int cur_int = 0;
			while(c != '\t' && c != '\n' && c != '\r' && cur < buflen) {
				cur_int *= 10;
				cur_int += (int)(c - '0');
				c = ra.readOrigBuf[cur++];
				if(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
					char cadd = intToPhred33(cur_int, solQuals_);
					cur_int = 0;
					assert_geq(cadd, 33);
					if(++nqual > this->trim5_) {
						r.qualBuf[qualoff++] = cadd;
					}
				}
			}
		} else {
			while(c != '\t' && c != '\n' && c != '\r') {
				if(c == ' ') {
					wrongQualityFormat(r.name);
					return false;
				}
				char cadd = charToPhred33(c, solQuals_, phred64Quals_);
				if(++nqual > this->trim5_) {
					r.qualBuf[qualoff++] = cadd;
				}
				if(cur >= buflen) break;
				c = ra.readOrigBuf[cur++];
			}
		}
		if(nchar > nqual) {
			tooFewQualities(r.name);
			return false;
		} else if(nqual > nchar) {
			tooManyQualities(r.name);
			return false;
		}
		r.qualBuf[seqan::length(r.patFw)] = '\0';
		_setBegin(r.qual, r.qualBuf);
		_setLength(r.qual, seqan::length(r.patFw));
		assert(c == '\t' || c == '\n' || c == '\r' || cur >= buflen);
	}
	ra.parsed = true;
	rb.parsed = paired;
	return true;
}

/**
 * Light-parse a batch of raw-format reads into given buffer.
 */
pair<bool, int> RawPatternSource::nextBatchFromFile(
	PerThreadReadBuf& pt,
	bool batch_a)
{
	int c = getc_unlocked(fp_);
	while(c >= 0 && (c == '\n' || c == '\r')) {
		c = getc_unlocked(fp_);
	}
	vector<Read>& readbuf = batch_a ? pt.bufa_ : pt.bufb_;
	size_t readi = 0;
	// Read until we run out of input or until we've filled the buffer
	for(; readi < pt.max_buf_ && c >= 0; readi++) {
		readbuf[readi].readOrigBufLen = 0;
		while(c >= 0 && c != '\n' && c != '\r') {
			readbuf[readi].readOrigBuf[readbuf[readi].readOrigBufLen++] = c;
			c = getc_unlocked(fp_);
		}
		while(c >= 0 && (c == '\n' || c == '\r')) {
			c = getc_unlocked(fp_);
		}
	}
	return make_pair(c < 0, readi);
}

/**
 * Finalize raw parsing outside critical section.
 */
bool RawPatternSource::parse(Read& r, Read& rb, TReadId rdid) const {
	assert(r.empty());
	assert_gt(r.readOrigBufLen, 0);
	size_t cur = 0;
	const size_t buflen = r.readOrigBufLen;

	// Parse sequence
	assert_eq(0, seqan::length(r.patFw));
	int nchar = 0, seqoff = 0;
	int c = r.readOrigBuf[cur++];

	if(color_ && asc2dnacat[c] > 0) {
		// First char is a DNA char (primer)
		if(asc2colcat[toupper(r.readOrigBuf[cur++])] <= 0) {
			// 2nd char isn't a color, so don't assume 'c' is primer
			cur -= 2;
		} else {
			// 'c' is primer
			r.primer = c;
		}
		c = r.readOrigBuf[cur++];
	}
	if(color_) {
		while(c != '\0') {
			assert(c != '\r' && c != '\n');
			if(c >= '0' && c < '4') {
				c = "ACGTN"[(int)c - '0'];
			}
			if(c == '.') {
				c = 'N';
			}
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(nchar++ >= this->trim5_) {
					assert_neq(0, asc2dnacat[c]);
					r.patBufFw[seqoff++] = charToDna5[c]; // ascii to int
				}
			}
			c = r.readOrigBuf[cur++];
		}
		r.color = true;
	} else {
		cur--;
		while(cur < buflen) {
			c = r.readOrigBuf[cur++];
			assert(c != '\r' && c != '\n');
			if(isalpha(c)) {
				assert_in(toupper(c), "ACGTN");
				if(nchar++ >= this->trim5_) {
					assert_neq(0, asc2dnacat[c]);
					r.patBufFw[seqoff++] = charToDna5[c];
				}
			}
		}
	}
	_setBegin(r.patFw, (Dna5*)r.patBufFw);
	// record amt trimmed from 5' end due to --trim5
	r.trimmed5 = (int)(nchar - seqoff);
	// record amt trimmed from 3' end due to --trim3
	int trim3 = (seqoff < this->trim3_) ? seqoff : this->trim3_;
	_setLength(r.patFw, seqoff - trim3);
	r.patBufFw[seqan::length(r.patFw)] = '\0';
	r.trimmed3 = trim3;
	
	// Give the name field a dummy value
	itoa10<TReadId>(rdid, r.nameBuf);
	_setBegin(r.name, r.nameBuf);
	_setLength(r.name, strlen(r.nameBuf));
	
	// Give the base qualities dummy values
	assert_eq(0, seqan::length(r.qual));
	const size_t len = seqan::length(r.patFw);
	for(size_t i = 0; i < len; i++) {
		r.qualBuf[i] = 'I';
	}
	_setBegin(r.qual, r.qualBuf);
	_setLength(r.qual, seqan::length(r.patFw));
	
	r.parsed = true;
	if(!rb.parsed && rb.readOrigBufLen > 0) {
		return parse(rb, r, rdid);
	}
	return true;
}

void wrongQualityFormat(const String<char>& read_name) {
	cerr << "Encountered a space parsing the quality string for read " << read_name << endl
	     << "If this is a FASTQ file with integer (non-ASCII-encoded) qualities, please" << endl
	     << "re-run Bowtie with the --integer-quals option." << endl;
	throw 1;
}

void tooFewQualities(const String<char>& read_name) {
	cerr << "Too few quality values for read: " << read_name << endl
		 << "\tare you sure this is a FASTQ-int file?" << endl;
	throw 1;
}

void tooManyQualities(const String<char>& read_name) {
	cerr << "Reads file contained a pattern with more than 1024 quality values." << endl
		 << "Please truncate reads and quality values and and re-run Bowtie" << endl;
	throw 1;
}

void tooManySeqChars(const String<char>& read_name) {
	cerr << "Reads file contained a pattern with more than 1024 sequence characters." << endl
		 << "Please truncate reads and quality values and and re-run Bowtie." << endl
		 << "Offending read: " << read_name << endl;
	throw 1;
}
