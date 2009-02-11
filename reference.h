#ifndef REFERENCE_H_
#define REFERENCE_H_

/**
 * Abstract parent class for a representation of the reference strings.
 */
class Reference {
public:
	virtual ~Reference() { }
	virtual void loadStretch(uint8_t *dest, uint32_t tidx, uint32_t toff, uint32_t count) = 0;
};

/**
 * Concrete reference representation that bulk-loads the reference from
 * the bit-pair-compacted binary file and stores it in memory also in
 * bit-pair-compacted format.
 */
class BitPairReference {

public:
	/**
	 * Load from Bowtie index.
	 */
	BitPairReference(const string& in,
	                 bool sanity = false,
	                 vector<string>* infiles = NULL,
	                 bool infilesSeq = false,
	                 bool useShmem = false) :
	sanity_(sanity),
	useShmem_(useShmem)
	{
		string s3 = in + ".3.ebwt";
		string s4 = in + ".4.ebwt";
		FILE *f3 = fopen(s3.c_str(), "rb");
		if(f3 == NULL) {
			cerr << "Could not open reference-string index file " << s3 << " for reading." << endl;
			exit(1);
		}
		// Read endianness sentinel, set 'swap'
		uint32_t one;
		bool swap = false;
		if(!fread(&one, 4, 1, f3)) {
			cerr << "Error reading first word from reference-structure index file " << s3 << endl;
			exit(1);
		}
		if(one != 1) {
			assert_eq(0x1000000, one);
			swap = true; // have to endian swap U32s
		}
		// Read # records
		uint32_t sz;
		if(!fread(&sz, 4, 1, f3)) {
			cerr << "Error reading # records from reference-structure index file " << s3 << endl;
			exit(1);
		}
		if(swap) sz = endianSwapU32(sz);
		// Read records
		nrefs_ = 0;
		// Cumulative count of all unambiguous characters on a per-
		// stretch 8-bit alignment (i.e. count of bytes we need to
		// allocate in buf_)
		uint32_t cumsz = 0;
		// Cumulative offset in reference string
		uint32_t cumoff = 0;
		for(uint32_t i = 0; i < sz; i++) {
			recs_.push_back(RefRecord(f3, swap));
			if(recs_.back().first) {
				// Remember that this is the first record for this
				// reference sequence
				refRecOffs_.push_back(recs_.size()-1);
				nrefs_++;
				cumoff = 0;
			}
			cumoff += recs_.back().off;
			refOffs_.push_back(cumoff);
			cumoff += recs_.back().len;
			cumsz += recs_.back().len;
		}
		refRecOffs_.push_back(recs_.size());
		assert_eq(sz, recs_.size());
		fclose(f3); // done with .3.ebwt file
		// Round cumsz up to nearest byte boundary
		if((cumsz & 3) != 0) {
			cumsz += (4 - (cumsz & 3));
		}
		assert_eq(0, cumsz & 3);
		FILE *f4 = fopen(s4.c_str(), "rb");
		buf_ = new uint8_t[cumsz >> 2];
		size_t ret = fread(buf_, 1, cumsz >> 2, f4);
		if(ret != (cumsz >> 2)) {
			cerr << "Only read " << ret << " bytes (out of " << (cumsz >> 2) << ") from reference index file " << s4 << endl;
			exit(1);
		}
		fclose(f4);
#ifndef NDEBUG
		if(sanity_) {
			// Compare the sequence we just read from the compact index
			// file to the true reference sequence.
			assert(infiles != NULL);
			std::vector<seqan::String<seqan::Dna5> > os; // for holding references
			if(infilesSeq) {
				for(size_t i = 0; i < infiles->size(); i++) {
					os.push_back(String<Dna5>((*infiles)[i]));
				}
			} else {
				readSequenceFiles<seqan::String<seqan::Dna5>, seqan::Fasta>(*infiles, os);
			}
			for(size_t i = 0; i < os.size(); i++) {
//				size_t olen = seqan::length(os[i]);
//				uint8_t *buf = new uint8_t[olen];
//				loadStretch(buf, i, 0, olen);
//				for(size_t j = 0; j < olen; j++) {
//					assert_eq(os[i][j], buf[j]);
//				}
//				delete[] buf;
			}
		}
#endif
	}

	virtual ~BitPairReference() {
		delete[] buf_;
	}

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 */
	virtual void loadStretch(uint8_t *dest,
	                         uint32_t tidx,
	                         uint32_t toff,
	                         uint32_t count)
	{
		uint32_t reci = refRecOffs_[tidx];
		uint32_t recf = refRecOffs_[tidx+1];
		assert_gt(recf, reci);
		uint32_t cur = 0;
		for(uint32_t i = reci; i < recf; i++) {
			if(i < recf-1) {
				assert_lt(refOffs_[i], refOffs_[i+1]);
			}
			// Holds the beginning of the stretch?
			if(refOffs_[i] <= toff && recs_[i].len + refOffs_[i] > toff) {

				//
				// TODO
				//

				uint8_t *ptr = buf_;
				// Calculate offset into this stretch to start capturing
				uint32_t off = toff - refOffs_[i];
				// Calculate how much of this stretch to capture
				uint32_t captureSz = min(count, recs_[i].len - off);
				uint32_t byOff = off >> 2;
				uint32_t bpOff = off & 3;
				uint32_t j = 0;
				// Do some base-pairs at the beginning
				//int beginDiff = 4 - (off & 3);
				for(; bpOff != 0; bpOff++) {
					dest[cur++] = ptr[byOff];

				}
				// Do base-pairs in the middle
				for(; j + 4 < captureSz; j += 4) {
					dest[cur++] =  ptr[byOff]       & 3;
					dest[cur++] = (ptr[byOff] >> 2) & 3;
					dest[cur++] = (ptr[byOff] >> 4) & 3;
					dest[cur++] = (ptr[byOff] >> 6) & 3;
					byOff++;
				}
				// Do
				for(; j + 4 < captureSz; j += 4) {

				}
				assert_leq(captureSz, count);
				count -= captureSz;
			}
		}
	}

protected:
	uint8_t **refs_;
	std::vector<RefRecord> recs_;
	std::vector<uint32_t>  refOffs_;
	std::vector<uint32_t>  refRecOffs_;
	uint8_t *buf_;   // the whole reference as a big bitpacked byte array
	uint32_t bufSz_; // size of buf_
	uint32_t nrefs_; //
	bool     sanity_;
	bool     useShmem_;
};

#endif
