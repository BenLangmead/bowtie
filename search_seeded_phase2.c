/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the second phase of the seeded, quality-
 * aware search routine.  It is implemented as a code fragment so that
 * it can be reused in both the half-index-in-memory and full-index-in-
 * memory situations.
 */
{
	if(!nofw) {
		// If we reach here, then cases 1R, 2R, and 3R have been
		// eliminated.  The next most likely cases are 1F, 2F and
		// 3F...
		params.setFw(true);  // looking at forward strand
		btf2.setReportExacts(false);
		btr2.setReportExacts(false);
		btf2.setQuery(patsrc->bufa());
		// Set up seed bounds
		if(qs < s) {
			btf2.setOffs(0, 0,
						(seedMms > 0)? qs5 : qs,
						(seedMms > 1)? qs5 : qs,
						(seedMms > 2)? qs5 : qs,
						(seedMms > 3)? qs5 : qs);
		} else {
			btf2.setOffs(0, 0,
						(seedMms > 0)? s5 : s,
						(seedMms > 1)? s5 : s,
						(seedMms > 2)? s5 : s,
						(seedMms > 3)? s5 : s);
		}
		// Do a 12/24 backtrack on the forward-strand read using
		// the mirror index.  This will find all case 1F, 2F
		// and 3F hits.
		if(btf2.backtrack()) {
			// The reverse complement hit, so we're done with this
			// read
			DONEMASK_SET(patid);
			continue;
		}

		if(sink->finishedWithStratum(0)) { // no more exact hits are possible
			DONEMASK_SET(patid);
			continue;
		}
	}

	// No need to collect partial alignments if we're not
	// allowing mismatches in the 5' seed half
	if(seedMms == 0) continue;

	if(!norc) {
		// If we reach here, then cases 1F, 2F, 3F, 1R, 2R, and 3R
		// have been eliminated, leaving us with cases 4F and 4R
		// (the cases with 1 mismatch in the 5' half of the seed)
		params.setFw(false);  // looking at reverse-comp strand
		// Set up seed bounds
		if(qs < s) {
			btr2.setOffs(0, 0,
						qs3,
						(seedMms > 1)? qs3 : qs,
						(seedMms > 2)? qs3 : qs,
						(seedMms > 3)? qs3 : qs);
		} else {
			btr2.setOffs(0, 0,
						s3,
						(seedMms > 1)? s3 : s,
						(seedMms > 2)? s3 : s,
						(seedMms > 3)? s3 : s);
		}
		btr2.setQuery(patsrc->bufa());
		btr2.setQlen(s); // just look at the seed
		// Find partial alignments for case 4R
		ASSERT_ONLY(bool done =) btr2.backtrack();
#ifndef NDEBUG
		vector<PartialAlignment> partials;
		assert(pamRc != NULL);
		pamRc->getPartials(patid, partials);
		if(done) assert_gt(partials.size(), 0);
		for(size_t i = 0; i < partials.size(); i++) {
			uint32_t pos0 = partials[i].entry.pos0;
			assert_lt(pos0, s5);
			uint8_t oldChar = (uint8_t)patRcRev[pos0];
			assert_neq(oldChar, partials[i].entry.char0);
			if(partials[i].entry.pos1 != 0xffff) {
				uint32_t pos1 = partials[i].entry.pos1;
				assert_lt(pos1, s5);
				oldChar = (uint8_t)patRcRev[pos1];
				assert_neq(oldChar, partials[i].entry.char1);
				if(partials[i].entry.pos2 != 0xffff) {
					uint32_t pos2 = partials[i].entry.pos2;
					assert_lt(pos2, s5);
					oldChar = (uint8_t)patRcRev[pos2];
					assert_neq(oldChar, partials[i].entry.char2);
				}
			}
		}
#endif
	}
}
