/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the third phase of the seeded, quality-
 * aware search routine.  It is implemented as a code fragment so that
 * it can be reused in both the half-index-in-memory and full-index-in-
 * memory situations.
 */
{
	if(!norc) {
		params.setFw(false);
		btr3.setReportExacts(true);

		btr3.setQuery(patsrc->bufa());
		// Get all partial alignments for this read's reverse
		// complement
		pals.clear();
		if(pamRc != NULL && pamRc->size() > 0) {
			// We can get away with an unsynchronized call because there
			// are no writers for pamRc in this phase
			pamRc->getPartialsUnsync(patid, pals);
			pamRc->clear(patid);
			assert_eq(0, pamRc->size());
		}
		bool done = false;
		if(pals.size() > 0) {
			// Partial alignments exist - extend them
			// Set up seed bounds
			if(qs < s) {
				btr3.setOffs(0, 0, qs, qs, qs, qs);
			} else {
				btr3.setOffs(0, 0, s, s, s, s);
			}
			for(size_t i = 0; i < pals.size(); i++) {
				seqan::clear(muts);
				uint8_t oldQuals =
					PartialAlignmentManager::toMutsString(
						pals[i], patRc, qualRev, muts, !noMaqRound);

				// Set the backtracking thresholds appropriately
				// Now begin the backtracking, treating the first
				// 24 bases as unrevisitable
				ASSERT_ONLY(String<Dna5> tmp = patRc);
				btr3.setMuts(&muts);
				done = btr3.backtrack(oldQuals);
				btr3.setMuts(NULL);
				assert_eq(tmp, patRc); // assert mutations were undone
				if(done) {
					// The reverse complement hit, so we're done with this
					// read
					DONEMASK_SET(patid);
					// Got a hit; stop processing partial
					// alignments
					break;
				}
			} // Loop over partial alignments
		}
		seqan::clear(muts);
		// Case 4R yielded a hit continue to next pattern
		if(done) continue;
		// If we're in two-mismatch mode, then now is the time to
		// try the final case that might apply to the reverse
		// complement pattern: 1 mismatch in each of the 3' and 5'
		// halves of the seed.
		bool gaveUp = false;
		if(seedMms >= 2) {
			btr23.setQuery(patsrc->bufa());
			// Set up special seed bounds
			if(qs < s) {
				btr23.setOffs(qs5, qs,
							  0,                         // unrevOff
							  (seedMms <= 2)? qs5 : 0,   // 1revOff
							  (seedMms < 3 )? qs  : qs5, // 2revOff
							  qs);                       // 3revOff
			} else {
				btr23.setOffs(s5, s,
							  0,                       // unrevOff
							  (seedMms <= 2)? s5 : 0,  // 1revOff
							  (seedMms < 3 )? s  : s5, // 2revOff
							  s);                      // 3revOff
			}
			done = btr23.backtrack();
			if(btr23.numBacktracks() == btr23.maxBacktracks()) {
				gaveUp = true;
			}
			if(done) {
				DONEMASK_SET(patid);
				btr23.resetNumBacktracks();
				continue;
			}
			btr23.resetNumBacktracks();
		}
	}

	if(nofw) { // no more 1-mm-in-seed hits are possible
		//DONEMASK_SET(patid);
		continue;
	}

	// If we reach here, then cases 1F, 2F, 3F, 1R, 2R, 3R and
	// 4R have been eliminated leaving only 4F.
	params.setFw(true);  // looking at forward strand
	btf3.setQuery(patsrc->bufa());
	btf3.setQlen(seedLen); // just look at the seed
	// Set up seed bounds
	if(qs < s) {
		btf3.setOffs(0, 0,
		             qs3,
		             (seedMms > 1)? qs3 : qs,
		             (seedMms > 2)? qs3 : qs,
		             (seedMms > 3)? qs3 : qs);
	} else {
		btf3.setOffs(0, 0,
		             s3,
		             (seedMms > 1)? s3 : s,
		             (seedMms > 2)? s3 : s,
		             (seedMms > 3)? s3 : s);
	}
	// Do a 12/24 seedling backtrack on the forward read
	// using the normal index.  This will find seedlings
	// for case 4F
	btf3.backtrack();
#ifndef NDEBUG
	vector<PartialAlignment> partials;
	pamFw->getPartials(patid, partials);
	for(size_t i = 0; i < partials.size(); i++) {
		uint32_t pos0 = partials[i].entry.pos0;
		assert_lt(pos0, s5);
		uint8_t oldChar = (uint8_t)patFw[pos0];
		assert_neq(oldChar, partials[i].entry.char0);
		if(partials[i].entry.pos1 != 0xffff) {
			uint32_t pos1 = partials[i].entry.pos1;
			assert_lt(pos1, s5);
			oldChar = (uint8_t)patFw[pos1];
			assert_neq(oldChar, partials[i].entry.char1);
			if(partials[i].entry.pos2 != 0xffff) {
				uint32_t pos2 = partials[i].entry.pos2;
				assert_lt(pos2, s5);
				oldChar = (uint8_t)patFw[pos2];
				assert_neq(oldChar, partials[i].entry.char2);
			}
		}
	}
#endif
}
