/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the fourth phase of the seeded, quality-
 * aware search routine.  It is implemented as a code fragment so that
 * it can be reused in both the half-index-in-memory and full-index-in-
 * memory situations.
 */
{
	if(!nofw) {
		params.setFw(true);
		btf4.setReportExacts(true);
		btf4.setQuery(patsrc->bufa());
		// Get all partial alignments for this read's reverse
		// complement
		pals.clear();
		if(pamFw != NULL && pamFw->size() > 0) {
			// We can get away with an unsynchronized call because there
			// are no writers for pamFw in this phase
			pamFw->getPartialsUnsync(patid, pals);
			pamFw->clear(patid);
			assert_eq(0, pamFw->size());
		}
		bool done = false;
		if(pals.size() > 0) {
			// Partial alignments exist - extend them
			// Set up seed bounds
			if(qs < s) {
				btf4.setOffs(0, 0, qs, qs, qs, qs);
			} else {
				btf4.setOffs(0, 0, s, s, s, s);
			}
			for(size_t i = 0; i < pals.size(); i++) {
				seqan::clear(muts);
				uint8_t oldQuals =
					PartialAlignmentManager::toMutsString(
						pals[i], patFwRev, qualRev, muts, !noMaqRound);
				assert_leq(oldQuals, qualThresh);
				// Set the backtracking thresholds appropriately
				// Now begin the backtracking, treating the first
				// 24 bases as unrevisitable
				ASSERT_ONLY(String<Dna5> tmp = patFwRev);
				btf4.setMuts(&muts);
				done = btf4.backtrack(oldQuals);
				btf4.setMuts(NULL);
				assert_eq(tmp, patFwRev); // assert mutations were undone
				if(done) {
					// Got a hit; stop processing partial
					// alignments
					break;
				}
			} // Loop over partial alignments
		}
		seqan::clear(muts);

		// Case 4F yielded a hit; continue to next pattern
		if(done) continue;

		if(sink->finishedWithStratum(1)) { // no more exact hits are possible
			continue;
		}

		// If we're in two-mismatch mode, then now is the time to
		// try the final case that might apply to the forward
		// pattern: 1 mismatch in each of the 3' and 5' halves of
		// the seed.
		bool gaveUp = false;
		if(seedMms >= 2) {
			btf24.setQuery(patsrc->bufa());
			// Set up seed bounds
			if(qs < s) {
				btf24.setOffs(qs5, qs,
							 0,                        // unrevOff
							 (seedMms <= 2)? qs5 : 0,  // 1revOff
							 (seedMms < 3)? qs : qs5,  // 2revOff
							 qs);                      // 3revOff
			} else {
				btf24.setOffs(s5, s,
							  0,                       // unrevOff
							  (seedMms <= 2)? s5 : 0,  // 1revOff
							  (seedMms < 3)?  s : s5,  // 2revOff
							  s);                      // 3revOff
			}
			done = btf24.backtrack();
			if(btf24.numBacktracks() == btf24.maxBacktracks()) {
				gaveUp = true;
			}
			if(done) {
				btf24.resetNumBacktracks();
				continue;
			}
			btf24.resetNumBacktracks();
		}
	}
}
