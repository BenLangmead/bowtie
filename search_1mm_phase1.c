/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the first phase of the 1-mismatch search
 * routine.  It is implemented as a code fragment so that it can be
 * reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	params.setFw(true);
	params.setEbwtFw(true);
	bt.setEbwt(&ebwtFw);
	bt.setReportExacts(true);

	assert_eq(0, sink->retainedHits().size());
	if(plen < 2) {
		cerr << "Error: Reads must be at least 2 characters long in 1-mismatch mode" << endl;
		exit(1);
	}

	// First, try exact hits for the forward-oriented read
	bt.setQuery(&patFw, &qualFw, &name);
	bt.setOffs(0, 0, s, s, s, s);
	if(bt.backtrack()) {
		sanityCheckExact(os, *sink, patFw, patid);
		DONEMASK_SET(patid);
		continue;
	}

	params.setFw(false);

	// Next, try exact hits for the reverse-complement read
	bt.setQuery(&patRc, &qualRc, &name);
	bt.setOffs(0, 0, s, s, s, s);
	if(bt.backtrack()) {
		sanityCheckExact(os, *sink, patRc, patid);
		DONEMASK_SET(patid);
		continue;
	}

	sink->finishedWithStratum(0); // no more exact hits are possible
	bt.setReportExacts(false);

	// Next, try hits with one mismatch on the 3' end for the reverse-complement read
	bt.setQuery(&patRc, &qualRc, &name);
	bt.setOffs(0, 0, s5, s, s, s); // 1 mismatch allowed in 3' half
	if(bt.backtrack()) {
		sanityCheckHits(patRc, *sink, patid, false /*fw*/, os,
		                false /*allowExact*/, false /*transpose*/);
		DONEMASK_SET(patid);
		continue;
	}

	params.setFw(true);

	// Next, try hits with one mismatch on the 3' end for the reverse-complement read
	bt.setQuery(&patFw, &qualFw, &name);
	bt.setOffs(0, 0, s5, s, s, s); // 1 mismatch allowed in 3' half
	if(bt.backtrack()) {
		sanityCheckHits(patFw, *sink, patid, true /*fw*/, os,
		                false /*allowExact*/, false /*transpose*/);
		DONEMASK_SET(patid);
		continue;
	}
}
