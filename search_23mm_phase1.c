/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the first phase of the 2/3-mismatch search
 * routine.  It is implemented as a code fragment so that it can be
 * reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	// If requested, check that this read has the same length
	// as all the previous ones
	btr1.setReportExacts(true);

	if(plen < 3 && two) {
		cerr << "Error: Read (" << name << ") is less than 3 characters long" << endl;
		throw 1;
	}
	else if(plen < 4) {
		cerr << "Error: Read (" << name << ") is less than 4 characters long" << endl;
		throw 1;
	}
	if(!nofw) {
		// Do an exact-match search on the forward pattern, just in
		// case we can pick it off early here
		params.setFw(true);
		btr1.setQuery(patsrc->bufa());
		btr1.setOffs(0, 0, plen, plen, plen, plen);
		if(btr1.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
	}
	if(!norc) {
		// Set up backtracker with reverse complement
		params.setFw(false);
		// Set up the revisitability of the halves
		btr1.setQuery(patsrc->bufa());
		btr1.setOffs(0, 0, s5, s5, two ? s : s5, s);
		if(btr1.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
	}
	if(nofw && sink->finishedWithStratum(0)) { // no more exact hits are possible
		DONEMASK_SET(patid);
		continue;
	}
}
