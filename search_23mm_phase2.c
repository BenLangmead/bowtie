/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the second phase of the 2/3-mismatch
 * search routine.  It is implemented as a code fragment so that it can
 * be reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	bt2.setReportExacts(false);

	if(!nofw) {
		// Set up the revisitability of the halves
		params.setFw(true);
		bt2.setQuery(patsrc->bufa());
		bt2.setOffs(0, 0, s5, s5, two? s : s5, s);
		if(bt2.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
		// if nofw is true, then we already did this
		if(sink->finishedWithStratum(0)) { // no more exact hits are possible
			DONEMASK_SET(patid);
			continue;
		}
	}

	if(!norc) {
		// Try 2/3 backtracks in the 3' half of the reverse complement read
		params.setFw(false);  // looking at reverse complement
		// Set up the revisitability of the halves
		bt2.setQuery(patsrc->bufa());
		bt2.setOffs(0, 0, s3, s3, two? s : s3, s);
		if(bt2.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
	}

	if(nofw && sink->finishedWithStratum(1)) {
		DONEMASK_SET(patid);
		continue;
	}
}
