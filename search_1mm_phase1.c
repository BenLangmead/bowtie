/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the first phase of the 1-mismatch search
 * routine.  It is implemented as a code fragment so that it can be
 * reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	bt.setEbwt(&ebwtFw);
	bt.setReportExacts(true);

	if(plen < 2) {
		cerr << "Error: Reads must be at least 2 characters long in 1-mismatch mode" << endl;
		throw 1;
	}

	if(!nofw) {
		// First, try exact hits for the forward-oriented read
		params.setFw(true);
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, s, s, s, s);
		if(bt.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
	}

	if(!norc) {
		params.setFw(false);

		// Next, try exact hits for the reverse-complement read
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, s, s, s, s);
		if(bt.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
	}

	if(sink->finishedWithStratum(0)) { // no more exact hits are possible
		// the sink tells us we needn't try 1-mismatch alignments
		DONEMASK_SET(patid);
		continue;
	}
	bt.setReportExacts(false);

	if(!norc) {
		// Next, try hits with one mismatch on the 3' end for the reverse-complement read
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, s5, s, s, s); // 1 mismatch allowed in 3' half
		if(bt.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
	}

	if(!nofw) {
		params.setFw(true);
		// Next, try hits with one mismatch on the 3' end for the reverse-complement read
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, s5, s, s, s); // 1 mismatch allowed in 3' half
		if(bt.backtrack()) {
			DONEMASK_SET(patid);
			continue;
		}
	}
}
