/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the third phase of the 2/3-mismatch
 * search routine.  It is implemented as a code fragment so that it can
 * be reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	if(!nofw) {
		// Try 2/3 backtracks in the 3' half of the forward read
		params.setFw(true);
		bt3.setReportExacts(false);
		bt3.setQuery(patsrc->bufa());
		bt3.setOffs(0, 0,
					s3,
					s3,
					two? s : s3,
					s);
		bool done = bt3.backtrack();
		if(done) continue;
		// no more 1-mismatch hits are possible after this point
		if(sink->finishedWithStratum(1)) {
			continue;
		}

		// Try a half-and-half on the forward read
		//bool gaveUp = false;
		bthh3.setQuery(patsrc->bufa());
		// Processing the forward pattern with the forward index;
		// s3 ("lo") half is on the right
		bthh3.setOffs(s3, s,
		              0,
		              two ? s3 : 0,
		              two ? s  : s3,
		              s);
		done = bthh3.backtrack();
		//if(bthh3.numBacktracks() == bthh3.maxBacktracks()) {
			//gaveUp = true;
		//}
		bthh3.resetNumBacktracks();
		if(done) {
			continue;
		}
	}

	if(!norc) {
		// Try a half-and-half on the reverse complement read
		//bool gaveUp = false;
		params.setFw(false);
		bthh3.setQuery(patsrc->bufa());
		// Processing the forward pattern with the forward index;
		// s5 ("hi") half is on the right
		bthh3.setOffs(s5, s,
					  0,
					  two ? s5 : 0,
					  two ? s  : s5,
					  s);
		bool done = bthh3.backtrack();
		//if(bthh3.numBacktracks() == bthh3.maxBacktracks()) {
		//	gaveUp = true;
		//}
		bthh3.resetNumBacktracks();
		if(done) {
			continue;
		}
	}
}
