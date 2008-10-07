/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the third phase of the 2/3-mismatch
 * search routine.  It is implemented as a code fragment so that it can
 * be reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	params.setFw(true);
	params.setEbwtFw(true);
	bt3.setReportExacts(false);

	// Try 2/3 backtracks in the 3' half of the forward read
	bt3.setQuery(&patFw, &qualFw, &name);
	bt3.setOffs(0, 0,
	            s3,
	            s3,
	            two? s : s3,
	            s);
	ASSERT_ONLY(uint64_t numHits = sink->numHits());
	bool done = bt3.backtrack();
	if(done) continue;

	// no more 1-mismatch hits are possible after this point
	sink->finishedWithStratum(1);

	// Try a half-and-half on the forward read
	bool gaveUp = false;
	bthh3.setQuery(&patFw, &qualFw, &name);
	// Processing the forward pattern with the forward index;
	// s3 ("lo") half is on the right
	bthh3.setOffs(s3, s,
	              0,
	              two ? s3 : 0,
	              two ? s  : s3,
	              s);
	ASSERT_ONLY(numHits = sink->numHits());
	done = bthh3.backtrack();
	if(bthh3.numBacktracks() == bthh3.maxBacktracks()) {
		gaveUp = true;
	}
	bthh3.resetNumBacktracks();
	if(done) {
		if(dumpHHHits != NULL) {
			(*dumpHHHits) << patFw << endl << qualFw << endl << "---" << endl;
		}
		continue;
	}

#ifndef NDEBUG
	// The forward version of the read doesn't hit
	// at all!  Check with the oracle to make sure it agrees.
	if(!gaveUp) {
		ASSERT_NO_HITS_FW(true);
	}
#endif
	// Try a half-and-half on the reverse complement read
	gaveUp = false;
	params.setFw(false);
	bthh3.setQuery(&patRc, &qualRc, &name);
	// Processing the forward pattern with the forward index;
	// s5 ("hi") half is on the right
	bthh3.setOffs(s5, s,
	              0,
	              two ? s5 : 0,
	              two ? s  : s5,
	              s);
	ASSERT_ONLY(numHits = sink->numHits());
	done = bthh3.backtrack();
	if(bthh3.numBacktracks() == bthh3.maxBacktracks()) {
		gaveUp = true;
	}
	bthh3.resetNumBacktracks();
	if(done) {
		if(dumpHHHits != NULL) {
			(*dumpHHHits) << patFw << endl << qualFw << endl << "---" << endl;
		}
		continue;
	}

#ifndef NDEBUG
	// The reverse-complement version of the read doesn't hit
	// at all!  Check with the oracle to make sure it agrees.
	if(!gaveUp) {
		ASSERT_NO_HITS_RC(true);
		if(dumpNoHits != NULL) {
			(*dumpNoHits) << patFw << endl << qualFw << endl << "---" << endl;
		}
	}
#endif
}
