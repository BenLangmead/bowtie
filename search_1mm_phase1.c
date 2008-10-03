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

	assert_eq(0, sink.retainedHits().size());
	if(plen < 2) {
		cerr << "Error: Reads must be at least 2 characters long in 1-mismatch mode" << endl;
		exit(1);
	}
	bool hit;
	ASSERT_ONLY(uint64_t numHits);

	// First, try exact hits for the forward-oriented read
	ASSERT_ONLY(numHits = sink.numHits());
	bt1.setQuery(&patFw, &qualFw, &name);
	bt1.setOffs(0, 0, s, s, s, s);
	hit = bt1.backtrack();
	assert(hit  || numHits == sink.numHits());
	assert(!hit || numHits <  sink.numHits());
	if(hit) {
		assert_eq(numHits+1, sink.numHits());
		sanityCheckExact(os, sink, patFw, patid);
		DONEMASK_SET(patid);
		continue;
	}

	params.setFw(false);

	// Next, try exact hits for the reverse-complement read
	bt1.setQuery(&patRc, &qualRc, &name);
	bt1.setOffs(0, 0, s, s, s, s);
	hit = bt1.backtrack();
	assert(hit  || numHits == sink.numHits());
	assert(!hit || numHits <  sink.numHits());
	if(hit) {
		assert_eq(numHits+1, sink.numHits());
		sanityCheckExact(os, sink, patRc, patid);
		DONEMASK_SET(patid);
		continue;
	}

	// Next, try hits with one mismatch on the 3' end for the reverse-complement read
	bt1.setQuery(&patRc, &qualRc, &name);
	bt1.setOffs(0, 0, s5, s, s, s); // 1 mismatch allowed in 3' half
	hit = bt1.backtrack();
	assert(hit  || numHits == sink.numHits());
	assert(!hit || numHits <  sink.numHits());
	if(hit) {
		assert_eq(numHits+1, sink.numHits());
		sanityCheckHits(patRc, sink, patid, false /*fw*/, os,
		                false /*allowExact*/, false /*transpose*/);
		DONEMASK_SET(patid);
		continue;
	}

	params.setFw(true);

	// Next, try hits with one mismatch on the 3' end for the reverse-complement read
	bt1.setQuery(&patFw, &qualFw, &name);
	bt1.setOffs(0, 0, s5, s, s, s); // 1 mismatch allowed in 3' half
	hit = bt1.backtrack();
	assert(hit  || numHits == sink.numHits());
	assert(!hit || numHits <  sink.numHits());
	if(hit) {
		assert_eq(numHits+1, sink.numHits());
		sanityCheckHits(patFw, sink, patid, true /*fw*/, os,
		                false /*allowExact*/, false /*transpose*/);
		DONEMASK_SET(patid);
		continue;
	}
}
