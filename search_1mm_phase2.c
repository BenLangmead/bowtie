/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the second phase of the 1-mismatch search
 * routine.  It is implemented as a code fragment so that it can be
 * reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	params.setFw(false);
	params.setEbwtFw(false);
	bool hit;
	ASSERT_ONLY(uint64_t numHits);

	// Next, try hits with one mismatch on the 3' end for the reverse-complement read
	ASSERT_ONLY(numHits = sink.numHits());
	bt2.setQuery(&patRc, &qualRc, &name);
	bt2.setOffs(0, 0, s3, s, s, s); // 1 mismatch allowed in 3' half
	hit = bt2.backtrack();
	assert(hit  || numHits == sink.numHits());
	assert(!hit || numHits <  sink.numHits());
	if(hit) {
		assert_eq(numHits+1, sink.numHits());
		sanityCheckHits(patRc, sink, patid, false /*fw*/, os,
		                false /*allowExact*/, true /*transpose*/);
		continue;
	} else {
		// No hits
		sanityCheckHits(patRc, sink, patid, false /*fw*/, os,
		                false /*allowExact*/, true /*transpose*/);
	}

	params.setFw(true);

	// Next, try hits with one mismatch on the 3' end for the reverse-complement read
	ASSERT_ONLY(numHits = sink.numHits());
	bt2.setQuery(&patFw, &qualFw, &name);
	bt2.setOffs(0, 0, s3, s, s, s); // 1 mismatch allowed in 3' half
	hit = bt2.backtrack();
	assert(hit  || numHits == sink.numHits());
	assert(!hit || numHits <  sink.numHits());
	if(hit) {
		assert_eq(numHits+1, sink.numHits());
		sanityCheckHits(patFw, sink, patid, true /*fw*/, os,
		                false /*allowExact*/, true /*transpose*/);
		continue;
	} else {
		// No hits
		sanityCheckHits(patFw, sink, patid, true /*fw*/, os,
		                false /*allowExact*/, true /*transpose*/);
	}
}
