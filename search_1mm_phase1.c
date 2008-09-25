/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the first phase of the 1-mismatch search
 * routine.  It is implemented as a code fragment so that it can be
 * reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
params.setFw(false);
params.setEbwtFw(true);
assert_eq(0, sink.retainedHits().size());
assert_eq(lastHits, sink.numHits());
uint32_t plen = length(patFw);
if(plen < 2) {
	cerr << "Error: Reads must be at least 2 characters long in 1-mismatch mode" << endl;
	exit(1);
}
// Create state for a search in the forward index
sfw.newQuery(&patRc, &name, &qualRc);
ebwtFw.search1MismatchOrBetter(sfw, params,
							   true,  // allow exact hits,
							   true); // inexact hits provisional
bool hit = sink.numHits() > lastHits;
// Set a bit indicating this pattern is done and needn't be
// considered by the 1-mismatch loop
if(sanity) sanityCheckHits(patRc, sink, patid, false, os, true, false);
assert_eq(0, sink.retainedHits().size());
if(hit) lastHits = sink.numHits();
if(oneHit && hit) {
	assert_eq(0, sink.numProvisionalHits());
	DONEMASK_SET(patid);
	continue;
}
params.setFw(true);
sfw.newQuery(&patFw, &name, &qualFw);
if(sink.numProvisionalHits() > 0) {
	// There is a provisional inexact match for the
	// reverse-complement read, so just try exact on the
	// forward-oriented read
	ebwtFw.search(sfw, params);
	if(sink.numHits() > lastHits) {
		// Got one or more exact hits from the reverse
		// complement; reject provisional hits
		sink.rejectProvisionalHits();
		if(sanity) sanityCheckHits(patFw, sink, patid, true, os, true, false);
	} else {
		// No exact hits from reverse complement; accept
		// provisional hits and finish with this read
		sink.acceptProvisionalHits();
		assert_gt(sink.numHits(), lastHits);
	}
	assert_eq(0, sink.numProvisionalHits());
	if(sink.numHits() > lastHits) {
		lastHits = sink.numHits();
		if(oneHit) {
			// Update doneMask
			DONEMASK_SET(patid);
			continue;
		}
	}
	assert_eq(0, sink.retainedHits().size());
} else {
	// There is no provisional inexact match for the
	// reverse-complement read, so try inexact on the
	// forward-oriented read
	ebwtFw.search1MismatchOrBetter(sfw, params,
								   true,   // allow exact hits
								   false); // no provisional hits
	bool hit = sink.numHits() > lastHits;
	// Set a bit indicating this pattern is done and needn't be
	// considered by the 1-mismatch loop
	if(sanity) sanityCheckHits(patFw, sink, patid, true, os, true, false);
	assert_eq(0, sink.retainedHits().size());
	if(hit) lastHits = sink.numHits();
	if(oneHit && hit) {
		// Update doneMask
		assert_eq(0, sink.numProvisionalHits());
		DONEMASK_SET(patid);
		continue;
	}
}
