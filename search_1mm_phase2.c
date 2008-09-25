/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the second phase of the 1-mismatch search
 * routine.  It is implemented as a code fragment so that it can be
 * reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
params.setFw(true);
params.setEbwtFw(false);
sbw.newQuery(&patFw, &name, &qualFw);
ebwtBw.search1MismatchOrBetter(sbw, params,
							   false,  // no exact hits
							   false); // no provisional hits
// Check all hits against a naive oracle
assert_eq(0, sink.numProvisionalHits());
if(sanity) sanityCheckHits(patFw, sink, patid, true, os, false, true);
assert_eq(0, sink.retainedHits().size());
// If the forward direction matched with one mismatch, ignore
// the reverse complement
if(oneHit && revcomp && sink.numHits() > lastHits) {
	lastHits = sink.numHits();
	continue;
}
if(!revcomp) continue;
params.setFw(false);
sbw.newQuery(&patRc, &name, &qualRc);
ebwtBw.search1MismatchOrBetter(sbw, params,
							   false,  // no exact hits
							   false); // no provisional hits
assert_eq(0, sink.numProvisionalHits());
if(sanity) sanityCheckHits(patRc, sink, patid, false, os, false, true);
assert_eq(0, sink.retainedHits().size());
params.setFw(true);
lastHits = sink.numHits();
