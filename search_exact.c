/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the exact-search routine.  It is
 * implemented as a code fragment so that it can be reused in both
 * paired and unpaired alignment.
 */
{
	uint32_t plen = length(patsrc->bufa().patFw);
	if(!nofw) {
		// Match against forward strand
		params.setFw(true);
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, plen, plen, plen, plen);
		// If we matched on the forward strand, ignore the reverse-
		// complement strand
		if(bt.backtrack()) {
			continue;
		}
	}
	if(!norc) {
		// Process reverse-complement read
		params.setFw(false);
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, plen, plen, plen, plen);
		bt.backtrack();
	}
}
