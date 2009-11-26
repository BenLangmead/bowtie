/*
 * This is a fragment, included from multiple places in ebwt_search.cpp.
 * It implements the logic of the second phase of the 1-mismatch search
 * routine.  It is implemented as a code fragment so that it can be
 * reused in both the half-index-in-memory and full-index-in-memory
 * situations.
 */
{
	bt.setEbwt(&ebwtBw);
	bt.setReportExacts(false);

	if(!norc) {
		params.setFw(false);
		// Next, try hits with one mismatch on the 3' end for the reverse-complement read
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, s3, s, s, s); // 1 mismatch allowed in 3' half
		if(bt.backtrack()) {
			continue;
		}
	}

	if(!nofw) {
		params.setFw(true);
		// Next, try hits with one mismatch on the 3' end for the reverse-complement read
		bt.setQuery(patsrc->bufa());
		bt.setOffs(0, 0, s3, s, s, s); // 1 mismatch allowed in 3' half
		if(bt.backtrack()) {
			continue;
		}
	}
}
