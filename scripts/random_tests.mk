# Run from base directory as "make -f scripts/random_tests.mk" (add -jN
# if multiple cores are available) 

SEED=0
ITERS=10000

OUTS = random_bsa_tests.out \
       random_qsort_tests.out \
       random_ebwt_tests.out \
       diff_sample_tests.out \
       random_diff_sample_tests.out

all: $(OUTS)

random_bsa_tests.out: blockwise_sa
	perl scripts/random_bsa_tester.pl $(SEED) $(ITERS) > $@ 2>&1

blockwise_sa:
	$(MAKE) $@

random_qsort_tests.out: multikey_qsort
	perl scripts/random_qsort_tester.pl $(SEED) $(ITERS) > $@ 2>&1

multikey_qsort:
	$(MAKE) $@

random_ebwt_tests.out: ebwt_build-with-asserts ebwt_search-with-asserts
	perl scripts/random_ebwt_tester.pl $(SEED) $(ITERS) > $@ 2>&1

ebwt_build-with-asserts:
	$(MAKE) $@

ebwt_search-with-asserts:
	$(MAKE) $@

diff_sample_tests.out: diff_sample-with-asserts
	sh scripts/diff_sample_tester.sh > $@ 2>&1

random_diff_sample_tests.out: diff_sample-with-asserts
	perl scripts/random_dc_tester.pl $(SEED) $(ITERS) > $@ 2>&1

diff_sample-with-asserts:
	$(MAKE) $@

.PHONY: clean
clean:
	rm -f $(OUTS)
