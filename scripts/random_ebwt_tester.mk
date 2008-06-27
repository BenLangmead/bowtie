all: random_ebwt_tester.1.out \
     random_ebwt_tester.2.out \
     random_ebwt_tester.3.out \
     random_ebwt_tester.4.out \
     random_ebwt_tester.5.out \
     random_ebwt_tester.6.out \
     random_ebwt_tester.7.out \
     random_ebwt_tester.8.out

random_ebwt_tester.%.out: ebwt_build-with-asserts ebwt_search-with-asserts
	perl scripts/random_ebwt_tester.pl > $@

ebwt_build-with-asserts:
	$(MAKE) ebwt_build-with-asserts

ebwt_search-with-asserts:
	$(MAKE) ebwt_search-with-asserts
