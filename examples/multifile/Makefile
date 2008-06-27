BUILD  = ../../ebwt_build-with-asserts
SEARCH = ../../ebwt_search-with-asserts

results: pat1.fa pat2.fa text.ebwt $(SEARCH)
	$(SEARCH) -s text.ebwt pat1.fa,pat2.fa > $@

text.ebwt: text1.fa text2.fa $(BUILD)
	$(BUILD) -s text1.fa,text2.fa $@

$(BUILD):
	$(MAKE) -C ../.. ebwt_build-with-asserts

$(SEARCH):
	$(MAKE) -C ../.. ebwt_search-with-asserts

.PHONY: clean
clean:
	rm -f text.ebwt results
