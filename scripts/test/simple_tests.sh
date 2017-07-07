#!/bin/sh

#  simple_tests.sh

make $* bowtie-align-s \
		bowtie-align-l \
		bowtie-align-s-debug \
		bowtie-align-l-debug \
		bowtie-build-s \
		bowtie-build-l \
		bowtie-build-s-debug \
		bowtie-build-l-debug && \
perl scripts/test/simple_tests.pl \
	--bowtie=./bowtie \
	--bowtie-build=./bowtie-build
