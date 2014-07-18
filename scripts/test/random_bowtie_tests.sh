#!/bin/bash

PREFIX=${0%/*}

SEED="$1"
if [ -z "$SEED" ] ; then
	SEED=77
else
	shift
fi

if [ "$1" == "-c" ] ; then
	make clean
	shift
fi

make allall

# Args: seed, outer, inner, tbase, trand, pbase, prand

echo "Long test emphasizing building..."
perl $PREFIX/random_bowtie_tests.pl $SEED 6000 10 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Short test emphasizing searching..."
perl $PREFIX/random_bowtie_tests.pl $SEED 1000 200 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Short test emphasizing searching with short patterns..."
perl $PREFIX/random_bowtie_tests.pl $SEED 1000 200 300 200 6 6
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Short test emphasizing searching with long patterns..."
perl $PREFIX/random_bowtie_tests.pl $SEED 1000 200 300 200 30 20
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi
echo "Short test emphasizing building..."
perl $PREFIX/random_bowtie_tests.pl $SEED 1200 10 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi

echo "Long test emphasizing searching..."
perl $PREFIX/random_bowtie_tests.pl $SEED 5000 200 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; exit 1 ; fi

