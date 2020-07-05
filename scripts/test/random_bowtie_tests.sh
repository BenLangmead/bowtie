#!/bin/sh

PREFIX=${0%/*}

SEED="$1"
if [ -z "$SEED" ] ; then
	SEED=77
else
	shift
fi

MAKE=make
gmake -v > /dev/null 2>&1
if [ $? -eq 0 ] ; then
    MAKE=gmake
fi

if [ "$1" == "-c" ] ; then
	$MAKE clean
	shift
fi

if [ -z "$BT2_PATH" -a "$USE_BT2_INDEX" == "1" ]; then
    git clone --recursive "https://github.com/BenLangmead/bowtie2.git"
    cd bowtie2 && $MAKE bowtie2-build-s-debug bowtie2-build-s \
                        bowtie2-build-l-debug bowtie2-build-l \
                        NO_TBB=1
    if [ $? -ne 0 ]; then
        echo "Unable to compile bowtie2 build binaries"
    fi
    cd ..
    BT2_PATH="`pwd`/bowtie2"
fi

$MAKE allall "$@"
MAKE_ARGS="$*"
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

