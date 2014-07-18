#!/bin/sh

# Spawn a bunch of parallel processes running random_bowtie_tests.sh
# and dumping results into text files.  Run from the bowtie dir.
PREFIX=${0%/*}

NUM=$1
if [ -z "$NUM" ] ; then
	NUM=4
fi

make allall
if [ "$?" != "0" ] ; then
	echo "Error during build"
	exit 1
fi

echo > .randpids
while [ $NUM -gt 0 ] ; do
	echo "Spawning: sh scripts/random_bowtie_tests.sh $NUM > .randstdout.$NUM 2> .randstderr.$NUM &"
	sh $PREFIX/random_bowtie_tests.sh $NUM > .randstdout.$NUM 2> .randstderr.$NUM &
	echo $! >> .randpids
	NUM=`expr $NUM - 1`
done

for p in `cat .randpids` ; do
	echo "Waiting for pid $p"
	wait $p
done
