#!/bin/sh

dir=`pwd`
while [ $# -ge 1 ]; do
	name=`$1 | sed 's/_.*//'`
	echo "Doing dataset $name in directory ../$1"
	cd ../$1
	sh summarize_all_top.sh > $dir/$name.results.txt
	cd $dir
	shift
done
