#!/bin/sh

dir=`pwd`
while [ $# -ge 1 ]; do
	name=`echo $1 | sed 's/_.*//'`
	echo "Doing dataset $name in directory ../$1"
	cd ../$1
	sh summarize_all_top.sh > $dir/$name.results.txt
	cd $dir
	awk '{print $1}' $name.results.txt > $name.results.names.txt
	# Wall clock time
	awk '{print $2}' $name.results.txt | cut -d, -f 2 > $name.results.times.txt
	# Max VM footprint
	awk '{print $2}' $name.results.txt | cut -d, -f 3 > $name.results.vmmax.txt
	# Max RES footprint
	awk '{print $2}' $name.results.txt | cut -d, -f 4 > $name.results.rsmax.txt
	shift
done

perl plot.pl chr22 chr2 whole
