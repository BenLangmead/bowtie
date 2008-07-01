#!/bin/sh

for d in `ls *.top` ; do
	echo $d
	e=`echo $d | sed 's/[^.]*[.]//'`
	e=`echo $e | sed 's/[.].*//'`
	perl summarize_top.pl $d $e | tail -4
done
