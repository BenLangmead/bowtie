#!/bin/sh

# Run summarize_top.pl on each .top file in the current directory,
# assuming that the executable string to grep from the .top files is
# the string between the first two dots in the filename
for d in `ls *.top` ; do
	echo -n "$d: "
	# e = everything after first .
	e=`echo $d | sed 's/[^.]*[.]//'`
	# e = everything before first . (everything between first two dots)
	e=`echo $e | sed 's/[.].*//'`
	perl summarize_top.pl $d $e | tail -1
done
