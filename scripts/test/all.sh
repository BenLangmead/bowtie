#!/bin/sh

# Run from bowtie directory
# Usually should do 'make allall' first

d=`dirname $0`

for i in `ls $d/*.pl` ; do
	if ! perl $i ; then
		echo "Error running $i; aborting..."
		exit 1
	fi
done
echo "ALL PASSED"
