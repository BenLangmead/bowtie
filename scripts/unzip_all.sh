#!/bin/sh

# Unzips all of the .gz files in the current directory

for d in `ls *.gz` ; do
	zcat < $d > `sed 's/\.gz$//'`
done
