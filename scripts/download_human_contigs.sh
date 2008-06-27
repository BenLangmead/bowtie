#!/bin/sh

# download_human_contigs.pl
#
# Usage: perl download_human_contigs.pl <dest_dir>
#
# Must have 'wget'; if using a Mac, you can get wget from fink
# (possibly also MacPorts - I'm not sure)
#

if [ ! -d $1 ] ; then
	echo "No such directory: $d"
	exit 1
fi
d=$1
if [ -z "$d" ] ; then
	d="."
fi

mkdir -p $d
pushd $d

ctg_site=ftp://ftp.ncbi.nih.gov/genomes/H_sapiens

for d in 1 2 3 4 5 6 7 8 9
do
    wget --passive-ftp ${ctg_site}/CHR_0$d/hs_ref_chr$d.mfa.gz
done

for d in 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    wget --passive-ftp ${ctg_site}/CHR_$d/hs_ref_chr$d.mfa.gz
done
popd
