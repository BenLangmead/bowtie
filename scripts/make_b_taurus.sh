#!/bin/sh

#
# Builds an index from UMD Freeze 2.0 of the Bos Taurus (cow) genome.
#

BOWTIE_BUILD_EXE=./bowtie-build

INPUTS="/fs/szasmg3/bos_taurus/UMD_Freeze2.0/Chr1.fa"
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 X U ; do
	INPUTS="${INPUTS},/fs/szasmg3/bos_taurus/UMD_Freeze2.0/Chr${chr}.fa"
done

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} b_taurus
${BOWTIE_BUILD_EXE} ${INPUTS} b_taurus

if [ "$?" = "0" ] ; then
	echo "b_taurus index built:"
	echo "   b_taurus.1.ebwt b_taurus.2.ebwt"
	echo "   b_taurus.3.ebwt b_taurus.4.ebwt"
	echo "   b_taurus.rev.1.ebwt b_taurus.rev.2.ebwt"
else
	echo "Index building failed; see error message"
fi
