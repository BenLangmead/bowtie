#!/bin/sh

#
# Downloads sequence for the mm9 version of M. musculus (mouse) from
# UCSC.  This script was used to build the mm9 Bowtie index.
#
# From README:
#
# This directory contains the Jul. 2007 assembly of the mouse genome
# (mm9, NCBI Build 37) in one gzip-compressed FASTA file per chromosome.
#
# This assembly was produced by the Mouse Genome Sequencing Consortium,
# and the National Center for Biotechnology Information (NCBI).
# See also: http://www.ncbi.nlm.nih.gov/mapview/map_search.cgi?taxid=10090
#

UCSC_MM9_BASE=ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

INPUTS=
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M; do
	if [ ! -f chr${c}.fa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${UCSC_MM9_BASE}/chr${c}.fa.gz
			fi
		else
			# Use wget
			wget ${UCSC_MM9_BASE}/chr${c}.fa.gz
		fi
		gunzip chr${c}.fa.gz
	fi
	
	if [ ! -f chr${c}.fa ] ; then
		echo "Could not find chr${c}.fa file!"
		exit 2
	fi
	if [ -n "$INPUTS" ] ; then
		INPUTS=$INPUTS,chr${c}.fa
	else
		INPUTS=chr${c}.fa
	fi
done

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} mm9
${BOWTIE_BUILD_EXE} ${INPUTS} mm9
if [ "$?" = "0" ] ; then
	echo "mm9 index built:"
	echo "   mm9.ebwt mm9.2.ebwt"
	echo "   mm9.3.ebwt mm9.4.ebwt"
	echo "   mm9.rev.1.ebwt mm9.rev.2.ebwt"
	echo "You may remove chrXX.fa"
else
	echo "Index building failed; see error message"
fi
