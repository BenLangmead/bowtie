#!/bin/sh

#
# Downloads sequence for C. elegans from wormbase.  This script
# was used to build the Bowtie index for C. elegans on September 6,
# 2008.
#

GENOMES_MIRROR=ftp://ftp.gramene.org/pub/wormbase/genomes

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

if [ ! -f current.dna.fa ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/c_elegans/sequences/dna/current.dna.fa.gz
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/c_elegans/sequences/dna/current.dna.fa.gz
	fi
	gunzip current.dna.fa.gz
fi

if [ ! -f current.dna.fa ] ; then
	echo "Could not find current.dna.fa file!"
	exit 2
fi

echo Running ${BOWTIE_BUILD_EXE} current.dna.fa c_elegans
${BOWTIE_BUILD_EXE} current.dna.fa c_elegans
