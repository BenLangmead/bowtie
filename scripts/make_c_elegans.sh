#!/bin/sh

#
# Downloads sequence for C. elegans from wormbase.  This script
# was used to build the Bowtie index for C. elegans.  WS190 was the
# latest freeze as of the build date.
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

if [ ! -f elegans.WS190.dna.fa ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/c_elegans/sequences/dna/elegans.WS190.dna.fa.gz
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/c_elegans/sequences/dna/elegans.WS190.dna.fa.gz
	fi
	gunzip elegans.WS190.dna.fa.gz
fi

if [ ! -f elegans.WS190.dna.fa ] ; then
	echo "Could not find elegans.WS190.dna.fa file!"
	exit 2
fi

echo Running ${BOWTIE_BUILD_EXE} elegans.WS190.dna.fa c_elegans
${BOWTIE_BUILD_EXE} elegans.WS190.dna.fa c_elegans
if [ "$?" = "0" ] ; then
	echo "c_elegans index built:"
	echo "   c_elegans.1.ebwt c_elegans.2.ebwt"
	echo "   c_elegans.3.ebwt c_elegans.4.ebwt"
	echo "   c_elegans.rev.1.ebwt c_elegans.rev.2.ebwt"
	echo "You may remove elegans.WS190.dna.fa"
else
	echo "Index building failed; see error message"
fi
