#!/bin/sh

#
# Downloads sequence for a D. melanogaster from wormbase.  This script
# was used to build the Bowtie index for D. melanogaster on August 14,
# 2008.
#

GENOMES_MIRROR=ftp://ftp.flybase.net/genomes/Drosophila_melanogaster

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

if [ ! -f dmel-all-chromosome-r5.9.fasta ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/dmel_r5.9_FB2008_06/fasta/dmel-all-chromosome-r5.9.fasta.gz -o dmel-all-chromosome-r5.9.fasta.gz
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/dmel_r5.9_FB2008_06/fasta/dmel-all-chromosome-r5.9.fasta.gz
	fi
	gunzip dmel-all-chromosome-r5.9.fasta.gz
fi

if [ ! -f dmel-all-chromosome-r5.9.fasta ] ; then
	echo "Could not find dmel-all-chromosome-r5.9.fasta file!"
	exit 2
fi

echo Running ${BOWTIE_BUILD_EXE} dmel-all-chromosome-r5.9.fasta d_melanogaster
${BOWTIE_BUILD_EXE} dmel-all-chromosome-r5.9.fasta d_melanogaster


