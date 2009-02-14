#!/bin/sh

#
# Downloads sequence for a S. cerevisiae from CYGD.  This script
# was used to build the Bowtie index for S. cerevisiae.
#

GENOMES_MIRROR=ftp://ftp.jgi-psf.org/pub/JGI_data

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

if [ ! -f fugu.041029.scaffolds.fasta ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/Fugu/v4.0/fugu.041029.scaffolds.fasta.gz -o fugu.041029.scaffolds.fasta.gz
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/Fugu/v4.0/fugu.041029.scaffolds.fasta.gz
	fi
	gunzip fugu.041029.scaffolds.fasta.gz
fi

if [ ! -f fugu.041029.scaffolds.fasta ] ; then
	echo "Could not find fugu.041029.scaffolds.fasta file!"
	exit 2
fi

echo Running ${BOWTIE_BUILD_EXE} fugu.041029.scaffolds.fasta t_rubripes
${BOWTIE_BUILD_EXE} fugu.041029.scaffolds.fasta t_rubripes
if [ "$?" = "0" ] ; then
	echo "t_rubripes index built:"
	echo "   t_rubripes.1.ebwt t_rubripes.2.ebwt"
	echo "   t_rubripes.3.ebwt t_rubripes.4.ebwt"
	echo "   t_rubripes.rev.1.ebwt t_rubripes.rev.2.ebwt"
	echo "You may remove fugu.041029.scaffolds.fasta"
else
	echo "Index building failed; see error message"
fi
