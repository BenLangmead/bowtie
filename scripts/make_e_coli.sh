#!/bin/sh

#
# Downloads the sequence for a strain of e. coli from NCBI and builds a
# Bowtie index for it
#

GENOMES_MIRROR=ftp://ftp.ncbi.nlm.nih.gov/genomes

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

if [ ! -f NC_008253.fna ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/Bacteria/Escherichia_coli_536/NC_008253.fna -o NC_008253.fna
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/Bacteria/Escherichia_coli_536/NC_008253.fna
	fi
fi

if [ ! -f NC_008253.fna ] ; then
	echo "Could not find NC_008253.fna file!"
	exit 2
fi

echo Running ${BOWTIE_BUILD_EXE} -t 8 NC_008253.fna e_coli
${BOWTIE_BUILD_EXE} -t 8 NC_008253.fna e_coli
if [ "$?" = "0" ] ; then
	echo "e_coli index built:"
	echo "   e_coli.1.ebwt e_coli.2.ebwt"
	echo "   e_coli.3.ebwt e_coli.4.ebwt"
	echo "   e_coli.rev.1.ebwt e_coli.rev.2.ebwt"
	echo "You may remove NC_008253.fna"
else
	echo "Index building failed; see error message"
fi
