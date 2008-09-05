#!/bin/sh

#
# Downloads sequence for C. elegans from wormbase.  This script
# was used to build the Bowtie index for C. elegans on August 14,
# 2008.
#

GENOMES_MIRROR=ftp://ftp.arabidopsis.org/home/tair

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

for c in 1 2 3 4 5 C M ; do
	if [ ! -f chr$c.fas ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/Sequences/whole_chromosomes/chr$c.fas -o chr$c.fas
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/Sequences/whole_chromosomes/chr$c.fas
		fi
	fi
	
	if [ ! -f chr$c.fas ] ; then
		echo "Could not find chr$c.fas file!"
		exit 2
	fi
done

echo Running ${BOWTIE_BUILD_EXE} chr1.fas,chr2.fas,chr3.fas,chr4.fas,chr5.fas,chrC.fas,chrM.fas  c_elegans
${BOWTIE_BUILD_EXE} chr1.fas,chr2.fas,chr3.fas,chr4.fas,chr5.fas,chrC.fas,chrM.fas  c_elegans
