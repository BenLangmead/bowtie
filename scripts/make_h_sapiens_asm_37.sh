#!/bin/sh

#
# Downloads sequence for H. sapiens (human) from NCBI.  This script was
# used to build the Bowtie index for H. sapiens 37.
#
# From README_CURRENT_BUILD:
# Organism: Homo sapiens (human)
# NCBI Build Number: 37    
# Version: 1
# Release date: 04 August 2009
#

GENOMES_MIRROR=ftp://ftp.ncbi.nih.gov/genomes

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
BASE_NAME=hs_ref_GRCh37_
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
	if [ ! -f ${BASE_NAME}chr$c.fa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/H_sapiens/Assembled_chromosomes/${BASE_NAME}chr$c.fa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/H_sapiens/Assembled_chromosomes/${BASE_NAME}chr$c.fa.gz
		fi
		gunzip ${BASE_NAME}chr$c.fa.gz
	fi
	
	if [ ! -f ${BASE_NAME}chr$c.fa ] ; then
		echo "Could not find ${BASE_NAME}chr$c.fa file!"
		exit 2
	fi
	if [ -z "$INPUTS" ] ; then
	    INPUTS=${BASE_NAME}chr$c.fa
	else
	    INPUTS=$INPUTS,${BASE_NAME}chr$c.fa
	fi
done

# Special case: get mitochondrial DNA from its home
if [ ! -f hs_ref_chrMT.fa ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/H_sapiens/CHR_MT/hs_ref_chrMT.fa.gz
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/H_sapiens/CHR_MT/hs_ref_chrMT.fa.gz
	fi
	gunzip hs_ref_chrMT.fa.gz
fi

INPUTS=$INPUTS,hs_ref_chrMT.fa

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} h_sapiens_37_asm
${BOWTIE_BUILD_EXE} ${INPUTS} h_sapiens_37_asm

if [ "$?" = "0" ] ; then
	echo "h_sapiens_37_asm index built:"
	echo "   h_sapiens_37_asm.1.ebwt h_sapiens_37_asm.2.ebwt"
	echo "   h_sapiens_37_asm.3.ebwt h_sapiens_37_asm.4.ebwt"
	echo "   h_sapiens_37_asm.rev.1.ebwt h_sapiens_37_asm.rev.2.ebwt"
	echo "You may remove hs_ref_chr*.fa"
else
	echo "Index building failed; see error message"
fi
