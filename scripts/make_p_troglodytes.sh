#!/bin/sh

#
# Downloads sequence for P. troglodytes (chimpanzee) from NCBI.  This
# script was used to build the Bowtie index for P. troglodytes.
#
# From README_CURRENT_BUILD:
#  Organism: Pan troglodytes (chimpanzee)
#  NCBI Build Number: 2    
#  Version: 1
#  Release date: 4 October 2006
#  Freeze date for component genomic sequences: March 2006
#  Freeze date for other genomic sequences: 21 April 2006
#  Freeze date for mRNAs/ESTs used for annotation: 20 July 2006
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

for c in 1 2A 2B 3 4 5 6 7 8 9 ; do
	if [ ! -f ptr_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/Pan_troglodytes/CHR_0$c/ptr_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/Pan_troglodytes/CHR_0$c/ptr_ref_chr$c.mfa.gz
		fi
		gunzip ptr_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f ptr_ref_chr$c.mfa ] ; then
		echo "Could not find ptr_ref_chr$c.mfa file!"
		exit 2
	fi
done

for c in 10 11 12 13 14 15 16 17 18 19 20 21 22 MT X Y  ; do
	if [ ! -f ptr_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/Pan_troglodytes/CHR_$c/ptr_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/Pan_troglodytes/CHR_$c/ptr_ref_chr$c.mfa.gz
		fi
		gunzip ptr_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f ptr_ref_chr$c.mfa ] ; then
		echo "Could not find ptr_ref_chr$c.mfa file!"
		exit 2
	fi
done

INPUTS=ptr_ref_chr1.mfa,ptr_ref_chr2A.mfa,ptr_ref_chr2B.mfa,ptr_ref_chr3.mfa,ptr_ref_chr4.mfa,ptr_ref_chr5.mfa,ptr_ref_chr6.mfa,ptr_ref_chr7.mfa,ptr_ref_chr8.mfa,ptr_ref_chr9.mfa,ptr_ref_chr10.mfa,ptr_ref_chr11.mfa,ptr_ref_chr12.mfa,ptr_ref_chr13.mfa,ptr_ref_chr14.mfa,ptr_ref_chr15.mfa,ptr_ref_chr16.mfa,ptr_ref_chr17.mfa,ptr_ref_chr18.mfa,ptr_ref_chr19.mfa,ptr_ref_chr20.mfa,ptr_ref_chr21.mfa,ptr_ref_chr22.mfa,ptr_ref_chrMT.mfa,ptr_ref_chrX.mfa,ptr_ref_chrY.mfa

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} p_troglodytes
${BOWTIE_BUILD_EXE} ${INPUTS} p_troglodytes
if [ "$?" = "0" ] ; then
	echo "p_troglodytes index built:"
	echo "   p_troglodytes.1.ebwt p_troglodytes.2.ebwt"
	echo "   p_troglodytes.3.ebwt p_troglodytes.4.ebwt"
	echo "   p_troglodytes.rev.1.ebwt p_troglodytes.rev.2.ebwt"
	echo "You may remove ptr_ref_chr*.mfa"
else
	echo "Index building failed; see error message"
fi
