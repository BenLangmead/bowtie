#!/bin/sh

#
# Downloads sequence for H. sapiens (human) from NCBI.  This script was
# used to build the Bowtie index for H. sapiens.
#
# From README_CURRENT_BUILD:
#  Organism: Homo sapiens (human)
#  NCBI Build Number: 36    
#  Version: 3
#  Release date: 24 March 2008
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

for c in 1 2 3 4 5 6 7 8 9 ; do
	if [ ! -f hs_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/H_sapiens/CHR_0$c/hs_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/H_sapiens/CHR_0$c/hs_ref_chr$c.mfa.gz
		fi
		gunzip hs_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f hs_ref_chr$c.mfa ] ; then
		echo "Could not find hs_ref_chr$c.mfa file!"
		exit 2
	fi
done

for c in 10 11 12 13 14 15 16 17 18 19 20 21 22 MT X Y  ; do
	if [ ! -f hs_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/H_sapiens/CHR_$c/hs_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/H_sapiens/CHR_$c/hs_ref_chr$c.mfa.gz
		fi
		gunzip hs_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f hs_ref_chr$c.mfa ] ; then
		echo "Could not find hs_ref_chr$c.mfa file!"
		exit 2
	fi
done

INPUTS=hs_ref_chr1.mfa,hs_ref_chr2.mfa,hs_ref_chr3.mfa,hs_ref_chr4.mfa,hs_ref_chr5.mfa,hs_ref_chr6.mfa,hs_ref_chr7.mfa,hs_ref_chr8.mfa,hs_ref_chr9.mfa,hs_ref_chr10.mfa,hs_ref_chr11.mfa,hs_ref_chr12.mfa,hs_ref_chr13.mfa,hs_ref_chr14.mfa,hs_ref_chr15.mfa,hs_ref_chr16.mfa,hs_ref_chr17.mfa,hs_ref_chr18.mfa,hs_ref_chr19.mfa,hs_ref_chr20.mfa,hs_ref_chr21.mfa,hs_ref_chr22.mfa,hs_ref_chrMT.mfa,hs_ref_chrX.mfa,hs_ref_chrY.mfa

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} h_sapiens
${BOWTIE_BUILD_EXE} ${INPUTS} h_sapiens
if [ "$?" = "0" ] ; then
	echo "h_sapiens index built:"
	echo "   h_sapiens.1.ebwt h_sapiens.2.ebwt"
	echo "   h_sapiens.3.ebwt h_sapiens.4.ebwt"
	echo "   h_sapiens.rev.1.ebwt h_sapiens.rev.2.ebwt"
	echo "You may remove hs_ref_chr*.mfa"
else
	echo "Index building failed; see error message"
fi
