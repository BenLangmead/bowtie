#!/bin/sh

#
# Downloads sequence for C. familiaris (dog) from NCBI.  This script
# was used to build the Bowtie index for C. familiaris.
#
# From README_CURRENT_BUILD:
#  Organism: Canis familiaris (dog)
#  NCBI Build Number: 2    
#  Version: 1
#  Release date: 8 September 2005
#  Freeze date for component genomic sequences: May 2005
#  Freeze date for other genomic sequences: 28 July 2005
#  Freeze date for mRNAs/ESTs used for annotation: 27 July 2005
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
	if [ ! -f cfa_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/Canis_familiaris/CHR_0$c/cfa_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/Canis_familiaris/CHR_0$c/cfa_ref_chr$c.mfa.gz
		fi
		gunzip cfa_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f cfa_ref_chr$c.mfa ] ; then
		echo "Could not find cfa_ref_chr$c.mfa file!"
		exit 2
	fi
done

# Note: the CHR_Y directory is empty as of today
for c in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 MT X  ; do
	if [ ! -f cfa_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/Canis_familiaris/CHR_$c/cfa_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/Canis_familiaris/CHR_$c/cfa_ref_chr$c.mfa.gz
		fi
		gunzip cfa_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f cfa_ref_chr$c.mfa ] ; then
		echo "Could not find cfa_ref_chr$c.mfa file!"
		exit 2
	fi
done

INPUTS=cfa_ref_chr1.mfa,cfa_ref_chr2.mfa,cfa_ref_chr3.mfa,cfa_ref_chr4.mfa,cfa_ref_chr5.mfa,cfa_ref_chr6.mfa,cfa_ref_chr7.mfa,cfa_ref_chr8.mfa,cfa_ref_chr9.mfa,cfa_ref_chr10.mfa,cfa_ref_chr11.mfa,cfa_ref_chr12.mfa,cfa_ref_chr13.mfa,cfa_ref_chr14.mfa,cfa_ref_chr15.mfa,cfa_ref_chr16.mfa,cfa_ref_chr17.mfa,cfa_ref_chr18.mfa,cfa_ref_chr19.mfa,cfa_ref_chr20.mfa,cfa_ref_chr21.mfa,cfa_ref_chr22.mfa,cfa_ref_chr23.mfa,cfa_ref_chr24.mfa,cfa_ref_chr25.mfa,cfa_ref_chr26.mfa,cfa_ref_chr27.mfa,cfa_ref_chr28.mfa,cfa_ref_chr29.mfa,cfa_ref_chr30.mfa,cfa_ref_chr31.mfa,cfa_ref_chr32.mfa,cfa_ref_chr33.mfa,cfa_ref_chr34.mfa,cfa_ref_chr35.mfa,cfa_ref_chr36.mfa,cfa_ref_chr37.mfa,cfa_ref_chr38.mfa,cfa_ref_chrMT.mfa,cfa_ref_chrX.mfa

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} c_familiaris
${BOWTIE_BUILD_EXE} ${INPUTS} c_familiaris
if [ "$?" = "0" ] ; then
	echo "c_familiaris index built:"
	echo "   c_familiaris.1.ebwt c_familiaris.2.ebwt"
	echo "   c_familiaris.3.ebwt c_familiaris.4.ebwt"
	echo "   c_familiaris.rev.1.ebwt c_familiaris.rev.2.ebwt"
	echo "You may remove cfa_ref_chr*.mfa"
else
	echo "Index building failed; see error message"
fi
