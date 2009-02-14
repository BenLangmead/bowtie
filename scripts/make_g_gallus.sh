#!/bin/sh

#
# Downloads sequence for G. gallus (chicken) from NCBI.  This script
# was used to build the Bowtie index for G. gallus.
#
# Organism: Gallus gallus (chicken)
# NCBI Build Number: 2    
# Version: 1
# Release date: 30 November 2006
# Freeze date for component genomic sequences: May 2006
# Freeze date for other genomic sequences: 19 October 2006
# Freeze date for mRNAs/ESTs used for annotation: 13 October 2006
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
	if [ ! -f gga_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/Gallus_gallus/CHR_0$c/gga_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/Gallus_gallus/CHR_0$c/gga_ref_chr$c.mfa.gz
		fi
		gunzip gga_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f gga_ref_chr$c.mfa ] ; then
		echo "Could not find gga_ref_chr$c.mfa file!"
		exit 2
	fi
done

# Note: the CHR_Y directory is empty as of today
for c in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 MT W Z ; do
	if [ ! -f gga_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/Gallus_gallus/CHR_$c/gga_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/Gallus_gallus/CHR_$c/gga_ref_chr$c.mfa.gz
		fi
		gunzip gga_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f gga_ref_chr$c.mfa ] ; then
		echo "Could not find gga_ref_chr$c.mfa file!"
		exit 2
	fi
done

INPUTS=gga_ref_chr1.mfa,gga_ref_chr2.mfa,gga_ref_chr3.mfa,gga_ref_chr4.mfa,gga_ref_chr5.mfa,gga_ref_chr6.mfa,gga_ref_chr7.mfa,gga_ref_chr8.mfa,gga_ref_chr9.mfa,gga_ref_chr10.mfa,gga_ref_chr11.mfa,gga_ref_chr12.mfa,gga_ref_chr13.mfa,gga_ref_chr14.mfa,gga_ref_chr15.mfa,gga_ref_chr16.mfa,gga_ref_chr17.mfa,gga_ref_chr18.mfa,gga_ref_chr19.mfa,gga_ref_chr20.mfa,gga_ref_chr21.mfa,gga_ref_chr22.mfa,gga_ref_chr23.mfa,gga_ref_chr24.mfa,gga_ref_chr25.mfa,gga_ref_chr26.mfa,gga_ref_chr27.mfa,gga_ref_chr28.mfa,gga_ref_chrMT.mfa,gga_ref_chrW.mfa,gga_ref_chrZ.mfa

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} g_gallus
${BOWTIE_BUILD_EXE} ${INPUTS} g_gallus
if [ "$?" = "0" ] ; then
	echo "g_gallus index built:"
	echo "   g_gallus.1.ebwt g_gallus.2.ebwt"
	echo "   g_gallus.3.ebwt g_gallus.4.ebwt"
	echo "   g_gallus.rev.1.ebwt g_gallus.rev.2.ebwt"
	echo "You may remove gga_ref_chr*.mfa"
else
	echo "Index building failed; see error message"
fi
