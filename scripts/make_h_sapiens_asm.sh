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

for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
	if [ ! -f hs_ref_chr$c.fa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/H_sapiens/Assembled_chromosomes/hs_ref_chr$c.fa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/H_sapiens/Assembled_chromosomes/hs_ref_chr$c.fa.gz
		fi
		gunzip hs_ref_chr$c.fa.gz
	fi
	
	if [ ! -f hs_ref_chr$c.fa ] ; then
		echo "Could not find hs_ref_chr$c.fa file!"
		exit 2
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

#echo Running ${BOWTIE_BUILD_EXE} --bmaxdivn 4 hs_ref_chr1.fa,hs_ref_chr2.fa,hs_ref_chr3.fa,hs_ref_chr4.fa,hs_ref_chr5.fa,hs_ref_chr6.fa,hs_ref_chr7.fa,hs_ref_chr8.fa,hs_ref_chr9.fa,hs_ref_chr10.fa,hs_ref_chr11.fa,hs_ref_chr12.fa,hs_ref_chr13.fa,hs_ref_chr14.fa,hs_ref_chr15.fa,hs_ref_chr16.fa,hs_ref_chr17.fa,hs_ref_chr18.fa,hs_ref_chr19.fa,hs_ref_chr20.fa,hs_ref_chr21.fa,hs_ref_chr22.fa,hs_ref_chrMT.fa,hs_ref_chrX.fa,hs_ref_chrY.fa h_sapiens_asm
#${BOWTIE_BUILD_EXE} --bmaxdivn 4 hs_ref_chr1.fa,hs_ref_chr2.fa,hs_ref_chr3.fa,hs_ref_chr4.fa,hs_ref_chr5.fa,hs_ref_chr6.fa,hs_ref_chr7.fa,hs_ref_chr8.fa,hs_ref_chr9.fa,hs_ref_chr10.fa,hs_ref_chr11.fa,hs_ref_chr12.fa,hs_ref_chr13.fa,hs_ref_chr14.fa,hs_ref_chr15.fa,hs_ref_chr16.fa,hs_ref_chr17.fa,hs_ref_chr18.fa,hs_ref_chr19.fa,hs_ref_chr20.fa,hs_ref_chr21.fa,hs_ref_chr22.fa,hs_ref_chrMT.fa,hs_ref_chrX.fa,hs_ref_chrY.fa h_sapiens_asm
cat hs_ref_chr1.fa hs_ref_chr2.fa hs_ref_chr3.fa hs_ref_chr4.fa hs_ref_chr5.fa hs_ref_chr6.fa hs_ref_chr7.fa hs_ref_chr8.fa hs_ref_chr9.fa hs_ref_chr10.fa hs_ref_chr11.fa hs_ref_chr12.fa hs_ref_chr13.fa hs_ref_chr14.fa hs_ref_chr15.fa hs_ref_chr16.fa hs_ref_chr17.fa hs_ref_chr18.fa hs_ref_chr19.fa hs_ref_chr20.fa hs_ref_chr21.fa hs_ref_chr22.fa hs_ref_chrMT.fa hs_ref_chrX.fa hs_ref_chrY.fa > hs_ref_all.fa
maq fasta2bfa hs_ref_all.fa h_sapiens_asm.bfa
if [ "$?" = "0" ] ; then
	echo "h_sapiens index built:"
	echo "   h_sapiens_asm.1.ebwt h_sapiens_asm.2.ebwt h_sapiens_asm.rev.1.ebwt h_sapiens_asm.rev.2.ebwt"
	echo "You may remove hs_ref_chr*.fa"
else
	echo "Index building failed; see error message"
fi
