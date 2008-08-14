#
# Downloads sequence for M. musculus (mouse) from NCBI.  This script
# was used to build the Bowtie index for M. musculus on August 14,
# 2008.
#
# From README_CURRENT_BUILD:
#  Organism: Mus musculus (mouse)
#  NCBI Build Number: 37
#  Version: 1
#  Release date: 05 July 2007
#
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
	if [ ! -f mm_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/M_musculus/CHR_0$c/mm_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/M_musculus/CHR_0$c/mm_ref_chr$c.mfa.gz
		fi
		gunzip mm_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f mm_ref_chr$c.mfa ] ; then
		echo "Could not find mm_ref_chr$c.mfa file!"
		exit 2
	fi
done

for c in 10 11 12 13 14 15 16 17 18 19 X Y  ; do
	if [ ! -f mm_ref_chr$c.mfa ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/M_musculus/CHR_$c/mm_ref_chr$c.mfa.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/M_musculus/CHR_$c/mm_ref_chr$c.mfa.gz
		fi
		gunzip mm_ref_chr$c.mfa.gz
	fi
	
	if [ ! -f mm_ref_chr$c.mfa ] ; then
		echo "Could not find mm_ref_chr$c.mfa file!"
		exit 2
	fi
done

if [ ! -f mm_ref_chrMT.fa ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/M_musculus/CHR_MT/mm_ref_chrMT.fa.gz
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/M_musculus/CHR_MT/mm_ref_chrMT.fa.gz
	fi
	gunzip mm_ref_chrMT.fa.gz
fi

if [ ! -f mm_ref_chrMT.fa ] ; then
	echo "Could not find mm_ref_chrMT.fa file!"
	exit 2
fi

echo Running ${BOWTIE_BUILD_EXE} mm_ref_chr1.mfa,mm_ref_chr2.mfa,mm_ref_chr3.mfa,mm_ref_chr4.mfa,mm_ref_chr5.mfa,mm_ref_chr6.mfa,mm_ref_chr7.mfa,mm_ref_chr8.mfa,mm_ref_chr9.mfa,mm_ref_chr10.mfa,mm_ref_chr11.mfa,mm_ref_chr12.mfa,mm_ref_chr13.mfa,mm_ref_chr14.mfa,mm_ref_chr15.mfa,mm_ref_chr16.mfa,mm_ref_chr17.mfa,mm_ref_chr18.mfa,mm_ref_chr19.mfa,mm_ref_chrMT.fa,mm_ref_chrX.mfa,mm_ref_chrY.mfa m_musculus
${BOWTIE_BUILD_EXE} mm_ref_chr1.mfa,mm_ref_chr2.mfa,mm_ref_chr3.mfa,mm_ref_chr4.mfa,mm_ref_chr5.mfa,mm_ref_chr6.mfa,mm_ref_chr7.mfa,mm_ref_chr8.mfa,mm_ref_chr9.mfa,mm_ref_chr10.mfa,mm_ref_chr11.mfa,mm_ref_chr12.mfa,mm_ref_chr13.mfa,mm_ref_chr14.mfa,mm_ref_chr15.mfa,mm_ref_chr16.mfa,mm_ref_chr17.mfa,mm_ref_chr18.mfa,mm_ref_chr19.mfa,mm_ref_chrMT.fa,mm_ref_chrX.mfa,mm_ref_chrY.mfa m_musculus
