#!/bin/sh

#
# Downloads sequence for the galGal3 version of G. gallus (chicken)
# from UCSC.
#

BASE_CHRS=chr1
# Add autosomes 2-28
i=2
while [ $i -lt 29 ] ; do
	BASE_CHRS="$BASE_CHRS chr$i"
	i=`expr $i + 1`
done
# Add autosome 32
BASE_CHRS="$BASE_CHRS chr32"
# Add sex chromosomes
BASE_CHRS="$BASE_CHRS chrW chrZ"
# Add mitochondrial
BASE_CHRS="$BASE_CHRS chrM"
# Add other unplaced, but autosome-associated
BASE_CHRS="$BASE_CHRS chr1_random chr2_random chr4_random chr5_random chr6_random chr7_random chr8_random chr10_random chr11_random chr12_random chr13_random chr16_random chr17_random chr18_random chr20_random chr22_random chr25_random chr28_random"
# Add other unplaced, but sex-chromosome-associated
BASE_CHRS="$BASE_CHRS chrW_random chrZ_random"
# Add totally unplaced
BASE_CHRS="$BASE_CHRS chrUn_random"

CHRS_TO_INDEX=$BASE_CHRS

FTP_BASE=ftp://hgdownload.cse.ucsc.edu/goldenPath/galGal3/chromosomes

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

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
for c in $CHRS_TO_INDEX ; do
	if [ ! -f ${c}.fa ] ; then
		F=${c}.fa.gz
		get ${FTP_BASE}/$F || (echo "Error getting $F" && exit 1)
		gunzip $F || (echo "Error unzipping $F" && exit 1)
	fi
	[ -n "$INPUTS" ] && INPUTS=$INPUTS,${c}.fa
	[ -z "$INPUTS" ] && INPUTS=${c}.fa
done

CMD="${BOWTIE_BUILD_EXE} $* ${INPUTS} galGal3"
echo Running $CMD
if $CMD ; then
	echo "galGal3 index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
