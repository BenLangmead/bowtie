#!/bin/sh

#
# Downloads sequence for a D. melanogaster from wormbase.  This script
# was used to build the Bowtie index for the D. melanogaster version
# numbered "dmel_r5.11_FB2008_08". 
#

GENOMES_MIRROR=ftp://ftp.flybase.net/genomes/Drosophila_melanogaster
FNAME=dmel-all-chromosome-r5.11.fasta

BOWTIE_BUILD_EXE=./bowtie-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
	if ! which bowtie-build ; then
		echo "Could not find bowtie-build in current directory or in PATH"
		exit 1
	else
		BOWTIE_BUILD_EXE=`which bowtie-build`
	fi
fi

if [ ! -f $FNAME ] ; then
	if ! which wget > /dev/null ; then
		echo wget not found, looking for curl...
		if ! which curl > /dev/null ; then
			echo curl not found either, aborting...
		else
			# Use curl
			curl ${GENOMES_MIRROR}/dmel_r5.11_FB2008_08/fasta/$FNAME.gz -o $FNAME.gz
		fi
	else
		# Use wget
		wget ${GENOMES_MIRROR}/dmel_r5.11_FB2008_08/fasta/$FNAME.gz
	fi
	gunzip $FNAME.gz
fi

if [ ! -f $FNAME ] ; then
	echo "Could not find $FNAME file!"
	exit 2
fi

echo Running ${BOWTIE_BUILD_EXE} $FNAME d_melanogaster
${BOWTIE_BUILD_EXE} $FNAME d_melanogaster
if [ "$?" = "0" ] ; then
	echo "d_melanogaster index built:"
	echo "   d_melanogaster.1.ebwt d_melanogaster.2.ebwt"
	echo "   d_melanogaster.3.ebwt d_melanogaster.4.ebwt"
	echo "   d_melanogaster.rev.1.ebwt d_melanogaster.rev.2.ebwt"
	echo "You may remove $FNAME"
else
	echo "Index building failed; see error message"
fi
