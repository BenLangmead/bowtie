#!/bin/sh

#
# Downloads sequence for A. mellifera (honey bee) from the Honeybee
# Genome Sequencing Consortium web site.
#
# Assembly version is: Amel_4.0 (March 10, 2006)
#
# The version downloaded is the "linear contigs" version, which, as far
# as I can tell, is the closest thing to an "assembled" honeybee genome.
#

GENOMES_MIRROR=ftp://ftp.hgsc.bcm.tmc.edu/pub/data/Amellifera/fasta/Amel20060310-freeze/linearScaffolds/

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
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 Un ; do
	name=Group${c}_20060310.fa
	if [ ! -f $name ] ; then
		if ! which wget > /dev/null ; then
			echo wget not found, looking for curl...
			if ! which curl > /dev/null ; then
				echo curl not found either, aborting...
			else
				# Use curl
				curl ${GENOMES_MIRROR}/$name.gz -o $name.gz
			fi
		else
			# Use wget
			wget ${GENOMES_MIRROR}/$name.gz
		fi
		gunzip $name.gz
	fi
	
	if [ ! -f $name ] ; then
		echo "Could not find name file!"
		exit 2
	fi
	if [ -n "$INPUTS" ] ; then
		INPUTS=$INPUTS,$name
	else
		INPUTS=$name
	fi
done

echo Running ${BOWTIE_BUILD_EXE} ${INPUTS} a_mellifera
${BOWTIE_BUILD_EXE} ${INPUTS} a_mellifera
if [ "$?" = "0" ] ; then
	echo "a_mellifera index built:"
	echo "   a_mellifera.1.ebwt a_mellifera.2.ebwt"
	echo "   a_mellifera.3.ebwt a_mellifera.4.ebwt"
	echo "   a_mellifera.rev.1.ebwt a_mellifera.rev.2.ebwt"
else
	echo "Index building failed; see error message"
fi
