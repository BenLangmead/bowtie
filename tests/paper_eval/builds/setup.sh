#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
BOWTIE_HOME=$HOME/workspace/bowtie
KG_READS=/fs/szasmg/langmead/reads/SRR001115/s_7_0000_0255
GENOMES_DIR=/fs/szasmg/langmead

# Make link for FASTA sequence
if [ ! -f hs_ref_${NAME}.mfa ] ; then
	ln -s -f ${GENOMES_DIR}/hs_ref_${NAME}.mfa hs_ref_${NAME}.mfa
fi

# Check that link isn't broken
if ! wc -c hs_ref_${NAME}.mfa 2> /dev/null > /dev/null ; then
    echo "Broken link for hs_ref_${NAME}.mfa; aborting..."
    exit 1
fi

# Build and copy ebwt_build and ebwt_build_packed to here
make -C ${BOWTIE_HOME} ebwt_build
cp ${BOWTIE_HOME}/ebwt_build .
if [ ! -f ebwt_build ] ; then
	echo "Could not build ebwt_build; aborting..."
	exit 1
fi
make -C ${BOWTIE_HOME} ebwt_build_packed
cp ${BOWTIE_HOME}/ebwt_build_packed .
if [ ! -f ebwt_build_packed ] ; then
	echo "Could not build ebwt_build_packed; aborting..."
	exit 1
fi

# Copy analysis scripts from bowtie dir
cp ${BOWTIE_HOME}/scripts/summarize_top.pl .
cp ${BOWTIE_HOME}/scripts/summarize_all_top.sh .
cp ${BOWTIE_HOME}/scripts/wrap.sh .
