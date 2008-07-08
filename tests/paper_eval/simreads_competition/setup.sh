#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}

# Make link for FASTA sequence
ln -s -f ../../hs_ref_${NAME}.mfa hs_ref_${NAME}.mfa
if ! wc -c hs_ref_${NAME}.mfa 2> /dev/null > /dev/null ; then
    echo "Broken link for hs_ref_${NAME}.mfa; aborting..."
    exit 1
fi

# Make link for Maq .bfa file
ln -s -f ../../hs_ref_${NAME}.bfa hs_ref_${NAME}.bfa
if ! wc -c hs_ref_${NAME}.bfa 2> /dev/null > /dev/null ; then
    echo "Broken link for hs_ref_${NAME}.bfa; aborting..."
    exit 1
fi

# Simulate 8M 35-bp reads if necessary
if [ ! -f ${NAME}_sim_8000000.fa ] ; then
    rm -f ${NAME}_sim*.bfq
    make -C ~/workspace/bowtie simreads
    cp ~/workspace/bowtie/simreads .
    ./simreads \
	-r 8000000 \
	-l 35 \
	hs_ref_${NAME}.mfa \
	${NAME}_sim_8000000.fa \
	${NAME}_sim_8000000.fq
fi

ln -s -f /fs/szasmg/langmead/ebwts/${NAME}.1.ebwt ${NAME}.1.ebwt
ln -s -f /fs/szasmg/langmead/ebwts/${NAME}.2.ebwt ${NAME}.2.ebwt
ln -s -f /fs/szasmg/langmead/ebwts/${NAME}.3.ebwt ${NAME}.3.ebwt
ln -s -f /fs/szasmg/langmead/ebwts/${NAME}.rev.1.ebwt ${NAME}.rev.1.ebwt
ln -s -f /fs/szasmg/langmead/ebwts/${NAME}.rev.2.ebwt ${NAME}.rev.2.ebwt
ln -s -f /fs/szasmg/langmead/ebwts/${NAME}.rev.3.ebwt ${NAME}.rev.3.ebwt

err=0
if ! wc -c ${NAME}.1.ebwt 2> /dev/null > /dev/null ; then err=1 ; fi
if ! wc -c ${NAME}.2.ebwt 2> /dev/null > /dev/null ; then err=1 ; fi
if ! wc -c ${NAME}.rev.1.ebwt 2> /dev/null > /dev/null ; then err=1 ; fi
if ! wc -c ${NAME}.rev.2.ebwt 2> /dev/null > /dev/null ; then err=1 ; fi

if [ "$err" = "0" ] ; then
    echo "One more more ebwt links are invalid; aborting..."
    exit 1
fi

# Convert fq file to a set of bfq files with 2M reads each
if [ ! -f ${NAME}_sim@6000001.bfq ] ; then
    maq fastq2bfq -n 2000000 ${NAME}_sim_8000000.fq ${NAME}_sim.bfq
fi

# Convert fq to one big bfq file
if [ ! -f ${NAME}_sim_8000000.bfq ] ; then
    maq fastq2bfq ${NAME}_sim_8000000.fq ${NAME}_sim_8000000.bfq
fi
