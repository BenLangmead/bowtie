#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
BOWTIE_HOME=$HOME/workspace/bowtie
KG_READS=/fs/szasmg/langmead/reads/SRR001115/s_7_0000_0255

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

# Make links for ebwt files
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
# (allow the .3.ebwt files to be absent for now)

if [ "$err" != "0" ] ; then
    echo "One more more ebwt links are invalid; aborting..."
    exit 1
fi

# Make link to 1000-Genomes reads
ln -s -f ${KG_READS}.fastq kg_reads.fq
ln -s -f ${KG_READS}.fa kg_reads.fa

# Convert fq file to a set of bfq files with 2M reads each
if [ ! -f kg_reads@6000001.bfq ] ; then
    maq fastq2bfq -n 2000000 kg_reads.fq kg_reads.bfq
fi

# Convert fq to one big bfq file
if [ ! -f kg_reads.bfq ] ; then
    maq fastq2bfq kg_reads.fq kg_reads.bfq
fi

# Copy analysis scripts from bowtie dir
cp ${BOWTIE_HOME}/scripts/summarize_top.pl .
cp ${BOWTIE_HOME}/scripts/summarize_all_top.sh .
