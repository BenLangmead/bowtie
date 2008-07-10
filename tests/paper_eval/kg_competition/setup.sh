#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
BOWTIE_HOME=$HOME/workspace/bowtie
KG_READS=/fs/szasmg/langmead/reads/SRR001115/s_7_0000_0255
MAKE_LINKS=1
GENOME_DIR=../..

if [ "$MAKE_LINKS" = "1" ] ; then
	# Make link for FASTA sequence
	ln -s -f ${GENOME_DIR}/hs_ref_${NAME}.mfa hs_ref_${NAME}.mfa
	if ! wc -c hs_ref_${NAME}.mfa 2> /dev/null > /dev/null ; then
	    echo "Broken link for hs_ref_${NAME}.mfa; aborting..."
	    exit 1
	fi
	
	# Make link for Maq .bfa file
	ln -s -f ${GENOME_DIR}/hs_ref_${NAME}.bfa hs_ref_${NAME}.bfa
	if ! wc -c hs_ref_${NAME}.bfa 2> /dev/null > /dev/null ; then
	    echo "Broken link for hs_ref_${NAME}.bfa; aborting..."
	    exit 1
	fi
else
	if [ ! -f hs_ref_${NAME}.mfa ] ; then
		cp ${GENOME_DIR}/hs_ref_${NAME}.mfa .
	fi
	if [ ! -f hs_ref_${NAME}.bfa ] ; then
		cp ${GENOME_DIR}/hs_ref_${NAME}.bfa .
	fi
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
if [ "$MAKE_LINKS" = "1" ] ; then
	ln -s -f ${KG_READS}.fastq kg_reads.fq
	ln -s -f ${KG_READS}.fa kg_reads.fa
else
	if [ ! -f kg_reads.fq ] ; then cp ${KG_READS}.fastq kg_reads.fq ; fi
	if [ ! -f kg_reads.fa ] ; then cp ${KG_READS}.fa kg_reads.fa ; fi
fi

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
