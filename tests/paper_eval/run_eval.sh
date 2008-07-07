#!/bin/bash

BOWTIE=../../
BUILD=$BOWTIE/ebwt_build
#PYTHON=/fs/sz-user-supported/Linux-i686/bin/python2.5
PYTHON=/usr/bin/python
AWKDIR=/Users/cole/scripts/Awk
ECHO="echo -e"

REF=ref.mfa

# Generate a reference in several contigs
$PYTHON rand_reads.py 1 50000 > $REF 2>/dev/null
$PYTHON rand_reads.py 1 25000 >> $REF 2>/dev/null
$PYTHON rand_reads.py 1 35000 >> $REF 2>/dev/null
$PYTHON rand_reads.py 1 1000 >> $REF 2>/dev/null
awk 'BEGIN{id=0}{if($1==">0") printf ">%d\n", id++; else print}' $REF > xxx
mv xxx $REF

awk 'BEGIN{id = 0} {if (substr ($1, 1, 1) == ">") { print id++, substr($1, 2) }}' ref.mfa > ref.ids

# Make an ebwt out of it
$BUILD --entireSA -d $REF ref

# Grab a single contig from the reference for simulated sequencing.
# Note that we can ~cole/cut-fasta.awk here if we want finer granularity 
# targeted sequencing
TARGET='>2'
awk -f $AWKDIR/extract-fasta.awk $TARGET $REF > target.fa

# Use a human-like SNP rate of 1/800 bp
MISMATCH_RATE=0.00125
DELETION_RATE=0.0
INSERTION_RATE=0.0
$PYTHON polymorph.py -m $MISMATCH_RATE -d $DELETION_RATE -i $INSERTION_RATE target.fa > target_snpped.fa 2>target_snpped.transcript

# Generate simulated reads for the target region
READ_LEN=30
# Desired number of reads is (Desired coverage) * (ref len) / (read_len)  
NUM_READS=5250
$BOWTIE/simreads -r $NUM_READS -l $READ_LEN target_snpped.fa reads.fa reads.fq

# Run maq (in easymaq) mode to obtain SNP calls for the simulated data
MAQ_OUT=maq_out
if [ -e $MAQ_OUT ] 
    then
    rm -rf $MAQ_OUT
fi

maq.pl easyrun -d $MAQ_OUT $REF reads.fq
cd $MAQ_OUT
python ../snp_eval.py ../target_snpped.transcript cns.snp > raw_snp_calls.txt
python ../snp_eval.py ../target_snpped.transcript cns.final.snp > final_snp_calls.txt
cd ..

# Now map the reads with Bowtie
BWT_OUT=bwt_out
if [ -e $BWT_OUT ] 
    then
    rm -rf $BWT_OUT
fi

mkdir bwt_out
cd $BWT_OUT
../../../ebwt_search  -1 -q -3 8 ../ref ../reads.fq > reads1@1.bwtmap
../$BOWTIE/bowtie_convert reads1@1.bwtmap reads1@1.map ../ref.ids
maq mapcheck -s ../maq_out/ref.bfa reads1@1.map >mapcheck.txt
maq assemble -s consensus.cns ../maq_out/ref.bfa reads1@1.map 2>assemble.log
maq cns2fq consensus.cns >cns.fq
maq cns2snp consensus.cns >cns.snp
maq.pl SNPfilter -d3 cns.snp > cns.filter.snp
awk '$5>=30' cns.filter.snp > cns.final.snp
python ../snp_eval.py ../target_snpped.transcript cns.snp > raw_snp_calls.txt
python ../snp_eval.py ../target_snpped.transcript cns.final.snp > final_snp_calls.txt

