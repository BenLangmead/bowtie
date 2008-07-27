#!/bin/sh

BOWTIE_HOME=$HOME/workspace/bowtie
GENOMES_HOME=../..

# Grab the latest and greatest simreads
make -C ${BOWTIE_HOME} simreads
cp ${BOWTIE_HOME}/simreads .

pid0=0
pid1=0
pid2=0
pid3=0
pid4=0

# Simulate all 5 sets of reads in parallel
 
if [ ! -f whole.sim_reads.fq ] ; then
   ./simreads -r 4000000 \
      ${GENOMES_HOME}/hs_ref_all.mfa whole.sim_reads.fa \
                         whole.sim_reads.fq &
   pid0=$!
fi
if [ ! -f whole.1435609516.sim_reads.fq ] ; then
   ./simreads -r 4000000 -c 1435609516 \
      ${GENOMES_HOME}/hs_ref_all.mfa whole.1435609516.sim_reads.fa \
                         whole.1435609516.sim_reads.fq &
   pid1=$!
fi
if [ ! -f whole.717804758.sim_reads.fq ] ; then
   ./simreads -r 4000000 -c 717804758 \
      ${GENOMES_HOME}/hs_ref_all.mfa whole.717804758.sim_reads.fa \
                         whole.717804758.sim_reads.fq &
   pid2=$!
fi
if [ ! -f whole.358902379.sim_reads.fq ] ; then
   ./simreads -r 4000000 -c 358902379 \
      ${GENOMES_HOME}/hs_ref_all.mfa whole.358902379.sim_reads.fa \
                         whole.358902379.sim_reads.fq &
   pid3=$!
fi
if [ ! -f whole.179451189.sim_reads.fq ] ; then
   ./simreads -r 4000000 -c 179451189 \
      ${GENOMES_HOME}/hs_ref_all.mfa whole.179451189.sim_reads.fa \
                         whole.179451189.sim_reads.fq &
   pid4=$!
fi

if [ "$pid0" != "0" ] ; then wait $pid0 ; fi
if [ "$pid1" != "0" ] ; then wait $pid1 ; fi
if [ "$pid2" != "0" ] ; then wait $pid2 ; fi
if [ "$pid3" != "0" ] ; then wait $pid3 ; fi
if [ "$pid4" != "0" ] ; then wait $pid4 ; fi
