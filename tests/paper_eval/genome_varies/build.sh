#!/bin/sh

BOWTIE_HOME=$HOME/workspace/bowtie
GENOMES_HOME=../..

# Grab the latest and greatest ebwt_build
make -C ${BOWTIE_HOME} ebwt_build
cp ${BOWTIE_HOME}/ebwt_build .

pid1=0
pid2=0
pid3=0
pid4=0

# Build all 4 non-whole-human indexes in parallel
 
if [ ! -f /fs/szasmg/langmead/ebwts/whole.1435609516.1.ebwt ] ; then
   ./ebwt_build --cutoff 1435609516 -d --bmaxDivN 6 \
      ${GENOMES_HOME}/hs_ref_all.mfa /fs/szasmg/langmead/ebwts/whole.1435609516 &
   pid1=$!
fi
if [ ! -f /fs/szasmg/langmead/ebwts/whole.717804758.1.ebwt ] ; then
   ./ebwt_build --cutoff 717804758 -d --bmaxDivN 6 \
      ${GENOMES_HOME}/hs_ref_all.mfa /fs/szasmg/langmead/ebwts/whole.717804758 &
   pid2=$!
fi
if [ ! -f /fs/szasmg/langmead/ebwts/whole.358902379.1.ebwt ] ; then
   ./ebwt_build --cutoff 358902379 -d --bmaxDivN 6 \
      ${GENOMES_HOME}/hs_ref_all.mfa /fs/szasmg/langmead/ebwts/whole.358902379 &
   pid3=$!
fi
if [ ! -f /fs/szasmg/langmead/ebwts/whole.179451189.1.ebwt ] ; then
   ./ebwt_build --cutoff 179451189 -d --bmaxDivN 6 \
      ${GENOMES_HOME}/hs_ref_all.mfa /fs/szasmg/langmead/ebwts/whole.179451189 &
   pid4=$!
fi

if [ "$pid1" != "0" ] ; then wait $pid1 ; fi
if [ "$pid2" != "0" ] ; then wait $pid2 ; fi
if [ "$pid3" != "0" ] ; then wait $pid3 ; fi
if [ "$pid4" != "0" ] ; then wait $pid4 ; fi
