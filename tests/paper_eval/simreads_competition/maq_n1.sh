#!/bin/sh

#
# Usage:   maq map [options] <out.map> <chr.bfa> <reads_1.bfq> [reads_2.bfq]
#
# Options: -1 INT      length of the first read (<64) [0]
#          -2 INT      length of the second read (<64) [0]
#          -m FLOAT    rate of difference between reads and references [0.001]
#          -e INT      maximum allowed sum of qualities of mismatches [70]
#          -d FILE     adapter sequence file [null]
#          -a INT      max distance between two paired reads [250]
#          -A INT      max distance between two RF paired reads [0]
#          -n INT      number of mismatches in the first 24bp [2]
#          -M c|g      methylation alignment mode [null]
#          -u FILE     dump unmapped and poorly aligned reads to FILE [null]
#          -H FILE     dump multiple/all 01-mismatch hits to FILE [null]
#          -C INT      max number of hits to output. >512 for all 01 hits. [250]
#          -s INT      seed for random number generator [random]
#          -N          record mismatch positions (max read length<=55)
#          -t          trim all reads (usually not recommended)
#          -c          match in the colorspace
#

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}

# This is version 0.6.6
MAQ=/fs/sz-user-supported/Linux-x86_64/bin/maq
MAQ_ARGS="-n 1"
EXTRA_EXT=".n1"
REF=hs_ref_${NAME}.bfa
READ_BASE=${NAME}_sim

# Maq on split-up read set where each unit has 2M reads, as per Heng Li's suggestion
if [ ! -f ${NAME}.maq.1${EXTRA_EXT}.map -a `wc -c whole.maq.1${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
   echo > ${NAME}.maq${EXTRA_EXT}.top
   sh wrap.sh ${NAME}.maq${EXTRA_EXT} \
      $MAQ map $MAQ_ARGS \
         ${NAME}.maq.1${EXTRA_EXT}.map \
         ${REF} \
         ${READ_BASE}\@1.bfq
   if [ ! -f ${NAME}.maq.1${EXTRA_EXT}.map -a `wc -c whole.maq.1${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
      echo "Failed to create legitimate map file: ${NAME}.maq.1${EXTRA_EXT}.map; aborting..."
      exit 1
   fi
fi
if [ ! -f ${NAME}.maq.2${EXTRA_EXT}.map -a `wc -c whole.maq.2${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
   sh wrap.sh ${NAME}.maq${EXTRA_EXT} \
      $MAQ map $MAQ_ARGS \
         ${NAME}.maq.2${EXTRA_EXT}.map \
         ${REF} \
         ${READ_BASE}\@2000001.bfq
   if [ ! -f ${NAME}.maq.2${EXTRA_EXT}.map -a `wc -c whole.maq.2${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
      echo "Failed to create legitimate map file: ${NAME}.maq.2${EXTRA_EXT}.map; aborting..."
      exit 1
   fi
fi
if [ ! -f ${NAME}.maq.3${EXTRA_EXT}.map -a `wc -c whole.maq.3${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
   sh wrap.sh ${NAME}.maq${EXTRA_EXT} \
      $MAQ map $MAQ_ARGS \
         ${NAME}.maq.3${EXTRA_EXT}.map \
         ${REF} \
         ${READ_BASE}\@4000001.bfq
   if [ ! -f ${NAME}.maq.3${EXTRA_EXT}.map -a `wc -c whole.maq.3${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
      echo "Failed to create legitimate map file: ${NAME}.maq.3${EXTRA_EXT}.map; aborting..."
      exit 1
   fi
fi
if [ ! -f ${NAME}.maq.4${EXTRA_EXT}.map -a `wc -c whole.maq.4${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
   sh wrap.sh ${NAME}.maq${EXTRA_EXT} \
      $MAQ map $MAQ_ARGS \
         ${NAME}.maq.4${EXTRA_EXT}.map \
         ${REF} \
         ${READ_BASE}\@6000001.bfq 
   if [ ! -f ${NAME}.maq.4${EXTRA_EXT}.map -a `wc -c whole.maq.4${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
      echo "Failed to create legitimate map file: ${NAME}.maq.4${EXTRA_EXT}.map; aborting..."
      exit 1
   fi
fi

# Merge all 2M-read ${EXTRA_EXT}.map files into one big map
if [ ! -f ${NAME}.maq${EXTRA_EXT}.map -a `wc -c whole.maq${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
   maq mapmerge ${NAME}.maq${EXTRA_EXT}.map \
                  ${NAME}.maq.1${EXTRA_EXT}.map \
                  ${NAME}.maq.2${EXTRA_EXT}.map \
                  ${NAME}.maq.3${EXTRA_EXT}.map \
                  ${NAME}.maq.4${EXTRA_EXT}.map
fi

# Maq, all reads in 1 shot, default parameters (2 mm allowed in first 24 bases)
if [ ! -f ${NAME}.maq.all1bfq${EXTRA_EXT}.map -a `wc -c whole.maq.all1bfq${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
   echo > ${NAME}.maq.all1bfq${EXTRA_EXT}.top
   sh wrap.sh ${NAME}.maq.all1bfq${EXTRA_EXT} \
      $MAQ map $MAQ_ARGS \
         ${NAME}.maq.all1bfq${EXTRA_EXT}.map \
         ${REF} \
         ${READ_BASE}.bfq
   if [ ! -f ${NAME}.maq.all1bfq${EXTRA_EXT}.map -a `wc -c whole.maq.all1bfq${EXTRA_EXT}.map | sed 's/ .*//'` -gt 10 ] ; then
      echo "Failed to create legitimate map file: ${NAME}.maq.all1bfq${EXTRA_EXT}.map; aborting..."
      exit 1
   fi
fi
