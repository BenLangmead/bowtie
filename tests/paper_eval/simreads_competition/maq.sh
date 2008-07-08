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

# Maq on split-up read set where each unit has 2M reads, as per Heng Li's suggestion
if [ ! -f ${NAME}.maq.1.map ] ; then
   echo > ${NAME}.maq.top
   sh wrap.sh ${NAME}.maq \
      $MAQ map \
         ${NAME}.maq.1.map \
         hs_ref_${NAME}.bfa \
         ${NAME}_sim\@1.bfq
fi
if [ ! -f ${NAME}.maq.1.map ] ; then
   sh wrap.sh ${NAME}.maq \
      $MAQ map \
         ${NAME}.maq.2.map \
         hs_ref_${NAME}.bfa \
         ${NAME}_sim\@2000001.bfq
fi
if [ ! -f ${NAME}.maq.1.map ] ; then
   sh wrap.sh ${NAME}.maq \
      $MAQ map \
         ${NAME}.maq.3.map \
         hs_ref_${NAME}.bfa \
         ${NAME}_sim\@4000001.bfq
fi
if [ ! -f ${NAME}.maq.1.map ] ; then
   sh wrap.sh ${NAME}.maq \
      $MAQ map \
         ${NAME}.maq.4.map \
         hs_ref_${NAME}.bfa \
         ${NAME}_sim\@6000001.bfq 
fi

if [ ! -f ${NAME}.maq.map ] ; then
   maq mapmerge ${NAME}.maq.map \
                ${NAME}.maq.1.map \
                ${NAME}.maq.2.map \
                ${NAME}.maq.3.map \
                ${NAME}.maq.4.map
fi

# Maq on split-up read set where each unit has 2M reads, as per Heng Li's suggestion
if [ ! -f ${NAME}.maq.n1.map ] ; then
    echo > ${NAME}.maq.n1.top
    sh wrap.sh ${NAME}.maq.n1 \
	$MAQ map -n 1 \
	  ${NAME}.maq.n1.1.map \
	  hs_ref_${NAME}.bfa \
	  ${NAME}_sim\@1.bfq
    sh wrap.sh ${NAME}.maq.n1 \
	$MAQ map -n 1 \
	  ${NAME}.maq.n1.2.map \
	  hs_ref_${NAME}.bfa \
	  ${NAME}_sim\@2000001.bfq
    sh wrap.sh ${NAME}.maq.n1 \
	$MAQ map -n 1 \
	  ${NAME}.maq.n1.3.map \
	  hs_ref_${NAME}.bfa \
	  ${NAME}_sim\@4000001.bfq
    sh wrap.sh ${NAME}.maq.n1 \
	$MAQ map -n 1 \
	  ${NAME}.maq.n1.4.map \
	  hs_ref_${NAME}.bfa \
	  ${NAME}_sim\@6000001.bfq 

    maq mapmerge ${NAME}.maq.n1.map ${NAME}.maq.n1.1.map \
	                            ${NAME}.maq.n1.2.map \
	                            ${NAME}.maq.n1.3.map \
	                            ${NAME}.maq.n1.4.map
fi

# Maq, all reads in 1 shot, default parameters (2 mm allowed in first 24 bases)
if [ ! -f ${NAME}.maq.all1bfq.map ] ; then
   echo > ${NAME}.maq.all1bfq.top
   sh wrap.sh ${NAME}.maq.all1bfq \
      $MAQ map \
         ${NAME}.maq.all1bfq.map \
         hs_ref_${NAME}.bfa \
         ${NAME}_sim_8000000.bfq
fi

# Maq, all reads in 1 shot, with -n 1 (only 1 mm allowed in first 24 bases)
if [ ! -f ${NAME}.maq.n1.all1bfq.top ] ; then
    echo > ${NAME}.maq.n1.all1bfq.top
    sh wrap.sh ${NAME}.maq.n1.all1bfq \
	$MAQ map -n 1 \
	  ${NAME}.maq.n1.all1bfq.map \
	  hs_ref_${NAME}.bfa \
	  ${NAME}_sim_8000000.bfq
fi
