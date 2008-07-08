#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
BOWTIE_HOME=$HOME/workspace/bowtie
READS=${NAME}_sim_8000000.fq

# Make ebwt_search
if [ ! -f ${BOWTIE_HOME}/ebwt_search ] ; then
   make -C ${BOWTIE_HOME} ebwt_search
   if [ ! -f ${BOWTIE_HOME}/ebwt_search ] ; then
      echo "Failed to build ebwt_search in ${BOWTIE_HOME}; aborting..."
      exit 1
   fi
fi

# Copy ebwt_search to here
cp ${BOWTIE_HOME}/ebwt_search .
./ebwt_search --version

# Run ebwt_search to produce hits
if [ ! -f ${NAME}.ebwt.hits ] ; then
   echo > ${NAME}.ebwt.top
   sh wrap.sh ${NAME}.ebwt \
     ./ebwt_search -1tqr ${NAME} ${READS} ${NAME}.ebwt.hits
fi
