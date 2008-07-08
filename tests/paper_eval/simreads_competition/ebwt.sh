#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}

if [ ! -f ~/workspace/bowtie/ebwt_search ] ; then
    make -C ~/workspace/bowtie ebwt_search
fi

cp ~/workspace/ebwt/ebwt_search .

echo > ${NAME}.ebwt.top
sh wrap.sh ${NAME}.ebwt \
  ./ebwt_search -1tqr ${NAME} ${NAME}_sim_8000000.fq ${NAME}.ebwt.hits
