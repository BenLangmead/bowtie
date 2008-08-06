#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
BOWTIE_HOME=$HOME/workspace/bowtie
READS=kg_reads.fq
DO_CVS_UPDATE=0
USE_FILTERED_READS=0

# Possibly switch to filtered read set
if [ "$USE_FILTERED_READS" = "1" ] ; then
        READS=kg_reads_filt.fq
        echo "Using filtered reads; READ_BASE is $READ_BASE"
fi

# Optionally do a cvs update in the Bowtie home
if [ "$DO_CVS_UPDATE" = "1" ] ; then
	pushd ${BOWTIE_HOME}
	cvs update -d
	popd
fi

# Make ebwt_search
make -C ${BOWTIE_HOME} ebwt_search
if [ ! -f ${BOWTIE_HOME}/ebwt_search ] ; then
   echo "Failed to build ebwt_search in ${BOWTIE_HOME}; aborting..."
   exit 1
fi

# Copy ebwt_search to here
cp ${BOWTIE_HOME}/ebwt_search .
./ebwt_search --version

if [ ! -f ${READS} ] ; then
	echo "Could not find reads file ${READS}!  Aborting..."
	exit 1
fi

if [ ! -f ${NAME}.1.ebwt ] ; then
	echo "Could not ebwt index file ${NAME}.1.ebwt!  Aborting..."
	exit 1
fi

# Run ebwt_search in Maq-like mode to produce hits
if [ ! -f ${NAME}.ebwt.hits ] ; then
   echo > ${NAME}.ebwt.top
   sh wrap.sh ${NAME}.ebwt \
     ./ebwt_search -tqr ${NAME} ${READS} ${NAME}.ebwt.hits
else
	echo "${NAME}.ebwt.hits already exists; skipping Maq-mode run"
fi

# Run ebwt_search in Maq -n 1 mode to produce hits
if [ ! -f ${NAME}.ebwt.n1.hits ] ; then
   echo > ${NAME}.ebwt.n1.top
   sh wrap.sh ${NAME}.ebwt.n1 \
     ./ebwt_search -n 1 -tqr ${NAME} ${READS} ${NAME}.ebwt.n1.hits
else
	echo "${NAME}.ebwt.n1.hits already exists; skipping Maq-mode -n 1 run"
fi
