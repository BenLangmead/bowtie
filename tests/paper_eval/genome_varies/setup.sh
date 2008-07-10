#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}

# Copy analysis scripts from bowtie dir
cp ${BOWTIE_HOME}/scripts/summarize_top.pl .
cp ${BOWTIE_HOME}/scripts/summarize_all_top.sh .
