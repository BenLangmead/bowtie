#!/bin/sh

if [ ! -f whole.numreads.both ] ; then
	comm -1 -2 whole.ebwt.reads.mapped.uniq whole.maq.reads.mapped.uniq | wc -l > whole.numreads.both
fi
if [ ! -f whole.numreads.ebwt ] ; then
	comm -2 -3 whole.ebwt.reads.mapped.uniq whole.maq.reads.mapped.uniq | wc -l > whole.numreads.ebwt
fi
if [ ! -f whole.numreads.maq ] ; then
	comm -1 -3 whole.ebwt.reads.mapped.uniq whole.maq.reads.mapped.uniq | wc -l > whole.numreads.maq
fi

inboth=`cat whole.numreads.both`
inebwt=`cat whole.numreads.ebwt`
inmaq=`cat whole.numreads.maq`
perl plot_reads.pl $inboth $inebwt $inmaq
