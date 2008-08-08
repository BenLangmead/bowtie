#!/bin/sh

DO_FILTERED=1

WORKSTATION=1
if [ `hostname` = "privet.umiacs.umd.edu" ] ; then
	WORKSTATION=0
fi
if [ `hostname` = "larch.umiacs.umd.edu" ] ; then
	WORKSTATION=0
fi

if [ ! -f whole.numreads.em.both ] ; then
	comm -1 -2 whole.ebwt.reads.mapped.uniq whole.maq.reads.mapped.uniq | wc -l > whole.numreads.em.both
fi
if [ ! -f whole.numreads.em.ebwt ] ; then
	comm -2 -3 whole.ebwt.reads.mapped.uniq whole.maq.reads.mapped.uniq | wc -l > whole.numreads.em.ebwt
fi
if [ ! -f whole.numreads.em.maq ] ; then
	comm -1 -3 whole.ebwt.reads.mapped.uniq whole.maq.reads.mapped.uniq | wc -l > whole.numreads.em.maq
fi

# Do the Bowtie/Maq comparison
inboth=`cat whole.numreads.em.both`
inebwt=`cat whole.numreads.em.ebwt`
inmaq=`cat whole.numreads.em.maq`
perl plot_reads.pl Maq maq_reads $inboth $inebwt $inmaq

if [ "$DO_FILTERED" = "1" ] ; then
	if [ ! -f whole.numreads.emf.both ] ; then
		comm -1 -2 whole.ebwt.filt.reads.mapped.uniq whole.maq.filt.reads.mapped.uniq | wc -l > whole.numreads.emf.both
	fi
	if [ ! -f whole.numreads.emf.ebwt ] ; then
		comm -2 -3 whole.ebwt.filt.reads.mapped.uniq whole.maq.filt.reads.mapped.uniq | wc -l > whole.numreads.emf.ebwt
	fi
	if [ ! -f whole.numreads.emf.maq ] ; then
		comm -1 -3 whole.ebwt.filt.reads.mapped.uniq whole.maq.filt.reads.mapped.uniq | wc -l > whole.numreads.emf.maq
	fi
	
	# Do the Bowtie/SOAP comparison
	inboth=`cat whole.numreads.emf.both`
	inebwt=`cat whole.numreads.emf.ebwt`
	inmaq=`cat whole.numreads.emf.maq`
	perl plot_reads.pl Maq maq_reads_filt $inboth $inebwt $inmaq
fi

if [ "$WORKSTATION" = "0" ] ; then
	if [ ! -f whole.numreads.es.both ] ; then
		comm -1 -2 whole.ebwt.2.reads.mapped.uniq whole.soap.v2.reads.mapped.uniq | wc -l > whole.numreads.es.both
	fi
	if [ ! -f whole.numreads.es.ebwt ] ; then
		comm -2 -3 whole.ebwt.2.reads.mapped.uniq whole.soap.v2.reads.mapped.uniq | wc -l > whole.numreads.es.ebwt
	fi
	if [ ! -f whole.numreads.es.soap ] ; then
		comm -1 -3 whole.ebwt.2.reads.mapped.uniq whole.soap.v2.reads.mapped.uniq | wc -l > whole.numreads.es.soap
	fi
	
	# Do the Bowtie/SOAP comparison
	inboth=`cat whole.numreads.es.both`
	inebwt=`cat whole.numreads.es.ebwt`
	insoap=`cat whole.numreads.es.soap`
	perl plot_reads.pl Soap soap_reads $inboth $inebwt $insoap
fi
