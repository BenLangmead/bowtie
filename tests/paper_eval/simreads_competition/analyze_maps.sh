#!/bin/sh

export TOT_READS=8000000

if [ -f whole.ebwt.hits ] ; then
	if [ ! -f whole.ebwt.reads.mapped ] ; then
		awk '{print $1}' whole.ebwt.hits > whole.ebwt.reads.mapped
	fi
	if [ ! -f whole.ebwt.reads.mapped.uniq ] ; then
		sort -u whole.ebwt.reads.mapped > whole.ebwt.reads.mapped.uniq
	fi
	num=`wc -l whole.ebwt.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l whole.ebwt.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Bowtie: Num hits: $num, num unique hits: $numuniq!"
		exit 1
	fi
	echo -n "Bowtie: % reads mapped: "
	perl -e "print $num * 100.0 / $TOT_READS"
	echo
fi

if [ -f whole.maq.map ] ; then
	if [ ! -f whole.maq.reads.mapped ] ; then
		maq mapview whole.maq.map | awk '{print $1}' > whole.maq.reads.mapped
	fi
	if [ ! -f whole.maq.reads.mapped.uniq ] ; then
		sort -u whole.maq.reads.mapped > whole.maq.reads.mapped.uniq
	fi
	num=`wc -l whole.maq.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l whole.maq.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Maq: Num hits: $num, num unique hits: $numuniq!"
		exit 1
	fi
	echo -n "Maq: % reads mapped: "
	perl -e "print $num * 100.0 / $TOT_READS"
	echo
fi

if [ -f whole.maq.n1.map ] ; then
	if [ ! -f whole.maq.n1.reads.mapped ] ; then
		maq mapview whole.maq.n1.map | awk '{print $1}' > whole.maq.n1.reads.mapped
	fi
	if [ ! -f whole.maq.n1.reads.mapped.uniq ] ; then
		sort -u whole.maq.n1.reads.mapped > whole.maq.n1.reads.mapped.uniq
	fi
	num=`wc -l whole.maq.n1.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l whole.maq.n1.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Maq -n 1: Num hits: $num, num unique hits: $numuniq!"
		exit 1
	fi
	echo -n "Maq -n 1: % reads mapped: "
	perl -e "print $num * 100.0 / $TOT_READS"
	echo
fi
