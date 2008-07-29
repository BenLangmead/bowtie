#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
export TOT_READS=8956597

if [ -f ${NAME}.ebwt.n1.hits ] ; then
	if [ ! -f ${NAME}.ebwt.n1.reads.mapped ] ; then
		awk '{print $1}' ${NAME}.ebwt.n1.hits > ${NAME}.ebwt.n1.reads.mapped
	fi
	if [ ! -f ${NAME}.ebwt.n1.reads.mapped.uniq ] ; then
		sort -u ${NAME}.ebwt.n1.reads.mapped > ${NAME}.ebwt.n1.reads.mapped.uniq
	fi
	num=`wc -l ${NAME}.ebwt.n1.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l ${NAME}.ebwt.n1.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Bowtie -n 1: Num hits: $num, num unique hits: $numuniq!"
		echo "  Will use $numuniq for % reads mapped calculation"
	fi
	echo -n "Bowtie -n 1: % reads mapped: "
	perl -e "print $numuniq * 100.0 / $TOT_READS"
	echo
else
	echo "Didn't find ${NAME}.ebwt.n1.hits"
fi

if [ -f ${NAME}.ebwt.hits ] ; then
	if [ ! -f ${NAME}.ebwt.reads.mapped ] ; then
		awk '{print $1}' ${NAME}.ebwt.hits > ${NAME}.ebwt.reads.mapped
	fi
	if [ ! -f ${NAME}.ebwt.reads.mapped.uniq ] ; then
		sort -u ${NAME}.ebwt.reads.mapped > ${NAME}.ebwt.reads.mapped.uniq
	fi
	num=`wc -l ${NAME}.ebwt.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l ${NAME}.ebwt.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Bowtie: Num hits: $num, num unique hits: $numuniq!"
		echo "  Will use $numuniq for % reads mapped calculation"
	fi
	echo -n "Bowtie: % reads mapped: "
	perl -e "print $numuniq * 100.0 / $TOT_READS"
	echo
else
	echo "Didn't find ${NAME}.ebwt.hits"
fi

if [ -f ${NAME}.maq.n1.map ] ; then
	if [ ! -f ${NAME}.maq.n1.reads.mapped ] ; then
		maq mapview ${NAME}.maq.n1.map | awk '{print $1}' > ${NAME}.maq.n1.reads.mapped
	fi
	if [ ! -f ${NAME}.maq.n1.reads.mapped.uniq ] ; then
		sort -u ${NAME}.maq.n1.reads.mapped > ${NAME}.maq.n1.reads.mapped.uniq
	fi
	num=`wc -l ${NAME}.maq.n1.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l ${NAME}.maq.n1.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Maq -n 1: Num hits: $num, num unique hits: $numuniq!"
		echo "  Will use $numuniq for % reads mapped calculation"
	fi
	echo -n "Maq -n 1: % reads mapped: "
	perl -e "print $numuniq * 100.0 / $TOT_READS"
	echo
else
	echo "Didn't find ${NAME}.maq.n1.map"
fi

if [ -f ${NAME}.maq.map ] ; then
	if [ ! -f ${NAME}.maq.reads.mapped ] ; then
		maq mapview ${NAME}.maq.map | awk '{print $1}' > ${NAME}.maq.reads.mapped
	fi
	if [ ! -f ${NAME}.maq.reads.mapped.uniq ] ; then
		sort -u ${NAME}.maq.reads.mapped > ${NAME}.maq.reads.mapped.uniq
	fi
	num=`wc -l ${NAME}.maq.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l ${NAME}.maq.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Maq: Num hits: $num, num unique hits: $numuniq!"
		echo "  Will use $numuniq for % reads mapped calculation"
	fi
	echo -n "Maq: % reads mapped: "
	perl -e "print $numuniq * 100.0 / $TOT_READS"
	echo
else
	echo "Didn't find ${NAME}.maq.map"
fi

if [ -f ${NAME}.soap.v1.map ] ; then
	if [ ! -f ${NAME}.soap.v1.reads.mapped ] ; then
		awk '{print $1}' ${NAME}.soap.v1.map > ${NAME}.soap.v1.reads.mapped
	fi
	if [ ! -f ${NAME}.soap.v1.reads.mapped.uniq ] ; then
		sort -u ${NAME}.soap.v1.reads.mapped > ${NAME}.soap.v1.reads.mapped.uniq
	fi
	num=`wc -l ${NAME}.soap.v1.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l ${NAME}.soap.v1.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Bowtie: Num hits: $num, num unique hits: $numuniq!"
		echo "  Will use $numuniq for % reads mapped calculation"
	fi
	echo -n "Soap -v 1: % reads mapped: "
	perl -e "print $numuniq * 100.0 / $TOT_READS"
	echo
else
	echo "Didn't find ${NAME}.soap.v1.map"
fi

if [ -f ${NAME}.soap.v2.map ] ; then
	if [ ! -f ${NAME}.soap.v2.reads.mapped ] ; then
		awk '{print $1}' ${NAME}.soap.v2.map > ${NAME}.soap.v2.reads.mapped
	fi
	if [ ! -f ${NAME}.soap.v2.reads.mapped.uniq ] ; then
		sort -u ${NAME}.soap.v2.reads.mapped > ${NAME}.soap.v2.reads.mapped.uniq
	fi
	num=`wc -l ${NAME}.soap.v2.reads.mapped | cut -d" " -f1`
	numuniq=`wc -l ${NAME}.soap.v2.reads.mapped.uniq | cut -d" " -f1`
	if [ $num -ne $numuniq ] ; then
		echo "Bowtie: Num hits: $num, num unique hits: $numuniq!"
		echo "  Will use $numuniq for % reads mapped calculation"
	fi
	echo -n "Soap -v 2: % reads mapped: "
	perl -e "print $numuniq * 100.0 / $TOT_READS"
	echo
else
	echo "Didn't find ${NAME}.soap.v2.map"
fi

