#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
#echo Using NAME: ${NAME}
export TOT_READS=8839010
export TOT_FILT_READS=8400865

do_one()
{
	NAME=$1
	MAP_BASE=$2
	MAP_EXT=$3
	NREADS=$4
	ISMAQ=$5
	MAP="${MAP_BASE}.${MAP_EXT}"
	RMAP="${MAP_BASE}.reads.mapped"
	RMAPU="${MAP_BASE}.reads.mapped.uniq"
	if [ -f ${MAP} ] ; then
		if [ ! -f ${MAP_BASE}.reads.mapped ] ; then
			if [ "$IS_MAQ" = "1" ] ; then
				awk '{print $1}' ${MAP} > ${RMAP}
			else
				maq mapview ${MAP} | awk '{print $1}' > ${RMAP}
			fi
		fi
		if [ ! -f ${RMAPU} ] ; then
			sort -u ${RMAP} > ${RMAPU}
		fi
		num=`wc -l ${RMAP} | cut -d" " -f1`
		numuniq=`wc -l ${RMAPU} | cut -d" " -f1`
		if [ $num -ne $numuniq ] ; then
			echo "${NAME}: Num hits: $num, num unique hits: $numuniq!"
			echo "  Will use $numuniq for % reads mapped calculation"
		fi
		echo -n "${NAME}: % reads mapped: "
		perl -e "print $numuniq * 100.0 / ${NREADS}"
		echo
	else
		echo "Didn't find ${MAP}"
	fi
}

# Bowtie
if [ 0 -gt 1 ] ; then
	do_one "Bowtie -n 1" "${NAME}.ebwt.n1" "hits" "${TOT_READS}" "0"
fi
do_one "Bowtie" "${NAME}.ebwt" "hits" "${TOT_READS}" "0"
do_one "Bowtie filtered" "${NAME}.filt.ebwt" "hits" "${TOT_FILT_READS}" "0"
do_one "Bowtie -v 2" "${NAME}.ebwt.2" "hits" "${TOT_READS}" "0"

# Maq
if [ 0 -gt 1 ] ; then
	do_one "Maq -n 1" "${NAME}.maq.n1" "map" "${TOT_READS}" "1"
fi
do_one "Maq" "${NAME}.maq" "map" "${TOT_READS}" "1"
if [ 0 -gt 1 ] ; then
	do_one "Maq -n 1 filtered" "${NAME}.maq.n1.filt" "map" "${TOT_FILT_READS}" "1"
fi
do_one "Maq filtered" "${NAME}.maq.filt" "map" "${TOT_FILT_READS}" "1"

# Soap
if [ 0 -gt 1 ] ; then
	do_one "Soap -v 1" "${NAME}.soap.v1" "map" "${TOT_READS}" "0"
fi
do_one "Soap" "${NAME}.soap.v2" "map" "${TOT_READS}" "0"
