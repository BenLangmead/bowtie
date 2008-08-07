#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
#echo Using NAME: ${NAME}
export TOT_READS=8839010
export TOT_FILT_READS=8400865
WORKSTATION=1
if [ `hostname` = "privet.umiacs.umd.edu" ] ; then
	WORKSTATION=0
fi
if [ `hostname` = "larch.umiacs.umd.edu" ] ; then
	WORKSTATION=0
fi
#echo "Don't forget to set WORKSTATION=? appropriately in plot.sh"
#echo -n "Currently set to: "
#if [ "$WORKSTATION" = "1" ] ; then
#	echo Workstation
#else
#	echo Server
#fi

do_one()
{
	RNAME=$1
	MAP_BASE=$2
	MAP_EXT=$3
	NREADS=$4
	ISMAQ=$5
	MAPFILE="${MAP_BASE}.${MAP_EXT}"
	RMAP="${MAP_BASE}.reads.mapped"
	RMAPU="${MAP_BASE}.reads.mapped.uniq"
	if [ -f "${MAPFILE}" ] ; then
		if [ ! -f ${MAP_BASE}.reads.mapped ] ; then
			if [ "${ISMAQ}" = "1" ] ; then
				maq mapview ${MAPFILE} | awk '{print $1}' > ${RMAP}
			else
				awk '{print $1}' ${MAPFILE} > ${RMAP}
			fi
		fi
		if [ ! -f ${RMAPU} ] ; then
			sort -u ${RMAP} > ${RMAPU}
		fi
		num=`wc -l ${RMAP} | cut -d" " -f1`
		numuniq=`wc -l ${RMAPU} | cut -d" " -f1`
		if [ $num -ne $numuniq ] ; then
			echo "${RNAME}: Num hits: $num, num unique hits: $numuniq!"
			echo "  Will use $numuniq for % reads mapped calculation"
		fi
		echo -n "${RNAME}: % reads mapped: "
		perl -e "print $numuniq * 100.0 / ${NREADS}"
		echo
	else
		echo "Didn't find ${MAPFILE}"
	fi
}

# Bowtie
if [ 0 -gt 1 ] ; then
	do_one "Bowtie -n 1" "${NAME}.ebwt.n1" "hits" "${TOT_READS}" "0"
fi
do_one "Bowtie" "${NAME}.ebwt" "hits" "${TOT_READS}" "0"
do_one "Bowtie filtered" "${NAME}.ebwt.filt" "hits" "${TOT_FILT_READS}" "0"
if [ "$WORKSTATION" = "0" ] ; then
	if [ 0 -gt 1 ] ; then
		do_one "Bowtie -v 1" "${NAME}.ebwt.1" "hits" "${TOT_READS}" "0"
	fi
	do_one "Bowtie -v 2" "${NAME}.ebwt.2" "hits" "${TOT_READS}" "0"
fi

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
if [ "$WORKSTATION" = "0" ] ; then
	if [ 0 -gt 1 ] ; then
		do_one "Soap -v 1" "${NAME}.soap.v1" "map" "${TOT_READS}" "0"
	fi
	do_one "Soap" "${NAME}.soap.v2" "map" "${TOT_READS}" "0"
fi
