#!/bin/sh

#
# Usage:  soap [options]
#        -a  <str>   query a file, *.fq or *.fa format
#        -d  <str>   reference sequences file, *.fa format
#        -o  <str>   output alignment file
#        -s  <int>   seed size, default=10. [read>18,s=8; read>22,s=10, read>26, s=12]
#        -v  <int>   maximum number of mismatches allowed on a read, <=5. default=2bp. For pair-ended alignment, this version will allow either 0 or 2 mismatches.
#        -g  <int>   maximum gap size allowed on a read, default=0bp
#        -w  <int>   maximum number of equal best hits to count, smaller will be faster, <=10000
#        -e  <int>   will not allow gap exist inside n-bp edge of a read, default=5bp
#        -z  <char>  initial quality, default=@ [Illumina is using '@', Sanger Institute is using '!']
#        -c  <int>   how to trim low-quality at 3-end?
#                    0:     don't trim;
#                    1-10:  trim n-bps at 3-end for all reads;
#                    11-20: trim first bp and (n-10)-bp at 3-end for all reads;
#                    21-30: trim (n-20)-bp at 3-end and redo alignment if the original read have no hit;
#                    31-40: trim first bp and (n-30)-bp at 3-end and redo alignment if the original read have no hit;
#                    41-50: iteratively trim (n-40)-bp at 3-end until getting hits;
#                    51-60: if no hit, trim first bp and iteratively trim (n-50)bp at 3-end until getting hits;
#                    default: 0
#        -f  <int>   filter low-quality reads containing >n Ns, default=5
#        -r  [0,1,2] how to report repeat hits, 0=none; 1=random one; 2=all, default=1
#        -t          read ID in output file, [name, order in input file], default: name
#        -n  <int>   do alignment on which reference chain? 0:both; 1:forward only; 2:reverse only. default=0
#        -p  <int>   number of processors to use, default=1

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
SOAP=`which soap.contig`
echo Using soap.contig: ${SOAP}
REF=hs_ref_${NAME}.mfa
READ_BASE=kg_reads

# Do 2-mismatch
if [ ! -f ${NAME}.soap.v2.map ] ; then
    echo Doing 2-mismatch...
    echo > ${NAME}.soap.v2.top
    sh wrap.sh ${NAME}.soap.v2 \
	$SOAP -v 2 -w 1 -f 35 -r 1 -n 0 \
	  -o ${NAME}.soap.v2.map \
	  -d ${REF} \
	  -a ${READ_BASE}.fq
	if [ ! -f ${NAME}.soap.v2.map ] ; then
		echo "soap (v2) did not produce legitimate .map file; aborting..."
		exit 1
	fi
fi

# Do 1-mismatch
if [ ! -f ${NAME}.soap.v1.map ] ; then
    echo Doing 1-mismatch...
    echo > ${NAME}.soap.v1.top
    sh wrap.sh ${NAME}.soap.v1 \
	$SOAP -v 1 -w 1 -f 35 -r 1 -n 0 \
	  -o ${NAME}.soap.v1.map \
	  -d ${REF} \
	  -a ${READ_BASE}.fq
	if [ ! -f ${NAME}.soap.v1.map ] ; then
		echo "soap (v1) did not produce legitimate .map file; aborting..."
		exit 1
	fi
fi
