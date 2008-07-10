#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
DO_SMALL_FOOTPRINT=0

./ebwt_build --version
./ebwt_build_packed --version

# ./ebwt_build -d --entireSA
if [ ! -f ${NAME}.esa.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.esa.top
	sh wrap.sh ${NAME}.ebwt_build.esa \
	    ./ebwt_build -d --entireSA -v hs_ref_${NAME}.mfa ${NAME}.esa \
	    2>&1 | tee ${NAME}.ebwt_build.esa.out
fi
# ./ebwt_build -d --dcv 256 --bmax 3000000000
if [ ! -f ${NAME}.blf.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.blf.top
	sh wrap.sh ${NAME}.ebwt_build.blf \
	    ./ebwt_build -d --dcv 256 --bmax 3000000000 -v hs_ref_${NAME}.mfa ${NAME}.blf \
	    2>&1 | tee ${NAME}.ebwt_build.blf.out
fi
# ./ebwt_build --d --bmaxDivN 7
if [ ! -f ${NAME}.bl7.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.bl7.top
	sh wrap.sh ${NAME}.ebwt_build.bl7 \
	    ./ebwt_build -d --bmaxDivN 7 -v hs_ref_${NAME}.mfa ${NAME}.bl7 \
	    2>&1 | tee ${NAME}.ebwt_build.bl7.out
fi
# ./ebwt_build_packed -d --dcv 256 --bmax 3000000000
if [ ! -f ${NAME}.pkf.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.pkf.top
	sh wrap.sh ${NAME}.ebwt_build.pkf \
	    ./ebwt_build_packed -d --dcv 256 --bmax 3000000000 -v hs_ref_${NAME}.mfa ${NAME}.pkf \
	    2>&1 | tee ${NAME}.ebwt_build.pkf.out
fi
# ./ebwt_build_packed -d --bmaxDivN 25 --dcv 4096
if [ "$DO_SMALL_FOOTPRINT" = "1" ] ; then
	if [ ! -f ${NAME}.pkt.1.ebwt ] ; then
		echo > ${NAME}.ebwt_build.pkt.top
		sh wrap.sh ${NAME}.ebwt_build.pkt \
		    ./ebwt_build_packed -d --bmaxDivN 25 --dcv 4096 -v hs_ref_${NAME}.mfa ${NAME}.pkt \
		    2>&1 | tee ${NAME}.ebwt_build.pkt.out
	fi
fi
