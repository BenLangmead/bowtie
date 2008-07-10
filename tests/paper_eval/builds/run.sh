#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}

./ebwt_build --version

if [ ! -f ${NAME}.esa.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.esa.top
	sh wrap.sh ${NAME}.ebwt_build.esa \
	    ./ebwt_build -d --entireSA -v hs_ref_${NAME}.mfa ${NAME}.esa \
	    2>&1 | tee ${NAME}.ebwt_build.esa.out
fi
if [ ! -f ${NAME}.blf.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.blf.top
	sh wrap.sh ${NAME}.ebwt_build.blf \
	    ./ebwt_build -d --bmax 3000000000 -v hs_ref_${NAME}.mfa ${NAME}.blf \
	    2>&1 | tee ${NAME}.ebwt_build.blf.out
fi
if [ ! -f ${NAME}.bl7.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.bl7.top
	sh wrap.sh ${NAME}.ebwt_build.bl7 \
	    ./ebwt_build -d --bmaxDivN 7 -v hs_ref_${NAME}.mfa ${NAME}.bl7 \
	    2>&1 | tee ${NAME}.ebwt_build.bl7.out
fi
if [ ! -f ${NAME}.pkf.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.pkf.top
	sh wrap.sh ${NAME}.ebwt_build.pkf \
	    ./ebwt_build_packed -d --bmax 3000000000 -v hs_ref_${NAME}.mfa ${NAME}.pkf \
	    2>&1 | tee ${NAME}.ebwt_build.pkf.out
fi
