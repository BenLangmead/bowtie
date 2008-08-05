#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
DO_SMALL_FOOTPRINT=1

./ebwt_build --version
./ebwt_build_packed --version

# Probably the fastest configuration that fits in a 16GB memory budget
# ./ebwt_build -d --dcv 256 --bmax 3000000000
if [ ! -f ${NAME}.blf.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.blf.top
	sh wrap.sh ${NAME}.ebwt_build.blf \
	    ./ebwt_build -d --dcv 256 --bmax 3000000000 -v hs_ref_${NAME}.mfa ${NAME}.blf \
	    2>&1 | tee ${NAME}.ebwt_build.blf.out
fi
# Probably the fastest configuration that fits in a 8GB memory budget
# ./ebwt_build --d --bmaxdivn 7
if [ ! -f ${NAME}.bl4.1.ebwt ] ; then
	echo > ${NAME}.ebwt_build.bl4.top
	sh wrap.sh ${NAME}.ebwt_build.bl4 \
	    ./ebwt_build -d --bmaxdivn 4 -v hs_ref_${NAME}.mfa ${NAME}.bl4 \
	    2>&1 | tee ${NAME}.ebwt_build.bl4.out
fi
# ./ebwt_build_packed -d --bmaxdivn 25 --dcv 4096
if [ "$DO_SMALL_FOOTPRINT" = "1" ] ; then
	# Probably the fastest configuration that fits in a 4GB memory budget
	if [ ! -f ${NAME}.pkl.1.ebwt ] ; then
		echo > ${NAME}.ebwt_build.pkl.top
		sh wrap.sh ${NAME}.ebwt_build.pkl \
		    ./ebwt_build -d --bmaxdivn 25 --dcv 4096 -v hs_ref_${NAME}.mfa ${NAME}.pkl \
		    2>&1 | tee ${NAME}.ebwt_build.pkl.out
	fi
	# Probably the fastest configuration that fits in a 2GB memory budget
	if [ ! -f ${NAME}.pkt.1.ebwt ] ; then
		echo > ${NAME}.ebwt_build.pkt.top
		sh wrap.sh ${NAME}.ebwt_build.pkt \
		    ./ebwt_build_packed -d --bmaxdivn 25 --dcv 4096 -v hs_ref_${NAME}.mfa ${NAME}.pkt \
		    2>&1 | tee ${NAME}.ebwt_build.pkt.out
	fi
fi
