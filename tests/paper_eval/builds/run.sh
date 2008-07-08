#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}

./ebwt_build --version

sh wrap.sh ${NAME}.esa ./ebwt_build -d --entireSA -v hs_ref_${NAME}.mfa ${NAME}.esa
sh wrap.sh ${NAME}.blf ./ebwt_build -d --bmax 3000000000 -v hs_ref_${NAME}.mfa ${NAME}.blf
sh wrap.sh ${NAME}.bl7 ./ebwt_build -d --bmaxDivN 7 -v hs_ref_${NAME}.mfa ${NAME}.bl7
