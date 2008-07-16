#!/bin/sh

# Goal:
#  Budget  Actual Max Virtual    Running Time
#    16GB              XX.XGB         YYm:ZZs
#     8GB               X.XGB         YYm:ZZs
#     4GB               X.XGB         YYm:ZZs
#     2GB               X.XGB         YYm:ZZs

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
if [ ! -f ${NAME}.results.txt ] ; then
	perl summarize_top.pl ${NAME}.ebwt_build.blf.top ebwt_build | tail -1 \
		| sed 's/.*: //' | sed 's/,[0-9]*$//' > ${NAME}.results.txt
	perl summarize_top.pl ${NAME}.ebwt_build.bl7.top ebwt_build | tail -1 \
		| sed 's/.*: //' | sed 's/,[0-9]*$//' >> ${NAME}.results.txt
	perl summarize_top.pl ${NAME}.ebwt_build.pkl.top ebwt_build | tail -1 \
		| sed 's/.*: //' | sed 's/,[0-9]*$//' >> ${NAME}.results.txt
	perl summarize_top.pl ${NAME}.ebwt_build.pkt.top ebwt_build | tail -1 \
		| sed 's/.*: //' | sed 's/,[0-9]*$//' >> ${NAME}.results.txt
fi

perl plot.pl
