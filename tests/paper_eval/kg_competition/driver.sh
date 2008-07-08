#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo "Using NAME: ${NAME}" > name.txt

sh ebwt.sh 2>&1 | tee ${NAME}.ebwt.out
sh maq.sh 2>&1  | tee ${NAME}.maq.out
sh maq_n1.sh 2>&1  | tee ${NAME}.maq_n1.out
sh soap.sh 2>&1 | tee ${NAME}.soap.out
