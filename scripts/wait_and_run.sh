#!/bin/sh

echo blah > .wait_and_run.1
echo blah > .wait_and_run.2
echo blah > .wait_and_run.3

while true ; do
	mv .wait_and_run.2 .wait_and_run.3
	mv .wait_and_run.1 .wait_and_run.2
	top -b -n 1 | perl $HOME/workspace/bowtie/scripts/top_get_cpu_col.pl -f 40 > .wait_and_run.1
	if [ `wc -c .wait_and_run.1 | awk '{print $1}'` == "0" ] ; then \
	if [ `wc -c .wait_and_run.2 | awk '{print $1}'` == "0" ] ; then \
	if [ `wc -c .wait_and_run.3 | awk '{print $1}'` == "0" ] ; then \
	    break;
	fi 
	fi 
	fi 
	sleep 1
done
echo All clear
