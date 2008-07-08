#!/bin/sh
str=$1
shift
top -b -d 5 >> $str.top &
pid=$!
$*
kill $pid
