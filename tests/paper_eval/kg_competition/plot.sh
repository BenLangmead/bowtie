#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
sh summarize_all_top.sh > $NAME.results.txt
sh analyze_maps.sh > $NAME.maps.txt

# Bowtie
echo -n "Bowtie," > $NAME.results.bowtie.txt
head -1 $NAME.results.txt | cut -d" " -f 2 | cut -d, -f 1 >> $NAME.results.bowtie.txt
echo -n "," >> $NAME.results.bowtie.txt
head -2 $NAME.maps.txt | tail -1 | cut -d: -f 2 >> $NAME.results.bowtie.txt

# Maq -n1
echo -n "Maq -n1, " > $NAME.results.maq.n1.txt
head -2 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1 >> $NAME.results.maq.n1.txt
echo -n " " >> $NAME.results.maq.n1.txt
head -4 $NAME.maps.txt | tail -1 | cut -d: -f 2 >> $NAME.results.maq.n1.txt

# Maq
echo -n "Maq, " > $NAME.results.maq.txt
head -3 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1 >> $NAME.results.maq.txt
echo -n " " >> $NAME.results.maq.txt
head -3 $NAME.maps.txt | tail -1 | cut -d: -f 2 >> $NAME.results.maq.txt

#perl plot.pl whole
