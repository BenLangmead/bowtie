#!/bin/sh

# Goal:
#
# WORKSTATION=1:
#                 Running Time        % Reads Mapped
#  Bowtie              YYm:ZZs                 XX.X%
#  Maq -n 1         Xh:YYm:ZZs                 XX.X%
#  Maq             XXh:YYm:ZZs                 XX.X%
#
# WORKSTATION=0:
#                 Running Time        % Reads Mapped
#  Bowtie              YYm:ZZs                 XX.X%
#  Maq -n 1         Xh:YYm:ZZs                 XX.X%
#  Maq             XXh:YYm:ZZs                 XX.X%
#  SOAP -v 1       XXh:YYm:ZZs                 XX.X%
#  SOAP            XXh:YYm:ZZs                 XX.X%

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
WORKSTATION=1
if [ ! -f $NAME.results.txt ] ; then
	perl summarize_top.pl $NAME.ebwt.n1.top ebwt_search | tail -1  > $NAME.results.txt
	perl summarize_top.pl $NAME.ebwt.top    ebwt_search | tail -1 >> $NAME.results.txt
	perl summarize_top.pl $NAME.maq.n1.top  maq         | tail -1 >> $NAME.results.txt
	perl summarize_top.pl $NAME.maq.top     maq         | tail -1 >> $NAME.results.txt
	perl summarize_top.pl $NAME.soap.v1.top soap        | tail -1 >> $NAME.results.txt
	perl summarize_top.pl $NAME.soap.v2.top soap        | tail -1 >> $NAME.results.txt
fi
if [ ! -f $NAME.maps.txt ] ; then
	sh analyze_maps.sh > $NAME.maps.txt
fi

# For a workstation:

# Bowtie -n 1
echo -n "Bowtie," > $NAME.results.bowtie.n1.txt
tmp=`head -1 $NAME.results.txt | cut -d" " -f 2 | cut -d, -f 1`
echo -n "$tmp," >> $NAME.results.bowtie.n1.txt 
tmp=`head -1 $NAME.results.txt | cut -d" " -f 2 | cut -d, -f 2`
echo -n "$tmp," >> $NAME.results.bowtie.n1.txt
head -2 $NAME.maps.txt | tail -1 | cut -d":" -f 3 >> $NAME.results.bowtie.n1.txt

# Bowtie
echo -n "Bowtie," > $NAME.results.bowtie.txt
tmp=`head -2 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1`
echo -n "$tmp," >> $NAME.results.bowtie.txt 
tmp=`head -2 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 2`
echo -n "$tmp," >> $NAME.results.bowtie.txt
head -3 $NAME.maps.txt | tail -1 | cut -d":" -f 3 >> $NAME.results.bowtie.txt

# Maq -n1
echo -n "Maq -n1," > $NAME.results.maq.n1.txt
tmp=`head -3 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1`
echo -n "$tmp," >> $NAME.results.maq.n1.txt
tmp=`head -3 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 2`
echo -n "$tmp," >> $NAME.results.maq.n1.txt
head -4 $NAME.maps.txt | tail -1 | cut -d":" -f 3 >> $NAME.results.maq.n1.txt

# Maq
echo -n "Maq," > $NAME.results.maq.txt
tmp=`head -4 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1`
echo -n "$tmp," >> $NAME.results.maq.txt
tmp=`head -4 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 2`
echo -n "$tmp," >> $NAME.results.maq.txt
head -5 $NAME.maps.txt | tail -1 | cut -d":" -f 3 >> $NAME.results.maq.txt

if [ "$WORKSTATION" = "0" ] ; then
	# SOAP -v 1
	echo -n "SOAP -v 1," > $NAME.results.soap.v1.txt
	tmp=`head -5 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1`
	echo -n "$tmp," >> $NAME.results.soap.v1.txt
	tmp=`head -5 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 2`
	echo -n "$tmp," >> $NAME.results.soap.v1.txt
	head -6 $NAME.maps.txt | tail -1 | cut -d":" -f 3 >> $NAME.results.soap.v1.txt

	# SOAP -v
	echo -n "SOAP," > $NAME.results.soap.txt
	tmp=`head -6 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1`
	echo -n "$tmp," >> $NAME.results.soap.txt
	tmp=`head -6 $NAME.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 2`
	echo -n "$tmp," >> $NAME.results.soap.txt
	head -7 $NAME.maps.txt | tail -1 | cut -d":" -f 3 >> $NAME.results.soap.txt
	perl plot.pl bowtie maq.n1 maq soap.v1 soap
else 
	perl plot.pl bowtie maq.n1 maq
fi
