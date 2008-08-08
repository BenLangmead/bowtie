#!/bin/sh

dir=`pwd`
NAME=`basename $dir | sed 's/_.*//'`
echo Using NAME: ${NAME}
WORKSTATION=1
if [ `hostname` = "privet.umiacs.umd.edu" ] ; then
	WORKSTATION=0
fi
if [ `hostname` = "larch.umiacs.umd.edu" ] ; then
	WORKSTATION=0
fi
echo "Don't forget to set WORKSTATION=? appropriately in plot.sh"
echo -n "Currently set to: "
if [ "$WORKSTATION" = "1" ] ; then
	echo Workstation
else
	echo Server
fi

do_one()
{
	RNAME=$1
	LINE=$2
	RESULTS=$3
	echo -n "${RNAME}," > ${NAME}.results.${RESULTS}.txt
	tmp=`head -${LINE} ${NAME}.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 1`
	echo -n "$tmp," >> ${NAME}.results.${RESULTS}.txt 
	tmp=`head -${LINE} ${NAME}.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 2`
	echo -n "$tmp," >> ${NAME}.results.${RESULTS}.txt
	tmp=`head -${LINE} ${NAME}.results.txt | tail -1 | cut -d" " -f 2 | cut -d, -f 3`
	echo -n "$tmp," >> ${NAME}.results.${RESULTS}.txt
	head -${LINE} $NAME.maps.txt | tail -1 | cut -d":" -f 3 >> $NAME.results.${RESULTS}.txt
}

if [ ! -f $NAME.results.txt ] ; then
	if [ 0 -gt 1 ] ; then
		perl summarize_top.pl $NAME.ebwt.top         ebwt_search | tail -1 >> $NAME.results.txt
		perl summarize_top.pl $NAME.ebwt.filt.top    ebwt_search | tail -1 >> $NAME.results.txt
		perl summarize_top.pl $NAME.ebwt.n1.top      ebwt_search | tail -1  > $NAME.results.txt
		perl summarize_top.pl $NAME.ebwt.n1.filt.top ebwt_search | tail -1  > $NAME.results.txt
		if [ "$WORKSTATION" = "0" ] ; then
			perl summarize_top.pl $NAME.ebwt.1.top       ebwt_search | tail -1 >> $NAME.results.txt
			perl summarize_top.pl $NAME.ebwt.2.top       ebwt_search | tail -1 >> $NAME.results.txt
		fi
		perl summarize_top.pl $NAME.maq.n1.top       maq         | tail -1 >> $NAME.results.txt
		perl summarize_top.pl $NAME.maq.top          maq         | tail -1 >> $NAME.results.txt
		perl summarize_top.pl $NAME.maq.n1.filt.top  maq         | tail -1 >> $NAME.results.txt
		perl summarize_top.pl $NAME.maq.filt.top     maq         | tail -1 >> $NAME.results.txt
		if [ "$WORKSTATION" = "0" ] ; then
			# Up the threshold on the SOAP runs; otherwise our figures
			# include multiple SOAP processes running in parallel, which
			# screws up the CPU Time/Wall clock time figures
			perl summarize_top.pl -t 900 $NAME.soap.v1.top soap | tail -1 >> $NAME.results.txt
			perl summarize_top.pl -t 900 $NAME.soap.v2.top soap | tail -1 >> $NAME.results.txt
		fi
	else
		perl summarize_top.pl $NAME.ebwt.top         ebwt_search | tail -1 >> $NAME.results.txt
		perl summarize_top.pl $NAME.ebwt.filt.top    ebwt_search | tail -1 >> $NAME.results.txt
		if [ "$WORKSTATION" = "0" ] ; then
			perl summarize_top.pl $NAME.ebwt.2.top       ebwt_search | tail -1 >> $NAME.results.txt
		fi
		perl summarize_top.pl $NAME.maq.top          maq         | tail -1 >> $NAME.results.txt
		perl summarize_top.pl $NAME.maq.filt.top     maq         | tail -1 >> $NAME.results.txt
		if [ "$WORKSTATION" = "0" ] ; then
			# Up the threshold on the SOAP runs; otherwise our figures
			# include multiple SOAP processes running in parallel, which
			# screws up the CPU Time/Wall clock time figures
			perl summarize_top.pl -t 900 $NAME.soap.v2.top soap | tail -1 >> $NAME.results.txt
		fi
	fi
fi

if [ ! -f $NAME.maps.txt ] ; then
	sh analyze_maps.sh > $NAME.maps.txt
fi

if [ 0 -gt 1 ] ; then
	# Bowtie
	do_one "Bowtie" "1" "bowtie"
	do_one "Bowtie filtered" "2" "bowtie.filt"
	do_one "Bowtie -n1" "3" "bowtie.n1" 
	do_one "Bowtie -n1 filtered" "4" "bowtie.n1.filt" 
	if [ "$WORKSTATION" = "0" ] ; then
		do_one "Bowtie -v 1" "5" "bowtie.v1"
		do_one "Bowtie -v 2" "6" "bowtie.v2"
		# Maq
		do_one "Maq -n 1" "7" "maq.n1"
		do_one "Maq -n 1 filtered" "8" "maq.n1.filt"
		do_one "Maq" "9" "maq"
		do_one "Maq filtered" "10" "maq.filt"
		# Soap
		do_one "Soap -v 1" "11" "soap.v1"
		do_one "Soap" "12" "soap.v2"
		perl plot.pl bowtie.n1 maq.n1 \
		             bowtie maq \
		             bowtie.n1.filt maq.n1.filt \
		             bowtie.filt maq.filt \
		             bowtie.v1 soap.v1 \
		             bowtie.v2 soap.v2
	else
		# Maq
		do_one "Maq -n 1" "5" "maq.n1"
		do_one "Maq -n 1 filtered" "6" "maq.n1.filt"
		do_one "Maq" "7" "maq"
		do_one "Maq filtered" "8" "maq.filt"
		# Soap
		perl plot.pl bowtie.n1 maq.n1 \
		             bowtie maq \
		             bowtie.n1.filt maq.n1.filt \
		             bowtie.filt maq.filt \
		             - - \
		             - -
	fi

else
	# No 1-mm

	# Bowtie
	do_one "Bowtie" "1" "bowtie"
	do_one "Bowtie filtered" "2" "bowtie.filt"
	if [ "$WORKSTATION" = "0" ] ; then
		do_one "Bowtie -v 2" "3" "bowtie.v2"
		# Maq
		do_one "Maq" "4" "maq"
		do_one "Maq filtered" "5" "maq.filt"
		# Soap
		do_one "Soap" "6" "soap.v2"
		perl plot.pl - - \
		             bowtie maq \
		             - - \
		             bowtie.filt maq.filt \
		             - - \
		             bowtie.v2 soap.v2
	else
		# Maq
		do_one "Maq" "3" "maq"
		do_one "Maq filtered" "4" "maq.filt"
		perl plot.pl - - \
		             bowtie maq \
		             - - \
		             bowtie.filt maq.filt \
		             - - \
		             - -
	fi
fi

echo Done
