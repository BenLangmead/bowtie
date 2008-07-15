#!/bin/sh

# 6 output files: 3 trim lengths x 2 mismatch levels

for t in 0 5 10 ; do
	# Truncate output files
	echo > trim$t_0mm.dat
	echo > trim$t_1mm.dat
	for s in whole.179451189 \
	         whole.358902379 \
	         whole.717804758 \
	         whole.1435609516 \
	         whole ; do
		perl summarize_top.pl $s.sim_reads.t$t.top ebwt_search \
		  | tail -1 \
		  | cut -d" " -f 2 \
		  | cut -d"," -f 1 >> trim$t_0mm.dat
		perl summarize_top.pl $s.sim_reads.1.t$t.top ebwt_search \
		  | tail -1 \
		  | cut -d" " -f 2 \
		  | cut -d"," -f 1 >> trim$t_1mm.dat
	done
done

# Run gplot, which will generate GNUplot output plot.gnuplot and
# render it to postscript in /tmp/gplot.eps.  Usually, you'll want to
# tinker with plot.gnuplot and run gnuplot manually to get the best
# output.
perl ~/software/gplot-1.3/gplot.pl \
  -plotcmds plot.gpl \
  -xlabel "Reference genome length (bases)" \
  -ylabel "Running time (CPU seconds)" \
  -type eps \
  -showgnuplot \
  -title "Running time to map 2M simulated reads onto various genome lengths" \
  -thickness 2 -color black -point box          -pointsize 1.3 \
     -name "35bp reads, exact matching"      trim0_0mm.dat \
  -thickness 2 -color black -point circle       -pointsize 1.3 \
     -name "30bp reads, exact matching"      trim5_0mm.dat \
  -thickness 2 -color black -point downtriangle -pointsize 1.3 \
     -name "25bp reads, exact matching"      trim10_0mm.dat \
  -thickness 2 -color black -point box          -pointsize 1.3 \
     -name "35bp reads, 1-mismatch matching" trim0_1mm.dat \
  -thickness 2 -color black -point circle       -pointsize 1.3 \
     -name "30bp reads, 1-mismatch matching" trim5_1mm.dat \
  -thickness 2 -color black -point downtriangle -pointsize 1.3 \
     -name "25bp reads, 1-mismatch matching" trim10_1mm.dat \
  | grep -v "GUPLOT" > plot.gnuplot

if [ ! -f /tmp/gplot.eps ] ; then
	echo "Error: didn't find gnuplot output in /tmp/gplot.eps"
	exit 1
fi

cp /tmp/gplot.eps .
