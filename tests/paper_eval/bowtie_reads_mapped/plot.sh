#!/bin/sh

# Create plot data and plot it with the assistance of a gplot wrapper
# around gnuplot

perl plot.pl

# Run gplot, which will generate GNUplot output plot.gnuplot and
# render it to postscript in /tmp/gplot.eps.  Usually, you'll want to
# tinker with plot.gnuplot and run gnuplot manually to get the best
# output.
perl ~/software/gplot-1.3/gplot.pl \
  -plotcmds plot.gpl \
  -xlabel "Trimmed read length" \
  -ylabel "% trimmed reads mapped" \
  -type eps \
  -showgnuplot \
  -title "Percentage of trimmed reads mapped for one lane of\n1000-Genomes Project Illumina reads (Accession # SRR001115)" \
  -thickness 2 -color black -point box          -pointsize 1.3 \
     -name "Trimmed reads with one exact hit" exactu.dat \
  -thickness 2 -color black -point circle       -pointsize 1.3 \
     -name "Above, plus reads with many exact hits" exact.dat \
  -thickness 2 -color black -point downtriangle -pointsize 1.3 \
     -name "Above, plus reads with one 1-mismatch hit" inexactu.dat \
  -thickness 2 -color black -point uptriangle   -pointsize 1.3 \
     -name "Above, plus reads with many 1-mismatch hits" inexact.dat \
  -thickness 2 -color black -point diamond      -pointsize 1.3 \
     -name "Above, plus reads mapped with Maq policy with 1 5' 24-bp mismatch" maq1.dat \
  -thickness 2 -color black -point pentagon     -pointsize 1.3 \
     -name "Above, plus reads mapped with Maq policy with 2 5' 24-bp mismatches" maq2.dat \
  | grep -v "GUPLOT" > plot.gnuplot

if [ ! -f /tmp/gplot.eps ] ; then
	echo "Error: didn't find gnuplot output in /tmp/gplot.eps"
	exit 1
fi

cp /tmp/gplot.eps .
