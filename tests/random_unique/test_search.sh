#!/bin/bash

CSAMapper=../../
BUILD=$CSAMapper/ebwt_build
#PYTHON=/fs/sz-user-supported/Linux-i686/bin/python2.5
PYTHON=/usr/bin/python
ECHO="echo -e"

# # Generate a reference that's 50000bp long
# $PYTHON rand_reads.py 1 50000 > ref.fna 2>/dev/null
# # 
# # # Make an ebwt out of it
# $BUILD --entireSA ref.fna ref
# # 
# # 
# # Clear previous meryl data, if any
# if [ -e ref22.mcdat ]
#     then
#     rm ref22.*
# fi
# 
# # Identify unique 22-mers
# echo "Indentifying unique 22-mers"
# meryl -B -f -s ref.fna -U 1 -m 22 -o ref22
# meryl -Dt -s ref22 | awk 'BEGIN{i = 0}{ if ($1 == ">1"){t++; printf ">%d\n",t; }else print $1}' > ref22-unique-mers.fna
# 
# mummer -maxmatch ref.fna ref22-unique-mers.fna > ref22-unique-mer.out 2>/dev/null

# # Make exact reads from unique 22-mers
# $PYTHON reads_from_mers.py -q ref.fna ref22-unique-mer.out \
#     exact_reads.fq exact.ebwt_out 2>/dev/null
# 
# if [ -e exact.actual_out ] 
#     then
#     rm exact.actual_out
# fi
# 
# #Run ebwt search on the reads...
# $CSAMapper/ebwt_search -q -k 22 ref exact_reads.fq > exact.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w exact.ebwt_out exact.actual_out > /dev/null; then
#     $ECHO "Exact forward matching:\t\tOK"
# else
#     $ECHO "Exact forward matching:\t\tFAILED"
# fi
# 
# # Make exact (reverse complemented) reads from unique 22-mers
# $PYTHON reads_from_mers.py -r ref.fna ref22-unique-mer.out \
#     exact_reads_rc.fna exact_rc.ebwt_out 2>/dev/null
# 
# if [ -e exact_rc.actual_out ] 
#     then
#     rm exact_rc.actual_out
# fi
# 
# # Run ebwt search on the reads...
# $CSAMapper/ebwt_search -k 22 -r ref exact_reads_rc.fna | \
#     awk '{if ($1 ~ /-/)  print }' > exact_rc.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w exact_rc.ebwt_out exact_rc.actual_out > /dev/null; then
#     $ECHO "Exact reverse matching:\t\tOK"
# else
#     $ECHO "Exact reverse matching:\t\tFAILED"
# fi
# 
# # Make reads from unique 22-mers that have mismatches at the 3' end
# $PYTHON reads_from_mers.py -q 3 ref.fna ref22-unique-mer.out \
#     3mis_reads.fna 3mis_reads.ebwt_out 2>/dev/null
# 
# if [ -e 3mis_reads.actual_out ] 
#     then
#     rm 3mis_reads.actual_out
# fi
# 
# # Run ebwt search on the mismatching reads...
# $CSAMapper/ebwt_search -k 22 ref 3mis_reads.fna > 3mis_reads.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w 3mis_reads.ebwt_out 3mis_reads.actual_out > /dev/null; then
#     $ECHO "3' mismatch forward matching:\tOK"
# else
#     $ECHO "3' mismatch forward matching:\tFAILED"
# fi
# 
# #Make reverse complemented reads from unique 22-mers that 
# #have mismatches at the 3' end
# $PYTHON reads_from_mers.py -r -q 3 ref.fna ref22-unique-mer.out \
#     3mis_reads_rc.fna 3mis_reads_rc.ebwt_out 2>/dev/null
# 
# if [ -e 3mis_reads_rc.actual_out ] 
#     then
#     rm 3mis_reads_rc.actual_out
# fi
# 
# # Run ebwt search on the mismatching reads...
# $CSAMapper/ebwt_search -r -k 22 ref 3mis_reads_rc.fna  | awk '{if ($1 ~ /-/)  print }' > 3mis_reads_rc.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w 3mis_reads_rc.ebwt_out 3mis_reads_rc.actual_out > /dev/null; then
#     $ECHO "3' mismatch reverse matching:\tOK"
# else
#     $ECHO "3' mismatch reverse matching:\tFAILED"
# fi

# Make reads from unique 22-mers that have 1 mismatch in the 5' end
$PYTHON reads_from_mers.py  -m -1 ref.fna ref22-unique-mer.out \
    1_mis_reads.fq 1_mis_reads.ebwt_out 2>/dev/null

if [ -e 1_mis_reads.actual_out ] 
    then
    rm 1_mis_reads.actual_out
fi

# Run ebwt search on the mismatching reads...
$CSAMapper/ebwt_search -1 -q ref 1_mis_reads.fq | sort -n -k4 > 1_mis_reads.actual_out
sort -n -k4 1_mis_reads.ebwt_out > xxx
mv xxx 1_mis_reads.ebwt_out

# and then check that they all map to where they are supposed to
if diff -w 1_mis_reads.ebwt_out 1_mis_reads.actual_out > /dev/null; then
    $ECHO "1-mismatch in 5' mer forward matching:\tOK"
else
    $ECHO "1-mismatch in 5' mer forward matching:\tFAILED"
fi

# Make reads from unique 22-mers that have 1 mismatch in the 22-mer and 
# several mismatches at the 3' end
$PYTHON reads_from_mers.py -q 3 -m -1 ref.fna ref22-unique-mer.out \
    1_mis_with_3mis_reads.fq 1_mis_with_3mis_reads.ebwt_out 2>/dev/null

if [ -e 1_mis_with_3mis_reads.actual_out ] 
    then
    rm 1_mis_with_3mis_reads.actual_out
fi

# Run ebwt search on the mismatching reads...
$CSAMapper/ebwt_search -1 -k 22 -q ref 1_mis_with_3mis_reads.fq | sort -n -k4 > 1_mis_with_3mis_reads.actual_out
sort -n -k4 1_mis_with_3mis_reads.ebwt_out > xxx
mv xxx 1_mis_with_3mis_reads.ebwt_out

# and then check that they all map to where they are supposed to
if diff -w 1_mis_with_3mis_reads.ebwt_out 1_mis_with_3mis_reads.actual_out > /dev/null; then
    $ECHO "1-mismatch in 5' mer + 3' mismatch forward matching:\tOK"
else
    $ECHO "1-mismatch in 5' mer + 3' mismatch forward matching:\tFAILED"
fi

