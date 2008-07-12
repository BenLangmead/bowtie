#!/bin/bash

CSAMapper=../../
BUILD=$CSAMapper/ebwt_build
#PYTHON=/fs/sz-user-supported/Linux-i686/bin/python2.5
PYTHON=/usr/bin/python
ECHO="echo -e"
export LC_ALL=C

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

####################################################################

# # Make exact reads from unique 22-mers
# $PYTHON reads_from_mers.py -q ref.fna ref22-unique-mer.out \
#     exact_e2e_reads.fq exact_e2e.ebwt_out 
# 
# if [ -e exact_e2e.actual_out ] 
#     then
#     rm exact_e2e.actual_out
# fi
# 
# #Run *end-to-end* ebwt search on the reads...
# $CSAMapper/ebwt_search -q ref exact_e2e_reads.fq > exact_e2e.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w exact_e2e.ebwt_out exact_e2e.actual_out > /dev/null; then
#     $ECHO "Exact e2e forward matching:\t\tOK"
# else
#     $ECHO "Exact e2e forward matching:\t\tFAILED"
# fi
# 
# ######################################################################
# 
# # Make exact (reverse complemented) reads from unique 22-mers
# $PYTHON reads_from_mers.py -q -r ref.fna ref22-unique-mer.out \
#     exact_e2e_reads_rc.fq exact_e2e_rc.ebwt_out 
# 
# if [ -e exact_e2e_rc.actual_out ] 
#     then
#     rm exact_e2e_rc.actual_out
# fi
# 
# # Run *end-to-end* ebwt search on the reads...
# $CSAMapper/ebwt_search -r ref -q exact_e2e_reads_rc.fq > exact_e2e_rc.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w exact_e2e_rc.ebwt_out exact_e2e_rc.actual_out > /dev/null; then
#     $ECHO "Exact e2e reverse matching:\t\tOK"
# else
#     $ECHO "Exact e2e reverse matching:\t\tFAILED"
# fi
# 
# ######################################################################
# 
# # Make exact reads from unique 22-mers
# $PYTHON reads_from_mers.py -q ref.fna ref22-unique-mer.out \
#     exact_ext_reads.fq exact_ext.ebwt_out 
# 
# if [ -e exact_ext.actual_out ] 
#     then
#     rm exact_ext.actual_out
# fi
# 
# # Run *extension* ebwt search on the reads...
# $CSAMapper/ebwt_search -q -k 22 ref exact_ext_reads.fq > exact_ext.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w exact_ext.ebwt_out exact_ext.actual_out > /dev/null; then
#     $ECHO "Exact ext. forward matching:\t\tOK"
# else
#     $ECHO "Exact ext. forward matching:\t\tFAILED"
# fi
# 
# ######################################################################
# 
# # Make exact (reverse complemented) reads from unique 22-mers
# $PYTHON reads_from_mers.py -q -r ref.fna ref22-unique-mer.out \
#     exact_ext_reads_rc.fq exact_ext_rc.ebwt_out 
# 
# if [ -e exact_ext_rc.actual_out ] 
#     then
#     rm exact_ext_rc.actual_out
# fi
# 
# # Run *extension* ebwt search on the reads...
# $CSAMapper/ebwt_search -k 22 -r ref -q exact_ext_reads_rc.fq > exact_ext_rc.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w exact_ext_rc.ebwt_out exact_ext_rc.actual_out > /dev/null; then
#     $ECHO "Exact ext. reverse matching:\t\tOK"
# else
#     $ECHO "Exact ext. reverse matching:\t\tFAILED"
# fi
# 
# ######################################################################
# 
# # Make reads from unique 22-mers that have mismatches at the 3' end
# $PYTHON reads_from_mers.py -q -3 3 ref.fna ref22-unique-mer.out \
#     3mis_reads.fq 3mis_reads.ebwt_out 
# 
# if [ -e 3mis_reads.actual_out ] 
#     then
#     rm 3mis_reads.actual_out
# fi
# 
# # Run ebwt search on the mismatching reads...
# $CSAMapper/ebwt_search -k 22 ref -q 3mis_reads.fq > 3mis_reads.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w 3mis_reads.ebwt_out 3mis_reads.actual_out > /dev/null; then
#     $ECHO "3' mismatch forward matching:\tOK"
# else
#     $ECHO "3' mismatch forward matching:\tFAILED"
# fi
# 
# ######################################################################
# 
# # Make reverse complemented reads from unique 22-mers that 
# # have mismatches at the 3' end
# $PYTHON reads_from_mers.py -r -3 3 -q ref.fna ref22-unique-mer.out \
#     3mis_reads_rc.fq 3mis_reads_rc.ebwt_out
# 
# if [ -e 3mis_reads_rc.actual_out ] 
#     then
#     rm 3mis_reads_rc.actual_out
# fi
# 
# # Run ebwt search on the mismatching reads...
# $CSAMapper/ebwt_search -q -r -k 22 ref 3mis_reads_rc.fq  > 3mis_reads_rc.actual_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w 3mis_reads_rc.ebwt_out 3mis_reads_rc.actual_out > /dev/null; then
#     $ECHO "3' mismatch reverse matching:\tOK"
# else
#     $ECHO "3' mismatch reverse matching:\tFAILED"
# fi
# 
# #####################################################################
# 
# # Make reads from unique 22-mers that have 1 mismatch in the 5' end
# $PYTHON reads_from_mers.py  -q -1 ref.fna ref22-unique-mer.out \
#     1_mis_reads.fq 1_mis_reads.ebwt_out
# 
# if [ -e 1_mis_reads.actual_out ] 
#     then
#     rm 1_mis_reads.actual_out
# fi
# 
# # Run ebwt search on the mismatching reads...
# $CSAMapper/ebwt_search -1 -q ref 1_mis_reads.fq | sort -n -k4 > 1_mis_reads.actual_out
# sort -n -k4 1_mis_reads.ebwt_out > xxx
# mv xxx 1_mis_reads.ebwt_out
# 
# # and then check that they all map to where they are supposed to
# if diff -w 1_mis_reads.ebwt_out 1_mis_reads.actual_out > /dev/null; then
#     $ECHO "1-mismatch in 5' mer forward matching:\tOK"
# else
#     $ECHO "1-mismatch in 5' mer forward matching:\tFAILED"
# fi

####################################################################

# Make reverse-complemented reads from unique 22-mers that have 1 mismatch 
# in the 5' end
$PYTHON reads_from_mers.py  -r -q -1 ref.fna ref22-unique-mer.out \
    1_mis_reads_rc.fq 1_mis_reads_rc.ebwt_out 

if [ -e 1_mis_reads.actual_out ] 
    then
    rm 1_mis_reads.actual_out
fi

# Run ebwt search on the mismatching reads...
$CSAMapper/ebwt_search -r -1 -q ref 1_mis_reads_rc.fq | sort -n -k4 > 1_mis_reads_rc.actual_out
sort -n -k4 1_mis_reads_rc.ebwt_out > xxx
mv xxx 1_mis_reads_rc.ebwt_out

# and then check that they all map to where they are supposed to
if diff -w 1_mis_reads_rc.ebwt_out 1_mis_reads_rc.actual_out > /dev/null; then
    $ECHO "1-mismatch in 5' mer reverse matching:\tOK"
else
    $ECHO "1-mismatch in 5' mer reverse matching:\tFAILED"
fi

# Run extension ebwt search on the mismatching reads...
$CSAMapper/ebwt_search -k 22 -r -1 -q ref 1_mis_reads_rc.fq | sort -n -k4 > 1_mis_reads_ext_rc.actual_out

# and then check that they all map to where they are supposed to
if diff -w 1_mis_reads_rc.ebwt_out 1_mis_reads_ext_rc.actual_out > /dev/null; then
    $ECHO "1-mismatch in 5' mer reverse extension matching:\tOK"
else
    $ECHO "1-mismatch in 5' mer reverse extension matching:\tFAILED"
fi



######################################################################

# Make reads from unique 22-mers that have 1 mismatch in the 22-mer and 
# several mismatches at the 3' end
$PYTHON reads_from_mers.py -3 3 -q -1 ref.fna ref22-unique-mer.out \
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


# Make reverse-complemented reads from unique 22-mers that have 1 mismatch 
# in the 22-mer and several mismatches at the 3' end
$PYTHON reads_from_mers.py -3 3 -q -1 -r ref.fna ref22-unique-mer.out \
    1_mis_with_3mis_reads_rc.fq 1_mis_with_3mis_reads_rc.ebwt_out 2>/dev/null

if [ -e 1_mis_with_3mis_reads_rc.actual_out ] 
    then
    rm 1_mis_with_3mis_reads_rc.actual_out
fi

# Run ebwt search on the mismatching reads...
$CSAMapper/ebwt_search -1 -k 22 -q -r ref 1_mis_with_3mis_reads_rc.fq | sort -n -k4 > 1_mis_with_3mis_reads_rc.actual_out
sort -n -k4 1_mis_with_3mis_reads_rc.ebwt_out > xxx
mv xxx 1_mis_with_3mis_reads_rc.ebwt_out

# and then check that they all map to where they are supposed to
if diff -w 1_mis_with_3mis_reads_rc.ebwt_out 1_mis_with_3mis_reads_rc.actual_out > /dev/null; then
    $ECHO "1-mismatch in 5' mer + 3' mismatch reverse matching:\tOK"
else
    $ECHO "1-mismatch in 5' mer + 3' mismatch reverse matching:\tFAILED"
fi
