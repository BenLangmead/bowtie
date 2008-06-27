#!/usr/bin/perl -w

#
# 
#
# Usage: perl do_test_fasta.pl <fasta_file> <ebwt>
#

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("s",\%options);

# Get fasta input file
defined($ARGV[0]) || die "Must specify fasta file as first argument";
my $fasta = $ARGV[0];

# Get EBWT input file
defined($ARGV[1]) || die "Must specify EBWT file as second argument";
my $ebwt = $ARGV[1];

if(!(defined $options{s})) {
	system("./ebwt_search-with-asserts -f $ebwt $fasta > .out.tmp") == 0 || die "ebwt_search command failed";
}

open(OUT, ".out.tmp") || die "Could not open .out.tmp";
open(FASTA, "grep '^>' $fasta |") || die "Could not open $fasta";
my $i = 0;
while(<OUT>) {
	my $fline = substr(<FASTA>, 1);
	chomp($fline);
	/$i[+]:.*<.+,$fline>/ || die "Did not find $fline ($i)";
	$i++;
}
