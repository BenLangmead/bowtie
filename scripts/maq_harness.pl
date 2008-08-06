#!/usr/bin/perl -w

use strict;
use warnings;

my $mapfile = ".tmp.map";
my $fafile  = ".tmp.fa";
my $bfafile = ".tmp.bfa";
my $fqfile  = ".tmp.fq";
my $bfqfile = ".tmp.bfq";

my $verbose = 1;
$#ARGV >= 1 || die "Takes at least two arguments: text, pat1";

my $text = shift @ARGV;
print "Text is: $text\n" if $verbose;

# Write text to a temporary fa file
open(FA, ">$fafile");
print FA ">tmp\n$text\n";
close(FA);

# Write reads to a temporary fq file
open(FQ, ">$fqfile");
my $i = 0;
for my $q (@ARGV) {
	my @qs = split(/_/, $q);
	$#qs >= 1 || die "Read should have two parts, separated by _: $q";
	length($qs[0]) == length($qs[1]) || die "Two parts of read should have same size: $q";
	print FQ "\@tmp$i\n$qs[0]\n+\n$qs[1]\n";
	$i++;
}
close(FQ);

system("maq fasta2bfa .tmp.fa .tmp.bfa") == 0 || die "maq fasta2bfa failed";
system("maq fastq2bfq .tmp.fq .tmp.bfq") == 0 || die "maq fastq2bfq failed";
system("maq map $mapfile $bfafile $bfqfile") == 0 || die "maq map failed";
system("maq mapview $mapfile");
