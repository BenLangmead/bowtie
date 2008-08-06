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

system("maq fasta2bfa .tmp.fa .tmp.bfa > /dev/null") == 0 || die "maq fasta2bfa failed";
system("maq fastq2bfq .tmp.fq .tmp.bfq > /dev/null") == 0 || die "maq fastq2bfq failed";
system("maq map $mapfile $bfafile $bfqfile > /dev/null") == 0 || die "maq map failed";

# fprintf(fpout, "%s\t%s\t%d\t%c\t%d\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
#       m1->name,               // read name
#		m->ref_name[m1->seqid], // reference name
#		(m1->pos>>1) + 1,       // reference offset
#		(m1->pos&1)? '-' : '+', // orientation of alignment
#		m1->dist,               // offset of mate
#		m1->flag,               // status of pair
#		m1->map_qual,           // "final mapping quality"
#		(signed char)m1->seq[MAX_READLEN-1], // single-end mapping quality
#		m1->alt_qual,           // the lower quality of the two ends
#		m1->info1&0xf,          // # mismatches (in the best hit?)
#		m1->info2,              // sum of errors of best hit
#		m1->c[0],               // # of exact hits
#		m1->c[1],               // # of 1-mismatch hits
#		m1->size);              // length of alignment

print "read_nm\tref_nm\trefoff\torient\tmate_off\tpair_stat\n";
system("maq mapview $mapfile");
