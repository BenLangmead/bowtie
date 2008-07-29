#!/usr/bin/perl -w

defined($ARGV[0])
	|| die "Must specify number of reads mapped by both Maq and Bowtie as first argument";
defined($ARGV[1])
	|| die "Must specify number of reads mapped by Bowtie but not Maq as second argument";
defined($ARGV[2])
	|| die "Must specify number of reads mapped by Maq but not Bowtie as third argument";

my $inboth = int($ARGV[0]);
my $inebwt = int($ARGV[1]);
my $inmaq  = int($ARGV[2]);
my $tot = $inboth + $inebwt + $inmaq;

my $bothpct = sprintf("%2.1f%%", $inboth * 100.0 / $tot);
my $ebwtpct = sprintf("%2.1f%%", $inebwt * 100.0 / $tot);
my $maqpct  = sprintf("%2.1f%%", $inmaq  * 100.0 / $tot);

print "  In either: $tot\n";
print "    In noth: $inboth ($bothpct)\n";
print "Bowtie only: $inebwt ($ebwtpct)\n";
print "   Maq only: $inmaq ($maqpct)\n";
