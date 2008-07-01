#!/usr/bin/perl -w

#
# Don't pass -t (or anything else that causes non-hit output) to
# ebwt_search if you're piping its output to this script
#

use warnings;
use strict;
use Getopt::Std;

my %options=();
getopts("sr:",\%options);
my $summarize = 0; $summarize = $options{s} if defined $options{s};
defined($options{r}) || die "Must specify number of reads with -r";
my $reads = $options{r};

my %cnt0 = ();
my %cnt1 = ();
my %both = ();

while(<>) {
	next unless /^[0-9]+/;
	chomp();
	my @cs = split(/:</);
	my $q = $cs[0];
	my $fw = ((chop($q) eq "+")? 1 : 0);
	$q = int($q);
	my @hs = split(/>,</, $cs[1]);
	foreach my $h (@hs) {
		my @cs = split(/,/, $h);
		my $diff = int($cs[1]) - int($cs[0]);
		my $mm = 0;
		chop($cs[2]) if length($cs[2]) > 1; 
		length($cs[2]) == 1 || die "mm string too long: $cs[2]";
		$mm = int($cs[2]);
		$mm == 0 || $mm == 1 || die "mm must be 0 or 1: $mm";
		my $cnt = (($mm == 0)? \%cnt0 : \%cnt1);
		$both{$q} = 0 unless defined($both{$q});
		$cnt->{$q} = 0 unless defined($cnt->{$q});
		$cnt->{$q} += $diff;
	}
}

my $qsThatHit = 0;
my $hit0s = 0;
my $hit1s = 0;
my $unique0s = 0;
my $gt1_0s = 0;
my $gt10_0s = 0;
my $gt100_0s = 0;
my $gt1000_0s = 0;
my $gt10000_0s = 0;
my $gt100000_0s = 0;
my $gt1000000_0s = 0;
my $unique1s = 0;
my $gt1_1s = 0;
my $gt10_1s = 0;
my $gt100_1s = 0;
my $gt1000_1s = 0;
my $gt10000_1s = 0;
my $gt100000_1s = 0;
my $gt1000000_1s = 0;

for my $k (keys %both) {
	$qsThatHit++;
	if(defined($cnt0{$k})) {
		$hit0s++;
		print "$k:<$cnt0{$k},0>\n" unless $summarize;
		$unique0s++    if $cnt0{$k} == 1;
		$gt1_0s++      if $cnt0{$k} > 1;
		$gt10_0s++     if $cnt0{$k} > 10;
		$gt100_0s++    if $cnt0{$k} > 100;
		$gt1000_0s++   if $cnt0{$k} > 1000;
		$gt10000_0s++  if $cnt0{$k} > 10000;
		$gt100000_0s++ if $cnt0{$k} > 100000;
		$gt1000000_0s++ if $cnt0{$k} > 1000000;
	} else {
		$hit1s++;
		print "$k:<$cnt1{$k},1>\n" unless $summarize;
		$unique1s++    if $cnt1{$k} == 1;
		$gt1_1s++      if $cnt1{$k} > 1;
		$gt10_1s++     if $cnt1{$k} > 10;
		$gt100_1s++    if $cnt1{$k} > 100;
		$gt1000_1s++   if $cnt1{$k} > 1000;
		$gt10000_1s++  if $cnt1{$k} > 10000;
		$gt100000_1s++ if $cnt1{$k} > 100000;
		$gt1000000_1s++ if $cnt1{$k} > 1000000;
	}
}

sub pctout($) {
	my $n = shift;
	print "$n (";
	printf "%2.1f%%", ($n*100.0)/$reads;
	print ")\n";
}

if($summarize) {
	print "Queries that hit: "; pctout($qsThatHit);
	print "   with at least 1 exact hit: "; pctout($hit0s);
	print "      1 exact hit: "; pctout($unique0s);
	print "      >1 exact hits: "; pctout($gt1_0s);
	print "        >10 exact hits: "; pctout($gt10_0s);
	print "          >100 exact hits: "; pctout($gt100_0s);
	print "            >1000 exact hits: "; pctout($gt1000_0s);
	print "              >10000 exact hits: "; pctout($gt10000_0s);
	print "                >100000 exact hits: "; pctout($gt100000_0s);
	print "                  >1000000 exact hits: "; pctout($gt1000000_0s);
	print "   with at least 1 1-mismatch hit (and no exact hits): "; pctout($hit1s);
	print "      1 1-mismatch hit: "; pctout($unique1s);
	print "      >1 1-mismatch hits: "; pctout($gt1_1s);
	print "        >10 1-mismatch hits: "; pctout($gt10_1s);
	print "          >100 1-mismatch hits: "; pctout($gt100_1s);
	print "            >1000 1-mismatch hits: "; pctout($gt1000_1s);
	print "              >10000 1-mismatch hits: "; pctout($gt10000_1s);
	print "                >100000 1-mismatch hits: "; pctout($gt100000_1s);
	print "                  >1000000 1-mismatch hits: "; pctout($gt1000000_1s);
}
