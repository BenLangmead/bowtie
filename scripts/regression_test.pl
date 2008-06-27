#/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("vas",\%options);

my $asserts = 0; # Enable -with-asserts version?
$asserts = $options{a} if defined($options{a});
my $sanity = 0;  # Enable Ebwt sanity checks?
$sanity = $options{s} if defined($options{s});
my $verbose = 0;
$verbose = $options{v} if defined($options{v});

my $vmatch = 1;

defined($ARGV[0]) || die "No reference text given";
my $ref = $ARGV[0];

defined ($ARGV[1]) || die "No query file given";
my $qry = $ARGV[1];

my $ebwt_args = "-fa";
my $vmatch_args = "-d -complete";

# Make version of query file with comments stripped
system("egrep -v '^#' $qry > .qry.fa");
system("egrep -v '^#' $ref > .ref.fa");

# Open query file
open(QRYS, "$qry");
{
	# Read arguments for ebwt_search and vmatch
	my $s = <QRYS>;
	chomp($s);
	my @ss = split (/:/, $s);
	$ebwt_args = $ss[1] if defined $ss[1];
	$s = <QRYS>;
	chomp($s);
	@ss = split (/:/, $s);
	$vmatch_args = $ss[1] if defined $ss[1];
	close(QRYS);
}

# Possible add sanity-checking arguments to ebwt_search
$ebwt_args .= " -s --orig .ref.fa" if $sanity;

print "ebwt_search args: $ebwt_args\n" if $verbose;
print "vmatch args: $vmatch_args\n" if $verbose;

# Build the Ebwt index
my $ebwt_build = "./ebwt_build";
$ebwt_build .= "-with-asserts" if $asserts;
my $cmd = "$ebwt_build -d .ref.fa --entireSA .regr";
print "$cmd\n" if $verbose;
system($cmd);

# Search the Ebwt index
my $ebwt_search = "./ebwt_search";
$ebwt_search .= "-with-asserts" if $asserts;
$cmd = "$ebwt_search $ebwt_args .regr $qry .regr.hits";
print "$cmd\n" if $verbose;
system($cmd);

# Write an alphabet file for vmatch that treats all wildcards as the
# same as As
open(ALPHA, ">.regr.alpha");
print ALPHA "aAnsywrkvbdhmNSYWRKVBDHM\n";
print ALPHA "cC\n";
print ALPHA "gG\n";
print ALPHA "tTuU\n";
print ALPHA "-\n";
close(ALPHA);

# Build a Vmatch index from the reference
$cmd = "mkvtree -db .ref.fa -smap .regr.alpha -indexname .regr -allout -pl";
print "$cmd\n" if $verbose;
system($cmd) == 0 || die "Error running mkvtree";

# Search the Vmatch index
$cmd = "vmatch -q .qry.fa $vmatch_args -v .regr > .regr.vmatch.hits";
print "$cmd\n" if $verbose;
system($cmd) == 0 || die "Error running vmatch";

my %expecteds = ();

open(QRYS, ".qry.fa");
while(<QRYS>) {
	next unless /^>[0-9]+[+-]:/;
	chomp;
	my @colon = split (/:</);
	my $q = substr($colon[0], 1);
	my @comma = split (/,</, $colon[1]);
	for my $h (@comma) {
		$expecteds{"$q:<$h"} = 0;
	}
}
close(QRYS);

my $matchups = 0;
my $unexpected = 0;
my $expectedNotFound = 0;

my $vmatchups = 0;
my $vunexpected = 0;
my $vexpectedNotFound = 0;

# Make a copy of the expecteds map
my %vexpecteds = %expecteds;

open(HITS, ".regr.hits");
while(<HITS>) {
	chomp;
	my @colon = split (/:</);
	my $q = $colon[0];
	my @comma = split (/,</, $colon[1]);
	for my $h (@comma) {
		my $k = "$q:<$h";
		if(defined($expecteds{$k})) {
			$matchups++;
			delete $expecteds{$k};
		} else {
			$unexpected++;
			print "Unexpected hit: $k\n";
		}
	}
}
close(HITS);

for my $k (keys %expecteds)  {
	print "Expected hit did not occur in Ebwt output: $k\n";
	$expectedNotFound++;
}

open(VHITS, ".regr.vmatch.hits");
while(<VHITS>) {
	# E.g.:
	# len seq soff o len  pat poff mm   pct id
	# 26    0 24   D 26     0   0   0   100.00
	next if /^#/;
	chomp;
	my @vl = split;
	my $o = ($vl[3] eq "D")? "+":"-";
	my $mm = ($vl[7] eq "-1")? "1":"0";
	my $s = "$vl[5]$o:<$vl[1],$vl[2],$mm>";
	if(defined($vexpecteds{$s})) {
		$vmatchups++;
		delete $vexpecteds{$s};
	} else {
		$vunexpected++;
		print "Unexpected Vmatch hit: $s\n";
	}
}
close(VHITS);

for my $k (keys %vexpecteds)  {
	print "Expected hit did not occur in Vmatch output: $k\n";
	$vexpectedNotFound++;
}

print "Matched up: $matchups\n";
print "Unexpected hits: $unexpected\n";
print "Expected hits not found: $expectedNotFound\n";

print "Matched up (Vmatch): $vmatchups\n";
print "Unexpected hits (Vmatch): $vunexpected\n";
print "Expected hits not found (Vmatch): $vexpectedNotFound\n";

if($unexpected + $expectedNotFound + $vunexpected + $vexpectedNotFound == 0) {
	print "PASSED\n";
}
