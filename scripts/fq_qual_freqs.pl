#!/usr/bin/perl -w

#
# Calculate and print "average" quality values for each position in a
# Solexa FASTQ file.
# 

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("c:",\%options);

my $cutoff = 0xffffffff; $cutoff = $options{c} if defined($options{c});

open(FQ, $ARGV[0]) || die "Could not open .fq";
my $lines = 0;
my @tots = (
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
);
my @mins = (
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
);
my @maxs = (
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
);
while(<FQ>) {
	next if /^@/;
	next if /^\+$/;
	next if /^[ACGTNacgtn]+$/;
	for(my $i = 0; $i < length($_) && $i < 40; $i++) {
		my $oi = (ord(substr($_, $i, 1)) - 32);
		$tots[$i] += $oi;
		if($oi > $mins[$i]) {
			$mins[$i] = $oi;
		}
		if($oi < $maxs[$i]) {
			$maxs[$i] = $oi;
		}
	}
	$lines++;
	last if $lines >= $cutoff;
}
# Print averages
print "Average: ";
for(my $i = 0; $i < 40; $i++) {
	last if $tots[$i] == 0;
	my $q = $tots[$i] * 1.0 / $lines;
	my $rq = int($q + 0.5);
	print chr($rq+32);
}
print "\n";
# Print mins and maxs
#for(my $i = 0; $i < 40; $i++) {
#	last if $tots[$i] == 0;
#	print "$i min/max: ".chr($mins[$i]+32)."/".chr($maxs[$i]+32)."\n";
#}
