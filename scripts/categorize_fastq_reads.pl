#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("vc:",\%options);
my $cutoff = 0; $cutoff = $options{c} if defined($options{c});
my $verbose = 0; $verbose = $options{v} if defined($options{v});

my $thresh = 70;

my $allGtThresh = 0;
my $allGtThreshDiv2 = 0;
my $allBut1GtThreshDiv2 = 0;
my $allBut2GtThreshDiv2 = 0;
my $allBut3GtThreshDiv2 = 0;
my $allBut4GtThreshDiv2 = 0;
my $allGtThreshDiv3 = 0;
#my $allBut1GtThreshDiv3 = 0;
#my $allBut2GtThreshDiv3 = 0;
my $other = 0;

my $pats = 0;

while(<>) {
	my $l1 = $_;
	my $l2 = <>;
	my $l3 = <>;
	my $l4 = <>;
	chomp($l4);
	my $gtThresh = 0;
	my $gtThreshDiv2 = 0;
	my $gtThreshDiv3 = 0;
	for(my $i = 0; $i < length($l4); $i++) {
		my $q = ord(substr($l4, $i, 1)) - 33;
		print "$q " if $verbose;
		$gtThresh++ if $q > $thresh;
		$gtThreshDiv2++ if $q > ($thresh/2);
		#$gtThreshDiv3++ if $q > ($thresh/3);
	}
	print "\n" if $verbose;
	if($gtThresh == length($l4)) {
		$allGtThresh++;
	} elsif($gtThreshDiv2 == length($l4)) {
		$allGtThreshDiv2++;
	} elsif($gtThreshDiv2 == length($l4)-1) {
		$allBut1GtThreshDiv2++;
	} elsif($gtThreshDiv2 == length($l4)-2) {
		$allBut2GtThreshDiv2++;
	} elsif($gtThreshDiv2 == length($l4)-3) {
		$allBut3GtThreshDiv2++;
	} elsif($gtThreshDiv2 == length($l4)-4) {
		$allBut4GtThreshDiv2++;
#	} elsif($gtThreshDiv3 == length($l4)) {
#		$allGtThreshDiv2++;
#	} elsif($gtThreshDiv3 == length($l4)-1) {
#		$allBut1GtThreshDiv3++;
#	} elsif($gtThreshDiv3 == length($l4)-2) {
#		$allBut2GtThreshDiv3++;
	} else {
		$other++;
	}
	$pats++;
	last if $pats == $cutoff;
}

print "Threshold: $thresh\n";
printf "   All > thresh: $allGtThresh (%.1f%%)\n", $allGtThresh * 100.0 / $pats;
printf "   All > thresh/2: $allGtThreshDiv2 (%.1f%%)\n", $allGtThreshDiv2 * 100.0 / $pats;
printf "   All but 1 > thresh/2: $allBut1GtThreshDiv2 (%.1f%%)\n", $allBut1GtThreshDiv2 * 100.0 / $pats;
printf "   All but 2 > thresh/2: $allBut2GtThreshDiv2 (%.1f%%)\n", $allBut2GtThreshDiv2 * 100.0 / $pats;
printf "   All but 3 > thresh/2: $allBut3GtThreshDiv2 (%.1f%%)\n", $allBut3GtThreshDiv2 * 100.0 / $pats;
printf "   All but 4 > thresh/2: $allBut4GtThreshDiv2 (%.1f%%)\n", $allBut4GtThreshDiv2 * 100.0 / $pats;
#printf "   All > thresh/3: $allGtThreshDiv3 (%.1f%%)\n", $allGtThreshDiv3 * 100.0 / $pats;
#printf "   All but 1 > thresh/3: $allBut1GtThreshDiv3 (%.1f%%)\n", $allBut1GtThreshDiv3 * 100.0 / $pats;
#printf "   All but 2 > thresh/3: $allBut2GtThreshDiv3 (%.1f%%)\n", $allBut2GtThreshDiv3 * 100.0 / $pats;
printf "   Other: $other (%.1f%%)\n", $other * 100.0 / $pats;
