#!/usr/bin/perl -w

use warnings;
use strict;

my %reads = ();

while(<STDIN>) {
	chomp;
	defined($reads{$_}) && die "Already added to set: $_";
	$reads{$_} = 1;
}

open(MAQ, "maq mapview $ARGV[0] |");
while(<MAQ>) {
	my @s = split;
	if(defined($reads{$s[0]})) {
		print $_;
	}
}
