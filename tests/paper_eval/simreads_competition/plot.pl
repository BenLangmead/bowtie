#!/usr/bin/perl -w

use strict;
use warnings;

defined($ARGV[0]) || die "Must specify run names";
my @runnames = $ARGV; # -> column names

sub readlines {
	my $f = shift;
	my @ret;
	open(FILE, $f) || die "Could not open $f";
	while(<FILE>) {
		chomp;
		push(@ret, $_)
	}
	return @ret;
}

foreach my $n (@runnames)  {
	my @names = readlines("$n.results.names.txt");
	my @times = readlines("$n.results.times.txt");
	my @vmmax = readlines("$n.results.vmmax.txt");
	my @rsmax = readlines("$n.results.rsmax.txt");
}
