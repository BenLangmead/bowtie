#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("f:",\%options);
my $filter = 0;
$filter = $options{f} if defined($options{f});

my $cpuCol = -1;
while(<>) {
	if(/PID.*USER.*PR/) {
	    # Grab column headers
	    my @line = split;
	    for(my $i = 0; $i <= $#line; $i++) {
	    	if($line[$i] eq "%CPU") {
	    		$cpuCol = $i;
	    		last;
	    	}
	    }
	}
	next unless /^[ ]*[0-9]+/;
	my @line = split;
	next if $line[$cpuCol] <= $filter;
	print "$line[$cpuCol]\n";
}

