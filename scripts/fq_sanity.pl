#!/usr/bin/perl -w

#
# Sanity check a given FASTQ file and abort with a descriptive message
# it it seems malformed.  Malformed files can be cleaned up with
# fq_filter.pl.
# 

use strict;
use warnings;

my $lines = 0;
my $seqlen = 0;
while(<>) {
	chomp;
	my $mod4 = $lines % 4;
	my $ulen = length($_);
	if($mod4 == 0) {
		/^@/ || die "Bad 0th line: $_";
	}
	elsif($mod4 == 1) {
		if($lines == 1) {
			$seqlen = length($_);
		} else {
			$ulen == $seqlen || die "Bad sequence length: $ulen (expected $seqlen)";
		}
	}
	elsif($mod4 == 2) {
		/^\+/ || die "Bad 2nd line: $_";
	}
	else {
		$ulen == $seqlen || die "Bad quality length: $ulen (expected $seqlen)";
	}
	$lines++;
}
print "PASSED\n";
