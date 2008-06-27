#!/usr/bin/perl -w

#
# Extract some random substrings of the sequence in the given fasta
# file and search for them in the ebwt.  We assume all substrings
# should be found, and we abort with an error message if any aren't.
#
# Usage: perl random_ebwt_tester2.pl <fasta_file> <ebwt>
#

use strict;
use warnings;

# Get fasta input file
defined($ARGV[0]) || die "Must specify fasta file as first argument";
my $fasta = $ARGV[0];

# Get EBWT input file
defined($ARGV[1]) || die "Must specify EBWT file as second argument";
my $ebwt = $ARGV[1];

my $seqno = -1;
my $off = 0;
open(FASTA, $fasta) || die "Could not open $fasta";
while(<FASTA>) {
	chomp;
	if(/^>/) {
		$seqno++;
		next;
	}
	my $oldoff = $off;
	$off += length($_);
	next if /[^actgACTG]/; # Skip lines with non-ATCG codes
	#my @ebwt_out = `./ebwt_search-with-asserts -c $ebwt $_`;
	#$#ebwt_out == 0 || die "output had more than one line";
	#$ebwt_out[0] =~ /0[+]:.*<.+,$oldoff>/ || die "output did not contain $oldoff";
	#print $ebwt_out[0];
	print ">$oldoff\n$_\n";
}
