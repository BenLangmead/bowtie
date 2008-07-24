#!/usr/bin/perl -w

#
# Filter out malformed reads in a FASTQ file and print filtered output
# to STDOUT.  Informative messages are printed to STDERR.
# 

use strict;
use warnings;

my $read = 0;
my $seqlen = 0;
my $lastwasbad = 0;
while(<>) {
	my $bad = 0;
	my $l1 = $_; chomp($l1);
	if($lastwasbad) {
		for(my $i = 0; $i < 5; $i++) {
			if($l1 =~ /^@/) {
				last;
			} else {
				$l1 = <>;
			}
		}
		unless($l1 =~ /^@/) {
			die "catastrophic failure!";
		}
		$lastwasbad = 0;
	}
	my $l2 = <>; chomp($l2);
	my $l3 = <>; chomp($l3);
	my $l4 = <>; chomp($l4);
	
	my $ulen2 = length($l2);
	my $ulen4 = length($l4);
	unless($l1 =~ /^@/) {
		print STDERR "Bad 1st line for read $read: $_";
		$bad = 1;
	}
	if($read == 0) {
		$seqlen = $ulen2;
	} else {
		if($ulen2 != $seqlen) {
			print STDERR "Bad sequence length for read $read: $ulen2 (expected $seqlen)\n";
			$bad = 1;
		}
	}
	unless($l2 =~ /^[ACTGNactgn]*$/) {
		print STDERR "Bad sequence content for read $read: $l2\n";
		$bad = 1;
	}
	unless($l3 =~ /^\+/) {
		print STDERR "Bad 3rd line for read $read: $_\n";
		$bad = 1;
	} 
	if($ulen4 != $seqlen) {
		print STDERR "Bad quality length for read $read: $ulen4 (expected $seqlen)\n";
		$bad = 1;
	}
	
	if(!$bad) {
		print "$l1\n$l2\n$l3\n$l4\n";
	} else {
		# Separate complaints
		print STDERR "---\n";
		$lastwasbad = 1;
	}
	$read++;
}
print STDERR "DONE\n";
