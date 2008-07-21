#!/usr/bin/perl -w

use strict;
use warnings;

#
# Quick-and-dirty way to sanity-check that a pattern set was read
# correctly.
#
# Example usage:
# $ /ebwt_search-with-asserts \
#      --dumpPats pats.txt -1frt -u 10000 \
#      $EBWTS/whole $READS/s_7_0000_0255.fa /dev/null
# $ perl scripts/check_fasta_patdump.pl < pats.txt
# $
#

my $line = 0;
my $seqno = -1;
my $lastName = "";
my $lastSeq = "";
my $turnovers = 1;
my $groupSize = 10000;
for(my $i = 0; $i < $turnovers; $i++) {
	my @seqs = ();
	while(<>) {
		my $name = $lastName;
		my $seq = "";
		my $headerLine = 0;
		if(/^>(.*)$/) {
			$seqno++; $name = $1;
			$headerLine = 1;
		} else {
			chomp;
			$seq = $_;
			push(@seqs, $seq);
		}
		if($headerLine) {
			if(($seqno % 2) == 0) {
				$name ne $lastName || die "New name $name ($seqno) should not match last name $lastName";
				$lastName = $name;
			} else {
				$name eq $lastName || die "New name $name ($seqno) should match last name $lastName";
			}
		}
		else {
			if(($seqno % 2) == 1 && length($seq) > 0) {
				$seq = reverse($seq);
				$seq =~ tr/aAcCgGtT/tTgGcCaA/;
				$seq eq $lastSeq || die "Seq should have been revcomp of lastSeq $lastSeq";
			}
			if(length($seq) > 0) {
				$lastSeq = $seq;
			}
		}
		$line++;
		last if (!$headerLine && $seqno == $groupSize*2-1);
	}
	$line == $groupSize*4 || die "Expected to leave first loop on line ".($groupSize*4)." but was $line";
	$seqno = -1;
	while(<>) {
		my $name = $lastName;
		my $seq = "";
		my $headerLine = 0;
		if(/^>(.*)$/) {
			$seqno++; $name = $1;
			$headerLine = 1;
		} else {
			chomp;
			$seq = $_;
			$seqs[$seqno] eq reverse($seq) || die "reverse-mismatch";
		}
		if($headerLine) {
			if(($seqno % 2) == 0) {
				$name ne $lastName || die "New name $name ($seqno) should not match last name $lastName";
				$lastName = $name;
			} else {
				$name eq $lastName || die "New name $name ($seqno) should match last name $lastName";
			}
		}
		else {
			if(($seqno % 2) == 1 && length($seq) > 0) {
				$seq = reverse($seq);
				$seq =~ tr/aAcCgGtT/tTgGcCaA/;
				$seq eq $lastSeq || die "Seq should have been revcomp of lastSeq $lastSeq";
			}
			if(length($seq) > 0) {
				$lastSeq = $seq;
			}
		}
		last if (!$headerLine && $seqno == $groupSize*2);
		$line++;
	}
}
