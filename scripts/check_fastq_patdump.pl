#!/usr/bin/perl -w

#
# Quick-and-dirty way to sanity-check that a pattern set was read
# correctly.
#
# Example usage:
# $ /ebwt_search-with-asserts \
#      --dumpPats pats.txt -1qrt -u 10000 \
#      $EBWTS/whole $READS/s_7_0000_0255.fastq /dev/null
# $ perl scripts/check_fastq_patdump.pl < pats.txt
# $
#

use strict;
use warnings;

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
		my $seqLine = 0;
		if(($line % 4) == 0) {
			/^@(.*)$/ || die "Expected a \@ line";
			$seqno++; $name = $1;
			$headerLine = 1;
		} elsif(($line % 4) == 1) {
			chomp;
			$seq = $_;
			push(@seqs, $seq);
			$seqLine = 1;
		}
		if($headerLine) {
			if(($seqno % 2) == 0) {
				$name ne $lastName || die "New name $name ($seqno) should not match last name $lastName";
				$lastName = $name;
			} else {
				$name eq $lastName || die "New name $name ($seqno) should match last name $lastName";
			}
		}
		if($seqLine) {
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
		if ($seqLine && $seqno == $groupSize*2-1) {
			<>; $line++; <>; $line++; # skip over quals
			last;
		}
	}
	$line == $groupSize*8 || die "Expected to leave first loop on line ".($groupSize*8)." but was $line";
	$seqno = -1;
	while(<>) {
		my $name = $lastName;
		my $seq = "";
		my $headerLine = 0;
		my $seqLine = 0;
		if(($line % 4) == 0) {
			/^@(.*)$/ || die "Expected a \@ line";
			$seqno++; $name = $1;
			$headerLine = 1;
		} elsif(($line % 4) == 1) {
			chomp;
			$seq = $_;
			$seqs[$seqno] eq reverse($seq) || die "reverse-mismatch";
			$seqLine = 1;
		}
		if($headerLine) {
			if(($seqno % 2) == 0) {
				$name ne $lastName || die "New name $name ($seqno) should not match last name $lastName";
				$lastName = $name;
			} else {
				$name eq $lastName || die "New name $name ($seqno) should match last name $lastName";
			}
		}
		if($seqLine) {
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
	}
}
