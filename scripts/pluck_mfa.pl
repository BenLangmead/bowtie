#/usr/bin/perl -w

#
# Pluck a sequence fragment out of a multi-fasta file
#

use strict;
use warnings;

my $qid = 0;   # query id 
my $qoff = 0;  # query offset
my $qlen = 40; # query length

defined($ARGV[0]) || die "Must specify fasta file as first arg";
my $fin = $ARGV[0];
$qid  = $ARGV[1] if defined($ARGV[1]);
$qoff = $ARGV[2] if defined($ARGV[2]);
$qlen = $ARGV[3] if defined($ARGV[3]);

my $l = -1;  # line number within fasta file
my $id = -1; # index off fasta entry
my $off = 0; # offset within fasta entry

open(FASTA, $fin) || die "Could not open fasta file $fin";

while(<FASTA>) {
	my $line = $_;
	$l++;
	if(/^>/) {
		$off = 0;
		$id++;
		next;
	}
	# If this is the query entry... 
	if($id == $qid) {
		$line =~ s/[^a-zA-Z]//; # remove non-alphebetical chars
		$line =~ tr/a-z/A-Z/;   # make everytihng uppercase
		my $len = length($line);
		if($qoff >= $len) {
			$qoff -= $len;
		} else {
			for(my $i = 0; $qlen > 0 && $qoff + $i < $len; $i++) {
				print substr($line, $qoff + $i, 1);
				$qlen--;
			}
			$qoff = 0;
		}
		if($qlen == 0) {
			print "\n"; last;
		}
	} elsif($id > $qid) {
		last;
	}
}
