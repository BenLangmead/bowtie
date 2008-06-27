#/usr/bin/perl -w

#
# Split a large fasta file up into a set of smaller files of a
# specified maximum size.
#

use strict;
use warnings;

defined($ARGV[0]) || die "Must specify fasta input file as first argument";
my $infa = $ARGV[0];
my $inbase = $infa;
$inbase =~ s/\.m?fa//;
defined($ARGV[1]) || die "Must specify max reads per file as second argument";
my $max = int($ARGV[1]);
print "input file: $infa\n";
print "max reads: $max\n";
open(IN, $infa) || die("could not open $infa");
my $reads = 0;
my $curfile = 0;
open(CUR, ">".$inbase."_".$curfile.".fa")
	|| die("could not open ".$inbase."_".$curfile.".fa");
print "Writing ".$inbase."_".$curfile.".fa\n";
while(<IN>) {
	$reads++ if /^>/;
	if($reads > $max) {
		$curfile++;
		$reads = 1;
		open(CUR, ">".$inbase."_".$curfile.".fa")
			|| die("could not open ".$inbase."_".$curfile.".fa");
		print "Writing ".$inbase."_".$curfile.".fa\n";
	}
	print CUR $_;
}
