#/usr/bin/perl -w

#
# Split a large fastq file up into a set of smaller files of a
# specified maximum size.
#

use strict;
use warnings;

defined($ARGV[0]) || die "Must specify fqsta input file as first argument";
my $infq = $ARGV[0];
my $inbase = $infq;
$inbase =~ s/\.f(ast)?q//;
defined($ARGV[1]) || die "Must specify max reads per file as second argument";
my $max = int($ARGV[1]);
print "input file: $infq\n";
print "max reads: $max\n";
open(IN, $infq) || die("could not open $infq");
my $reads = 0;
my $line = 0;
my $curfile = 0;
open(CUR, ">".$inbase."_".$curfile.".fq")
	|| die("could not open ".$inbase."_".$curfile.".fq");
print "Writing ".$inbase."_".$curfile.".fq\n";
while(<IN>) {
	$reads++ if (/^[@]/ && ($line % 4 == 0));
	if($reads > $max) {
		$curfile++;
		$reads = 1;
		open(CUR, ">".$inbase."_".$curfile.".fq")
			|| die("could not open ".$inbase."_".$curfile.".fq");
		print "Writing ".$inbase."_".$curfile.".fq\n";
	}
	$line++;
	print CUR $_;
}
