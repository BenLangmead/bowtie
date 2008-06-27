#/usr/bin/perl -w

#
# Usage: perl summarize_solexa.pl <foo_seq.txt[.gz]> [<read_limit>]
#
# Will print various statistics for occurrence rates of dots and
# characters in (possibly gzipped) Solexa _seq.txt file.  Will stop
# after read_limit reads.
#

use strict;
use warnings;

my $cutoff = -1;
$cutoff = $ARGV[1] || defined($ARGV[1]);

# Handles up to 64 positions
my @positions = (
	{}, {}, {}, {}, {}, {}, {}, {},
	{}, {}, {}, {}, {}, {}, {}, {},
	{}, {}, {}, {}, {}, {}, {}, {},
	{}, {}, {}, {}, {}, {}, {}, {},
	{}, {}, {}, {}, {}, {}, {}, {},
	{}, {}, {}, {}, {}, {}, {}, {},
	{}, {}, {}, {}, {}, {}, {}, {},
	{}, {}, {}, {}, {}, {}, {}, {},
);
my $with_dots = 0;
my $with_8_dots = 0;
my %lengths = ();
my %numDots = ();

sub print_hash {
	my ($hs) = @_;
	my $tot = 0;
	my $k;
	for $k (keys %{$hs}) { $tot += $hs->{$k} };
	for $k ("A", "C", "G", "T", ".") {
		my $val = 0;
		$val = $hs->{$k} if defined($hs->{$k});
		my $frac = $val * 100.0 / $tot;
		printf "%2.1f\t", $frac;
	}
	print "\n";
}

my $i = 0;
my $o = 0;
if($ARGV[0] =~ /\.gz$/) {
	open(SOL, "zcat $ARGV[0] |") || die "Couldn't open gzipped Solexa file";
} else {
	open(SOL, "cat $ARGV[0] |") || die "Couldn't open text Solexa file";
}
while(<SOL>) {
	$o++;
	my @line = split;
	my $seq = $line[4];
	my $dots = 0;
	for(my $j = 0; $j < length($seq); $j++) {
		if(substr($seq, $j, 1) eq ".") { $dots++; }
	}
	if($dots > 0) {
		$with_dots++;
		$with_8_dots++ if $dots > 7;
	} else {
		$i++;
	}
	if($dots <= 9) { $dots = "0$dots"; }
	$numDots{$dots} = 0 unless defined($numDots{$dots});
	$numDots{$dots}++;
	$lengths{length($seq)} = 0 unless defined($lengths{length($seq)});
	$lengths{length($seq)}++;
	for(my $j = 0; $j < length($seq); $j++) {
		my $c = substr($seq, $j, 1);
		$positions[$j]->{$c} = 0 unless defined($positions[$j]->{$c});
		$positions[$j]->{$c}++;
	}
	last if $i == $cutoff;
}

print "# Reads: $o\n";
printf "  # Good reads: $i (%2.1f%% of total)\n", ($i * 100.0 / $o);
printf "  # Reads with dots: $with_dots (%2.1f%% of total)\n", ($with_dots * 100.0 / $o);
printf "    # Reads with 8 or more dots: $with_8_dots (%2.1f%% of total)\n", ($with_8_dots * 100.0 / $o);
print "Number of reads (RHS) with given length (LHS):\n";
for my $k (sort keys %lengths) {
	print "  $k: $lengths{$k}\n";
}
print "Number of reads (RHS) with given number of dots (LHS):\n";
for my $k (sort keys %numDots) {
	print "  $k: $numDots{$k}\n";
}
print "Per-position occurrence rates:\n";
print "  Pos\tA\tC\tG\tT\t.\n";
print "  ---\t----\t----\t----\t----\t----\n";
for(my $j = 0; $j < $#positions; $j++) {
	my @ks = keys %{$positions[$j]};
	last if $#ks < 0; # no bases at this position
	print "  $j\t";
	print_hash($positions[$j]);
}
