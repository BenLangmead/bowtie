#/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("r:",\%options);

my $rate1 = 100;
my $rate2 = 200;

my @dnaMap = ("A", "T", "C", "G");

while(<>) {
	if(/^>/) {
		print $_;
		next;
	} 
	# sequence line
	my $r1 = int(rand($rate1));
	my $r2 = int(rand($rate2));
	my $s = $_;
	if($r1 < length($_)-1) {
		substr($s, $r1, 1, $dnaMap[int(rand(4))]);
		length($s) == length($_) || die "Length changed";
	}
	if($r2 < length($_)-1) {
		substr($s, $r2, 1, $dnaMap[int(rand(4))]);
		length($s) == length($_) || die "Length changed";
	}
	length($s) == length($_) || die "Length changed";
	print $s;
}
