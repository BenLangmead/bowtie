#/usr/bin/perl -w

use strict;
use warnings;

my %vhitsFw = ();
my %vhitsRc = ();

my $lastS5 = "-1";
my $lastD = "D";
while(<>) {
	my @s = split;
	next if $#s < 7;
	next if /^#/;
	my $vhits = ($s[3] eq "D") ? \%vhitsFw : \%vhitsRc;
	$vhits->{$s[5]} = [] unless defined ($vhits->{$s[5]});
	push(@{$vhits->{$s[5]}}, (int($s[1]), int($s[2])));
	if($s[5] ne $lastS5 && $lastS5 ne "-1") {
		print "$lastS5";
		if($s[3] eq "D") {
			print "+";
		} else {
			print "-";
		}
		print ": @{$vhits->{$lastS5}}\n";
	}
	$lastS5 = $s[5];
	$lastD = $s[3];
}
print "$lastS5";
if($lastD eq "D") {
	print "+";
} else {
	print "-";
}
my $vhits = ($lastD eq "D") ? \%vhitsFw : \%vhitsRc;
print ": @{$vhits->{$lastS5}}\n";
