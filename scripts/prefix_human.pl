#/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;

my @prefixes = ( 100000000,
	             250000000,
	             500000000,
	            1000000000,
	            2000000000 );
my %handles = ();
my $whole = new FileHandle ">hs_ref_whole.mfa";

for my $f (@prefixes) {
	$handles{$f} = new FileHandle ">hs_ref_first$f.mfa";
	print { $handles{$f} } ">$f\n";
}
print $whole ">whole\n";

my $sofar = 0;
while(<>) {
	next if /^>/;
	chomp;
	for my $f (@prefixes) {
		if($sofar < $f) {
			next if($f == 2000000000 && $sofar < 1000000000);
			print { $handles{$f} } "$_";
		}
	}
	if($sofar > 2000000000) {
		print $whole "$_";
	}
	$sofar += length($_);
}

for my $f (@prefixes) {
	print { $handles{$f} } "\n";
	$handles{$f}->close();
}
print $whole "\n";
$whole->close();
