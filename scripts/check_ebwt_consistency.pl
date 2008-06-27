#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %options=();
getopts("c",\%options);
my $canonicalize = 0; $canonicalize = $options{c} if defined $options{c};

foreach my $n ("1", "2") {
	# Mapping from XXX from first observed XXX_*.ebwt to full filename of
	# observed file.
	my %map = ();
	open(LS, "ls -a *.$n.ebwt|");
	while(<LS>) {
		chomp;
		my @d = split /\./;
		$d[$#d]   eq "ebwt" || die;
		$d[$#d-1] eq $n || die;
		my @s = split /_/;
		if(defined $map{$s[0]}) {
			# We've already seen another file with the same prefix, so
			# diff this file against that one, asserting that they should
			# be the same
			print "Comparing $_ to $map{$s[0]}\n";
			system("diff $_ $map{$s[0]}") == 0 || die "$_ didn't match $map{$s[0]}";
			system("rm -f $_") if $canonicalize;
		} else {
			# Keep this for comparison with others later
			$map{$s[0]} = $_;
		}
	}
	if($canonicalize) {
		foreach my $f (keys %map) {
			# Make canonical .ebwt
			system("mv $map{$f} $f.$n.ebwt");
		}
	}
}
