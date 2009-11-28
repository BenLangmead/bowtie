#!/usr/bin/perl -w

# Used to convert MANUAL.markdown to MANUAL.

use strict;
use warnings;

my $lastBlank = 0;

while(<>) {
	# Skip comments
	next if /^\s*<!--/;
	next if /^\s*!/;
	next if /^\s*-->/;
	# Skip internal links
	next if /\[.*\]: #/;
	# Skip HTML
	next if /^\s?\s?\s?<.*>\s*$/;
	# Strip [`...`]
	s/\[`/`/g;
	s/`\]/`/g;
	# Strip [#...]
	s/\[#[^\]]*\]//g;
	# Turn hashes into spaces
	s/^####/   /;
	s/^###/ /;
	if(/^\s*$/) {
		next if $lastBlank;
		$lastBlank = 1;
	} else {
		$lastBlank = 0;
	}
	print $_;
}
