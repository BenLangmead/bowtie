#!/usr/bin/perl -w

#
# Throw lots of random but interesting test cases at the extended
# Burrows-Wheeler transform builder.
#
# Usage: perl random_tester.pl [rand seed] [# outer iters] [# inner iters]
#

use List::Util qw[min max];
use Getopt::Std;

my %options=();
getopts("mhn",\%options);

if(defined $options{h}) {
	print "Usage: perl random_bowtie_tests.pl.pl seed outer inner tbase trand pbase prand\n";
	exit 0;
}

unless(defined $options{n}) {
	system("make bowtie-debug bowtie-build-debug bowtie-build-packed-debug") == 0 || die "Error building";
}

my @policies = (
	"-n 2",
	"-n 1",
	"-n 0",
	"-v 2",
	"-v 1",
	"-v 0"
);

sub pickPolicy {
	my $r = int(rand($#policies + 1));
	return $policies[$r];
}

my $seed = 0;
$seed = int $ARGV[0] if defined($ARGV[0]);
srand $seed;

my $outer = 5000;
$outer = int $ARGV[1] if defined($ARGV[1]);
my $limit = $outer;

my $inner = 100;
$inner = int $ARGV[2] if defined($ARGV[2]);

my $tbase = 100;
$tbase = int $ARGV[3] if defined($ARGV[3]);
my $trand = 300;
$trand = int $ARGV[4] if defined($ARGV[4]);

my $pbase = 10;
$pbase = int $ARGV[5] if defined($ARGV[5]);
my $prand = 30;
$prand = int $ARGV[6] if defined($ARGV[6]);

my $verbose = 0;
my $exitOnFail = 1;
my @dnaMap = ('A', 'T', 'C', 'G',
              'N',
              'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'X');

sub randGap() {
	my $or = int(rand(4));
	my $gap = "";
	if($or == 0) {
		my $ir = int(rand(100))+1;
		if(($ir & 1) == 1) {
			for(my $i = 0; $i < $ir; $i++) {
				$gap .= 'N';
			}
		} else {
			for(my $i = 0; $i < $ir; $i++) {
				$gap .= '-';
			}
		}
	}
	return $gap;
}

# Utility function generates a random DNA string of the given length
sub randDna($) {
	my $num = shift;
	my $i;
	my $t = '';
	for($i = 0; $i < $num; $i++) {
		my $or = int(rand(50));
		if($or == 0) {
			# Add a random, possibly ambiguous character
			$t .= $dnaMap[int(rand($#dnaMap+1))];
		} elsif($or == 1) {
			# Add a random-length streak of Ns (max: 20)
			my $streak = int(rand(20))+1;
			for(my $j = $i; $j < $num && $j < $streak; $j++) {
				$t .= 'N';
			}
		} else {
			# Add a random non-ambiguous character
			$t .= $dnaMap[int(rand(4))];
		}
	}
	return $t;
}

# Utility function that returns the reverse complement of its argument
sub reverseComp($) {
	my $r = shift;
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

# Add some random quality values to encourage excercising the
# backtracking code
sub addQual($) {
	my $r = shift;
	my $len = length($r);
	$r .= ":";
	for(my $i = 0; $i < $len; $i++) {
		my $c = "-";
		while(not $c =~ /[0-9A-Z]/) {
			$c = chr(33 + int(rand(40)));
		}
		$r .= $c;
	}
	return $r;
}

# Trim whitespace from a string argument
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

# Build an Ebwt based on given arguments
sub build {
	my($t, $lineRate, $linesPerSide, $offRate, $ftabChars, $chunkRate, $endian) = @_;
	my $ret = 0;
	
	my $file1 = "-c";
	my $file2 = "\"$t\"";
	if(substr($t, 0, 1) eq '-') {
		# Add backslash to escape first dash
		$file2 = "\"\\$t\"";
	}
	my $fasta = int(rand(2)) == 0;
	if($fasta) {
		open FA, ">.randtmp.fa" || die "Could not open temporary fasta file";
		my @seqs = split(/,/, $t);
		for(my $i = 0; $i <= $#seqs; $i++) {
			print FA ">ref$i\n";
			print FA "$seqs[$i]\n";
		}
		close(FA);
		$file1 = "-f";
		$file2 = ".randtmp.fa";
	}
	
	# Do unpacked version
	my $cmd = "./bowtie-build-debug -s $file1 --offrate $offRate --ftabchars $ftabChars --chunkrate $chunkRate $endian $file2 .tmp";
	print "$cmd\n";
	my $out = trim(`$cmd 2>&1`);
	if($out eq "") {
		$ret++;
	} else {
		print "$out\n";
		if($exitOnFail) {
			exit 1;
		}
	}

	# Do packed version and assert that it matches unpacked version
	# (sometimes, but not all the time because it takes a while)
	if(int(rand(4)) == 5) {
		$cmd = "./bowtie-build-packed-debug -s $file1 --offrate $offRate --ftabchars $ftabChars --chunkrate $chunkRate $endian $file2 .tmp.packed";
		print "$cmd\n";
		$out = trim(`$cmd 2>&1`);
		if($out eq "") {
			if(system("diff .tmp.1.ebwt .tmp.packed.1.ebwt") != 0) {
				if($exitOnFail) {
					exit 1;
				}
			} elsif(system("diff .tmp.2.ebwt .tmp.packed.2.ebwt") != 0) {
				if($exitOnFail) {
					exit 1;
				}
			} else {
				$ret++;
			}
		} else {
			print "$out\n";
			if($exitOnFail) {
				exit 1;
			}
		}
	}
	
	return $ret;
}

# Search for a pattern in an existing Ebwt
sub search {
	my($t, $p, $policy, $oneHit, $requireResult) = @_;
	if($oneHit || 1) {
		$oneHit = "";
	} else {
		$oneHit = "-a";
	}
	my $cmd = "./bowtie-debug $policy --concise --orig \"$t\" $oneHit -s -c .tmp \"$p\"";
	print "$cmd\n";
	my $out = trim(`$cmd 2>.tmp.stderr`);
	
	# Bad exitlevel?
	if($? != 0) {
		print "Exitlevel: $?\n";
		if($exitOnFail) {
			my $err = trim(`cat .tmp.stderr 2> /dev/null`);
			print "Stdout:\n$out\nStderr:\n$err\n";
			exit 1;
		}
		return 0;
	}
	my $err = trim(`cat .tmp.stderr 2> /dev/null`);
	
	# Yielded no results when we were expecting some?
	if($out eq "No results" && $requireResult) {
		print "Expected results but got \"No Results\"\n";
		if($exitOnFail) {
			print "Stdout:\n$out\nStderr:\n$err\n";
			exit 1;
		}
		return 0;
	} elsif($out eq "No results") {
		# Yielded no results, but that's OK
		return 1;
	}
	
	# No output?
	if($out eq "") {
		print "Expected some output but got none\n";
		exit 1 if($exitOnFail);
		return 0;
	}
	
	# Parse output to see if any of it is bad
	my @outlines = split('\n', $out);
	foreach(@outlines) {
		chomp;
		print "$_\n";
		# Result should look like "4+:<4,231,0>,<7,111,0>,<7,112,1>,<4,234,0>"
		unless(m/^[01-9]+[+-]?[:](?:<[01-9]+,[01-9]+,[01-9]+>[,]?)+\s*$/) {
			print "Results malformed\n";
			print "$out\n";
			if($exitOnFail) {
				print "Stdout:\n$out\nStderr:\n$err\n";
				exit 1;
			}
			return 0;
		}
	}
	
	# Success
	return 1;
}

my $pass = 0;
my $tests = 0;
my $fail = 0;

for(; $outer > 0; $outer--) {

	# Generate random parameters
	my $lineRate = 4 + int(rand(7));     # Must be >= 4
	my $linesPerSide = 1 + int(rand(3));
	my $offRate = int(rand(16));         # Can be anything
	my $ftabChars = 1 + int(rand(8));    # Must be >= 1
	my $chunkRate = 1 + int(rand(10));   # Must be >= 1
	my $big = int(rand(2));
	my $endian = '';
	if($big) {
		$endian = "--big";
	} else {
		$endian = "--little";
	}

	# Generate random text(s)
	my $nt = int(rand(10)) + 1;
	my $t = '';
	my $tt = '';
	for(my $i = 0; $i < $nt; $i++) {
		my $tlen = $tbase + int(rand($trand));
		$tt = randDna($tlen);                # add text meat
		$t .= (randGap() . $tt . randGap()); # add random padding
		if($i < $nt-1) { $t .= ","; }        # add comma separator
	}
	
	# Run the command to build the Ebwt from the random text
	$pass += build($t, $lineRate, $linesPerSide, $offRate, $ftabChars, $chunkRate, $endian);
	last if(++$tests > $limit);

	my $in = $inner;
	for(; $in >= 0; $in--) {
		# Generate random pattern(s) based on text
		my $pfinal = '';
		my $np = int(rand(10)) + 1;
		for(my $i = 0; $i < $np; $i++) {
			my $pl = int(rand(length($tt))) - 10;
			$pl = max($pl, 0);
			$pl = min($pl, length($tt));
			my $plen = int(rand($prand)) + $pbase;
			my $pr = min($pl + $plen, length($tt));
			my $p = substr $tt, $pl, $pr - $pl;
			if(length($p) == 0 || index($p, ",") != -1) {
				$i--; next;
			}
			if(0) {
				# Add some random padding on either side
				my $lpad = randDna(max(0, int(rand(20)) - 10));
				my $rpad = randDna(max(0, int(rand(20)) - 10));
				$p = $lpad . $p . $rpad;
			}
			if((int(rand(2)) == 0)) {
				$p = reverseComp($p);
			}
			$p = addQual($p);
			$pfinal .= $p;
			if($i < $np-1) {
				$pfinal .= ","
			}
		}
		
		# Run the command to search for the pattern from the Ebwt
		my $oneHit = (int(rand(3)) == 0);
		my $policy = pickPolicy();
		my $expectResult = 1;
		for(my $i = 0; $i < length($pfinal); $i++) {
			my $c = substr($pfinal, $i, 1);
			if($c ne 'A' && $c ne 'C' && $c ne 'G' && $c ne 'T') {
				$expectResult = 0;
				last;
			}
		}
		$pass += search($t, $pfinal, $policy, $oneHit, $expectResult); # require 1 or more results
		last if(++$tests > $limit);
	}

	$in = $inner;
	for(; $in >= 0; $in--) {
		# Generate random pattern *not* based on text
		my $pfinal = '';
		my $np = int(rand(10)) + 1;
		for(my $i = 0; $i < $np; $i++) {
			my $plen = int(rand($prand)) + $pbase;
			my $p = randDna($plen);
			if(int(rand(2)) == 0) {
				$p = reverseComp($p);
			}
			$p = addQual($p);
			$pfinal .= $p;
			if($i < $np-1) {
				$pfinal .= ","
			}
		}

		# Run the command to search for the pattern from the Ebwt
		my $oneHit = (int(rand(3)) == 0);
		my $policy = pickPolicy();
		$pass += search($t, $pfinal, $policy, $oneHit, 0); # do not require any results
		last if(++$tests > $limit);
	}

	#system("rm -f .tmp.?.ebwt .tmp.packed.?.ebwt");
}

print "$pass tests passed, $fail failed\n";
exit 1 if $fail > 0;
