#!/usr/bin/perl -w

#
# Throw lots of random test cases, including both suffix sorts and
# regular sorts, at the multikey quicksorter.
#
# Usage: perl random_qsort_tester.pl [rand seed] [# outer iters]
#

use List::Util qw[min max];

my $seed = 0;
$seed = int $ARGV[0] if defined($ARGV[0]);
srand $seed;

my $outer = 10000;
$outer = int $ARGV[1] if defined($ARGV[1]);

my $verbose = 0;
my $exitOnFail = 1;
my @dnaMap = ('A', 'T', 'C', 'G');

# Utility function generates a random DNA string of the given length
sub randDna($) {
	my $num = shift;
	my $i;
	my $t = '';
	for($i = 0; $i < $num; $i++) {
		$t .= $dnaMap[int(rand(4))];
	}
	return $t;
}

# Trim whitespace from a string argument
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

# Build an Ebwt based on given arguments
sub qsort {
	my($t, $suf, $d) = @_;
	if($suf) {
		$suf = "-x";
		if($d > 0) {
			$d = 1 << $d;
			$suf .= " -d $d";
		}
	} else {
		$suf = "";
	}
	my $cmd = "./multikey_qsort -s -q -c $suf $t";
	print "$cmd\n";
	my $out = trim(`$cmd 2>&1`);
	if($out eq "") {
		return 1;
	} else {
		print "$out\n";
		if($exitOnFail) {
			exit 1;
		}
		return 0;
	}
}

my $pass = 0;

for(my $o = 0; $o < $outer; $o++) {
	my $suf = ((int rand(2)) == 0);
	my $t = '';
	my $d = 0;
	if(!$suf) {
		# Generate random texts and sort them
		my $nt = int(rand(100)) + 1;
		for(my $i = 0; $i < $nt; $i++) {
			my $tlen = 1 + int(rand(30));
			if($tlen > 15) {
				$tlen = 1 + int(rand(60));
				if($tlen > 30) {
					$tlen = 1 + int(rand(120));
					if($tlen > 60) {
						$tlen = 1 + int(rand(10000));
					}
				}
			}
			$t .= randDna($tlen);
			if($i < $nt-1) {
				$t .= ",";
			}
		}
	} else {
		# Generate a random text and sort a sample of its suffixes
		my $nt = int(rand(1009)) + 1;
		$t = randDna($nt);
		# Don't pick a v (thereby disabling difference cover) 50% of the time
		if(int(rand(2))) {
			if(int(rand(2))) {
				# Pick larger v's < 50% of the time
				$d = int(rand(10)) + 2;
			} else {
				# Pick a small v at least 50% of the time
				$d = int(rand(3)) + 2;
			}
		}
	}
	
	# Run the command to build the Ebwt from the random text
	$pass += qsort($t, $suf, $d);
}

print "$pass out of $outer tests passed\n";
exit 1 if $pass < $outer;
