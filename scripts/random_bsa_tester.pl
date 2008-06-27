#!/usr/bin/perl -w

#
# Throw lots of random DNA strings and difference-cover/bmax parameters
# at the blockwise_sa and die if it dies. 
#
# Usage: perl random_bsa_tester.pl [rand seed] [# outer iters]
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
sub bsa {
	my($t, $d, $b) = @_;
	my $cmd = "./blockwise_sa -q -s -b $b -d $d $t";
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
	# Generate random text(s)
	my $nt = int(rand(2100)) + 1;
	my $t = randDna($nt);
	my $d = int rand(10) + 2; # min: 2, max: 11
	$d = 1 << $d;
	my $b = int rand($nt) + 1;
	
	# Run the command to build the Ebwt from the random text
	$pass += bsa($t, $d, $b);
}

print "$pass out of $outer tests passed\n";
exit 1 if $pass < $outer;
