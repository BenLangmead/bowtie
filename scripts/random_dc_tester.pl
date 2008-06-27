#!/usr/bin/perl -w

#
# Throw lots of random test cases at the difference-sample builder.
#
# Usage: perl random_dc_tester.pl [rand seed] [# outer iters] [textcap] [vcap]
#

use List::Util qw[min max];

my $seed = 0;
$seed = int $ARGV[0] if defined($ARGV[0]);
srand $seed;

my $outer = 10000;
$outer = int $ARGV[1] if defined($ARGV[1]);

my $tcap = 4100;
$tcap = int $ARGV[2] if defined($ARGV[2]);

my $vcap = 12;
$vcap = int $ARGV[3] if defined($ARGV[3]);

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
sub dc {
	my($v, $t) = @_;
	my $cmd = "./diff_sample-with-asserts -c -q -s $v $v $t";
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
	
	# No output?
	if($out ne "") {
		print "Expected no output but got some:\n$out\n";
		exit 1 if($exitOnFail);
		return 0;
	}
	
	# Success
	return 1;
}

my $pass = 0;
my $fail = 0;

for(; $outer > 0; $outer--) {

	# Generate random parameters
	my $tlen = int(rand($tcap))+1;
	my $t = randDna($tlen);

	# Generate random parameters
	my $v = int(rand($vcap-2+1))+2;
	$v = 1 << $v;
	
	# Run the command to build the Ebwt from the random text
	$pass += dc($v, $t);
}

print "$pass tests passed, $fail failed\n";
exit 1 if $fail > 0;
