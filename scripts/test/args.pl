#!/usr/bin/perl -w

##
# args.pl
#
# Basic tests to ensure that bad combinations of arguments are rejected
# and good ones are accepted.
#

my $debug = "-debug";
my $bowtie = "./bowtie$debug";
if(system("$bowtie --version") != 0) {
	$bowtie = `which bowtie$debug`;
	chomp($bowtie);
	if(system("$bowtie --version") != 0) {
		die "Could not find bowtie in current directory or in PATH\n";
	}
}

my @bad = (
	"-n 6".
	"-v 6",
	"-n 4",
	"-v 4",
	"-v 2 -n 4",
	"-v -1",
	"-n -10",
	"-3 -3",
	"-5 -1",
	"-e -1",
	"-l 4",
	"-l 0"
);

my @good = (
	"-n 0",
	"-n 1",
	"-n 2",
	"-n 3",
	"-v 0",
	"-v 1",
	"-v 2",
	"-v 3",
	"-v 3 -n 3"
);

for my $a (@bad) {
	system("$bowtie $a e_coli reads/e_coli_1000.fq /dev/null") != 0 || die;
	print "PASSED: bad args \"$a\"\n";
}
for my $a (@good) {
	system("$bowtie $a e_coli reads/e_coli_1000.fq /dev/null") == 0 || die;
	print "PASSED: good args \"$a\"\n";
}
