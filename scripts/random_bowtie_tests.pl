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
	system("make bowtie-debug bowtie-build-debug bowtie-build-packed-debug bowtie-maptool-debug bowtie-inspect-debug") == 0 || die "Error building";
}

my @policies = (
	"-n 3",
	"-n 2",
	"-n 1",
	"-n 0",
	"-v 3",
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

sub nonACGTtoN {
	my $t = shift;
	$t =~ tr/-NnMmRrWwSsYyKkVvHhDdBbXx/N/;
	$t =~ /[ACGTN]+/ || die "Bad N-ized DNA string: $t";
	return $t;
}

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
	my $noAmbigs = int(rand(3)) == 0;
	for($i = 0; $i < $num; $i++) {
		my $or = int(rand(50));
		if($or == 0 && !$noAmbigs) {
			# Add a random, possibly ambiguous character
			$t .= $dnaMap[int(rand($#dnaMap+1))];
		} elsif($or == 1 && !$noAmbigs) {
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
		while(not $c =~ /[0-9A-Z\/=@%]/) {
			$c = chr(33 + int(rand(41)));
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
	my($t, $lineRate, $linesPerSide, $offRate, $isaRate, $ftabChars, $endian) = @_;
	my $ret = 0;
	
	my $file1 = "-c";
	my $file2 = "\"$t\"";
	if(substr($t, 0, 1) eq '-') {
		# Add backslash to escape first dash
		$file2 = "\"\\$t\"";
	}
	
	# Write reference sequences to a FASTA file
	open FA, ">.randtmp$seed.fa" || die "Could not open temporary fasta file";
	my @seqs = split(/,/, $t);
	for(my $i = 0; $i <= $#seqs; $i++) {
		print FA ">$i\n";
		print FA "$seqs[$i]\n";
	}
	close(FA);
	
	# Make a version of the FASTA file where all non-A/C/G/T characters
	# are Ns.
	open FAN, ">.randtmp$seed.ns.fa" || die "Could not open temporary fasta file";
	for(my $i = 0; $i <= $#seqs; $i++) {
		print FAN ">$i\n";
		my $t = nonACGTtoN($seqs[$i]);
		print FAN "$t\n";
	}
	close(FAN);
	
	my $fasta = int(rand(2)) == 0;
	if($fasta) {
		# Use the FASTA file as input
		$file1 = "-f";
		$file2 = ".randtmp$seed.fa";
	}
	
	my $bucketArg = "";
	my $bucketRand = int(rand(3));
	if($bucketRand == 0) {
		$bucketArg = "--bmaxdivn ";
		$bucketArg .= (int(rand(30))+1);
	} elsif($bucketRand == 1) {
		$bucketArg = "-a ";
	}
	
	my $bsearch_arg = "";
	if(int(rand(2)) == 0) {
		$bsearch_arg = "--oldpmap";
	}
	
	my $isaArg = "";
	if($isaRate >= 0) {
		$isaArg = "--isarate $isaRate";
	}

	my $args = "-q $isaArg --sanity $file1 --offrate $offRate --ftabchars $ftabChars $bsearch_arg $bucketArg $endian $file2";
	
	# Do unpacked version
	my $cmd = "./bowtie-build-debug $args .tmp$seed";
	system("echo \"$cmd\" > .tmp$seed.cmd");
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
	
	# Use bowtie-inspect to compare the output of bowtie-build to the
	# original reference sequences
	$cmd = "./bowtie-inspect-debug -a -1 .tmp$seed > .tmp$seed.inspect.ref";
	print "$cmd\n";
	system($cmd) == 0 || die "$cmd - failed";
	$cmd = "diff .randtmp$seed.ns.fa .tmp$seed.inspect.ref";
	print "$cmd\n";
	system($cmd) == 0 || die "$cmd - failed";

	# Do packed version and assert that it matches unpacked version
	# (sometimes, but not all the time because it takes a while)
	if(int(rand(4)) == 0) {
		$cmd = "./bowtie-build-packed-debug $args .tmp$seed.packed";
		print "$cmd\n";
		$out = trim(`$cmd 2>&1`);
		if($out eq "") {
			if(system("diff .tmp$seed.1.ebwt .tmp$seed.packed.1.ebwt") != 0) {
				if($exitOnFail) {
					exit 1;
				}
			} elsif(system("diff .tmp$seed.2.ebwt .tmp$seed.packed.2.ebwt") != 0) {
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
	my($t, $p, $policy, $oneHit, $requireResult, $offRate, $isaRate) = @_;
	
	my $patarg = "-c";
	my $patstr = "\"$p\"";
	my $outformat = int(rand(3));
	my $maptool_cmd = "";
	my $outfile = ".tmp$seed.out";
	if($outformat == 0) {
		$outformat = ""; # default (verbose)
		$maptool_cmd = " && ./bowtie-maptool-debug -d -C .tmp$seed.out"; # default
	} elsif($outformat == 1) {
		$outformat = "-b"; # binary
		$maptool_cmd = " && ./bowtie-maptool-debug -b -C .tmp$seed.out"; # binary
	} else {
		$outformat = "--concise"; # concise
	}
	my $format = int(rand(5));
	if($format == 0) {
		# FASTA
		open FA, ">.randread$seed.fa" || die "Could not open temporary fasta file";
		my @cs = split /[,]/, $p;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			print FA ">\n$cms[0]\n";
		}
		close FA;
		$patarg = "-f";
		$patstr = ".randread$seed.fa";
	} elsif($format == 1) {
		# FASTQ with ASCII qualities
		open FQ, ">.randread$seed.fq" || die "Could not open temporary fastq file";
		my @cs = split /[,]/, $p;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			print FQ "@\n$cms[0]\n+\n$cms[1]\n";
		}
		close FQ;
		$patarg = "-q";
		$patstr = ".randread$seed.fq";
	} elsif($format == 2) {
		# FASTQ with integer qualities
		open FQ, ">.randread$seed.integer.fq" || die "Could not open temporary solexa fastq file";
		my @cs = split /[,]/, $p;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			print FQ "@\n$cms[0]\n+\n";
			for(my $i = 0; $i < length($cms[1]); $i++) {
				my $q = substr($cms[1], $i, 1);
				$q = ord($q) - 33;
				print FQ "$q ";
			}
			print FQ "\n";
		}
		close FQ;
		$patarg = "-q --integer-quals";
		$patstr = ".randread$seed.integer.fq";
	} elsif($format == 3) {
		# Raw
		open RAW, ">.randread$seed.raw" || die "Could not open temporary raw file";
		my @cs = split /[,]/, $p;
		foreach my $c (@cs) {
			my @cms = split /[:]/, $c;
			print RAW "$cms[0]\n";
		}
		close RAW;
		$patarg = "-r";
		$patstr = ".randread$seed.raw";
	}
	
	my $isaArg = "";
	if($isaRate >= 0) {
		$isaArg = "--isarate $isaRate";
	}
	
	my $phased = "";
	if(int(rand(2)) == 0) {
		$phased = "-z";
	}
	
	my $khits = "-k 1";
	my $mhits = 0;
	if($phased eq "-z") {
		# A phased search may optionally be a non-stratified all-hits
		# search
		if(int(rand(2)) == 0) {
			$khits = "-a --nostrata";
		}
	} else {
		if(int(rand(2)) == 0) {
			$khits = "-a";
		} else {
			$khits = "-k " . (int(rand(20))+2);
		}
		if(int(rand(2)) == 0) {
			$requireResult = 0;
			$mhits = (int(rand(20))+2);
		}
		if(int(rand(2)) == 0) {
			$khits .= " --best";
		}
		if(int(rand(2)) == 0) {
			$khits .= " --nostrata";
		}
	}
	if($mhits > 0) {
		$khits .= " -m $mhits";
	}
	
	if($oneHit || 1) {
		$oneHit = "";
	} else {
		$oneHit = "-a";
	}
	my $offRateStr = "";
	if(int(rand(3)) == 0) {
		$offRateStr = "--offrate " . ($offRate + 1 + int(rand(4)));
	}
	my $cmd = "./bowtie-debug $policy $khits $outformat $isaArg $offRateStr --orig \"$t\" $phased $oneHit --sanity $patarg .tmp$seed $patstr $outfile $maptool_cmd";
	print "$cmd\n";
	my $out = trim(`$cmd 2>.tmp$seed.stderr`);
	
	# Bad exitlevel?
	if($? != 0) {
		print "Exitlevel: $?\n";
		if($exitOnFail) {
			my $err = trim(`cat .tmp$seed.stderr 2> /dev/null`);
			print "Stdout:\n$out\nStderr:\n$err\n";
			exit 1;
		}
		return 0;
	}
	my $err = trim(`cat .tmp$seed.stderr 2> /dev/null`);
	
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
	my %outhash = ();
	my %readcount = ();
	my $lastread = "";
	foreach(@outlines) {
		chomp;
		print "$_\n";
		!defined($outhash{$_}) || die "Result $_ appears in output twice";
		$outhash{$_} = 1;
		next if /^Reported/;
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
		/^([01-9]+)[+-]?[:]/;
		my $read = $1;
		if(($read ne $lastread) && ($phased eq "")) {
			die "Read $read appears multiple times non-consecutively" if defined($readcount{$read});
		}
		$lastread = $read;
		$readcount{$read}++ if defined($readcount{$read});
		$readcount{$read} = 1 unless defined($readcount{$read});
		if($mhits > 0) {
			$readcount{$read} <= $mhits || die "Read $read matched more than $mhits times";
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
	my $isaRate = int(rand(8))-1;        # Can be anything
	my $ftabChars = 1 + int(rand(8));    # Must be >= 1
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
	my @tts;
	for(my $i = 0; $i < $nt; $i++) {
		my $tlen = $tbase + int(rand($trand));
		my $tt = randDna($tlen);             # add text meat
		push(@tts, $tt);
		$t .= (randGap() . $tt . randGap()); # add random padding
		if($i < $nt-1) { $t .= ","; }        # add comma separator
	}
	
	# Run the command to build the Ebwt from the random text
	$pass += build($t, $lineRate, $linesPerSide, $offRate, $isaRate, $ftabChars, $endian);
	last if(++$tests > $limit);

	my $in = $inner;
	for(; $in >= 0; $in--) {
		# Generate random pattern(s) based on text
		my $pfinal = '';
		my $np = int(rand(30)) + 1;
		for(my $i = 0; $i < $np; $i++) {
			my $tt = $tts[int(rand($#tts))];
			my $pl = int(rand(length($tt))) - 10;
			$pl = max($pl, 4);
			$pl = min($pl, length($tt));
			my $plen = int(rand($prand)) + $pbase;
			my $pr = min($pl + $plen, length($tt));
			my $p = substr $tt, $pl, $pr - $pl;
			# Check for empty patter or pattern that spans a comma
			if(length($p) == 0 || index($p, ",") != -1) {
				$i--; next;
			}
			# Optionally add mutations to pattern (but not the first)
			if($i > 0) {
				my $nummms = int(rand(4));
				for(my $j = 0; $j < $nummms; $j++) {
					substr($p, int(rand(length($p))), 1, randDna(1));
				}
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
			last if($c eq ',');
			if($c ne 'A' && $c ne 'C' && $c ne 'G' && $c ne 'T') {
				$expectResult = 0;
				last;
			}
		}
		$pass += search($t, $pfinal, $policy, $oneHit, $expectResult, $offRate, $isaRate); # require 1 or more results
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
		$pass += search($t, $pfinal, $policy, $oneHit, 0, $offRate, $isaRate); # do not require any results
		last if(++$tests > $limit);
	}

	#system("rm -f .tmp.?.ebwt .tmp.packed.?.ebwt");
}

print "$pass tests passed, $fail failed\n";
exit 1 if $fail > 0;
