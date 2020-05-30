##
# inspect.pl
#
# Basic tests to ensure bowtie-inspect is behaving as expected.
#
# OLD_BOWTIE=~/software/bowtie-0.12.5 && \
# perl scripts/test/inspect.pl \
#    -bowtie-build2=$OLD_BOWTIE/bowtie-build \
#    -bowtie-inspect2=$OLD_BOWTIE/bowtie-inspect
#

use strict;
use warnings;
use Getopt::Long;

my $bowtie = "./bowtie";
my $bowtie2 = "./bowtie";
my $bowtie_build = "./bowtie-build";
my $bowtie_build2 = "./bowtie-build";
my $bowtie_inspect = "./bowtie-inspect";
my $bowtie_inspect2 = "./bowtie-inspect";
my $ref = "";
my $seed = 0;
my $debug = 0;
my $debug_old = 0;
my $nodie = 0;

GetOptions (
	"ref:s"                => \$ref,
	"bowtie:s"             => \$bowtie,
	"bowtie2:s"            => \$bowtie2,
	"bowtie-build:s"       => \$bowtie_build,
	"bowtie-build2:s"      => \$bowtie_build2,
	"bowtie-inspect:s"     => \$bowtie_inspect,
	"bowtie-inspect2:s"    => \$bowtie_inspect2,
	"seed:i"               => \$seed,
	"nodie"                => \$nodie,
	"debug"                => \$debug,
	"debug2"               => \$debug_old);

my $bowtie_d = "${bowtie} --debug";
my $bowtie2_d = "${bowtie2} --debug";
my $bowtie_build_d = "${bowtie_build} --debug";
my $bowtie_build2_d = "${bowtie_build2} --debug";
my $bowtie_inspect_d = "${bowtie_inspect} --debug";
my $bowtie_inspect2_d = "${bowtie_inspect2} --debug";

my @cases = (
	">\nN\n>0\nATCTAG\n>\nN\n",
	">-\nNA\n",
	">-\nNNA\n",
	">-\nNAA\n",
	">\nN\n>0\nATCTAG\n>\nN\n",
	">a\nN\n>b\nATCTAG\n>c\nN\n>d\nATCTAG\n",
	">a\nN\n>b\nATCTAG\n>c\nN\n",
	">\nN\n>0\nATCTAG\n>\nN\n",
	">a\nN\n>a\nNNNATCTAGNN\n>b\nNNN\n>c\nNNNNNN\n>d\nNNNNNA\n",
	">a\nNA\n>zzz\nNNNATCTAGNN\n>b\nNNN\n>c\nNNNNNN\n>d\nNA\n",
	">a\nAN\n>a\nNNNATCTAGNNA\n>b\nNNN\n>a\nNNNNNN\n>d\nAN\n"
);

sub mydie($) {
	my $msg = shift;
	print STDERR "$msg\n";
	if($nodie) {
		print "Press enter to continue...";
		my $tmp = <STDIN>;
	} else {
		die;
	}
}

##
# Given a fasta file string, strip away all sequences that consist
# entirely of gaps.
#
sub stripAllGaps($$) {
	my ($fa, $col) = @_;
	my @l = split(/[\r\n]+/, $fa);
	my $ret = "";
	while(1) {
		my $name = shift @l;
		last unless defined($name);
		substr($name, 0, 1) eq ">" || mydie("Unexpected name line:\n$name\nin fasta:\n$_[0]");
		my $seq = shift @l;
		if($col) {
			if($seq =~ /[ACGT][ACGT]/i) {
				$ret .= "$name\n$seq\n";
			}
		} else {
			if($seq =~ /[ACGT]/i) {
				$ret .= "$name\n$seq\n";
			}
		}
	}
	return $ret;
}

##
# Given a string in nucleotide space, convert to colorspace.
#
sub colorize($$) {
	my ($s, $nucs) = @_;
	defined($s) || die;
	my %cmap = (
		"AA" => "0", "CC" => "0", "GG" => "0", "TT" => "0",
		"AC" => "1", "CA" => "1", "GT" => "1", "TG" => "1",
		"AG" => "2", "GA" => "2", "CT" => "2", "TC" => "2",
		"AT" => "3", "TA" => "3", "CG" => "3", "GC" => "3",
		"NA" => ".", "NC" => ".", "NG" => ".", "NT" => ".",
		"AN" => ".", "CN" => ".", "GN" => ".", "TN" => ".",
		"NN" => "."
	);
	my %nmap = ("0" => "A", "1" => "C", "2" => "G", "3" => "T", "." => "N");
	my $ret = "";
	for(my $i = 0; $i < length($s)-1; $i++) {
		my $di = uc substr($s, $i, 2);
		$di =~ tr/-NnMmRrWwSsYyKkVvHhDdBbXx/N/;
		defined($cmap{$di}) || mydie("Bad dinuc: $di\n");
		$ret .= ($nucs ? $nmap{$cmap{$di}} : $cmap{$di});
	}
	return $ret;
}

##
# Return version of argument with whitespace trimmed off either end.
#
sub trim($) {
	my $ret = $_[0];
	$ret =~ s/^\s+//; $ret =~ s/\s+$//;
	return $ret;
}

##
# Return version of argument with given character trimmed off either end.
#
sub trimc($$) {
	my $ret = $_[0];
	my $c = $_[1];
	$ret =~ s/^$c+//; $ret =~ s/$c+$//;
	return $ret;
}

##
# Given a fasta string of references, parse out the sequences and
# return them joined with commas.
#
sub refsToReads($$) {
	my ($fa, $col) = @_;
	my @l = split(/[\r\n]+/, $fa);
	my $fn = ".inspect.pl.refsToReads.fa";
	open(OUT, ">$fn") || die;
	while(1) {
		my $name = shift @l;
		last unless defined($name);
		my $seq = shift @l;
		$seq = trimc($seq, 'N');
		$seq = colorize($seq, 1) if $col;
		$seq =~ s/N.*$//i;
		$seq !~ /N/i || die;
		if(length($seq) > 0) {
			print OUT "$name\n$seq\n";
		}
	}
	close(OUT);
	return $fn;
}

##
# Given a fasta string or references and a hash ref, fill the hash ref
# with with reference sequences keyed by name.
#
sub refmap($$$) {
	my ($fa, $map, $col) = @_;
	my @l = split(/[\r\n]+/, $fa);
	my $namelessref = 0;
	while(1) {
		my $name = shift @l;
		last unless defined($name);
		$name = substr($name, 1);
		next if scalar(@l) == 0 || $l[0] =~ /^>/;
		my $seq = shift @l;
		next if !defined($seq) || length($seq) == 0;
		if($name eq "") {
			$name = $namelessref;
			$namelessref++;
		}
		my $trimseq = trimc($seq, 'N');
		$trimseq =~ s/N.*$//i;
		$trimseq !~ /N/i || die;
		if($col) {
			$trimseq = "" unless $trimseq =~ /[acgt][acgt]/i;
		}
		$map->{seq}->{$name}->{$seq}++ if $seq ne "";
		$map->{trimseq}->{$name}->{$trimseq}++ if $trimseq ne "";
		$map->{trimseq}->{$name}{str} = $trimseq if $trimseq ne "";
	}
}

##
# Given output from running bowtie (first arg) and the fasta string of
# references (second), make sure that every reference is represented
# at least once in the bowtie output.  (at least b/c some references
# might be substrings of others.)
#
sub reconcileAlsWithRefs($$$) {
	my ($als, $fa, $col) = @_;
	my %rm = ();
	refmap($fa, \%rm, $col);
	return if $col;
	# Make a map from ref strs to the alignments that mapped
	my @l = split(/[\r\n]+/, $als);
	my %hits = ();
	for my $al (@l) {
		my @ls = split(/\t/, $al);
		my ($rd, $fw, $ref, $off, $seq) = ($ls[0], $ls[1], $ls[2], $ls[3], $ls[4]);
		$hits{$ref} = $seq if ($rd eq $ref && $fw eq '+');
	}
	for my $i (keys %{$rm{trimseq}}) {
		defined($hits{$i}) || mydie("No hit for reference $i:\n$als");
		defined($rm{trimseq}{$i}{$hits{$i}}) ||
			mydie("Hit sequence:\n\"$hits{$i}\"\n doesn't match ref ".
			      "sequence:\n\"".$rm{trimseq}{$i}{str}."\"\n$als");
	}
}

##
# Check that two strings are the same and, if they're not, print the
# provided error message along with the output of a 'diff' between
# them.
#
sub match($$$) {
	if($_[0] ne $_[1]) {
		open(D1, ">.inspect.pl.$seed.d1") || die;
		open(D2, ">.inspect.pl.$seed.d2") || die;
		print D1 "$_[0]\n";
		print D2 "$_[1]\n";
		close(D1);
		close(D2);
		system("diff -uw .inspect.pl.$seed.d1 .inspect.pl.$seed.d2");
		mydie("$_[2]");
	}
}

##
# Given a fasta file string, strip away all sequences that consist
# entirely of gaps.
#
sub colorizeFasta($) {
	my @l = split(/[\r\n]+/, $_[0]);
	my $ret = "";
	while(1) {
		my $name = shift @l;
		last unless defined($name);
		substr($name, 0, 1) eq ">" || die;
		my $seq = shift @l;
		my $cseq = colorize($seq, 1);
		$ret .= "$name\n$cseq\n" if $cseq ne "";
	}
	return $ret;
}

if($ref ne "") {
	open(REF, $ref) || mydie("Could not open -ref $ref");
	@cases = ();
	push @cases, "";
	while(<REF>) {
		chomp;
		unless(/^>/) { s/[^ACGT]/N/gi; }
		$cases[0] .= "$_\n";
	}
	close(REF);
}

my $fn = ".inspect.pl.tmp.$seed.fa";
for my $ca (@cases) {
	for my $col (0,) {
		my ($c, $e) = split(/:/, $ca);
		if(stripAllGaps($c, $col) eq "") {
			print "Skipping test case because it had no unambiguous stretches\n";
			next;
		}
		open(TMP, ">$fn") || die;
		print TMP $c;
		close(TMP);
		my $bb = $debug ? $bowtie_build_d : $bowtie_build;
		$bb .= " -C" if $col;
		my $bi = $debug ? $bowtie_inspect_d : $bowtie_inspect;
		$bi .= " -a -1";
		if($bowtie_build2 ne "") {
			my $cmdEnd = "$fn $fn >/dev/null && $bowtie_inspect -s --extra $fn | grep -v '^Sequence' | grep -v 'refnames.size' | grep -v 'Reverse'";
			my $bbo = $debug_old ? $bowtie_build2_d : $bowtie_build2;
			$bbo .= " -C" if $col;
			my $cmd = "$bbo $cmdEnd";
			print "$cmd\n";
			my $bb2out = `$cmd`;
			print "Old bowtie-build:\n$bb2out\n";
			$cmd = "$bb $cmdEnd";
			print "$cmd\n";
			my $bbout = `$cmd`;
			print "New bowtie-build:\n$bbout\n";
			match($bbout, $bb2out, "RefRecords from new and old bowtie-build (above) didn't match");
		}
		my $cmd = "$bb $fn $fn >/dev/null";
		print "$cmd\n";
		system($cmd) == 0 || mydie("Exitlevel $? from command '$cmd'");
		$cmd = "$bi $fn";
		print "$cmd\n";
		my $io = trim(`$cmd`);
		$? == 0 || mydie("Exitlevel $? from command '$cmd'");
		my $msg = "Output from bowtie-inspect:\n$io\ndidn't match input:\n";
		if(defined($e)) {
			my $e2 = trim(stripAllGaps($e, $col));
			match($io, $e2, "$msg$e2");
		} else {
			my $c2 = trim(stripAllGaps($c, $col));
			match($io, $c2, "$msg$c2");
		}
		print "$io\n";
		$cmd = "$bi -e $fn";
		print "$cmd\n";
		$io = trim(`$cmd`);
		$? == 0 || mydie("Exitlevel $? from command '$cmd'");
		$msg = "Output from bowtie-inspect -e:\n$io\ndidn't match input:\n";
		if(defined($e)) {
			# Colorspace quandry: strip gaps then colorize, or vice versa?
			my $e2 = ($col ? colorizeFasta($e) : $e);
			$e2 = trim(stripAllGaps($e2, 0));
			match($io, $e2, "$msg$e2");
		} else {
			my $c2 = ($col ? colorizeFasta($c) : $c);
			$c2 = trim(stripAllGaps($c2, 0));
			match($io, $c2, "$msg$c2");
		}
		print "$io\n";
		if($bowtie_inspect2 ne "") {
			my $bio = $debug_old ? $bowtie_inspect2_d : $bowtie_inspect2;
			$bio .= " -a -1";
			$cmd = "$bio $fn";
			print "$cmd\n";
			$io = trim(`$cmd`);
			$? == 0 || mydie("Exitlevel $? from command '$cmd'");
			$msg = "Output from bowtie-inspect:\n$io\ndidn't match input:\n";
			if(defined($e)) {
				my $e2 = trim(stripAllGaps($e, $col));
				match($io, $e2, "$msg$e2");
			} else {
				my $c2 = trim(stripAllGaps($c, $col));
				match($io, $c2, "$msg$c2");
			}
			print "$io\n";
		}
		
		#
		# Now query the index using all of the reference strings as
		# queries and using '-a -v 0'.  This will test whether bowtie
		# (and potentially also an older version of bowtie) agrees with
		# bowtie-inspect and bowtie-build about what's in the index.
		#
		my $ba = $debug ? $bowtie_d : $bowtie;
		$ba .= " -C --col-keepends" if $col;
		$cmd = "$ba -f -a -v 0 $fn ".refsToReads($c, $col);
		print "$cmd\n";
		my $bao = trim(`$cmd`);
		reconcileAlsWithRefs($bao, $c, $col);
		if($bowtie2 ne "") {
			my $ba2 = $debug ? $bowtie2_d : $bowtie2;
			$ba2 .= " -C --col-keepends" if $col;
			$cmd = "$ba2 -f -a -v 0 $fn ".refsToReads($c, $col);
			print "$cmd\n";
			my $bao2 = trim(`$cmd`);
			reconcileAlsWithRefs($bao2, $c, $col);
		}
		
		print "PASSED\n";
	}
}
