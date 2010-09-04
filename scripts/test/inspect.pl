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
my $bowtie_build = "./bowtie-build";
my $bowtie_build2 = "./bowtie-build";
my $bowtie_inspect = "./bowtie-inspect";
my $bowtie_inspect2 = "./bowtie-inspect";
my $ref = "";
my $debug = 0;
my $debug_old = 0;

GetOptions (
	"ref:s"                => \$ref,
	"bowtie:s"             => \$bowtie,
	"bowtie-build:s"       => \$bowtie_build,
	"bowtie-build2:s"      => \$bowtie_build2,
	"bowtie-inspect:s"     => \$bowtie_inspect,
	"bowtie-inspect2:s"    => \$bowtie_inspect2,
	"debug"                => \$debug,
	"debug2"               => \$debug_old);

my $bowtie_d = "${bowtie}-debug";
my $bowtie_build_d = "${bowtie_build}-debug";
my $bowtie_build2_d = "${bowtie_build2}-debug";
my $bowtie_inspect_d = "${bowtie_inspect}-debug";
my $bowtie_inspect2_d = "${bowtie_inspect2}-debug";

my @cases = (
	">-\nNA\n",
	">-\nNNA\n",
	">-\nNAA\n",
	">\nN\n>0\nATCTAG\n>\nN\n",
	">\nN\n>0\nATCTAG\n>\nN\n",
	">a\nN\n>b\nATCTAG\n>c\nN\n>d\nATCTAG\n",
	">a\nN\n>b\nATCTAG\n>c\nN\n",
	">\nN\n>0\nATCTAG\n>\nN\n",
	">a\nN\n>a\nNNNATCTAGNN\n>b\nNNN\n>c\nNNNNNN\n>d\nNNNNNA\n",
	">a\nNA\n>zzz\nNNNATCTAGNN\n>b\nNNN\n>c\nNNNNNN\n>d\nNA\n",
	">a\nAN\n>a\nNNNATCTAGNNA\n>b\nNNN\n>a\nNNNNNN\n>d\nAN\n"
);

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
		substr($name, 0, 1) eq ">" || die "Unexpected name line:\n$name\nin fasta:\n$_[0]";
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
		defined($cmap{$di}) || die "Bad dinuc: $di\n";
		$ret .= ($nucs ? $nmap{$cmap{$di}} : $cmap{$di});
	}
	return $ret;
}

sub trim($) {
	my $ret = $_[0];
	$ret =~ s/^\s+//; $ret =~ s/\s+$//;
	return $ret;
}

##
#
#
sub match($$$) {
	if($_[0] ne $_[1]) {
		open(D1, ">.inspect.pl.d1") || die;
		open(D2, ">.inspect.pl.d2") || die;
		print D1 "$_[0]\n";
		print D2 "$_[1]\n";
		close(D1);
		close(D2);
		system("diff -uw .inspect.pl.d1 .inspect.pl.d2");
		die "$_[2]";
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
	open(REF, $ref) || die "Could not open -ref $ref";
	@cases = ();
	push @cases, "";
	while(<REF>) {
		chomp;
		unless(/^>/) { s/[^ACGT]/N/gi; }
		$cases[0] .= "$_\n";
	}
	close(REF);
}

my $fn = ".inspect.pl.tmp.fa";
for my $ca (@cases) {
	for my $col (0, 1) {
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
			my $cmdEnd = "$fn $fn >/dev/null && $bowtie_inspect -s --extra $fn | grep -v '^Sequence'";
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
		system($cmd) == 0 || die "Exitlevel $? from command '$cmd'\n";
		$cmd = "$bi $fn";
		print "$cmd\n";
		my $io = trim(`$cmd`);
		$? == 0 || die "Exitlevel $? from command '$cmd'\n";
		my $msg = "Output from bowtie-inspect:\n$io\ndidn't match input:\n";
		if(defined($e)) {
			my $e2 = trim(stripAllGaps($e, $col));
			match($io, $e2, "$msg$e2");
		} else {
			my $c2 = trim(stripAllGaps($c, $col));
			match($io, $c2, "$msg$c2");
		}
		print $io;
		$cmd = "$bi -e $fn";
		print "$cmd\n";
		$io = trim(`$cmd`);
		$? == 0 || die "Exitlevel $? from command '$cmd'\n";
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
		print $io;
		if($bowtie_inspect2 ne "") {
			my $bio = $debug_old ? $bowtie_inspect2_d : $bowtie_inspect2;
			$bio .= " -a -1";
			$cmd = "$bio $fn";
			print "$cmd\n";
			$io = trim(`$cmd`);
			$? == 0 || die "Exitlevel $? from command '$cmd'\n";
			$msg = "Output from bowtie-inspect:\n$io\ndidn't match input:\n";
			if(defined($e)) {
				my $e2 = trim(stripAllGaps($e, $col));
				match($io, $e2, "$msg$e2");
			} else {
				my $c2 = trim(stripAllGaps($c, $col));
				match($io, $c2, "$msg$c2");
			}
			print $io;
		}
		print "PASSED\n";
	}
}
