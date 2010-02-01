#!/usr/bin/perl -w

##
# cs_trim.pl
#
# Basic tests to ensure that colorspace primer trimming is working as
# expected.
#

use strict;
use warnings;

my $debug = "-debug";
my $bowtie = "./bowtie$debug";
if(system("$bowtie --version") != 0) {
	$bowtie = `which bowtie$debug`;
	chomp($bowtie);
	if(system("$bowtie --version") != 0) {
		die "Could not find bowtie in current directory or in PATH\n";
	}
}

if(! -f "e_coli_c.1.ebwt") {
	print STDERR "Making colorspace e_coli index\n";
	system("make bowtie-build") && die;
	system("bowtie-build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

sub readToFastq {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPFQ, ">$fname" || die "Could not open $fname for writing\n";
	print TMPFQ "\@r\n$r[0]\n+\n$r[1]\n";
	close(TMPFQ);
}

sub readToFasta {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPFA, ">$fname" || die "Could not open $fname for writing\n";
	print TMPFA ">r\n$r[0]\n";
	close(TMPFA);
}

sub readToRaw {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPR, ">$fname" || die "Could not open $fname for writing\n";
	print TMPR "$r[0]\n";
	close(TMPR);
}

sub readToTabbed {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	system("rm -f $fname");
	open TMPT, ">$fname" || die "Could not open $fname for writing\n";
	print TMPT "r\t$r[0]\t$r[1]\n";
	close(TMPT);
}

sub readToTabbed2 {
	my ($rstr1, $rstr2, $fname) = @_;
	my @r1 = split(/[:]/, $rstr1);
	my @r2 = split(/[:]/, $rstr2);
	system("rm -f $fname");
	open TMPT, ">$fname" || die "Could not open $fname for writing\n";
	print TMPT "r\t$r1[0]\t$r1[1]\t$r2[0]\t$r2[1]\n";
	close(TMPT);
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
		$di =~ s/[URYMKWSBDHV]/N/gi;
		defined($cmap{$di}) || die "Bad dinuc: $di\n";
		$ret .= ($nucs ? $nmap{$cmap{$di}} : $cmap{$di});
	}
	return $ret;
}

# Utility function that returns the reverse complement of its argument
sub reverseComp($) {
	my $r = shift;
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

sub btrun {
	my ($name, $args, $num) = @_;
	my $cmd = "$bowtie $args";
	print "$cmd\n";
	open BTIE, "$cmd |" || die;
	while(<BTIE>) {
		my @s = split;
		my ($seq, $quals) = ($s[4], $s[5]);
		my $sl = length($seq);
		my $el = $num;
		$sl == $el || die "Expected seq length $el, got $sl\n";
		my $ql = length($quals);
		$ql == $el || die "Expected qual length $el, got $ql\n";
	}
	close(BTIE);
	$? == 0 || die "bowtie returned non-zero status $?\n";
	print "PASSED $name\n";
}

my $m1 = colorize("CTGACTGCAACGGGCAATATGTCTCTGTGTGGA", 0);
my $m2 = colorize("CATCACCATTACCACAGGTAACGGTGCGGGCTG", 0);

my $m1n = colorize("CTGACTGCAACGGGCAATATGTCTCTGTGTGGA", 1);
my $m2n = colorize("CATCACCATTACCACAGGTAACGGTGCGGGCTG", 1);

my $n1 = "CTGACTGCAACGGGCAATATGTCTCTGTGTGGA";
my $n2 = reverseComp("CATCACCATTACCACAGGTAACGGTGCGGGCTG");

my $q = "ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ";

my @reads = (
	# Trim me
	"T0$m1:$q",
	# Trim me
	"A1$m2:$q",
	# Don't trim me
	"$n1:".(substr($q, 0, -2)),
	# Don't trim me
	"$n2:".(substr($q, 0, -2)),
	# Don't trim me
	"TT$m1:$q",
	# Don't trim me
	"CG$m2:$q",
	# Don't trim me
	"31$m1:$q",
	# Don't trim me
	"20$m2:$q",
	# Don't trim me
	"AC$m1n:$q",
	# Don't trim me
	"CC$m2n:$q"
);

btrun("trim/-c/paired",   "-C -c e_coli_c -1 $reads[0] -2 $reads[1]", length($m1)-1);
btrun("trim/-c/unpaired", "-C -c e_coli_c $reads[0],$reads[1]", length($m1)-1);

readToFastq($reads[0], ".tmp1.fq");
readToFastq($reads[1], ".tmp2.fq");
btrun("no-trim/fastq/paired",   "-C e_coli_c -q -1 .tmp1.fq -2 .tmp2.fq", length($m1)-1);
btrun("no-trim/fastq/unpaired", "-C e_coli_c -q .tmp1.fq,.tmp2.fq", length($m1)-1);
btrun("trim5/fastq/unpaired", "-5 3 -C e_coli_c -q .tmp1.fq,.tmp2.fq", length($m1)-1-3);
btrun("trim35/fastq/unpaired", "-5 5 -3 3 -C e_coli_c -q .tmp1.fq,.tmp2.fq", length($m1)-1-8);

readToFasta($reads[0], ".tmp1.fa");
readToFasta($reads[1], ".tmp2.fa");
btrun("no-trim/fasta/paired",   "-C e_coli_c -f -1 .tmp1.fa -2 .tmp2.fa", length($m1)-1);
btrun("no-trim/fasta/unpaired", "-C e_coli_c -f .tmp1.fa,.tmp2.fa", length($m1)-1);
btrun("trim5/fasta/unpaired", "-5 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa", length($m1)-1-3);
btrun("trim35/fasta/unpaired", "-5 5 -3 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa", length($m1)-1-8);

readToRaw($reads[0], ".tmp1.raw");
readToRaw($reads[1], ".tmp2.raw");
btrun("no-trim/raw/paired",   "-C e_coli_c -r -1 .tmp1.raw -2 .tmp2.raw", length($m1)-1);
btrun("no-trim/raw/unpaired", "-C e_coli_c -r .tmp1.raw,.tmp2.raw", length($m1)-1);
btrun("trim5/raw/unpaired", "-5 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw", length($m1)-1-3);
btrun("trim35/raw/unpaired", "-5 5 -3 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw", length($m1)-1-8);

btrun("trim/-c/paired",   "-c e_coli -1 $reads[2] -2 $reads[3]", length($n1));
btrun("trim/-c/unpaired", "-c e_coli $reads[2],$reads[3]", length($n1));

for(my $i = 4; $i <= $#reads; $i += 2) {
	btrun("no-trim/-c/paired",   "-C -c e_coli_c -1 $reads[$i] -2 $reads[$i+1]", length($m1)+1);
	btrun("no-trim/-c/unpaired", "-C -c e_coli_c $reads[$i],$reads[$i+1]", length($m1)+1);
	btrun("trim5/-c/unpaired", "-5 3 -C -c e_coli_c $reads[$i],$reads[$i+1]", length($m1)+1-3);
	btrun("trim35/-c/unpaired", "-5 5 -3 3 -C -c e_coli_c $reads[$i],$reads[$i+1]", length($m1)+1-8);
	
	readToFastq($reads[$i], ".tmp1.fq");
	readToFastq($reads[$i+1], ".tmp2.fq");
	btrun("no-trim/fastq/paired",   "-C e_coli_c -q -1 .tmp1.fq -2 .tmp2.fq", length($m1)+1);
	btrun("no-trim/fastq/unpaired", "-C e_coli_c -q .tmp1.fq,.tmp2.fq", length($m1)+1);
	btrun("trim5/fastq/unpaired", "-5 3 -C e_coli_c -q .tmp1.fq,.tmp2.fq", length($m1)+1-3);
	btrun("trim35/fastq/unpaired", "-5 5 -3 3 -C e_coli_c -q .tmp1.fq,.tmp2.fq", length($m1)+1-8);

	readToFasta($reads[$i], ".tmp1.fa");
	readToFasta($reads[$i+1], ".tmp2.fa");
	btrun("no-trim/fasta/paired",   "-C e_coli_c -f -1 .tmp1.fa -2 .tmp2.fa", length($m1)+1);
	btrun("no-trim/fasta/unpaired", "-C e_coli_c -f .tmp1.fa,.tmp2.fa", length($m1)+1);
	btrun("trim5/fasta/unpaired", "-5 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa", length($m1)+1-3);
	btrun("trim35/fasta/unpaired", "-5 5 -3 3 -C e_coli_c -f .tmp1.fa,.tmp2.fa", length($m1)+1-8);

	readToRaw($reads[$i], ".tmp1.raw");
	readToRaw($reads[$i+1], ".tmp2.raw");
	btrun("no-trim/raw/paired",   "-C e_coli_c -r -1 .tmp1.raw -2 .tmp2.raw", length($m1)+1);
	btrun("no-trim/raw/unpaired", "-C e_coli_c -r .tmp1.raw,.tmp2.raw", length($m1)+1);
	btrun("trim5/raw/unpaired", "-5 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw", length($m1)+1-3);
	btrun("trim35/raw/unpaired", "-5 5 -3 3 -C e_coli_c -r .tmp1.raw,.tmp2.raw", length($m1)+1-8);

	readToTabbed($reads[$i], ".tmp1.tabbed");
	readToTabbed($reads[$i+1], ".tmp2.tabbed");
	readToTabbed2($reads[$i], $reads[$i+1], ".tmp.tabbed");
	btrun("no-trim/tabbed/paired",   "-C e_coli_c --12 .tmp.tabbed", length($m1)+1);
	btrun("no-trim/tabbed/unpaired", "-C e_coli_c --12 .tmp1.tabbed,.tmp2.tabbed", length($m1)+1);
	btrun("trim5/tabbed/unpaired", "-5 3 -C e_coli_c --12 .tmp1.tabbed,.tmp2.tabbed", length($m1)+1-3);
	btrun("trim35/tabbed/unpaired", "-5 5 -3 3 -C e_coli_c --12 .tmp1.tabbed,.tmp2.tabbed", length($m1)+1-8);
}
