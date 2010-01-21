#!/usr/bin/perl -w

##
# cs_trim.pl
#
# Basic tests to ensure that colorspace primer trimming is working as
# expected.
#

use strict;
use warnings;

my $bowtie = "./bowtie";
if(system("$bowtie --version") != 0) {
	$bowtie = `which bowtie`;
	chomp($bowtie);
	if(system("$bowtie --version") != 0) {
		die "Could not find bowtie in current directory or in PATH\n";
	}
}

sub readToFastq {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	open TMP, ">$fname" || die "Could not open $fname for writing\n";
	print TMP "@r\n$r[0]\n+\n$r[1]\n";
	close(TMP);
}

sub readToFasta {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	open TMP, ">$fname" || die "Could not open $fname for writing\n";
	print TMP ">r\n$r[0]\n";
	close(TMP);
}

sub readToRaw {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	open TMP, ">$fname" || die "Could not open $fname for writing\n";
	print TMP "$r[0]\n";
	close(TMP);
}

sub readToTabbed {
	my ($rstr, $fname) = @_;
	my @r = split(/[:]/, $rstr);
	open TMP, ">$fname" || die "Could not open $fname for writing\n";
	print TMP "$r[0]\t$r[1]\n";
	close(TMP);
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

my $m1 = colorize("CTGACTGCAACGGGCAATATGTCTCTGTGTGGA", 0);
my $m2 = colorize("CATCACCATTACCACAGGTAACGGTGCGGGCTG", 0);

my $m1n = colorize("CTGACTGCAACGGGCAATATGTCTCTGTGTGGA", 1);
my $m2n = colorize("CATCACCATTACCACAGGTAACGGTGCGGGCTG", 1);

my @reads = (
	# Trim me
	"T0$m1:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
	# Trim me
	"A1$m2:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
	# Don't trim me
	"TT$m1:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
	# Don't trim me
	"CG$m2:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
	# Don't trim me
	"31$m1:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
	# Don't trim me
	"20$m2:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
	# Don't trim me
	"AC$m1n:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
	# Don't trim me
	"CC$m2n:".
	"ABCDEFGHIJKLMNOPQRSTUVWXYZZZZZZZZZ",
);

{
	my @r1 = split(/[:]/, $reads[0]);
	my @r2 = split(/[:]/, $reads[1]);
	my $cmd = "$bowtie --ff -C -c e_coli_c -1 $reads[0] -2 $reads[1]";
	print "$cmd\n";
	open BTIE, "$cmd |" || die;
	while(<BTIE>) {
		my @s = split;
		my ($seq, $quals) = ($s[4], $s[5]);
		my $sl = length($seq);
		my $el = length($m1)-1;
		$sl == $el || die "Expected seq length $el, got $sl\n";
		my $ql = length($quals);
		$ql == $el || die "Expected qual length $el, got $ql\n";
	}
	close(BTIE);
}

for(my $i = 2; $i <= $#reads; $i += 2) {
	my @r1 = split(/[:]/, $reads[$i]);
	my @r2 = split(/[:]/, $reads[$i+1]);
	open BTIE, "$bowtie --ff -C -c e_coli_c -1 $reads[$i] -2 $reads[$i+1] |" || die;
	while(<BTIE>) {
		my @s = split;
		my ($seq, $quals) = ($s[4], $s[5]);
		length($seq) == length($m1+1) || die;
		length($quals) == length($m1+1) || die;
	}
	close(BTIE);
	
	readToFastq($reads[$i], ".tmp1.fq");
	readToFastq($reads[$i+1], ".tmp2.fq");
	open BTIE, "$bowtie --ff -C -q e_coli_c -1 .tmp1.fq -2 .tmp2.fq |" || die;
	while(<BTIE>) {
		my @s = split;
		my ($seq, $quals) = ($s[4], $s[5]);
		length($seq) == length($m1+1) || die;
		length($quals) == length($m1+1) || die;
	}
	close(BTIE);

	readToFasta($reads[$i], ".tmp1.fa");
	readToFasta($reads[$i+1], ".tmp2.fa");
	open BTIE, "$bowtie --ff -C -f e_coli_c -1 .tmp1.fa -2 .tmp2.fa |" || die;
	while(<BTIE>) {
		my @s = split;
		my ($seq, $quals) = ($s[4], $s[5]);
		length($seq) == length($m1+1) || die;
		length($quals) == length($m1+1) || die;
	}
	close(BTIE);

	readToRaw($reads[$i], ".tmp1.raw");
	readToRaw($reads[$i+1], ".tmp2.raw");
	open BTIE, "$bowtie --ff -C -r e_coli_c -1 .tmp1.raw -2 .tmp2.raw |" || die;
	while(<BTIE>) {
		my @s = split;
		my ($seq, $quals) = ($s[4], $s[5]);
		length($seq) == length($m1+1) || die;
		length($quals) == length($m1+1) || die;
	}
	close(BTIE);

	readToTabbed($reads[$i], ".tmp1.tabbed");
	readToTabbed($reads[$i+1], ".tmp2.tabbed");
	open BTIE, "$bowtie --ff -C --12 e_coli_c -1 .tmp1.tabbed -2 .tmp2.tabbed |" || die;
	while(<BTIE>) {
		my @s = split;
		my ($seq, $quals) = ($s[4], $s[5]);
		length($seq) == length($m1+1) || die;
		length($quals) == length($m1+1) || die;
	}
	close(BTIE);
}
