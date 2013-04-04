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

my $bowtie_d = "./bowtie-debug";
if(system("$bowtie_d --version") != 0) {
	$bowtie_d = `which bowtie-debug`;
	chomp($bowtie_d);
	if(system("$bowtie_d --version") != 0) {
		die "Could not find bowtie-debug in current directory or in PATH\n";
	}
}

my $samtools = "samtools";
if(! -x $samtools) {
	$samtools = `which samtools`;
	chomp($samtools);
	if($samtools eq "" || ! -x $samtools) {
		die "Could not find samtools in current directory or in PATH\n";
	}
}

if(! -f "e_coli_c.1.ebwt") {
	print STDERR "Making colorspace e_coli index\n";
	my $bowtie_build = "./bowtie-build";
	if(system("$bowtie_build --version") != 0) {
		print STDERR "Could not execute ./bowtie-build; looking in PATH...\n";
		$bowtie_build = `which $bowtie_build`;
		chomp($bowtie_build);
		if(system("$bowtie_build --version") != 0) {
			die "Could not find bowtie-build in current directory or in PATH\n";
		}
	}
	system("$bowtie_build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

sub run($) {
	my $cmd = shift;
	print "$cmd\n";
	return system($cmd);
}

system("rm -f .samtools.pl.*");
run("$bowtie_d -S e_coli reads/e_coli_10000snp.fq .samtools.pl.sam") && die;
run("$samtools view -bS -o .samtools.pl.bam .samtools.pl.sam") && die;
run("$samtools sort .samtools.pl.bam .samtools.pl.sorted") && die;
open SAM, "$samtools pileup -cv -f genomes/NC_008253.fna .samtools.pl.sorted.bam |" || die;
my $snps = 0;
while(<SAM>) {
	print $_;
	$snps++;
}
close(SAM);
$? == 0 || die "samtools pileup quit with exitlevel $?\n";
$snps == 10 || die "Wrong number of SNPs output by samtools; expected 10, got $snps\n";
print "PASSED\n";
system("rm -f .samtools.pl.*");

run("$bowtie_d -S -C -f e_coli_c reads/e_coli_10000snp.csfasta .samtools.pl.sam") && die;
run("$samtools view -bS -o .samtools.pl.bam .samtools.pl.sam") && die;
run("$samtools sort .samtools.pl.bam .samtools.pl.sorted") && die;
open SAM, "$samtools pileup -cv -f genomes/NC_008253.fna .samtools.pl.sorted.bam |" || die;
$snps = 0;
while(<SAM>) {
	print $_;
	$snps++;
}
close(SAM);
$? == 0 || die "samtools pileup quit with exitlevel $?\n";
$snps == 10 || die "Wrong number of SNPs output by samtools; expected 10, got $snps\n";
print "PASSED\n";
