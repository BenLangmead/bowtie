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
	system("make bowtie-build") && die;
	system("bowtie-build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

sub run($) {
	my $cmd = shift;
	print "$cmd\n";
	return system($cmd);
}

run("$bowtie_d -S e_coli reads/e_coli_10000snp.fq .samtools.pl.sam") && die;
run("$samtools view -bS -o .samtools.pl.bam .samtools.pl.sam") && die;
run("$samtools sort .samtools.pl.bam .samtools.pl.sorted") && die;
run("$samtools pileup -cv -f genomes/NC_008253.fna .samtools.pl.sorted.bam") && die;
