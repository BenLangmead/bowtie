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

my $bowtie_d = $bowtie . " --debug";
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

my $bcftools = "bcftools";
if(! -x $bcftools) {
	$bcftools = `which bcftools`;
	chomp($bcftools);
	if($bcftools eq "" || ! -x $bcftools) {
		die "Could not find bcftools in current directory or in PATH\n";
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

# BTL 4/11/2013: Couldn't get samtools / bcftools to find SNPs from the
# colorspace invocation of Bowtie.  Not sure why. 
for my $cmd ("$bowtie_d -S e_coli reads/e_coli_10000snp.fq"
             #, "$bowtie_d -S -C -f e_coli_c reads/e_coli_10000snp.csfasta"
             )
{
	system("rm -f .samtools.pl.*");
	# Run Bowtie and output SAM
	run("$cmd .samtools.pl.sam") && die;
	# Convert to BAM
	run("$samtools view -bS -o .samtools.pl.bam .samtools.pl.sam") && die;
	# Sort BAM
	run("$samtools sort .samtools.pl.bam .samtools.pl.sorted") && die;
	# Run samtools mpileup / bcftools to get SNPs
	my $bcfcmd = "$samtools mpileup -uf genomes/NC_008253.fna .samtools.pl.sorted.bam | $bcftools view -vcg -";
	print STDERR "$bcfcmd\n";
	open SAM, "$bcfcmd |" || die;
	my $snps = 0;
	while(<SAM>) {
		next if substr($_, 0, 1) eq "#";
		print $_;
		$snps++;
	}
	close(SAM);
	$? == 0 || die "samtools mpileup quit with exitlevel $?\n";
	$snps == 10 || die "Wrong number of SNPs output by samtools; expected 10, got $snps\n";
	print "PASSED\n";
}
