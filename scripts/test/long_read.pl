#!/usr/bin/perl -w

##
# long_read.pl
#
# Basic tests to ensure that long reads are handled properly.
#

use strict;
use warnings;

my $bowtie   = "./bowtie";
my $bowtie_d = "./bowtie-debug";
if(system("$bowtie --version") != 0) {
	$bowtie = `which bowtie`;
	chomp($bowtie);
	if(system("$bowtie --version") != 0) {
		die "Could not find bowtie in current directory or in PATH\n";
	}
}
if(system("$bowtie_d --version") != 0) {
	$bowtie_d = `which bowtie-debug`;
	chomp($bowtie_d);
	if(system("$bowtie_d --version") != 0) {
		die "Could not find bowtie-debug in current directory or in PATH\n";
	}
}

if(! -f "e_coli_c.1.ebwt") {
	print STDERR "Making colorspace e_coli index\n";
	system("make bowtie-build") && die;
	system("bowtie-build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

my @reads = (
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG".
	"CAACGGTAACCGGTAGCGTGGAGTCATGCCCCTGGACCGTCATGCCCAGTGAATCGCCTACCAGCATGAC:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDDDDDDDDDDDDDDDAAADDDDDDD;;;;DDDEFFFDEGGGGDDDDDDDDDCCCDDDDDD",
	""
);

sub btrun {
	my ($read, $color) = @_;
	my $args = $color ? "-C" : "";
	my $cmd = "$bowtie $args -c e_coli $read";
	print "$cmd\n";
	system($cmd) && die;
	$cmd = "$bowtie_d $args -c e_coli $read";
	print "$cmd\n";
	system($cmd) && die;
}

for my $r (@reads) {
	btrun($r, 0);
}
