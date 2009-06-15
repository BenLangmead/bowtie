#!/usr/bin/perl -w

# 
# reconcile_alignments.pl
#
#  Author: Ben Langmead
#    Date: 6/14/2009
#
# Reconcile and sanity-check an input read file, an output hits file
# (the default, verbose kind), and an output unaligned-read file
# (generated with the --un option) to be sure they're consistent with
# each other.  If the --max option was used during alignment, then the
# --max file should also be specified.  If the --al option was also
# used, specify the --al file as the last argument to include it in the
# sanity check.
#
# Usage: perl reconcile_alignments.pl \
#            [-k <int>] [-a] [-m <int>] \
#            [-u <int>] [-f|-q|-r] \
#            <input read file> \
#            <hits file> \
#            <--un file> \
#            <--max file> \
#            <--al file>
#

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("rfqak:m:u:",\%options);

my $khits = 1;
$khits = int($options{k}) if defined($options{k});
$khits = 999999 if $options{a};
my $maxhits = 999999;
$maxhits = int($options{m}) if defined($options{m});
my $num_reads = -1;
$num_reads = $options{u} if defined($options{u});

my $format = "fastq";
$format = "fasta" if defined($options{f});
$format = "raw" if defined($options{r});

# Utility function that returns the reverse complement of its argument
sub reverseComp($) {
	my $r = shift;
	$r = reverse($r);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	return $r;
}

defined($ARGV[0]) || die "Must specify input read file as first arg";
defined($ARGV[1]) || die "Must specify alignment file as second arg";
defined($ARGV[2]) || die "Must specify unaligned-read file as third arg";

my $read_files  = $ARGV[0];
my $algn_file   = $ARGV[1];
my $un_file     = $ARGV[2];
my $max_file    = "";
$max_file = $ARGV[3] if defined($ARGV[3]);
my $al_file     = "";
$al_file = $ARGV[4] if defined($ARGV[4]);

my %hits_hash = ();
my %al_hash = ();
my %un_hash = ();
my %max_hash = ();

# Go through Bowtie-produced alignments
my $hits = 0;
my $distinctHits = 0;
open(ALGN, $algn_file);
{
	my %num_hash = (); # for checking -k/-m
	while(<ALGN>) {
		my @s = split /\t/;
		my $name = $s[0];
		my $add_read = 1;
		$hits++;
		if($khits > 1) {
			if(defined($num_hash{$name})) {
				$num_hash{$name}++;
				$add_read = 0; # already added this one
			} else {
				$num_hash{$name} = 1;
			$distinctHits++;
			}
			if($num_hash{$name} > $khits) {
				die "Number of alignments for read $name exceeds -k limit of $khits";
			} elsif($num_hash{$name} > $maxhits) {
				die "Number of alignments for read $name exceeds -m limit of $maxhits";
			}
		} else {
			defined($hits_hash{$name}{seq}) && 
				die "Read w/ name $name appears twice in alignment file";
			$distinctHits++;
		}
		if($add_read) {
			my $fw = ($s[1] eq "+");
			my $seq = $s[4];
			my $qual = $s[5];
			if(!$fw) {
				# Reverse / reverse-comp
				$seq = reverseComp($seq);
				$qual = reverse $qual;
			}
			$hits_hash{$name}{seq} = $seq;
			$hits_hash{$name}{qual} = (($format eq "fasta")? "" : $qual);
		}
	}
}
close(ALGN);

sub get_read($) {
	my $fh = shift;
	my ($name, $seq, $qual) = ("", "", "");
	if($format eq "fastq") {
		$name = <$fh>;
		return ("", "", "") unless defined($name);
		chomp($name);
		$name = substr($name, 1);
		$seq = <$fh>;
		chomp($seq);
		my $tmp = <$fh>;
		$qual = <$fh>;
		chomp($qual);
	} elsif($format eq "fasta") {
		$name = <$fh>;
		return ("", "", "") unless defined($name);
		chomp($name);
		$name = substr($name, 1);
		$seq = <$fh>;
		chomp($seq);
	} else {
		$format eq "raw" || die;
		die "Raw format not supported in reconcile_alignment.pl; read names required";
	}
	return ($name, $seq, $qual);
}

# Go through entries of the FASTQ file for the unaligned reads
my $uns = 0;
my $UN;
open $UN, $un_file;
while(1) {
	my ($name, $seq, $qual) = get_read($UN);
	last if $name eq "";
	$uns++;
	defined($hits_hash{$name}) &&
		die "Read $name appears both in hits file $algn_file and in --un file $un_file";
	defined($un_hash{$name}) &&
		die "Read $name appears more than once in --un file $un_file";
	$un_hash{$name}{seq} = $seq;
	$un_hash{$name}{qual} = $qual;
}
close($UN);

my $maxs = 0;
if($max_file ne "") {
	my $MAX;
	open $MAX, $max_file;
	# Go through entries of the MAX file for the unaligned reads
	while(1) {
		my ($name, $seq, $qual) = get_read($MAX);
		last if $name eq "";
		$maxs++;
		defined($hits_hash{$name}) &&
			die "Read $name appears both in hits file $algn_file and in --max file $max_file";
		defined($un_hash{$name})   &&
			die "Read $name appears both in --un file $un_file and in --max file $max_file";
		defined($max_hash{$name})  &&
			die "Read $name appears in --max file $max_file more than once";
		$max_hash{$name}{seq} = $seq;
		$max_hash{$name}{qual} = $qual;
	}
	close($MAX);
}

my $als = 0;
if($al_file ne "") {
	my $AL;
	open $AL, $al_file;
	# Go through entries of the MAX file for the unaligned reads
	while(1) {
		my ($name, $seq, $qual) = get_read($AL);
		last if $name eq "";
		$als++;
		defined($hits_hash{$name}) ||
			die "Read $name appears --al file $al_file but not in hits file $algn_file";
		defined($un_hash{$name})   &&
			die "Read $name appears both in --un file $un_file and in --al file $al_file";
		defined($max_hash{$name})  &&
			die "Read $name appears both in --max file $max_file and in --al file $al_file";
		defined($al_hash{$name})  &&
			die "Read $name appears in --al file $al_file more than once";
		$al_hash{$name}{seq} = $seq;
		$al_hash{$name}{qual} = $qual;
	}
	close($AL);
}

my @read_list = split(/,/, $read_files);
my $reads = 0;
for my $read_file (@read_list) {
	# Go through entries of the FASTQ file for the input reads and make
	# sure that each entry is mirrored by an entry either in the alignment
	# file or in the unaligned-read file.
	my $READ;
	open $READ, $read_file;
	my $patid = 0;
	while(1) {
		my ($name, $seq, $qual) = get_read($READ);
		last if $name eq "";
		$reads++;
		if(defined($hits_hash{$name})) {
			$hits_hash{$name}{seq} eq $seq ||
				die "Read $name in hits file $algn_file has different sequence ".
				    "from input read.\nHit: $hits_hash{$name}{seq}\nInput: $seq";
			# Qualities can be legitimately different
			#$hits_hash{$name}{qual} eq $qual ||
			#	die "Read $name in hits file $algn_file has different sequence ".
			#	    "from input read.\nHit: $hits_hash{$name}{qual}\nInput: $qual";
		}
		elsif(defined($un_hash{$name})) {
			$un_hash{$name}{seq} eq $seq ||
				die "Read $name in --un file $un_file has different sequence ".
				    "from input read.\nHit: $un_hash{$name}{seq}\nInput: $seq";
			$un_hash{$name}{qual} eq $qual ||
				die "Read $name in --un file $un_file has different sequence ".
				    "from input read.\nHit: $un_hash{$name}{qual}\nInput: $qual";
		}
		elsif(defined($max_hash{$name})) {
			$max_hash{$name}{seq} eq $seq ||
				die "Read $name in --max file $max_file has different sequence ".
				    "from input read.\nHit: $max_hash{$name}{seq}\nInput: $seq";
			$max_hash{$name}{qual} eq $qual ||
				die "Read $name in --max file $max_file has different sequence ".
				    "from input read.\nHit: $max_hash{$name}{qual}\nInput: $qual";
		}
		else {
			die "Read with name $name appears in input, but not in any of the output files";
		}
		if(defined($al_hash{$name})) {
			$al_hash{$name}{seq} eq $seq ||
				die "Read $name in --al file $al_file has different sequence ".
				    "from input read.\nHit: $al_hash{$name}{seq}\nInput: $seq";
			$al_hash{$name}{qual} eq $qual ||
				die "Read $name in --al file $al_file has different sequence ".
				    "from input read.\nHit: $al_hash{$name}{qual}\nInput: $qual";
		}
		$patid++;
		last if $patid == $num_reads;
	}
	close($READ);
}

if($al_file ne "") {
	$als == $distinctHits || die "number of --al $als does not match number of distinct hits $distinctHits";
}
$distinctHits + $uns + $maxs == $reads ||
	die "distinct hits $distinctHits, --un $uns, --max $maxs doesn't add to reads $reads";

print "PASSED; processed $hits hits, $reads reads, $uns --un, $maxs --max, $als --al\n";
