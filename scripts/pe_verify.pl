#!/usr/bin/perl -w

#
# Verifies that Bowtie's paired-end mode gives alignments that are
# consistent with the alignments produced in single-end mode with -a
# and --nostrata options.
#
# Run this from the Bowtie directory.
#
# Usage: perl ../scripts/pe_verify.pl [mode] [index] [mates1] [nates2]
#
# Defaults are:
#  [mode] = -n 2
#  [index] = e_coli
#  [mates1] = reads/e_coli_1000_1.fq
#  [mates2] = reads/e_coli_1000_2.fq
#
# E.g.: perl ../scripts/pe_verify.pl -v 0
#       perl ../scripts/pe_verify.pl -v 0 e_coli reads/e_coli_1000_1.fq reads/e_coli_100000_e1_2.fq
#

use warnings;
use strict;
use Getopt::Std;

my %options=();
getopts("dv:n:e:t",\%options);

my $debug = $options{d} if defined($options{d});
my $match_mode = "-n 2";
$match_mode = "-v " . $options{v} if defined($options{v});
$match_mode = "-n " . $options{n} if defined($options{n});

my $extra_args = "";
$extra_args = $options{e} if defined($options{e});

print "Using match mode: $match_mode\n";

my $bowtie_dir = ".";
my $bowtie_exe = "bowtie";
$bowtie_exe .= "-debug" if $debug;

my $index  = "e_coli";
$index = $ARGV[0] if defined($ARGV[0]);
my $reads1 = "reads/e_coli_1000_1.fq";
$reads1 = $ARGV[1] if (defined($ARGV[1]));
my $reads2 = "reads/e_coli_1000_2.fq";
$reads2 = $ARGV[2] if (defined($ARGV[2]));

# Infer input type so we can provide Bowtie with appropriate option
if($reads1 =~ /\.fa/) {
	$reads2 =~ /\.fa/ || die "Reads files $reads1 and $reads2 have different extensions";
	$extra_args .= " -f ";
} elsif($reads1 =~ /\.raw/) {
	$reads2 =~ /\.raw/ || die "Reads files $reads1 and $reads2 have different extensions";
	$extra_args .= " -r ";
} elsif(!($reads1 =~ /\.fq/)) {
	(!($reads2 =~ /\.fq/)) || die "Reads files $reads1 and $reads2 have different extensions";
	$extra_args .= " -c ";
}

# Inner and outer differences to use when determining whether two
# single-end alignments satisfy the mating constraint
my $inner = 1;
my $outer = 250;

my $sesMustHavePes = 0; # force se pairs to have corresponding pe pairs 

# Defaults are for Illumina short-insert library
my $m1fw = 1;
my $m2fw = 0;

system("make -C $bowtie_dir $bowtie_exe") == 0 || die;

# Run Bowtie not in paired-end mode for first mate file
my $bowtie_se_cmd1 = "$bowtie_dir/$bowtie_exe $match_mode -y $extra_args -a --refidx --nostrata $index $reads1";

# Run Bowtie not in paired-end mode for second mate file
my $bowtie_se_cmd2 = "$bowtie_dir/$bowtie_exe $match_mode -y $extra_args -a --refidx --nostrata $index $reads2";

# Run Bowtie in paired-end mode
my $bowtie_pe_cmd = "$bowtie_dir/$bowtie_exe $match_mode -y $extra_args --refidx $index -1 $reads1 -2 $reads2";
print "$bowtie_pe_cmd\n";
open BOWTIE_PE, "$bowtie_pe_cmd |";

my $pes = 0;
my $pesFw = 0;
my $pesRc = 0;
my %peHash = ();
while(<BOWTIE_PE>) {
	next if /^Reported/;
	my $line1 = $_;
	my $line2 = <BOWTIE_PE>;
	my @l1s = split(/[\t]/, $line1);
	my @l2s = split(/[\t]/, $line2);
	$#l1s >= 5 || die "Paired-end alignment not formatted correctly: $line1";
	$#l2s >= 5 || die "Paired-end alignment not formatted correctly: $line2";
	my @l1rs = split(/\//, $l1s[0]);
	my @l2rs = split(/\//, $l2s[0]);
	$#l1rs >= 1 || die "Read not formatted correctly: $l1s[0]";
	$#l2rs >= 1 || die "Read not formatted correctly: $l2s[0]";
	$l1rs[0] eq $l2rs[0] || die "Before-/ parts of read names don't match: $l1rs[0], $l2rs[0]";
	my $mate1 = ($l1rs[$#l1rs] eq "1");
	my $mate1str = $mate1 ? "1" : "0";
	my $basename = $l1rs[0];
	my $loff = int($l1s[3]);
	my $roff = int($l2s[3]);
	$loff < $roff || die "Left/right mate seem to be switched for $basename";
	$roff - $loff >= $inner || die "Paired-end alignment seems to violate lower bound:\n$line1\n$line2";
	$roff - $loff <= $outer || die "Paired-end alignment seems to violate upper bound:\n$line1\n$line2";
	my $read1Short = $l1s[4];
	my $qual1Short = $l1s[5];
	my $read2Short = $l2s[4];
	my $qual2Short = $l2s[5];
	my $read1Mms = "";
	$read1Mms = $l1s[7] if defined($l1s[7]);
	my $read2Mms = "";
	$read2Mms = $l2s[7] if defined($l2s[7]);
	my $val = "$basename $mate1str $l1s[2] $l1s[3] $l2s[3] $read1Short:$qual1Short $read2Short:$qual2Short $read1Mms $read2Mms";
	die if defined ($peHash{$basename});
	$peHash{$basename} = $val;
	if($mate1) { $pesFw++; } else { $pesRc++; }
	$pes++;
}
close(BOWTIE_PE);

my $ses = 0;
my %seHash = ();
my %seMatedHash = ();
my %unmatchedSe = ();
my $unmatchedSes = 0;
my $unmatchedSeReads = 0;

print "$bowtie_se_cmd1\n";
open BOWTIE_SE1, "$bowtie_se_cmd1 |";
while(<BOWTIE_SE1>) {
	next if /^Reported/;
	$ses++;
	my @ls = split(/[\t]/, $_);
	my @lrs = split(/\//, $ls[0]);
	my $basename = $lrs[0];
	my $ref = $ls[2];
	my $off = int($ls[3]);
	my $len = length($ls[4]);
	my $readShort = $ls[4];
	my $qualShort = $ls[5];
	my $mms = "";
	$mms = $ls[7] if defined($ls[7]);
	my $key = "$ref $ls[1] $off $len $readShort $qualShort $mms";
	push @{ $seHash{$basename}{$ref} }, $key;
}
close(BOWTIE_SE1);

print "$bowtie_se_cmd2\n";
open BOWTIE_SE2, "$bowtie_se_cmd2 |";
open UNMATCHED_SE, ">.unmatched.se";
while(<BOWTIE_SE2>) {
	next if /^Reported/;
	$ses++;
	my @ls = split(/[\t]/, $_);
	my @lrs = split(/\//, $ls[0]);
	my $basename = $lrs[0];
	my $ref = $ls[2];
	my $off = int($ls[3]);
	my $len = length($ls[4]);
	my $fw = ($ls[1] eq "+");
	my $readShort = $ls[4];
	my $qualShort = $ls[5];
	my $mms = "";
	$mms = $ls[7] if defined($ls[7]);
	my $key = "$ref $ls[1] $off $len $readShort $qualShort $mms";
	# Is the other mate already aligned?
	if(defined($seHash{$basename}{$ref})) {
		# Get all of the alignments for the mate
		for my $om (@{ $seHash{$basename}{$ref} }) {
			my @oms = split(/ /, $om);
			$#oms == 6 || die "Wrong number of elements for oms: $#oms";
			my $oref = $oms[0];
			my $ofw = ($oms[1] eq "+");
			my $ooff = int($oms[2]);
			my $olen = int($oms[3]);
			my $oreadShort = $oms[4];
			my $oqualShort = $oms[5];
			my $omms = "";
			$omms = $oms[6] if defined($oms[6]);
			$oref eq $ref || die "Refs don't match: $oref, $ref";
			my $diff;
			my $peKey = "$basename ";
			if($ooff > $off) {
				# The #1 mate is on the right
				$diff = $ooff - $off + $olen;
				# upstream mate contains downstream one?
				#next if $off + $len >= $ooff + $olen;
				# mates are at the same position?
				next if $ooff == $off;
				next if ($diff < $inner || $diff > $outer);
				next if $ofw == $m1fw;
				next if $fw == $m2fw;
				$peKey .= "0 $ref $off $ooff $readShort:$qualShort $oreadShort:$oqualShort $mms $omms";
			} else {
				# The #1 mate is on the left
				$diff = $off - $ooff + $len;
				# upstream mate contains downstream one?
				#next if $ooff + $olen >= $off + $len;
				# mates are at the same position?
				next if $ooff == $off;
				next if ($diff < $inner || $diff > $outer);
				next if $fw != $m2fw;
				next if $ofw != $m1fw;
				$peKey .= "1 $ref $ooff $off $oreadShort:$oqualShort $readShort:$qualShort $omms $mms";
			}
			# Found a legitimate paired-end alignment using a pair of
			# single-end alignments
			next if $seMatedHash{$basename}; # already found corresponding paired-end
			if($sesMustHavePes) {
				defined($peHash{$basename}) ||
					die "Found single-end alignment for $basename, but no paired-end";
			} else {
				if(!defined($peHash{$basename})) {
					if(!defined($unmatchedSe{$basename})) {
						$unmatchedSe{$basename} = 0;
						$unmatchedSeReads++;
					}
					$unmatchedSe{$basename}++;
					$unmatchedSes++;
					print UNMATCHED_SE "$om\n$key\n";
				}
			}
			if(defined($peHash{$basename}) && $peHash{$basename} eq $peKey) {
				delete $peHash{$basename};
				$seMatedHash{$basename} = 1;
			} else {
				#print "No matchup:\n$peHash{$basename}\n$peKey\n";
			}
			#print "Found alignment for mate $otherMate of $ls[0]; diff: $diff\n";
		}
	}
}
close(BOWTIE_SE2);
close(UNMATCHED_SE);

my $die = 0;
for my $peKey (keys %peHash) {
	print "Paired-end $peKey has a paired-end alignment without single-end support\n";
	$die++;
}
$die && die "Found $die paired-end reads with no corresponding single-end mates";

if($unmatchedSes > 0) {
	print "Total of $unmatchedSes unmatched single-end alignments found for $unmatchedSes distinct pairs\n";
	print "Ref orientation off len seq quals mms\n";
	system("cat .unmatched.se");
	die;
}

print "PASSED; analyzed $pes paired-end ($pesFw fw, $pesRc rc) and $ses single-end alignments\n";
