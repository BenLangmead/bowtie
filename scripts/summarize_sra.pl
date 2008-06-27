#!/usr/bin/perl -w

#
# Summarize the contents of the SRA database.
# Usage:
#   perl summarize_sra.pl
#

use strict;
use warnings;
use Getopt::Std;

my %options=();
my %sraCount=();
my %sraInst=();
my %sraStud=();
my %sraInd=();
my %sraCent=();
my %sraPair=();
my %sraDesign=();

my %indCount=();
my %indInst=();
my %indStud=();
my %indInd=();
my %indCent=();
my %indPair=();
my %indDesign=();

system("curl 'http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=table&f=run&m=data&s=run' -o .sra.html 2> /dev/null") == 0 || die "curl failed";
open(SRA, "cat .sra.html | sed 's/<[/]tr>/\\n/' |") || die "Failed to open .sra.html";

my $sra;
my $inst;
my $stud;
my $ind;
my $cent;
my $design;
my $last = "";
while(<SRA>) {
	$inst   = $1 if /Instrument:\s*<span>\s*([^<]*)/;
	$stud   = $1 if /Study:\s*<span>\s*([^<]*)/;
	$cent   = $1 if /Center:\s*<span>\s*([^<]*)/;
	$design = $1 if /Design:\s*<span>\s*([^<]*)/;
	$ind    = $1 if /individual ([a-zA-Z]+[0-9]+)/;
	if(/Accession:\s*<span>\s*SRA([0-9][0-9][0-9][0-9][0-9][0-9])/) {
		unless(defined($design)) {
			print "WARNING: Undefined design!\n"; $last = $_; next;
		}
		$sra = $1 . '!' . $design;
		$sraCount{$sra}++;
		$sraInst{$sra}   = $inst if defined($inst);
		$sraStud{$sra}   = $stud if defined($stud);
		$sraCent{$sra}   = $cent if defined($cent);
		$sraInd{$sra}    = $ind if defined($ind);
		$sraDesign{$sra} = $design if defined($design);
		$sraPair{$sra}   = "yes" if defined($design) && $design =~ /pair.*library/i;
		
		if(defined($ind)) {
			my $indKey = $ind . '!' . $inst . '!';
			if(defined($design) && $design =~ /pair.*library/i) {
				$indKey .= "yes";
			}
			$indCount{$indKey}++;
			$indInst{$indKey}   = $inst if defined($inst);
			$indStud{$indKey}   = $stud if defined($stud);
			$indCent{$indKey}   = $cent if defined($cent);
			$indInd{$indKey}    = $ind if defined($ind);
			$indDesign{$indKey} = $design if defined($design);
			$indPair{$indKey}   = "yes" if defined($design) && $design =~ /pair.*library/i;
		}
		
		undef $inst   if defined($inst);
		undef $stud   if defined($stud);
		undef $cent   if defined($cent);
		undef $design if defined($design);
		undef $ind    if defined($ind);
	}
	$last = $_;
}

print "By SRA+Design:\n";
print "  Submission accession,# Runs,Instrument,Paired-end?,Study,Center,[Individual]\n";
for my $k (sort keys %sraCount) {
	if($sraCount{$k} >= 1) {
		my @sraAndDesign = split(/!/, $k);
		print "  $sraAndDesign[0],";
		print "$sraCount{$k},";
		print "$sraInst{$k},";
		if(defined($sraPair{$k})) { print "yes,"; } else { print "no,"; }
		print "$sraStud{$k},";
		print "$sraCent{$k}";
		print ",$sraInd{$k}" if defined($sraInd{$k});
		print "\n";
		print "  $sraDesign{$k}\n";
	}
}

print "By Individual+Instrument:\n";
print "  Individual,# Runs,Instrument,Paired-end?,Study,Center\n";
for my $k (sort keys %indCount) {
	if($indCount{$k} >= 1) {
		my @indInst = split(/!/, $k);
		print "  $indInst[0],";
		print "$indCount{$k},";
		print "$indInst{$k},";
		if(defined($indPair{$k})) { print "yes,"; } else { print "no,"; }
		print "$indStud{$k},";
		print "$indCent{$k}";
		print ",$indInd{$k}" if defined($sraInd{$k});
		print "\n";
		#print "  $indDesign{$k}\n";
	}
}
