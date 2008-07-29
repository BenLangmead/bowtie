#!/usr/bin/perl -w

use strict;
use warnings;

defined($ARGV[0]) || die "Must specify run names";
my @runnames = @ARGV; # -> column names

my @names = ("Bowtie -n 1",
             "Bowtie",
             "Maq -n 1",
             "Maq",
             "Soap -v 1",
             "Soap");

system("cp headerinc.tex kg.tex") == 0 || die ("Must have headerinc.tex");
open(KG, ">>kg.tex") || die "Could not open >>kg.tex";
print KG "\\begin{document}\n";
print KG "\\begin{table}[tp]\n";
print KG "\\scriptsize\n";
print KG "\\begin{tabular}{lrrrrr}";
print KG "\\toprule\n";
print KG " & \\multirow{2}{*}{CPU Time} & Wall clock & Bowtie  & \\multicolumn{2}{c}{Reads mapped} \\\\\n";
print KG " &                            & time       & Speedup & Overall    & w/r/t Bowtie \\\\[3pt]\n";
print KG "\\toprule\n";

my $bowtieSecs = 0;
my $bowtiePct = 0;
for(my $ni = 0; $ni <= $#runnames; $ni++) {
	my $n = $runnames[$ni];
	next if $n =~ /^-/;
	my $l = readfline("whole.results.$n.txt", 0);
	my @ls = split(/,/, $l);
	my $rt = toMinsSecsHrs($ls[1]);
	my $wrt = toMinsSecsHrs($ls[2]);
	my $pct = trim($ls[3]);
	$bowtieSecs = $ls[2] if $ni == 0;
	$bowtiePct = $pct if $ni == 0;
	my $speedup = ($ni == 0 ? "-" : sprintf("%2.1fx", $ls[2] * 1.0 / $bowtieSecs));
	my $moreReads = ($ni == 0 ? "-" : sprintf("%2.1f\\%%", abs($pct - $bowtiePct) * 100.0 / $bowtiePct));
	my $moreReadsSign = ($pct >= $bowtiePct)? "+" : "-";
	$moreReadsSign = "" if $ni == 0;
	print KG "$names[$ni] & $rt & $wrt & $speedup & ";
	printf KG "%2.1f\\%%", $pct;
	print KG " & $moreReadsSign$moreReads \\\\";
	print KG "\\midrule" if $ni < $#runnames;
	print KG "\n";
}

print KG "\\bottomrule\n";
print KG "\\end{tabular}\n";
print KG "\\caption{".
	"CPU time ".
	"for mapping 8.96M 35bp ".
	"Illumina/Solexa reads (1 lane's worth) against the whole human ".
	"genome on a workstation with a 2.40GHz Intel Core 2 Q6600 and 2 GB of RAM. ".
	"Reads were originally extracted as part of the 1000-Genomes project pilot. ".
	"They were downloaded from the NCBI Short Read archive, accession \\#SRR001115. ".
	"Reference sequences were the ".
	"contigs of Genbank human genome build 36.3. ".
	"Soap was not run against the whole-human reference because its ".
	"memory footprint exceeds physical RAM. ".  
	"For the Maq runs, the ".
	"reads were first divided into chunks of 2M reads each, ".
	"as per the Maq Manual.".
	"}\n";
print KG "\\end{table}\n";
print KG "\\end{document}\n";
close(KG);

sub trim {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub readlines {
	my $f = shift;
	my @ret;
	open(FILE, $f) || die "Could not open $f";
	while(<FILE>) {
		chomp;
		push(@ret, $_)
	}
	return @ret;
}

sub readfline {
	my $f = shift;
	my $l = shift;
	my @ret;
	open(FILE, $f) || die "Could not open $f";
	while(<FILE>) {
		chomp;
		push(@ret, $_)
	}
	return $ret[$l] if $l <= $#ret;
	return "";
}

sub toMinsSecsHrs {
	my $s = shift;
	my $hrs  = int($s / 60 / 60);
	my $mins = int(($s / 60) % 60);
	my $secs = int($s % 60);
	while(length($secs) < 2) { $secs = "0".$secs; }
	if($hrs > 0) {
		while(length($mins) < 2) { $mins = "0".$mins; }
		return $hrs."h:".$mins."m:".$secs."s";
	} else {
		return $mins."m:".$secs."s";
	}
}

sub commaize {
	my $s = shift;
	return $s if length($s) <= 3;
	my $t = commaize(substr($s, 0, length($s)-3)) . "," . substr($s, length($s)-3, 3);
	return $t;
}
