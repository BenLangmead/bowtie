#!/usr/bin/perl -w

use strict;
use warnings;

my @names  = ("16 GB", "8 GB", "4 GB", "2 GB");
my @exts   = ("blf", "bl4", "pkl", "pkt");
my @packed = ("no", "no", "no", "yes");  # TODO: get this automatically

system("cp headerinc.tex builds.tex") == 0 || die ("Must have headerinc.tex");
open(BUILD, ">>builds.tex") || die "Could not open >>builds.tex";
print BUILD "\\begin{document}\n";
print BUILD "\\begin{table}[tp]\n";
print BUILD "\\scriptsize\n";
print BUILD "\\begin{tabular}{rrrrcr}";
print BUILD "\\toprule\n";
print BUILD "Target peak & Actual peak & \\# Suffix   & Difference   & 2-bit-per-base & \\multirow{2}{*}{Build time} \\\\\n";
print BUILD "footprint   & footprint   & array blocks & cover period & references     & \\\\ \n";
print BUILD "\\toprule\n";

for(my $ni = 0; $ni <= $#names; $ni++) {
	my $n = $names[$ni];
	my $l = readfline("whole.results.txt", $ni);
	my @ls = split(/,/, $l);
	my $rt = toMinsHrs($ls[1]);
	my $gb = $ls[2] / 1024 / 1024;
	$gb = sprintf("%.2f", $gb);
	my $blocksFw = `grep 'Getting block' whole.ebwt_build.$exts[$ni].out | head -1 | sed 's/.* of //'`;
	# Not showing the fw and bw block counts separately because it's
	# too difficult to explain in the column header
	#my $blocksBw = `grep 'Getting block' whole.ebwt_build.$exts[$ni].out | tail -1 | sed 's/.* of //'`;
	my $dcv = `grep 'Difference-cover' whole.ebwt_build.$exts[$ni].out | sed 's/.*: //'`;
	print BUILD "$n & $gb GB & $blocksFw & $dcv & $packed[$ni] & $rt \\\\";
	print BUILD "\\midrule" if $ni < $#names;
	print BUILD "\n";
}

print BUILD "\\bottomrule\n";
print BUILD "\\end{tabular}\n";
print BUILD "\\caption{".
	"Maximum virtual memory footprint and wall-clock running times ".
	"for various configurations of the Bowtie index-building tool ".
	"run on a single CPU of a server with a 1.8 Ghz AMD Opteron 875 ".
	"processor and 32 GB of RAM. ".
	"In all cases, the index produced is about 2.2 GB on disk. ".
	"TODO: RECALCULATE ON PRIVET OR LARCH FOR CONSISTENCY.".
	"}\n";
print BUILD "\\end{table}\n";
print BUILD "\\end{document}\n";
close(BUILD);

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

sub toMinsHrs {
	my $s = shift;
	my $hrs  = int($s / 60 / 60);
	my $mins = int(($s / 60) % 60);
	if($hrs > 0) {
		while(length($mins) < 2) { $mins = "0".$mins; }
		return $hrs."h:".$mins."m";
	} else {
		return $mins."m";
	}
}

sub commaize {
	my $s = shift;
	return $s if length($s) <= 3;
	my $t = commaize(substr($s, 0, length($s)-3)) . "," . substr($s, length($s)-3, 3);
	return $t;
}
