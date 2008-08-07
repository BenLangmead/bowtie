#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("w",\%options);
my $workstation = 0; $workstation = $options{w} if defined($options{w});

defined($ARGV[0]) || die "Must specify run names";
my @runnames = @ARGV; # -> column names
print "Run names: @runnames\n";

my @names = ("Bowtie -n 1",
             "Bowtie",
             "Bowtie filtered",
             "Bowtie -v 2",
             "Maq -n 1",
             "Maq -n 1 filtered",
             "Maq",
             "Maq filtered",
             "Soap -v 1",
             "Soap");

system("cp headerinc.tex kg.tex") == 0 || die ("Must have headerinc.tex");
open(KG, ">>kg.tex") || die "Could not open >>kg.tex";
print KG "\\begin{document}\n";
print KG "\\begin{table}[tp]\n";
print KG "\\scriptsize\n";
#print KG "\\begin{tabular}{lrrrrrr}";
print KG "\\begin{tabular}{lrrrrr}\n";
#if($workstation) {
#	print KG "\\multicolumn{6}{c}{\\small{2.4 GHz Intel Core 2 workstation with 2 GB of RAM}}\\\\[3pt]\n";
#} else {
#	print KG "\\multicolumn{6}{c}{\\small{2.4 GHz AMD Opteron 850 server with 32 GB of RAM}}\\\\[3pt]\n";
#}
print KG "\\toprule\n";
#print KG " & \\multirow{2}{*}{CPU Time} & Wall clock & Bowtie  & \\multicolumn{2}{c}{Reads mapped} \\\\\n";
#print KG " &                            & time       & Speedup & Overall    & w/r/t Bowtie \\\\[3pt]\n";
print KG " & \\multirow{2}{*}{CPU Time} & Wall clock & Bowtie  & Peak virtual & Reads  \\\\\n";
print KG " &                            & time       & Speedup & memory usage & mapped \\\\[3pt]\n";
print KG "\\toprule\n";

my $bowtieSecs = 0;
my $bowtiePct = 0;
my $rows = 0;
my $addedExtraHrule = 0;
for(my $ni = 0; $ni <= $#runnames; $ni++) {
	my $n = $runnames[$ni];
	if($n =~ /^-/) {
		# No results for this tool; don't print a row
		print "Skipping results for $n due to dash\n";
		next;
	}
	if(!$addedExtraHrule && !($n =~ /1$/)) {
		# Add an extra horizontal rule separating the -n/-v 1 results
		# from the standard results
		print KG "\\midrule\n";
		$addedExtraHrule = 1;
	}
	# Print \midrule before all rows except the first
	if($rows > 0) {
		print KG "\\midrule\n";
	}
	$rows++;
	my $l = readfline("whole.results.$n.txt", 0);
	my @ls = split(/,/, $l);
	my $rt = toMinsSecsHrs($ls[1]);
	my $wrt = toMinsSecsHrs($ls[2]);
	my $pct = trim($ls[4]);
	my $vm = int((trim($ls[3]) + 512) / 1024);
	$vm = commaize($vm);
	my $isBowtie = ($n =~ /bowtie/i);
	$bowtieSecs = $ls[2] if $isBowtie;
	$bowtiePct = $pct if $isBowtie;
	my $speedup = ($isBowtie ? "-" : sprintf("%2.1fx", $ls[2] * 1.0 / $bowtieSecs));
	my $moreReads = ($isBowtie ? "-" : sprintf("%2.1f\\%%", abs($pct - $bowtiePct) * 100.0 / $bowtiePct));
	my $moreReadsSign = ($pct >= $bowtiePct)? "+" : "-";
	$moreReadsSign = "" if $isBowtie;
	print KG "$names[$ni] & $rt & $wrt & $speedup & $vm MB & ";
	printf KG "%2.1f\\%%", $pct;
	#print KG " & $moreReadsSign$moreReads \\\\";
	print KG "\\\\";
}
print KG "\n";

print KG "\\bottomrule\n";
print KG "\\end{tabular}\n";
print KG "\\caption{";
print KG
	"Performance measurements for mapping 8.96M 35bp Illumina/Solexa ".
	"reads against the whole human genome on a single CPU of a ";
if($workstation) {
	print KG "workstation with a 2.40GHz Intel Core 2 Q6600 processor and 2 GB of RAM. ";
} else {
	print KG "server with a 2.4 GHz AMD Opteron 850 processor and 32 GB of RAM. ";
}
print KG "Bowtie speedup is calculated with respect to wall clock time. ";
print KG
	"Both CPU time and wall clock times are included to demonstrate ".
	"that no one tool suffers disproportionately from I/O pauses or ".
	"contention with other processes on the system. ";
if($workstation) {
	print KG
		"Note that Maq indexes the reads as it maps them".
	    ", whereas Bowtie ".
	    "requires that an index of the genome be pre-built.  The cost ".
	    "of building the Bowtie index is not included in these ".
	    "timings since we expect that in practice that cost will ".
	    "be rapidly amortized across multiple mapping jobs, or ".
	    "that the researcher will simply download a pre-built ".
	    "index from a shared repository in much less time than is ".
	    "required to build from scratch. "
} else {
	print KG
		"Note that Maq (resp. Soap) indexes the reads (resp. ".
		"genome) as it maps, whereas the Bowtie mapper ".
	    "requires a pre-built index of the genome.  The cost ".
	    "of building the Bowtie index is not included in these ".
	    "timings since we expect that in practice that cost will ".
	    "be rapidly amortized across multiple mapping jobs, or ".
	    "that the researcher will simply download a pre-built ".
	    "index from a shared repository in much less time than is ".
	    "required to build from scratch. "
}
print KG
	"Reads are taken from the 1000-Genomes project pilot ".
	"via the NCBI Short Read archive, accession ".
	"\\#SRR001115 and trimmed to 35bps. ".
	"Reference sequences were the ".
	"contigs of Genbank human genome build 36.3. ";
if($workstation) {
	print KG
		"Soap was not run because its ".
		"memory footprint would have exceeded the physical RAM of the ".
		"workstation. ";  
}
print KG
	"For the Maq runs, the ".
	"reads were first divided into chunks of 2M reads each, ".
	"as per the Maq Manual. ";
if($workstation) {
	print KG "Maq v0.6.6 was used. ";
} else {
	print KG "Soap v1.10 and Maq v0.6.6 were used. ";
}
print KG "}\n";
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
