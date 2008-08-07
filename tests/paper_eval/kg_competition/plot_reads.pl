#!/usr/bin/perl -w

sub commaize {
	my $s = shift;
	return $s if length($s) <= 3;
	my $t = commaize(substr($s, 0, length($s)-3)) . "," . substr($s, length($s)-3, 3);
	return $t;
}

defined($ARGV[0])
	|| die "Must specify alternative program as first argument";
defined($ARGV[1])
	|| die "Must specify output file prefix as second argument";
my $alt = $ARGV[0];
my $pfx = $ARGV[1];

defined($ARGV[2])
	|| die "Must specify number of reads mapped by both $alt and Bowtie as first argument";
defined($ARGV[3])
	|| die "Must specify number of reads mapped by Bowtie but not $alt as second argument";
defined($ARGV[4])
	|| die "Must specify number of reads mapped by $alt but not Bowtie as third argument";

my $inboth = int($ARGV[2]);
my $inebwt = int($ARGV[3]);
my $inalt  = int($ARGV[4]);
my $tot = $inboth + $inebwt + $inalt;

my $inbothc = commaize($inboth);
my $inebwtc = commaize($inebwt);
my $inaltc  = commaize($inalt);
my $totc    = commaize($tot);

my $bothpct = sprintf("%2.1f", $inboth * 100.0 / $tot);
my $ebwtpct = sprintf("%2.1f", $inebwt * 100.0 / $tot);
my $maqpct  = sprintf("%2.1f", $inalt  * 100.0 / $tot);

#
# Print plain-text summary
#

print "     In either: $tot\n";
print "       In both: $inboth ($bothpct%)\n";
print "In Bowtie only: $inebwt ($ebwtpct%)\n";
print "   In $alt only: $inalt ($maqpct%)\n";

#
# Print TeX
#

system("cp headerinc.tex $pfx.tex") == 0 || die ("Must have headerinc.tex");
open(READS, ">>$pfx.tex") || die "Could not open >>$pfx.tex";
print READS "\\begin{document}\n";
print READS "\\begin{table}[tp]\n";
print READS "\\scriptsize\n";
print READS "\\begin{tabular}{cccc}\n";

print READS "\\multicolumn{4}{c}{Reads mapped} \\\\[3pt] \n";
print READS "By $alt or     & By both $alt & By Bowtie    & By $alt         \\\\ \n";
print READS "Bowtie or both & and Bowtie   & but not $alt & but not Bowtie \\\\ \n";
print READS "\\toprule\n";
print READS "$totc & $inbothc ($bothpct\\%) & $inebwtc ($ebwtpct\\%) & $inaltc ($maqpct\\%) \\\\ \n";
print READS "\\bottomrule\n";

print READS "\\end{tabular}\n";
print READS "\\caption{".
	"Blah".
	"}\n";
print READS "\\end{table}\n";
print READS "\\end{document}\n";

close(READS);
