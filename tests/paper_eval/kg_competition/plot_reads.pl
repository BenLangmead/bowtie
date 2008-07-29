#!/usr/bin/perl -w

sub commaize {
	my $s = shift;
	return $s if length($s) <= 3;
	my $t = commaize(substr($s, 0, length($s)-3)) . "," . substr($s, length($s)-3, 3);
	return $t;
}

defined($ARGV[0])
	|| die "Must specify number of reads mapped by both Maq and Bowtie as first argument";
defined($ARGV[1])
	|| die "Must specify number of reads mapped by Bowtie but not Maq as second argument";
defined($ARGV[2])
	|| die "Must specify number of reads mapped by Maq but not Bowtie as third argument";

my $inboth = int($ARGV[0]);
my $inebwt = int($ARGV[1]);
my $inmaq  = int($ARGV[2]);
my $tot = $inboth + $inebwt + $inmaq;

my $inbothc = commaize($inboth);
my $inebwtc = commaize($inebwt);
my $inmaqc  = commaize($inmaq);
my $totc    = commaize($tot);

my $bothpct = sprintf("%2.1f", $inboth * 100.0 / $tot);
my $ebwtpct = sprintf("%2.1f", $inebwt * 100.0 / $tot);
my $maqpct  = sprintf("%2.1f", $inmaq  * 100.0 / $tot);

#
# Print plain-text summary
#

print "     In either: $tot\n";
print "       In both: $inboth ($bothpct%)\n";
print "In Bowtie only: $inebwt ($ebwtpct%)\n";
print "   In Maq only: $inmaq ($maqpct%)\n";

#
# Print TeX
#

system("cp headerinc.tex reads.tex") == 0 || die ("Must have headerinc.tex");
open(READS, ">>reads.tex") || die "Could not open >>reads.tex";
print READS "\\begin{document}\n";
print READS "\\begin{table}[tp]\n";
print READS "\\scriptsize\n";
print READS "\\begin{tabular}{cccc}\n";

print READS "\\multicolumn{4}{c}{Reads mapped} \\\\[3pt] \n";
print READS "By Maq or      & By both Maq & By Bowtie   & By Maq         \\\\ \n";
print READS "Bowtie or both & and Bowtie  & but not Maq & but not Bowtie \\\\ \n";
print READS "\\toprule\n";
print READS "$totc & $inbothc ($bothpct\\%) & $inebwtc ($ebwtpct\\%) & $inmaqc ($maqpct\\%) \\\\ \n";
print READS "\\bottomrule\n";

print READS "\\end{tabular}\n";
print READS "\\caption{".
	"Blah".
	"}\n";
print READS "\\end{table}\n";
print READS "\\end{document}\n";

close(READS);
