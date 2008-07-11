#!/usr/bin/perl -w

use strict;
use warnings;

defined($ARGV[0]) || die "Must specify run names";
my @runnames = @ARGV; # -> column names

open(RUNTIME, ">runtime.tex") || die "Could not open >runtime.tex";
print RUNTIME "\\documentclass{article}\n";
print RUNTIME "\\begin{document}\n";
print RUNTIME "\\begin{tabular}{ | l || ";
for(my $i = 0; $i <= $#runnames; $i++) {
	print RUNTIME "r | ";
}
print RUNTIME "}\n";
print RUNTIME "\\hline\n";
print RUNTIME " & Human Chromosome 22 & Human Chromosome 2 & Whole Human Genome \\\\ \\hline \n";
print RUNTIME "\\end{tabular}\n";

open(MEMORY, ">memory.tex") || die "Could not open >memory.tex";
print MEMORY "\\documentclass{article}\n";
print MEMORY "\\begin{document}\n";
print MEMORY "\\begin{tabular}{ | l || ";
for(my $i = 0; $i <= $#runnames; $i++) {
	print MEMORY "r | ";
}
print MEMORY "}\n";
print MEMORY "\\hline\n";
print MEMORY " & Human Chromosome 22 & Human Chromosome 2 & Whole Human Genome \\\\ \\hline \n";
print MEMORY "\\end{tabular}\n";

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

foreach my $n (@runnames)  {
	my @names = readlines("$n.results.names.txt");
	my @times = readlines("$n.results.times.txt");
	my @vmmax = readlines("$n.results.vmmax.txt");
	my @rsmax = readlines("$n.results.rsmax.txt");
}

close(RUNTIME);
close(MEMORY);
