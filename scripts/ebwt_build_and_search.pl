#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("v",\%options);
my $verbose = $options{v} if defined($options{v});

defined($ARGV[0]) || die "Must define texts as first argument";
defined($ARGV[1]) || die "Must define patterns as second argument";

my $text = shift @ARGV;
my $pats = shift @ARGV;
my $args = "";
$args = join(" ", @ARGV) if defined(@ARGV);

my $cmd;
$cmd = "./ebwt_build-with-asserts -c -d -s $text .tmp2";
print "$cmd\n" if $verbose;
system($cmd);

$cmd = "./ebwt_search-with-asserts $args -c -1 --concise -a -s --orig $text .tmp2 $pats";
print "$cmd\n" if $verbose;
system($cmd);
