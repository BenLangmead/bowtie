#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("v",\%options);
my $verbose = $options{v} if defined($options{v});

defined($ARGV[0]) || die "Must define texts as first argument";
defined($ARGV[1]) || die "Must define patterns as second argument";

my $cmd;
$cmd = "./ebwt_build-with-asserts  -c -d -s $ARGV[0] .tmp2";
print "$cmd\n" if $verbose;
system($cmd);

$cmd = "./ebwt_search-with-asserts -c -1 --concise -a -s --orig $ARGV[0] .tmp2 $ARGV[1]";
print "$cmd\n" if $verbose;
system($cmd);
