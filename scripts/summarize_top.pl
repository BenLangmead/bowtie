#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("v",\%options);

my $verbose = 0;
$verbose = $options{v} if defined $options{v};

defined($ARGV[0]) || die "Must specify top output file as first argument";
my $topFile = $ARGV[0];

defined($ARGV[1]) || die "Must specify search string as second argument";
my $appName = $ARGV[1];

my %vmmax = (); my %vmtot = ();
my %rsmax = (); my %rstot = ();
my %secstot = ();
my %secslast = ();
my %lines = ();

# Trim whitespace from a string argument
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub parseSz {
	my ($s, $name) = @_;
	if($s =~ /([0-9\.]+)m/) {
		print "Parsing as MB: $s; $1\n" if $verbose;
		$s = int($1 * 1024 * 1024);
	} elsif($s =~ /([0-9\.]+)g/) {
		print "Parsing as GB: $s; $1\n" if $verbose;
		$s = int($1 * 1024 * 1024 * 1024);
	} else {
		$s =~ /^([0-9\.]+)$/ || die "Bad format for $name: $s";
		print "Parsing as B: $s; $1\n" if $verbose;
		$s = int($1);
	}
	print "Got $name: $s\n" if $verbose;
	return $s;
}

open(TOP, $topFile);
while(<TOP>) {
	$_ = trim($_);
	next if not /$appName/;
	my @line = split;
	my $pid = int($line[0]);
	
	# Update line count
	$lines{$pid} = 0 unless defined($lines{$pid});
	$lines{$pid}++;
	
	# Parse VIRT
	my $vm = parseSz($line[4], "VIRT");
	$vmmax{$pid} = 0 unless defined($vmmax{$pid});
	$vmtot{$pid} = 0 unless defined($vmtot{$pid});
	if($vm > $vmmax{$pid}) {
		$vmmax{$pid} = $vm;
	}
	$vmtot{$pid} += $vm;

	# Parse RES
	my $rs = parseSz($line[5], "RES");
	$rsmax{$pid} = 0 unless defined($rsmax{$pid});
	$rstot{$pid} = 0 unless defined($rstot{$pid});
	if($rs > $rsmax{$pid}) {
		$rsmax{$pid} = $rs;
	}
	$rstot{$pid} += $rs;
	
	my @ts = split(/:/, $line[10]);
	$ts[1] =~ s/\..*//; # chop everything starting at the dot
	my $secsacc = int($ts[1]);
	$secsacc += (int($ts[0]) * 60);
	$secslast{$pid} = 0 unless defined($secslast{$pid});
	$secstot{$pid} = 0 unless defined($secstot{$pid});
	if($secsacc < $secslast{$pid}) {
		$secstot{$pid} += $secslast{$pid};
	}
	$secslast{$pid} = $secsacc;
}

# Report
for my $k (keys %secstot) {
	$secstot{$k} += $secslast{$k} if defined($secslast{$k});

	my $vmavg = int($vmtot{$k} / $lines{$k});
	my $rsavg = int($rstot{$k} / $lines{$k});
	
	my $vmmaxK = int($vmmax{$k} / 1024);
	my $vmmaxM = int($vmmaxK / 1024);
	
	my $vmavgK = int($vmavg / 1024);
	my $vmavgM = int($vmavgK / 1024);
	
	my $rsmaxK = int($rsmax{$k} / 1024);
	my $rsmaxM = int($rsmaxK / 1024);
	
	my $rsavgK = int($rsavg / 1024);
	my $rsavgM = int($rsavgK / 1024);
	
	my $min = int($secstot{$k}/60);
	while(length($min) < 2) { $min = "0".$min; }
	my $secs = ($secstot{$k} % 60);
	while(length($secs) < 2) { $secs = "0".$secs; }
	
	print "For pid: $k\n";
	print "  VIRT: max=$vmmaxM MB, avg=$vmavgM MB\n";
	print "  RES: max=$rsmaxM MB, avg=$rsavgM MB\n";
	print "  Total time: $min:$secs\n";
}
