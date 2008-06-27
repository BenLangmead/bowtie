#!/usr/bin/perl -w

use Getopt::Std;

my %options=();
getopts("v",\%options);

defined($ARGV[0]) || die "Must specify top output file as first argument";
my $topFile = $ARGV[0];

my $appName = "ebwt_build";
$appName = $ARGV[1] if defined($ARGV[1]);

my $vmmax = 0; my $vmtot = 0;
my $rsmax = 0; my $rstot = 0;
my $secstot = 0;
my $secsacc = 0;
my $secslast = 0;
my $lines = 0;
my $verbose = 0;
$verbose = $options{v} if defined $options{v};

# Trim whitespace from a string argument
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub parseSz {
	my ($s, $name) = @_;
	if($s =~ /([0-9]+)m/) {
		print "Parsing as MB: $s; $1\n" if $verbose;
		$s = int($1) * 1024 * 1024;
	} elsif($s =~ /([0-9]+)g/) {
		print "Parsing as GB: $s; $1\n" if $verbose;
		$s = int($1) * 1024 * 1024 * 1024;
	} else {
		$s =~ /^([0-9]+)$/ || die "Bad format for $name: $s";
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
	$lines++;
	
	# Parse VIRT
	my $vm = parseSz($line[4], "VIRT");
	if($vm > $vmmax) {
		$vmmax = $vm;
	}
	$vmtot += $vm;

	# Parse RES
	my $rs = parseSz($line[5], "RES");
	if($rs > $rsmax) {
		$rsmax = $rs;
	}
	$rstot += $rs;
	
	my @ts = split(/:/, $line[10]);
	$ts[1] =~ s/\..*//; # chop everything starting at the dot
	$secsacc = int($ts[1]);
	$secsacc += (int($ts[0]) * 60);
	if($secsacc < $secslast) {
		$secstot += $secslast;
	}
	$secslast = $secsacc;
}
$secstot += $secslast;

my $vmavg = int($vmtot / $lines);
my $rsavg = int($rstot / $lines);

my $vmmaxK = $vmmax / 1024; my $vmmaxM = $vmmaxK / 1024;
my $vmavgK = $vmavg / 1024; my $vmavgM = $vmavgK / 1024;

my $rsmaxK = $rsmax / 1024; my $rsmaxM = $rsmaxK / 1024;
my $rsavgK = $rsavg / 1024; my $rsavgM = $rsavgK / 1024;

#print "VIRT: max=$vmmax ($vmmaxK KB, $vmmaxM MB), avg=$vmavg ($vmavgK KB, $vmavgM MB)\n";
#print "RES: max=$rsmax ($rsmaxK KB, $rsmaxM MB), avg=$rsavg ($rsavgK KB, $rsavgM MB)\n";
print "VIRT: max=$vmmaxM MB, avg=$vmavgM MB\n";
print "RES: max=$rsmaxM MB, avg=$rsavgM MB\n";
print "Total time: ".($secstot/60)."m:".($secstot%60)."s\n";
