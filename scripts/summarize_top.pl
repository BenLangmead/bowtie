#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("vt:",\%options);

my $verbose = 0;
$verbose = $options{v} if defined $options{v};

defined($ARGV[0]) || die "Must specify top output file as first argument";
my $topFile = $ARGV[0];

defined($ARGV[1]) || die "Must specify search string as second argument";
my $appName = $ARGV[1];

# Any process that lasts fewer than this many seconds doesn't count
my $threshold = 120;
$threshold = int($options{t}) if defined $options{t};

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

# Trim whitespace from a string argument
sub minSecToSeconds {
	my $time = shift;
	my @ts = split(/:/, $time);
	$#ts == 1 || die "m:s time string couldn't be parsed: $time";
	$ts[1] =~ s/\..*//; # chop everything starting at the dot
	my $secsacc = int($ts[1]);
	$secsacc += (int($ts[0]) * 60);
	return $secsacc;
}

# Trim whitespace from a string argument
sub minSecHrToSeconds {
	my $time = shift;
	my @ts = split(/:/, $time);
	$#ts == 2 || die "h:m:s time string couldn't be parsed: $time";
	#$ts[1] =~ s/\..*//; # chop everything starting at the dot
	my $secsacc = int($ts[2]);
	$secsacc += (int($ts[1]) * 60);
	$secsacc += (int($ts[0]) * 60 * 60);
	return $secsacc;
}

sub parseSz {
	my ($s, $name) = @_;
	if($s =~ /([0-9\.]+)m/) {
		$s = int($1 * 1024 * 1024);
	} elsif($s =~ /([0-9\.]+)g/) {
		$s = int($1 * 1024 * 1024 * 1024);
	} else {
		$s =~ /^([0-9\.]+)/ || die "Bad format for $name: $s";
		$s = int($1);
	}
	return $s;
}

# Aware of two top output formats:
#   PID USER      PR  NI %CPU    TIME+  %MEM  VIRT  RES  SHR S COMMAND
#   PID USER      PR  NI  VIRT  RES  SHR S %CPU %MEM    TIME+  COMMAND

my $virtCol = 4;  # or 7
my $resCol  = 5;  # or 8
my $timeCol = 10; # or 5

my $lastTopTime = 0;
my $timeTurnovers = 0;
open(TOP, $topFile);
while(<TOP>) {
	$_ = trim($_);
	if(/PID.*USER.*PR/) {
	    # Grab column headers
	    my @line = split;
	    if($line[7] eq "VIRT") {
			$virtCol = 7;
			$resCol  = 8;
			$timeCol = 5;
	    }
	}
	if(/^top/) {
		my @s = split;
		$s[2] =~ /[0-9]+:[0-9]+:[0-9]+/ || die "Malformed time: $s[2]";
		my $toptime = minSecHrToSeconds($s[2]);
		$timeTurnovers++ if $toptime < $lastTopTime;
		$lastTopTime = $toptime;
		next;
	}
	next if not /$appName/;
	my @line = split;
	my $pid = int($line[0]);
	
	# Update line count
	$lines{$pid} = 0 unless defined($lines{$pid});
	$lines{$pid}++;
	
	# Parse VIRT
	my $vm = parseSz($line[$virtCol], "VIRT");
	$vmmax{$pid} = 0 unless defined($vmmax{$pid});
	$vmtot{$pid} = 0 unless defined($vmtot{$pid});
	if($vm > $vmmax{$pid}) {
		$vmmax{$pid} = $vm;
	}
	$vmtot{$pid} += $vm;

	# Parse RES
	my $rs = parseSz($line[$resCol], "RES");
	$rsmax{$pid} = 0 unless defined($rsmax{$pid});
	$rstot{$pid} = 0 unless defined($rstot{$pid});
	if($rs > $rsmax{$pid}) {
		$rsmax{$pid} = $rs;
	}
	$rstot{$pid} += $rs;
	
	my $secsacc = minSecToSeconds($line[$timeCol]);
	$secslast{$pid} = 0 unless defined($secslast{$pid});
	$secstot{$pid} = 0 unless defined($secstot{$pid});
	if($secsacc < $secslast{$pid}) {
		$secstot{$pid} += $secslast{$pid};
	}
	$secslast{$pid} = $secsacc;
}
close(TOP);

# Report
my $max_vmmax = 0;
my $max_rsmax = 0;
my $secstottot = 0;
for my $k (keys %secstot) {
	$secstot{$k} += $secslast{$k} if defined($secslast{$k});
	if($verbose) {
		print "Skipping pid $k b/c secstot $secstot{$k} is less than threshold $threshold\n";
	}
	next if $secstot{$k} < $threshold;

	my $vmavg = int($vmtot{$k} / $lines{$k});
	my $rsavg = int($rstot{$k} / $lines{$k});
	
	$max_vmmax = $vmmax{$k} if $vmmax{$k} > $max_vmmax;
	$max_rsmax = $rsmax{$k} if $rsmax{$k} > $max_rsmax;
	
	my $vmmaxK = int($vmmax{$k} / 1024);
	my $vmmaxM = int($vmmaxK / 1024);
	
	my $vmavgK = int($vmavg / 1024);
	my $vmavgM = int($vmavgK / 1024);
	
	my $rsmaxK = int($rsmax{$k} / 1024);
	my $rsmaxM = int($rsmaxK / 1024);
	
	my $rsavgK = int($rsavg / 1024);
	my $rsavgM = int($rsavgK / 1024);
	
	$secstottot += $secstot{$k};
	my $min = int($secstot{$k}/60);
	while(length($min) < 2) { $min = "0".$min; }
	my $secs = ($secstot{$k} % 60);
	while(length($secs) < 2) { $secs = "0".$secs; }
	
	print "For pid: $k\n";
	if($vmmaxM > 0 && $vmavgM > 0) {
		print "  VIRT: max=$vmmaxM MB, avg=$vmavgM MB\n";
	} else {
		print "  VIRT: max=$vmmaxK KB, avg=$vmavgK KB\n";
	}
	if($rsmaxM > 0 && $rsavgM > 0) {
		print "  RES: max=$rsmaxM MB, avg=$rsavgM MB\n";
	} else {
		print "  RES: max=$rsmaxK KB, avg=$rsavgK KB\n";
	}
	print "  Total time: $min:$secs\n";
}

my $max_vmmaxK = int($max_vmmax / 1024);
my $max_vmmaxM = int($max_vmmaxK / 1024);
my $max_rsmaxK = int($max_rsmax / 1024);
my $max_rsmaxM = int($max_rsmaxK / 1024);
my $min = int($secstottot/60);
while(length($min) < 2) { $min = "0".$min; }
my $secs = ($secstottot % 60);
while(length($secs) < 2) { $secs = "0".$secs; }

print "Overall:\n";
if($max_vmmaxM > 0) {
	print "  VIRT: max=$max_vmmaxM MB\n";
} else {
	print "  VIRT: max=$max_vmmaxK KB\n";
}
if($max_rsmaxM > 0) {
	print "  RES: max=$max_rsmaxM MB\n";
} else {
	print "  RES: max=$max_rsmaxK KB\n";
}

# Overall wall-clock time
my $wallipre = trim(`egrep '^top' $topFile | head -1 | cut -d' ' -f 3`);
$wallipre =~ /[0-9]+:[0-9]+:[0-9]+/ || die "Malformed time: $wallipre; topFile: $topFile";
my $wallfpre = trim(`egrep '^top' $topFile | tail -1 | cut -d' ' -f 3`);
$wallfpre =~ /[0-9]+:[0-9]+:[0-9]+/ || die "Malformed time: $wallfpre; topFile: $topFile";
my $walli = minSecHrToSeconds($wallipre);
my $wallf = minSecHrToSeconds($wallfpre);
$wallf += ($timeTurnovers * 24 * 60 * 60);
$walli >= 0 || die "Bad initial wall-clock time: $walli";
$wallf >= 0 || die "Bad final wall-clock time: $wallf";
$wallf >= $walli || die "Final wall-clock time $wallf does not exceed initial $walli";
my $walltotsecs = $wallf-$walli;

my $wmin = int(($wallf-$walli)/60);
while(length($wmin) < 2) { $wmin = "0".$wmin; }
my $wsecs = (($wallf-$walli) % 60);
while(length($wsecs) < 2) { $wsecs = "0".$wsecs; }

print "  Total CPU time: $min:$secs\n";
print "  Total wall-clock time: $wmin:$wsecs\n";
print "$secstottot,$walltotsecs,$max_vmmaxK,$max_rsmaxK\n";
