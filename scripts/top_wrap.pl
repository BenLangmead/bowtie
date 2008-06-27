#!/usr/bin/perl -w

use strict;

if($#ARGV < 2) {
	print "Usage: perl top_wrap.pl <run_name> <outdir> <cmd...>\n";
	exit 1;
}

my $runname = shift @ARGV;
my $outdir = shift @ARGV;

die "Empty runname argument" if $runname eq "";
die "Empty outdir argument" if $outdir eq "";

$SIG{CHLD} = "IGNORE";
my $child = fork();
defined($child) || die "Could not fork!";
if($child == 0) {
	# I'm the child
	my $cmd = "time ";
	$cmd .= join " ", @ARGV;
	print "Spawning child with command: $cmd\n";
	system(@ARGV);
} else {
	# Give other process a little time to catch up
	sleep 1; 
	# Try batch-mode, specific-PID-mode, 1-iteration-mode 'top'
	print "Writing top output for pid $child to $outdir/top.$runname.out\n";
	system("echo > $outdir/top.$runname.out"); # Truncate file
	while(1) {
		sleep 1;
		if(system("ps $child | grep $child") == 0) {
			system("top -n 1 -b | grep $ARGV[0] 1>> $outdir/top.$runname.out 2>&1");
		} else {
			print "Child w/ pid $child died; top child exiting\n";
			exit 1;
		}
	}
}
