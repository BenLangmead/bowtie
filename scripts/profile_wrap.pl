#!/usr/bin/perl -w

use strict;

if($#ARGV < 2) {
	print "Usage: perl profile_wrap.pl <run_name> <outdir> <cmd...>\n";
	exit 1;
}

my $runname = shift @ARGV;
my $outdir = shift @ARGV;

die "Empty runname argument" if $runname eq "";
die "Empty outdir argument" if $outdir eq "";

my $child = fork();
defined($child) || die "Could not fork!";
if($child == 0) {
	# I'm the child
	my $cmd = "time ";
	$cmd .= join " ", @ARGV;
	print "Spawning child with command: $cmd\n";
	exec(@ARGV);
}
else {
	# I'm the parent; see what processes I can spawn
	
	# Try 'heap', which, as far as I can tell, is only on OS X
	if(system("which heap > /dev/null") == 0) {
		print "Found heap; writing heap output to $outdir/heap.$runname.out\n";
		my $heapChild = fork();
		defined($heapChild) || die "Could not fork heap child!";
		if($heapChild == 0) {
			system("echo > $outdir/heap.$runname.out");
			while(1) {
				sleep 1;
				if(system("ps $child | grep $child > /dev/null") == 0) {
					system("heap $child 1>> $outdir/heap.$runname.out 2>&1");
				} else {
					print "Child w/ pid $child died; heap child exiting\n";
					exit 1;
				}
			}
		}
	}
	
	# Try 'vmmap'
	if(system("which vmmap > /dev/null") == 0) {
		print "Found vmmap; writing vmmap output to $outdir/vmmap.$runname.out\n";
		my $vmmapChild = fork();
		defined($vmmapChild) || die "Could not fork vmmap child!";
		if($vmmapChild == 0) {
			system("echo > $outdir/vmmap.$runname.out");
			while(1) {
				sleep 1;
				if(system("ps $child | grep $child > /dev/null") == 0) {
					system("vmmap $child 1>> $outdir/vmmap.$runname.out 2>&1");
				} else {
					print "Child w/ pid $child died; vmmap child exiting\n";
					exit 1;
				}
			}
		}
	}
	
	# Try batch-mode, specific-PID-mode, 1-iteration-mode 'top'
	if(system("top -n 1 -b -p $$ > /dev/null 2>&1") == 0) {
		print "Found top; writing top output to $outdir/top.$runname.out\n";
		my $topChild = fork();
		defined($topChild) || die "Could not fork top child!";
		if($topChild == 0) {
			system("echo > $outdir/top.$runname.out");
			while(1) {
				sleep 1;
				if(system("ps $child | grep $child > /dev/null") == 0) {
					system("top -n 1 -b -p $child 1>> $outdir/top.$runname.out 2>&1");
				} else {
					print "Child w/ pid $child died; top child exiting\n";
					exit 1;
				}
			}
		}
	}
}
