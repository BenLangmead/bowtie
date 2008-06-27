#!/usr/bin/perl -w
#######################################################################
# This software has been created by Genome Research Limited (GRL).    # 
#                                                                     #
# GRL hereby grants permission to use, copy, modify and distribute    # 
# this software and its documentation for non-commercial purposes     # 
# without fee at the user's own risk on the basis set out below.      #
#                                                                     #
# GRL neither undertakes nor accepts any duty whether contractual or  # 
# otherwise in connection with the software, its use or the use of    # 
# any derivative, and makes no representations or warranties, express #
# or implied, concerning the software, its suitability, fitness for   #
# a particular purpose or non-infringement.                           #
#                                                                     #
# In no event shall the authors of the software or GRL be responsible # 
# or liable for any loss or damage whatsoever arising in any way      # 
# directly or indirectly out of the use of this software or its       # 
# derivatives, even if advised of the possibility of such damage.     #
#                                                                     #
# Our software can be freely distributed under the conditions set out # 
# above, and must contain this copyright notice.                      #
#######################################################################

# Author: James Bonfield, March 2006
#
# Reads _seq.txt files and _prb.txt/_qcal.txt files in the main
# Bustard directory to produce a fastq output file. Sequence names are
# generated based on the input filename and their position in the
# file. Optionally an additional 'command' may be supplied to run as a
# filter to modify the name, eg to make the names unique across runs
#
# FH 5.21.08 Modified 'directory/filename hackery' to remove underscores 

use strict;
use Cwd;
use Getopt::Long;

# Argument parsing
my %opts = ("trim" => 0, "trimr" => 0, "calibrated" => 0, "quiet" => 0);

exit unless GetOptions(\%opts, "trim=i", "trimr=i", "name_sub=s", "calibrated|c", "quiet|q");

# Compute log-odds to fastq-scale mapping
my @rescale;
for (my $i = -100; $i < 100; $i++) {
    $rescale[$i+100] = chr (041+int(10*log(1+10**($i/10))/log(10)+.499));
}
my %hmap = ('A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, '.' => 0);

foreach my $f (@ARGV) {
    process_file($f);
}

# Processes a single _seq.txt (the filename) and it's corresponding _prb.txt
# quality file.
sub process_file {
    my ($in_fn) = @_;
    local $"="";

    print STDERR "processing $in_fn\n" unless $opts{quiet};

    # Directory/filename hackery
    $in_fn = getcwd() . "/$in_fn" if ($in_fn !~ /^\//);
    my ($dir, $fn) = ($in_fn =~ m:(.*)/(.*):);

    if ($fn !~ /(.*)seq.txt$/) {
	print STDERR "Filename '$fn' is not a *seq.txt file\n";
	return;
    }

    my $base = $1;

    # Default 'op' is based on directory name, or $opts{name_sub}
    my $old_dir = getcwd();
    chdir($dir);
    my @dirs = split('/', getcwd());

    chdir($old_dir);
    my ($machine, $run) = ("unknown", "unknown");
    if (exists($dirs[-4]) &&
	$dirs[-4] =~ m:^[0-9]+_slxa-([^-_]+).*[-_]([0-9]+)$:) {
	($machine, $run) = ($1, $2+0);
    }    
    my $op;
    if (exists($opts{name_sub})) {
	$op = $opts{name_sub};
    } else {
	#$op = "\"IL${machine}_${run}_\${lane}_\${tile}_\${x}_\${y}\"";
	$op = "\"IL\${lane}_\${tile}_\${x}_\${y}\"";
    }

    # Open the files
    my $qfn;
    if ($opts{calibrated}) {
	$qfn = "$opts{gerald_dir}/${base}_qcal.txt";
    } else {
	$qfn = "${base}prb.txt";
    }
    my $count = 0;

    open (my $qualfh,  "<", "$dir/$qfn")
	|| die "Couldn't open prb file '$qfn': $@";
    open (my $seqfh,   "<", "$dir/$fn") 
	|| die "Couldn't open seq file '$fn': $@";

    # Format of seq is <lane> <tile> <xpos> <ypos> <sequence>
    while (<$seqfh>) {
	my $q = (<$qualfh>);

    # Decode quality into 1 per called base
	my ($lane,$tile,$x,$y,$seq) = split(/\s+/,$_);
	my $i = 0;
	my @qual;
	if ($opts{calibrated}) {
	    chomp($q);
	    @qual = map {ord($_) - 64} split("", $q);
	} else {
	    my @qa = split('\t', $q);
	    @qual = map {[split()]->[$hmap{substr($seq, $i++, 1)}]} @qa;
	}

    # Map from log-odds to Phred scale
	@qual = map { $rescale[$_+100] } @qual;

	++$count;
	$_ = eval $op;
	die if $@;

	$seq =~ s/[^ACGT]/N/gi;
	print "\@$_\n" . substr($seq, $opts{trim}, length($seq) - $opts{trim} - $opts{trimr}) .
	    "\n+\n" . substr("@qual", $opts{trim}, length($seq) - $opts{trim} - $opts{trimr}) . "\n";
    }

    close($seqfh);
    close($qualfh);
}

__END__

=head1 NAME

solexa2fastq - converts Bustard seq and prb files into fastq format

=head1 SYNPOSIS

 solexa2fastq [options] s_1_1_seq.txt ...
 solexa2fastq [options] bustard_dir ...

=head1 DESCRIPTION

This converts a Bustard output files (I<*_seq.txt> and I<*_prb.txt>,
and optionally GERALD*/*_qcal.txt) into Fastq format. Only the
I<seq.txt> files need to be specified on the command line as the
program will automatically read either the associated I<prb.txt> or
I<qcal.txt> file. Alternatively specifying a directory acts as a
synonym for reading all I<*_seq.txt> files from within that
directory. Any number of files or directories may be specified, with
all translations being concatenated and written to stdout.

=head1 OPTIONS

=over

=item --trim <integer>

This specifies how many bases to trim from the start of the
sequence. This now defaults this is "0", but older runs will likely
require "1". (Should this be specifying use-bases and an
nYYYYY.. string instead to allow trimming off the other end too?)

=item --name_sub <perlop>

This is a perl expression applied to the generate the reading name.
It may use the perl variables $machine, $run, $lane, $tile, $x, $y and
$count with the the first two being computed from the current working
directory and the latter being an auto-incremented counter. By default
it uses "IL\${machine}_\${run}_\${lane}_\${tile}_\${x}_\${y}".

For example to use the older solexa2fastq default of
s_<lane>_<tile>_<counter> you'd specified:

 --name_sub '"s_${lane}_${tile}_${count}"'

A more complex example could be:

 --name_sub 'sprintf("s_%d_%04d_%X_%X", $run, $tie>

This specifies the GERALD output directory to use with the --calibrated
