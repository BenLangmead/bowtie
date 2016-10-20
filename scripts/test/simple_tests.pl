#!/usr/bin/perl -w

##
# Give simple tests with known results to bowtie.
#

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin); 
use lib $Bin;
use List::Util qw(max min);
use Data::Dumper;
use DNA;
use Clone qw(clone);
use Test::Deep;

my $bowtie = "";
my $bowtie_build = "";

GetOptions(
	"bowtie=s"       => \$bowtie,
	"bowtie-build=s" => \$bowtie_build) || die "Bad options";

if(! -x $bowtie || ! -x $bowtie_build) {
	my $bowtie_dir = `dirname $bowtie`;
	my $bowtie_exe = `basename $bowtie`;
	my $bowtie_build_exe = `basename $bowtie_build`;
	chomp($bowtie_dir);
	chomp($bowtie_exe);
	chomp($bowtie_build_exe);
	system("make -C $bowtie_dir $bowtie_exe $bowtie_build_exe") && die;
}

(-x $bowtie)       || die "Cannot run '$bowtie'";
(-x $bowtie_build) || die "Cannot run '$bowtie_build'";

my %prog_pairs = ($bowtie => $bowtie_build, $bowtie." --large-index " => $bowtie_build." --large-index ");
 
my @cases = (

	# Check paired-end exclusions

	{ name     => "Paired-end 1",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG      TAGATGG
	#                  ^2            ^16
	#                                CCATCTA
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "CCATCTA" ],
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => [{ "2,16" => 1 }] },

	{ name     => "Paired-end 2",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG TTTTATA
	#                  ^2       ^11
	#                           TATAAAA
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "TATAAAA" ],
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => [{ "2,11" => 1 }] },

	{ name     => "Paired-end 3",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                       AAGCTTT
	#                  ^2   ^7
	#                       AAAGCTT
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "AAAGCTT" ],
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => [{ "2,7" => 1 }] },

	{ name     => "Paired-end 4, containment excluded",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                   ACGAAAG
	#                  ^2   ^7
	#                   CTTTCGT
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "CTTTCGT" ],
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => [{ "*,*" => 1 }] },

	{ name     => "Paired-end 5, allow contain",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                   ACGAAAG
	#                  ^2   ^7
	#                   CTTTCGT
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "CTTTCGT" ],
	  args     => [ "--allow-contain -v 0",
	                "--allow-contain -n 0" ],
	  pairhits => [ { "2,3" => 1 } ] },

	{ name     => "Paired-end 6, allow contain",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "CTTTCGTT" ],
	  args     => [ "--allow-contain -v 0",
	                "--allow-contain -n 0" ],
	  pairhits => [ { "2,2" => 1 } ] },

	{ name     => "Paired-end 7, allow contain",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "CTTTCGTT" ],
	  args     => [ "--allow-contain -v 0 -m 1",
	                "--allow-contain -n 0 -m 1" ],
	  pairhits => [ { "2,2" => 1 } ] },

	{ name     => "Paired-end 8, allow contain",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "AACGAAAG" ],
	  args     => [ "--allow-contain --ff -v 0 -m 1",
	                "--allow-contain --ff -n 0 -m 1" ],
	  pairhits => [ { "2,2" => 1 } ] },

	{ name     => "Paired-end 9, allow contain",
	  ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   => [ "AACGAAAG" ],
	  mate2s   => [ "AACGATAG" ],
	  args     => [ "--allow-contain --ff -v 1 -m 1",
	                "--allow-contain --ff -n 1 -m 1" ],
	  pairhits => [ { "2,2" => 1 } ] },

	# Check basic -m funtionality

	{ name     => "Checking -m, 1",
	  ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => [ "-v 0",
	              "-n 0" ],
	  report =>   "-m 2 -a",
	  hits   => [ { 0 => 1, 8 => 1 } ] },

	{ name     => "Checking -m, 2",
	  ref    => [ "TTGTTCGTTTGTTCGTTTGTTCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => [ "-v 0",
	              "-n 0" ],
	  report =>   "-m 2 -a",
	  hits   => [ { } ] },
	
	# Check basic aedits field functionality

	{ name     => "Checking edits 1",
	  ref    => [ "TTGCCCGT" ],
	  reads  => [ "TTGTTCGT" ],
	  args   => [ "-v 2", "-n 2" ],
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "3:C>T,4:C>T" ],
	  orient =>   "+" },

	{ name     => "Checking edits 2",
	  ref    => [ "TTGTTCGT" ],
	  reads  => [ "ACGGGCAA" ],
	  args   => [ "-v 2",
	              "-n 2" ],
	  hits   => [ { 0 => 1 } ],
	  edits  => [ "3:T>C,4:T>C" ],
	  orient =>   "-" },

	{ name     => "Checking edits 3",
	  ref   => [ "ACGTTCGT" ],
	  reads => [ "GTTC" ],
	  args  => [ "-v 0",
	             "-n 0" ],
	  hits  => [ { 2 => 1 } ] },
	
	{ name  => "Colorspace 1",
	  ref   => [ "AAACGAAAGCTTTTATAGATGGGG" ],
	  reads => [    "132002320003332231" ],
	  args  => [ "-C -v 0",
	             "-C -n 0",
				 "-C -n 1",
				 "-C -v 1",
				 "-C -v 2",
				 "-C -n 2" ],
	  hits  => [ { 3 => 1 } ],
	  color => 1 },

	{ name  => "Colorspace 1",
	  ref   => [ "AAACGAAAGCTTTTATAGATGGGG" ],
	  reads => [    "132002320003332231" ],
	  args  => [ "-C -v 0",
	             "-C -n 0",
				 "-C -n 1",
				 "-C -v 1",
				 "-C -v 2",
				 "-C -n 2" ],
	  hits  => [ { 3 => 1 } ],
	  color => 1 },

	{ name  => "Colorspace FASTA",
	  ref   => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	  #           0123456789012345678901234
	  #                     1         2
	  #               CGAAAGCTTGTATAGAT
	  fasta =>   ">r0\n132002320113332231\n",
	  args  => [ "-C -v 2",
				 "-C -n 2" ],
	  hits  => [ { 4 => 1 } ],
	  color => 1 },

	{ name  => "Colorspace FASTQ",
	  ref   => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	  fastq =>   "\@r0\n132002320113332231\n+\nIIIIIIIIIIIIIIIIII\n",
	  args  => [ "-C -v 2",
				 "-C -n 2" ],
	  hits  => [ { 4 => 1 } ],
	  color => 1 },

	{ name   => "Colorspace Tabbed",
	  ref    => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	  tabbed =>   "r0\t132002320113332231\tIIIIIIIIIIIIIIIIII\n",
	  args   => [ "-C -v 2",
	              "-C -n 2" ],
	  hits   => [ { 4 => 1 } ],
	  color  => 1 },

	{ name   => "Colorspace raw",
	  ref    => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	  raw    =>   "132002320113332231\n",
	  args   => [ "-C -v 2",
	              "-C -n 2" ],
	  hits   => [ { 4 => 1 } ],
	  color  => 1 },

	{ name  => "Colorspace 3",
	  ref   => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	  reads => [     "132002320113332231" ],
	  args  => [ "-C -v 1",
				 "-C -n 1" ],
	  hits =>  [ { } ],
	  color => 1 },

	{ name  => "Colorspace 4",
	  ref   => [ "ATATATGTCGACATATATATATATAT" ],
	  reads => [     "3311232113333333" ],
	  args  => [ "-C -v 0",
				 "-C -n 0" ],
	  hits  => [ { 4 => 1 } ],
	  color => 1 },
);

##
# Take a list of reference sequences and write them to a temporary
# FASTA file of the given name.
#
sub writeFasta($$) {
	my ($l, $fa) = @_;
	open(FA, ">$fa") || die "Could not open $fa for writing";
	my $idx = 0;
	for(my $i = 0; $i < scalar(@$l); $i++) {
		print FA ">$idx\n".$l->[$i]."\n";
		$idx++;
	}
	close(FA);
}

##
# Take a lists of named reads/mates and write them to appropriate
# files.
#
sub writeReads($$$$$$$$$) {
	my (
		$reads,
		$quals,
		$mate1s,
		$qual1s,
		$mate2s,
		$qual2s,
		$names,
		$fq1,
		$fq2) = @_;

	open(FQ1, ">$fq1") || die "Could not open '$fq1' for writing";
	open(FQ2, ">$fq2") || die "Could not open '$fq2' for writing";
	my $pe = (defined($mate1s) && $mate1s ne "");
	if($pe) {
		for (0..scalar(@$mate1s)-1) {
			my $m1 = $mate1s->[$_];
			my $m2 = $mate2s->[$_];
			my $q1 = $qual1s->[$_];
			my $q2 = $qual2s->[$_];
			my $nm = $names->[$_];
			defined($m1) || die;
			defined($m2) || die;
			$q1 = $q1 || ("I" x length($m1));
			$q2 = $q2 || ("I" x length($m2));
			$nm = $nm || "r$_";
			print FQ1 "\@$nm/1\n$m1\n+\n$q1\n";
			print FQ2 "\@$nm/2\n$m2\n+\n$q2\n";
		}
	} else {
		for (0..scalar(@$reads)-1) {
			my $read = $reads->[$_];
			defined($read) || die;
			my $qual = $quals->[$_];
			my $nm = $names->[$_];
			$qual = $qual || ("I" x length($read));
			$nm = $nm || "r$_";
			print FQ1 "\@$nm\n$read\n+\n$qual\n";
		}
	}
	close(FQ1);
	close(FQ2);
}

##
# Run bowtie with given arguments
#
sub runbowtie($$$$$$$$$$$$$$$$$$$$$$$) {

	my (
		$do_build,
		$large_idx,
		$color,
		$debug_mode,
		$args,
		$fa,
		$reportargs,       #5
		$read_file_format,
		$read_file,
		$mate1_file,
		$mate2_file,
		$reads,
		$quals,
		$mate1s,
		$qual1s,
		$mate2s,
		$qual2s,
		$names,
		$ls,
		$rawls,
		$header_ls,
		$raw_header_ls,
		$should_abort) = @_;

	my $idx_type = "";
	$args .= " --quiet";
	$reportargs = "-a" unless defined($reportargs);
	$args .= " $reportargs";
	if ($large_idx){
	    $idx_type = "--large-index";
	}
	# Write the reference to a fasta file
	print "References:\n";
	open(FA, $fa) || die;
	while(<FA>) { print $_; }
	close(FA);
	if($do_build) {
		my $build_args = "";
		$build_args .= " -C " if $color;
		my $cmd = "$bowtie_build $idx_type --quiet --sanity $build_args $fa .simple_tests.tmp";
		print "$cmd\n";
		system($cmd);
		($? == 0) || die "Bad exitlevel from bowtie-build: $?";
	}
	my $pe = (defined($mate1s) && $mate1s ne "");
	$pe = $pe || (defined($mate1_file));
	my $mate1arg;
	my $mate2arg;
	my $readarg;
	my $formatarg = "-c";
	my ($readstr, $m1str, $m2str) = (undef, undef, undef);
	$readstr = join(",", @$reads)  if defined($reads);
	$m1str   = join(",", @$mate1s) if defined($mate1s);
	$m2str   = join(",", @$mate2s) if defined($mate2s);
	if(defined($read_file) || defined($mate1_file)) {
		defined($read_file_format) || die;
		my $ext = "";
		if($read_file_format eq "fastq") {
			$formatarg = "-q";
			$ext = ".fq";
		} elsif($read_file_format eq "tabbed") {
			$formatarg = "--12";
			$ext = ".tab";
		} elsif($read_file_format eq "cline_reads") {
			$formatarg = "-c";
			$readarg = $read_file;
			$mate1arg = $mate1_file;
			$mate2arg = $mate2_file;
		} elsif($read_file_format eq "cont_fasta_reads") {
			$formatarg = "";
			$readarg = $read_file;
			$mate1arg = $mate1_file;
			$mate2arg = $mate2_file;
		} elsif($read_file_format eq "fasta") {
			$formatarg = "-f";
			$ext = ".fa";
		} elsif($read_file_format eq "qseq") {
			$formatarg = "--qseq";
			$ext = "_qseq.txt";
		} elsif($read_file_format eq "raw") {
			$formatarg = "-r";
			$ext = ".raw";
		} else {
			die "Bad format: $read_file_format";
		}
		if($formatarg ne "-c") {
			if(defined($read_file)) {
				# Unpaired
				open(RD, ">.simple_tests$ext") || die;
				print RD $read_file;
				close(RD);
				$readarg = ".simple_tests$ext";
			} else {
				defined($mate1_file) || die;
				defined($mate2_file) || die;
				# Paired
				open(M1, ">.simple_tests.1$ext") || die;
				print M1 $mate1_file;
				close(M1);
				open(M2, ">.simple_tests.2$ext") || die;
				print M2 $mate2_file;
				close(M2);
				$mate1arg = ".simple_tests.1$ext";
				$mate2arg = ".simple_tests.2$ext";
			}
		}
	} else {
		writeReads(
			$reads,
			$quals,
			$mate1s,
			$qual1s,
			$mate2s,
			$qual2s,
			$names,
			".simple_tests.1.fq",
			".simple_tests.2.fq");
		$mate1arg = ".simple_tests.1.fq";
		$mate2arg = ".simple_tests.2.fq";
		$formatarg = "-q";
		$readarg = $mate1arg;
	}
	# Possibly add debug mode string
	my $debug_arg = "";
	$debug_arg = "--debug" if $debug_mode;
	my $cmd;
	if($pe) {
		# Paired-end case
		$cmd = "$bowtie $debug_arg $idx_type $args .simple_tests.tmp $formatarg -1 $mate1arg -2 $mate2arg";
	} else {
		# Unpaired case
		$cmd = "$bowtie $debug_arg $idx_type $args .simple_tests.tmp $formatarg $readarg";
	}
	print "$cmd\n";
	open(BT, "$cmd |") || die "Could not open pipe '$cmd |'";
	while(<BT>) {
		print $_;
		chomp;
		if(substr($_, 0, 1) eq "@") {
			push @$header_ls, [ split(/\t/, $_, -1) ];
			push @$raw_header_ls, $_;
		} else {
			push @$ls, [ split(/\t/, $_, -1) ];
			push @$rawls, $_;
		}
	}
	close(BT);
	($? == 0 ||  $should_abort) || die "bowtie aborted with exitlevel $?\n";
	($? != 0 || !$should_abort) || die "bowtie failed to abort!\n";
}

##
# Compare a hash ref of expected SAM flags with a hash ref of observed SAM
# flags.
#
sub matchSamOptionalFlags($$) {
	my ($flags, $ex_flags) = @_;
	my %ex = ();
	for(keys %$ex_flags) {
		my ($nm, $ty, $vl) = split(/:/, $_);
		defined($vl) || die "Could not parse optional flag field \"$_\"";
		($ex{$nm}{ty}, $ex{$nm}{vl}) = ($ty, $vl);
	}
	for(keys %$flags) {
		my ($ex_ty, $ex_vl);
		if(defined($ex{$_})) {
			($ex_ty, $ex_vl) = ($ex{$_}{ty}, $ex{$_}{vl});
		} else {
			($ex_ty, $ex_vl) = ("i", "0");
		}
		defined($ex_ty) || die;
		defined($ex_vl) || die;
		my ($ty, $vl) = ($flags->{$_}{ty}, $flags->{$_}{vl});
		defined($ty) || die;
		defined($vl) || die;
		$ex_ty eq $ty ||
			die "Expected SAM optional flag $_ to have type $ex_ty, had $ty";
		$ex_vl eq $vl ||
			die "Expected SAM optional flag $_ to have value $ex_vl, had $vl";
	}
	return 1;
}

my $tmpfafn = ".simple_tests.pl.fa";
my $last_ref = undef;
foreach my $large_idx (undef,1) {
	foreach my $debug_mode (undef,1) {
		for (my $ci = 0; $ci < scalar(@cases); $ci++) {
			my $c = $cases[$ci];
			last unless defined($c);
			# If there's any skipping of cases to be done, do it here prior to the
			# eq_deeply check
			my $do_build = 0;
			unless(defined($last_ref) && eq_deeply($c->{ref}, $last_ref)) {
				writeFasta($c->{ref}, $tmpfafn);
				$do_build = 1;
			}
			for my $a (@{$c->{args}}) {
				$last_ref = $c->{ref};
				# For each set of arguments...
				my $case_args = $a;
				$case_args = "" unless defined($a);
				my $first = 1; # did we build the index yet?
				# Forward, then reverse-complemented
				my $fwlo = ($c->{nofw} ? 1 : 0);
				my $fwhi = ($c->{norc} ? 0 : 1);
				for(my $fwi = $fwlo; $fwi <= $fwhi; $fwi++) {
					my $fw = ($fwi == 0);
					my $sam = 1;

					my $reads      = $c->{reads};
					my $quals      = $c->{quals};
					my $m1s        = $c->{mate1s};
					my $q1s        = $c->{qual1s};
					my $m2s        = $c->{mate2s};
					my $q2s        = $c->{qual2s};
					my $color      = $c->{color} || 0;

					my $read_file  = undef;
					my $mate1_file = undef;
					my $mate2_file = undef;

					$read_file  = $c->{fastq}   if defined($c->{fastq});
					$read_file  = $c->{tabbed}  if defined($c->{tabbed});
					$read_file  = $c->{fasta}   if defined($c->{fasta});
					$read_file  = $c->{qseq}    if defined($c->{qseq});
					$read_file  = $c->{raw}     if defined($c->{raw});
					$read_file  = $c->{cline_reads} if defined($c->{cline_reads});
					$read_file  = $c->{cont_fasta_reads} if defined($c->{cont_fasta_reads});

					$mate1_file = $c->{fastq1}  if defined($c->{fastq1});
					$mate1_file = $c->{tabbed1} if defined($c->{tabbed1});
					$mate1_file = $c->{fasta1}  if defined($c->{fasta1});
					$mate1_file = $c->{qseq1}   if defined($c->{qseq1});
					$mate1_file = $c->{raw1}    if defined($c->{raw1});
					$mate1_file = $c->{cline_reads1} if defined($c->{cline_reads1});
					$mate1_file = $c->{cont_fasta_reads1} if defined($c->{cont_fasta_reads1});

					$mate2_file = $c->{fastq2}  if defined($c->{fastq2});
					$mate2_file = $c->{tabbed2} if defined($c->{tabbed2});
					$mate2_file = $c->{fasta2}  if defined($c->{fasta2});
					$mate2_file = $c->{qseq2}   if defined($c->{qseq2});
					$mate2_file = $c->{raw2}    if defined($c->{raw2});
					$mate2_file = $c->{cline_reads2} if defined($c->{cline_reads2});
					$mate2_file = $c->{cont_fasta_reads2} if defined($c->{cont_fasta_reads2});

					my $read_file_format = undef;
					if(!defined($reads) && !defined($m1s) && !defined($m2s)) {
						defined($read_file) || defined($mate1_file) || die;
						$read_file_format = "fastq"  if defined($c->{fastq})  || defined($c->{fastq1});
						$read_file_format = "tabbed" if defined($c->{tabbed}) || defined($c->{tabbed});
						$read_file_format = "fasta"  if defined($c->{fasta})  || defined($c->{fasta1});
						$read_file_format = "qseq"   if defined($c->{qseq})   || defined($c->{qseq1});
						$read_file_format = "raw"    if defined($c->{raw})    || defined($c->{raw1});
						$read_file_format = "cline_reads" if defined($c->{cline_reads}) || defined($c->{cline_reads1});
						$read_file_format = "cont_fasta_reads" if defined($c->{cont_fasta_reads}) || defined($c->{cont_fasta_reads1});
						next unless $fw;
					}
					# Run bowtie
					my @lines = ();
					my @rawlines = ();
					my @header_lines = ();
					my @header_rawlines = ();
					print $c->{name}." " if defined($c->{name});
					print "(fw:".($fw ? 1 : 0).", sam:$sam)\n";
					my $mate1fw = 1;
					my $mate2fw = 0;
					$mate1fw = $c->{mate1fw} if defined($c->{mate1fw});
					$mate2fw = $c->{mate2fw} if defined($c->{mate2fw});
					if(!$fw) {
						# Reverse-complement the reads
						my @s = (); @s = @$reads if defined($reads);
						my @q = (); @q = @$quals if defined($quals);
						# Reverse-complement mates and switch mate1 with mate2
						my @m1 = (); @m1 = @$m1s if defined($m1s);
						my @m2 = (); @m2 = @$m2s if defined($m2s);
						my @q1 = (); @q1 = @$q1s if defined($q1s);
						my @q2 = (); @q2 = @$q2s if defined($q2s);
						for(0..scalar(@s)-1) {
							if($color) {
								$s[$_] = reverse($s[$_]);
							} else {
								$s[$_] = DNA::revcomp($s[$_]);
							}
							$q[$_] = reverse $q[$_] if $_ < scalar(@q);
						}
						if($mate1fw == $mate2fw) {
							if($color) {
								for(0..$#m1) { $m1[$_] = reverse $m1[$_]; }
								for(0..$#m2) { $m2[$_] = reverse $m2[$_]; }
							} else {
								for(0..$#m1) { $m1[$_] = DNA::revcomp($m1[$_]); }
								for(0..$#m2) { $m2[$_] = DNA::revcomp($m2[$_]); }
							}
							for(0..$#q1) { $q1[$_] = reverse $q1[$_]; }
							for(0..$#q2) { $q2[$_] = reverse $q2[$_]; }
						}
						$reads = \@s if defined($reads);
						$quals = \@q if defined($quals);
						$m1s   = \@m2 if defined($m1s);
						$q1s   = \@q2 if defined($q1s);
						$m2s   = \@m1 if defined($m2s);
						$q2s   = \@q1 if defined($q2s);
					}
					#if(defined($m2s)) {
					#	$a .= " --";
					#	$a .= ($mate1fw ? "f" : "r");
					#	$a .= ($mate2fw ? "f" : "r");
					#}
					my $args = "$a";
					$args .= " -S" if $sam;
					runbowtie(
						$do_build && $first,
						$large_idx,
						$color,
						$debug_mode,
						$args,
						$tmpfafn,
						$c->{report},
						$read_file_format, # formate of read/mate files
						$read_file,        # read file
						$mate1_file,       # mate #1 file
						$mate2_file,       # mate #2 file
						$reads,            # read list
						$quals,            # quality list
						$m1s,              # mate #1 sequence list
						$q1s,              # mate #1 quality list
						$m2s,              # mate #2 sequence list
						$q2s,              # mate #2 quality list
						$c->{names},
						\@lines,
						\@rawlines,
						\@header_lines,
						\@header_rawlines,
						$c->{should_abort});
					$first = 0;
					my $pe = defined($c->{mate1s}) && $c->{mate1s} ne "";
					$pe = $pe || defined($mate1_file);
					$pe = $pe || $c->{paired};
					my ($lastchr, $lastoff, $lastoff_orig) = ("", -1, -1);
					# Keep temporary copies of hits and pairhits so that we can
					# restore for the next orientation
					my $hitstmp = [];
					$hitstmp = clone($c->{hits}) if defined($c->{hits});
					my $pairhitstmp = [];
					$pairhitstmp = clone($c->{pairhits}) if defined($c->{pairhits});
					my $pairhits_orig_tmp = [];
					$pairhits_orig_tmp = clone($c->{pairhits_orig}) if defined($c->{pairhits_orig});
					# Record map from already-seen read name, read sequence and
					# quality to the place on the reference where it's reported.
					# This allows us to check that the pseudo-random generator
					# isn't mistakenly yielding different alignments for identical
					# reads.
					my %seenNameSeqQual = ();
					if(defined($c->{lines})) {
						my $l = scalar(@lines);
						$l == $c->{lines} || die "Expected $c->{lines} lines, got $l";
					}
					for my $li (0 .. scalar(@lines)-1) {
						my $l = $lines[$li];
						my ($readname, $orient, $chr, $off_orig, $off, $seq, $qual, $mapq,
							$oms, $editstr, $flagstr, $samflags, $cigar, $rnext, $pnext,
							$tlen);
						my %samoptflags = ();
						if($sam) {
							scalar(@$l) >= 11 ||
								die "Bad number of fields; expected at least 11 got ".
									scalar(@$l).":\n$rawlines[$li]\n";
							($readname, $samflags, $chr, $off) = @$l[0..3];
							($seq, $qual) = @$l[9..10];
							$orient = ((($samflags >> 4) & 1) == 0) ? "+" : "-";
							$mapq  = $l->[4]; # mapping quality
							$cigar = $l->[5]; # CIGAR string
							$rnext = $l->[6]; # ref seq of next frag in template
							$pnext = $l->[7]; # position of next frag in template
							$tlen  = $l->[8]; # template length
							if($pnext == 0) { $pnext = "*"; } else { $pnext--; }
							for(my $m = 11; $m < scalar(@$l); $m++) {
								next if $l->[$m] eq "";
								my ($nm, $ty, $vl) = split(/:/, $l->[$m]);
								defined($vl) ||
									die "Could not parse optional flag field $m: ".
										"\"$l->[$m]\"";
								$samoptflags{$nm}{ty} = $ty;
								$samoptflags{$nm}{vl} = $vl;
							}
							if($off > 0) { $off--; }
							else { $off = "*"; }
							$off_orig = $off;
							$off = "*" if $cigar eq "*";
						} else {
							scalar(@$l) == 9 ||
								die "Bad number of fields; expected 9 got ".
									scalar(@$l).":\n$rawlines[$li]\n";
							($readname, $orient, $chr, $off, $seq,
							 $qual, $oms, $editstr, $flagstr) = @$l;
							$off_orig = $off;
						}
						$readname =~ s/^\s+//;
						if($c->{check_random}) {
							my $rsqKey = "$readname\t$orient\t$seq\t$qual";
							my $rsqVal = "$chr\t$off";
							if(defined($seenNameSeqQual{$rsqKey})) {
								$seenNameSeqQual{$rsqKey} eq $rsqVal ||
									die "Two hits for read/seq/qual:\n$rsqKey\n".
										"had different alignments:\n".
										"$seenNameSeqQual{$rsqKey}\n$rsqVal\n";
							}
							$seenNameSeqQual{$rsqKey} = $rsqVal;
						}
						$readname ne "" || die "readname was blank:\n".Dumper($c);
						my $rdi = $readname;
						$rdi = substr($rdi, 1) if substr($rdi, 0, 1) eq "r";
						my $mate = 0;
						if($readname =~ /\//) {
							($rdi, $mate) = split(/\//, $readname);
							defined($rdi) || die;
						}
						$rdi = $c->{idx_map}{$rdi} if defined($c->{idx_map}{$rdi});
						$rdi ne "" || die "rdi was blank:\nreadname=$readname\n".Dumper($c);
						if($rdi != int($rdi)) {
							# Read name has non-numeric characters.  Figure out
							# what number it is by scanning the names list.
							my $found = 0;
							for(my $i = 0; $i < scalar(@{$c->{names}}); $i++) {
								if($c->{names}->[$i] eq $readname) {
									$rdi = $i;
									$found = 1;
									last;
								}
							}
							$found || die "No specified name matched reported name $readname";
						}
						# Make simply-named copies of some portions of the test case
						# 'hits'
						my %hits = ();
						%hits = %{$c->{hits}->[$rdi]} if
							defined($c->{hits}->[$rdi]);
						# 'flags'
						my $flags = undef;
						$flags = $c->{flags}->[$rdi] if
							defined($c->{flags}->[$rdi]);
						# 'samflags'
						my $ex_samflags = undef;
						$ex_samflags = $c->{ex_samflags}->[$rdi] if
							defined($c->{ex_samflags}->[$rdi]);
						# 'samflags_map'
						my $ex_samflags_map = undef;
						$ex_samflags_map = $c->{samflags_map}->[$rdi] if
							defined($c->{samflags_map}->[$rdi]);
						# 'samoptflags'
						my $ex_samoptflags = undef;
						$ex_samoptflags = $c->{samoptflags}->[$rdi] if
							defined($c->{samoptflags}->[$rdi]);
						# 'cigar'
						my $ex_cigar = undef;
						$ex_cigar = $c->{cigar}->[$rdi] if
							defined($c->{cigar}->[$rdi]);
						# 'cigar_map'
						my $ex_cigar_map = undef;
						$ex_cigar_map = $c->{cigar_map}->[$rdi] if
							defined($c->{cigar_map}->[$rdi]);
						# 'mapq_hi' - boolean indicating whether mapq is hi/lo
						my $ex_mapq_hi = undef;
						$ex_mapq_hi = $c->{mapq_hi}->[$rdi] if
							defined($c->{mapq_hi}->[$rdi]);
						# 'mapq'
						my $ex_mapq = undef;
						$ex_mapq = $c->{mapq}->[$rdi] if
							defined($c->{mapq}->[$rdi]);
						# 'mapq_map'
						my $ex_mapq_map = undef;
						$ex_mapq_map = $c->{mapq_map}->[$rdi] if
							defined($c->{mapq_map}->[$rdi]);
						# 'rnext_map'
						my $ex_rnext_map = undef;
						$ex_rnext_map = $c->{rnext_map}->[$rdi] if
							defined($c->{rnext_map}) && defined($c->{rnext_map}->[$rdi]);
						# 'pnext_map'
						my $ex_pnext_map = undef;
						$ex_pnext_map = $c->{pnext_map}->[$rdi] if
							defined($c->{pnext_map}) && defined($c->{pnext_map}->[$rdi]);
						# 'tlen_map'
						my $ex_tlen_map = undef;
						$ex_tlen_map = $c->{tlen_map}->[$rdi] if
							defined($c->{tlen_map}) && defined($c->{tlen_map}->[$rdi]);
						# 'flags_fw'
						my $flags_fw = undef;
						$flags_fw = $c->{flags_fw}->[$rdi] if
							defined($c->{flags_fw}->[$rdi]);
						# 'flags_rc'
						my $flags_rc = undef;
						$flags_rc = $c->{flags_rc}->[$rdi] if
							defined($c->{flags_rc}->[$rdi]);
						# 'pairhits'
						my %pairhits = ();
						%pairhits = %{$c->{pairhits}->[$rdi]} if
							defined($c->{pairhits}->[$rdi]);
						# 'pairhits_orig'
						my %pairhits_orig = ();
						%pairhits_orig = %{$c->{pairhits_orig}->[$rdi]} if
							defined($c->{pairhits_orig}->[$rdi]);
						# 'pairflags'
						my %pairflags = ();
						%pairflags = %{$c->{pairflags}->[$rdi]} if
							defined($c->{pairflags}->[$rdi]);
						# 'hits_are_superset'
						my $hits_are_superset = 0;
						$hits_are_superset = $c->{hits_are_superset}->[$rdi] if
							defined($ci);
						# edits
						my $ex_edits = undef;
						$ex_edits = $c->{edits}->[$rdi] if
							defined($c->{edits}->[$rdi]);
						if(!$sam) {
							# Bowtie flags
							if(defined($flags)) {
								$flagstr eq $flags ||
									die "Expected flags=\"$flags\", got \"$flagstr\"";
							}
							if(defined($flags_fw) && $fw) {
								$flagstr eq $flags_fw ||
									die "Expected flags=\"$flags_fw\", got \"$flagstr\"";
							}
							if(defined($flags_rc) && !$fw) {
								$flagstr eq $flags_rc ||
									die "Expected flags=\"$flags_rc\", got \"$flagstr\"";
							}
							if(defined($c->{flag_map})) {
								if(defined($c->{flag_map}->[$rdi]->{$off})) {
									$flagstr eq $c->{flag_map}->[$rdi]->{$off} ||
										die "Expected flags=\"$c->{flag_map}->[$rdi]->{$off}\"".
											" at offset $off, got \"$flagstr\"";
								}
							}
						}
						if($sam) {
							# SAM flags
							if(defined($ex_samflags)) {
								$samflags eq $ex_samflags ||
									die "Expected flags $ex_samflags, got $samflags";
							}
							if(defined($ex_samflags_map)) {
								if(defined($c->{samflags_map}->[$rdi]->{$off})) {
									my $ex = $c->{samflags_map}->[$rdi]->{$off};
									$samflags eq $ex || die
										"Expected FLAGS value $ex at offset $off, got $samflags"
								} else {
									die "Expected to see alignment with offset $off parsing samflags_map";
								}
							}
							# CIGAR string
							if(defined($ex_cigar)) {
								$cigar eq $ex_cigar ||
									die "Expected CIGAR string $ex_cigar, got $cigar";
							}
							if(defined($ex_cigar_map)) {
								if(defined($c->{cigar_map}->[$rdi]->{$off})) {
									my $ex = $c->{cigar_map}->[$rdi]->{$off};
									$cigar eq $ex || die
										"Expected CIGAR string $ex at offset $off, got $cigar"
								} else {
									die "Expected to see alignment with offset $off parsing cigar_map";
								}
							}
							# MAPQ
							if(defined($ex_mapq)) {
								$mapq eq $ex_mapq ||
									die "Expected MAPQ $ex_mapq, got $mapq";
							}
							if(defined($ex_mapq_map)) {
								if(defined($c->{mapq_map}->[$rdi]->{$off})) {
									my $ex = $c->{mapq_map}->[$rdi]->{$off};
									$mapq eq $ex || die
										"Expected MAPQ string $ex at offset $off, got $mapq"
								} else {
									die "Expected to see alignment with offset $off parsing mapq_map";
								}
							}
							# MAPQ
							if(defined($ex_mapq_hi)) {
								if($ex_mapq_hi == 0) {
									$mapq < 20 || die "Expected MAPQ < 20, got $mapq";
								} else {
									$mapq >= 20 || die "Expected MAPQ >= 20, got $mapq";
								}
							}
							if(defined($ex_mapq_map)) {
								if(defined($c->{mapq_map}->[$rdi]->{$off})) {
									my $ex = $c->{mapq_map}->[$rdi]->{$off};
									$mapq eq $ex || die
										"Expected MAPQ string $ex at offset $off, got $mapq"
								} else {
									die "Expected to see alignment with offset $off parsing mapq_map";
								}
							}
							# SAM optional flags
							if(defined($ex_samoptflags)) {
								matchSamOptionalFlags(\%samoptflags, $ex_samoptflags);
							}
							if(defined($c->{samoptflags_map})) {
								if(defined($c->{samoptflags_map}->[$rdi]->{$off})) {
									matchSamOptionalFlags(
										\%samoptflags,
										$c->{samoptflags_map}->[$rdi]->{$off});
								} else {
									die "Expected to see alignment with offset $off parsing samoptflags_map";
								}
							}
							if(defined($c->{samoptflags_flagmap})) {
								if(defined($c->{samoptflags_flagmap}->[$rdi]->{$samflags})) {
									matchSamOptionalFlags(
										\%samoptflags,
										$c->{samoptflags_flagmap}->[$rdi]->{$samflags});
								} else {
									die "Expected to see alignment with flag $samflags parsing samoptflags_flagmap";
								}
							}
							# RNEXT map
							if(defined($c->{rnext_map})) {
								if(defined($c->{rnext_map}->[$rdi]->{$off})) {
									my $ex = $c->{rnext_map}->[$rdi]->{$off};
									$rnext eq $ex || die
										"Expected RNEXT '$ex' at offset $off, got '$rnext'"
								} else {
									die "Expected to see alignment with offset $off parsing rnext_map".Dumper($c);
								}
							}
							# PNEXT map
							if(defined($c->{pnext_map})) {
								if(defined($c->{pnext_map}->[$rdi]->{$off})) {
									my $ex = $c->{pnext_map}->[$rdi]->{$off};
									$pnext eq $ex || die
										"Expected PNEXT '$ex' at offset $off, got '$pnext'"
								} else {
									die "Expected to see alignment with offset $off parsing pnext_map";
								}
							}
							# TLEN map
							if(defined($c->{tlen_map})) {
								if(defined($c->{tlen_map}->[$rdi]->{$off})) {
									my $ex = $c->{tlen_map}->[$rdi]->{$off};
									$tlen eq $ex || die
										"Expected TLEN '$ex' at offset $off, got '$tlen'"
								} else {
									die "Expected to see alignment with offset $off parsing tlen_map";
								}
							}
						}
						if($pe && $lastchr ne "") {
							my $offkey_orig = $lastoff.",".$off_orig;
							$offkey_orig = $off_orig.",".$lastoff_orig if $off_orig eq "*";

							my $offkey = $lastoff.",".$off;
							$offkey = $off.",".$lastoff if $off eq "*";

							if($lastoff ne "*" && $off ne "*") {
								$offkey = min($lastoff, $off).",".max($lastoff, $off);
							}
							if(defined($c->{pairhits}->[$rdi])) {
								defined($pairhits{$offkey}) ||
									die "No such paired off as $offkey in pairhits list: ".Dumper(\%pairhits)."\n";
								$c->{pairhits}->[$rdi]->{$offkey}--;
								delete $c->{pairhits}->[$rdi]->{$offkey} if $c->{pairhits}->[$rdi]->{$offkey} == 0;
								%pairhits = %{$c->{pairhits}->[$rdi]};
							}
							if(defined($c->{pairhits_orig}->[$rdi])) {
								defined($pairhits_orig{$offkey_orig}) ||
									die "No such paired off as $offkey in pairhits_orig list: ".Dumper(\%pairhits_orig)."\n";
								$c->{pairhits_orig}->[$rdi]->{$offkey_orig}--;
								delete $c->{pairhits_orig}->[$rdi]->{$offkey_orig} if $c->{pairhits_orig}->[$rdi]->{$offkey_orig} == 0;
								%pairhits_orig = %{$c->{pairhits_orig}->[$rdi]};
							}
							($lastchr, $lastoff, $lastoff_orig) = ("", -1, -1);
						} elsif($pe) {
							# Found an unpaired alignment from aligning a pair
							my $foundSe =
								defined($c->{pairhits}->[$rdi]) &&
								$c->{pairhits}->[$rdi]->{$off};
							if($foundSe) {
								$c->{pairhits}->[$rdi]->{$off}--;
								delete $c->{pairhits}->[$rdi]->{$off}
									if $c->{pairhits}->[$rdi]->{$off} == 0;
								%pairhits = %{$c->{pairhits}->[$rdi]};
							} else {
								($lastchr, $lastoff) = ($chr, $off);
							}
							# Found an unpaired alignment from aligning a pair
							$foundSe =
								defined($c->{pairhits_orig}->[$rdi]) &&
								$c->{pairhits_orig}->[$rdi]->{$off_orig};
							if($foundSe) {
								$c->{pairhits_orig}->[$rdi]->{$off_orig}--;
								delete $c->{pairhits_orig}->[$rdi]->{$off_orig}
									if $c->{pairhits_orig}->[$rdi]->{$off_orig} == 0;
								%pairhits_orig = %{$c->{pairhits_orig}->[$rdi]};
							} else {
								($lastchr, $lastoff, $lastoff_orig) = ($chr, $off, $off_orig);
							}
						} else {
							if(defined($c->{hits}->[$rdi]) && $off ne "*") {
								defined($hits{$off}) ||
									die "No such off as $off in hits list: ".Dumper(\%hits)."\n";
								$c->{hits}->[$rdi]->{$off}--;
								delete $c->{hits}->[$rdi]->{$off} if $c->{hits}->[$rdi]->{$off} == 0;
								%hits = %{$c->{hits}->[$rdi]};
							}
						}
						if(!$sam && defined($ex_edits)) {
							my $eds = $l->[7];
							$eds eq $ex_edits ||
								die "For edit string, expected \"$ex_edits\" got \"$eds\"\n";
						}
					}
					# Go through all the per-read
					my $klim = 0;
					$klim = scalar(@{$c->{hits}}) if defined($c->{hits});
					$klim = max($klim, scalar(@{$c->{pairhits}})) if defined($c->{pairhits});
					for (my $k = 0; $k < $klim; $k++) {
						# For each read
						my %hits     = %{$c->{hits}->[$k]}     if defined($c->{hits}->[$k]);
						my %pairhits = %{$c->{pairhits}->[$k]} if defined($c->{pairhits}->[$k]);
						my %pairhits_orig = %{$c->{pairhits_orig}->[$k]} if defined($c->{pairhits_orig}->[$k]);
						my $hits_are_superset = $c->{hits_are_superset}->[$k];
						# Check if there are any hits left over
						my $hitsLeft = scalar(keys %hits);
						if($hitsLeft != 0 && !$hits_are_superset) {
							print Dumper(\%hits);
							die "Had $hitsLeft hit(s) left over at position $k";
						}
						my $pairhitsLeft = scalar(keys %pairhits);
						if($pairhitsLeft != 0 && !$hits_are_superset) {
							print Dumper(\%pairhits);
							die "Had $pairhitsLeft hit(s) left over at position $k";
						}
						my $pairhits_orig_Left = scalar(keys %pairhits_orig);
						if($pairhits_orig_Left != 0 && !$hits_are_superset) {
							print Dumper(\%pairhits_orig);
							die "Had $pairhits_orig_Left hit(s) left over at position $k";
						}
					}

					$c->{hits} = $hitstmp;
					$c->{pairhits} = $pairhitstmp;
					$c->{pairhits_orig} = $pairhits_orig_tmp;
				}
				$last_ref = undef if $first;
			}
		}
    }
}
print "PASSED\n";
