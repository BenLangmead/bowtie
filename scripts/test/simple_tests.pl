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

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG      TAGATGG
	#                  ^2            ^16
	#                                CCATCTA
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "CCATCTA",
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => { "2,16" => 1 } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG TTTTATA
	#                  ^2       ^11
	#                           TATAAAA
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "TATAAAA",
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => { "2,11" => 1 } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                       AAGCTTT
	#                  ^2   ^7
	#                       AAAGCTT
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "AAAGCTT",
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => { "2,7" => 1 } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                   ACGAAAG
	#                  ^2   ^7
	#                   CTTTCGT
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "CTTTCGT",
	  args     => [ "-v 0",
	                "-n 0" ],
	  pairhits => { } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                   ACGAAAG
	#                  ^2   ^7
	#                   CTTTCGT
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "CTTTCGT",
	  args     => [ "--allow-contain -v 0",
	                "--allow-contain -n 0" ],
	  pairhits => { "2,3" => 1 } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "CTTTCGTT",
	  args     => [ "--allow-contain -v 0",
	                "--allow-contain -n 0" ],
	  pairhits => { "2,2" => 1 } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "CTTTCGTT",
	  args     => [ "--allow-contain -v 0 -m 1",
	                "--allow-contain -n 0 -m 1" ],
	  pairhits => { "2,2" => 1 } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "AACGAAAG",
	  args     => [ "--allow-contain --ff -v 0 -m 1",
	                "--allow-contain --ff -n 0 -m 1" ],
	  pairhits => { "2,2" => 1 } },

	{ ref      => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	#                  AACGAAAG 
	#                  AACGAAAG
	#                  ^2   ^7
	#                  CTTTCGTT
	  mate1s   =>   "AACGAAAG",
	  mate2s   =>   "AACGATAG",
	  args     => [ "--allow-contain --ff -v 1 -m 1",
	                "--allow-contain --ff -n 1 -m 1" ],
	  pairhits => { "2,2" => 1 } },

	# Check basic -m funtionality

	{ ref    => [ "TTGTTCGTTTGTTCGT" ],
	  reads  =>   "TTGTTCGT",
	  args   => [ "-v 0",
	              "-n 0" ],
	  report =>   "-m 2 -a",
	  hits   => { 0 => 1, 8 => 1 } },

	{ ref    => [ "TTGTTCGTTTGTTCGTTTGTTCGT" ],
	  reads  =>   "TTGTTCGT",
	  args   => [ "-v 0",
	              "-n 0" ],
	  report =>   "-m 2 -a",
	  hits   => { } },
	
	# Check basic aedits field functionality

	{ ref    => [ "TTGCCCGT" ],
	  reads  =>   "TTGTTCGT",
	  args   => [ "-v 2",
	              "-n 2" ],
	  hits   => { 0 => 1 },
	  edits  =>   "3:C>T,4:C>T",
	  orient =>   "+" },

	{ ref    => [ "TTGTTCGT" ],
	  reads  =>   "ACGGGCAA",
	  args   => [ "-v 2",
	              "-n 2" ],
	  hits   => { 0 => 1 },
	  edits  =>   "3:T>C,4:T>C",
	  orient =>   "-" },

	{ ref   => [ "ACGTTCGT" ],
	  reads =>   "GTTC",
	  args  => [ "-v 0",
	             "-n 0" ],
	  hits => { 2 => 1 } },
	
	{ ref   => [ "AAACGAAAGCTTTTATAGATGGGG" ],
	  reads =>      "132002320003332231",
	  args  => [ "-C -v 0",
	             "-C -n 0",
				 "-C -n 1",
				 "-C -v 1",
				 "-C -v 2",
				 "-C -n 2" ],
	  hits => { 3 => 1 },
	  color => 1 },
	{ ref   => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	  reads =>       "132002320113332231",
	  args  => [ "-C -v 2",
				 "-C -n 2" ],
	  hits => { 4 => 1 },
	  color => 1 },
	{ ref   => [ "AAAACGAAAGCTTTTATAGATGGGG" ],
	  reads =>       "132002320113332231",
	  args  => [ "-C -v 1",
				 "-C -n 1" ],
	  hits => { },
	  color => 1 },
	{ ref   => [ "ATATATGTCGACATATATATATATAT" ],
	  reads =>       "3311232113333333",
	  args  => [ "-C -v 0",
				 "-C -n 0" ],
	  hits => { 4 => 1 },
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
# Run bowtie with given arguments
#
sub runBowtie {
	my ($run_prog,
        $build_prog,
		$args,
		$color,
		$fa,
		$reportargs,
		$reads,
		$mate1s,
		$mate2s,
		$ls,
		$rawls) = @_;
	$args .= " --quiet";
	$reportargs = $reportargs || "-a";
	$args .= " $reportargs";
	# Write the reference to a fasta file
	my $build_args = ($color ? "-C" : "");
	my $cmd = "$build_prog --quiet $build_args $fa .simple_tests.tmp";
	print "$cmd\n";
	system($cmd);
	($? == 0) || die "Bad exitlevel from bowtie-build: $?";
	my $pe = (defined($mate1s) && $mate1s ne "");
	if($pe) {
		# Paired-end case
		$cmd = "$run_prog $args .simple_tests.tmp -1 $mate1s -2 $mate2s";
		print "$cmd\n";
		open(BT, "$cmd |") || die "Could not open pipe '$cmd |'";
		while(<BT>) {
			my $m1 = $_;
			my $m2 = <BT>;
			defined($m2) || die;
			print $m1;
			print $m2;
			chomp($m1);
			chomp($m2);
			push @$ls,    [ split(/\t/, $m1, -1) ];
			push @$ls,    [ split(/\t/, $m2, -1) ];
			push @$rawls, $m1;
			push @$rawls, $m2;
		}
		close(BT);
	} else {
		# Unpaired case
		$cmd = "$run_prog $args .simple_tests.tmp $reads";
		print "$cmd\n";
		open(BT, "$cmd |") || die "Could not open pipe '$cmd |'";
		while(<BT>) {
			print $_;
			chomp;
			push @$ls,    [ split(/\t/, $_, -1) ];
			push @$rawls, $_;
		}
		close(BT);
	}
	($? == 0) || die "bowtie exited with level $?\n";
}

my $tmpfafn = ".simple_tests.pl.fa";
for my $c (@cases) {
    while( my ($run_prg, $bld_prg) = each(%prog_pairs)){
	   writeFasta($c->{ref}, $tmpfafn);
	   # For each set of arguments...
	   for my $a (@{$c->{args}}) {
		   # Run bowtie
		   my @lines = ();
		   my @rawlines = ();
		   my %hits = ();
		   %hits = %{$c->{hits}} if defined($c->{hits});
		   my %pairhits = ();
		   %pairhits = %{$c->{pairhits}} if defined($c->{pairhits});
		   print $c->{name}."\n" if defined($c->{name});
		   my $color = 0;
		   $color = $c->{color} if defined($c->{color});
		   runBowtie(
               $run_prg,
               $bld_prg,
			   "$a -c",
			   $color,
			   $tmpfafn,
			   $c->{report},
			   $c->{reads},
			   $c->{mate1s},
			   $c->{mate2s},
			   \@lines,
			   \@rawlines);
		   my $pe = defined($c->{mate1s}) && $c->{mate1s} ne "";
		   my ($lastchr, $lastoff) = ("", -1);
		   for(my $li = 0; $li < scalar(@lines); $li++) {
			   my $l = $lines[$li];
			   scalar(@$l) == 8 || die "Bad number of fields; expected 8 got ".scalar(@$l).":\n$rawlines[$li]\n";
			   next if $l->[1] eq '*';
			   my ($chr, $off) = ($l->[0], $l->[3]);
			   if($pe && $lastchr ne "") {
				   my $offkey = min($lastoff, $off).",".max($lastoff, $off);
				   defined($pairhits{$offkey}) || die "No such paired off as $offkey in pairhits list: ".%{$c->{pairhits}}."\n";
				   $pairhits{$offkey}--;
				   delete $pairhits{$offkey} if $pairhits{$offkey} == 0;
				   ($lastchr, $lastoff) = ("", -1);
			   } elsif($pe) {
				   ($lastchr, $lastoff) = ($chr, $off);
			   } else {
				   defined($hits{$off}) || die "No such off as $off in hits list: ".%{$c->{hits}}."\n";
				   $hits{$off}--;
				   delete $hits{$off} if $hits{$off} == 0;
			   }
			   my $eds = $l->[-1];
			   !defined($c->{edits})  || $eds eq $c->{edits}  || die "For edit string, expected \"$c->{edits}\" got \"$eds\"\n";
		   }
		   my $hitsLeft = scalar(keys %hits);
		   $hitsLeft == 0 || die "Had $hitsLeft hit(s) left over";
		   my $pairhitsLeft = scalar(keys %pairhits);
		   $pairhitsLeft == 0 || die "Had $pairhitsLeft hit(s) left over";
	   }
   }
}
print "PASSED\n";
