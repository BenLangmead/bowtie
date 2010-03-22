#!/usr/bin/perl -w

##
# mapability.pl
#
# Calculate mapability of each reference position assuming reference is
# bisulfite treated.
#

use strict;
use warnings;
use Getopt::Long;

my $fa = "";
my $win = 50;
my $freq = 1;
my $bowtie = "";
my $bowtie_arg = "";
my $pol = "-v 3";
my $fwidx = "";
my $rcidx = "";
my $btargs = "-t --norc -S --sam-nohead -M 1 --mm";
my $debug = 0;

if(defined($ENV{BOWTIE_HOME})) {
	$bowtie = "$ENV{BOWTIE_HOME}/bowtie";
	unless(-x $bowtie) { $bowtie = "" };
}
if($bowtie eq "") {
	$bowtie = `which bowtie 2>/dev/null`;
	chomp($bowtie);
	unless(-x $bowtie) { $bowtie = "" };
}
$bowtie = "./bowtie" if ($bowtie eq "" && -x "./bowtie");

GetOptions(
	"fasta=s"     => \$fa,
	"window=i"    => \$win,
	"frequency=i" => \$freq,
	"bowtie=s"    => \$bowtie_arg,
	"policy=s"    => \$pol,
	"fwidx=s"     => \$fwidx,
	"rcidx=s"     => \$rcidx,
	"debug"       => \$debug
) || die;

print STDERR "Bowtie: found: $bowtie; given: $bowtie_arg\n";
print STDERR "Input fasta: $fa\n";
print STDERR "FW index: $fwidx\n";
print STDERR "RC index: $rcidx\n";
print STDERR "Alignment policy: $pol\n";
print STDERR "Window size: $win\n";
print STDERR "Frequency: $freq\n";

$fa ne "" || die "Must specify -fasta\n";
$fwidx ne "" || die "Must specify -fwidx\n";
$rcidx ne "" || die "Must specify -rcidx\n";
-f "$fwidx.1.ebwt" || die "Could not find -fwidx index file $fwidx.1.ebwt\n";
-f "$rcidx.1.ebwt" || die "Could not find -rcidx index file $rcidx.1.ebwt\n";

$bowtie = $bowtie_arg if $bowtie_arg ne "";
unless(-x $bowtie) {
	# No bowtie? die
	if($bowtie_arg ne "") {
		die "Specified -bowtie, \"$bowtie\" doesn't exist or isn't executable\n";
	} else {
		die "bowtie couldn't be found in BOWTIE_HOME, PATH, or current directory; please specify -bowtie\n";
	}
}

my $running = 0;
my $name = ""; # name of sequence currently being processed
my %lens = ();
my @names = ();
my $totlen = 0;

##
# Read lengths of all the entries in all the input fasta files.
#
sub readLens($) {
	my $ins = shift;
	my @is = split(/[,]/, $ins);
	for my $i (@is) {
		open IN, "$i" || die "Could not open $i\n";
		my $name = "";
		while(<IN>) {
			chomp;
			if(substr($_, 0, 1) eq '>') {
				next if /\?[0-9]*$/; # Skip >?50000 lines
				$name = substr($_, 4); # Chop off >FW:/>RC:
				my @ns = split(/\s+/, $name);
				$name = $ns[0]; # Get short name
				push @names, $name;
				print STDERR "Saw name $name\n";
			} else {
				$name ne "" || die;
				$lens{$name} += length($_); # Update length
				$totlen += length($_);
			}
		}
		close(IN);
	}
}
print STDERR "Reading fasta lengths\n";
readLens($fa);
print STDERR "  read ".scalar(keys %lens)." sequences with total length $totlen\n";

my @last;
for(my $i = 0; $i < $win; $i++) { push @last, 0 };
sub clearLast {
	for(my $i = 0; $i < $win; $i++) { $last[$i] = 0 };
}

print STDERR "Opening bowtie pipes\n";
my $fwcmd = "$bowtie -F $win,$freq $btargs $pol $fwidx $fa";
my $rccmd = "$bowtie -F $win,$freq $btargs $pol $rcidx $fa";
print STDERR "Forward command: $fwcmd\n";
print STDERR "Reverse comp command: $rccmd\n";
open BTFW, "$fwcmd |" || die "Couldn't open pipe '$fwcmd |'\n";
open BTRC, "$rccmd |" || die "Couldn't open pipe '$rccmd |'\n";

# 10_554  4       *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      XM:i:0
# 10_555  4       *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      XM:i:0
# 10_556  4       *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      XM:i:0
# 10_557  4       *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      XM:i:0
# 10_558  4       *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      XM:i:0
# 10_559  4       *       0       0       *       *       0       0       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      XM:i:0

print STDERR "Reading...\n";
my $ln = 0;
my $cur = 0;
my $lastc = "\n";
while(1) {
	my $fwl = <BTFW>;
	my $rcl = <BTRC>;
	
	last unless defined($fwl) && defined($rcl);
	$ln++;

	# TODO: should probably be a space instead of an underscore
	my @fws = split(/\t/, $fwl);
	my @fws1 = split(/_/, $fws[0]);
	my @rcs = split(/\t/, $rcl);
	my @rcs1 = split(/_/, $rcs[0]);
	
	my $cname = join("_", @fws1[0..$#fws1-1]);
	$cname eq join("_", @fws1[0..$#fws1-1]) || die "Name mismatch on line $ln:\n$fwl\n$rcl\n";
	my $off = $fws1[-1];
	$off eq $rcs1[-1] || die "Offset mismatch on line $ln:\n$fwl\n$rcl\n";
	my $mappedFw = ($fws[-1] =~ /XM:i/ ? 0 : 1);
	my $mappedRc = ($rcs[-1] =~ /XM:i/ ? 0 : 1);
	my $mapable = $mappedFw != $mappedRc;
	
	if($name ne $cname) {
		if($name ne "") {
			# Flush remaining characters from previous name
			if($debug) {
				print STDERR "Read $cur characters of sequence $name, with lenth ".$lens{$name}."\n";
				print STDERR "Flushing...\n";
				my $tmp = <STDIN>;
			}
			defined($lens{$name}) || die;
			for(; $cur < $lens{$name}; $cur++) {
				$running -= $last[$cur % $win];
				$last[$cur % $win] = 0;
				$lastc = chr($running + 64);
				print $lastc;
				if((($cur+1) % 60) == 0) {
					$lastc = "\n";
					print $lastc;
				}
			}
		}
		$name = $cname;
		defined($lens{$name}) || die;
		print "\n" unless $lastc eq "\n";
		print ">$name\n";
		$lastc = "\n";
		$cur = 0;
		$running = 0;
		clearLast();
	}
	
	if($debug) {
		#print STDERR "Read $cur characters of sequence $name, with lenth ".$lens{$name}."\n";
		#print STDERR "Flushing...\n";
		#my $tmp = <STDIN>;
	}
	
	$running -= $last[$cur % $win];
	$last[$cur % $win] = $mapable;
	$running += $last[$cur % $win];

	$running <= $win || die "running counter $running got higher than window size $win\n";
	$lastc = chr($running + 64);
	print $lastc;
	if((($cur+1) % 60) == 0) {
		$lastc = "\n";
		print $lastc;
	}

	$cur++;
}
defined($lens{$name}) || die;
for(; $cur < $lens{$name}; $cur++) {
	$running -= $last[$cur % $win];
	$last[$cur % $win] = 0;
	$lastc = chr($running + 64);
	print $lastc;
	if((($cur+1) % 60) == 0) {
		$lastc = "\n";
		print $lastc;
	}
}
close(BTFW);
$? == 0 || die "Bad exitlevel from forward bowtie: $?\n";
close(BTRC);
$? == 0 || die "Bad exitlevel from reverse-comp bowtie: $?\n";
