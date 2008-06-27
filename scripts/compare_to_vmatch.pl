#/usr/bin/perl -w

#
# Compare a (non-binary) Ebwt hit list against a Vmatch hit list and
# report any Ebwt hits that may be false positives or false negatives.
# Assumes that Ebwt hit list is the result of a "pick-1-random-hit"
# run, as opposed to a "report-all-hits" run.  The Vmatch output need
# not go past the "query sequence number" column.
#
# Usage: perl compare_to_vmatch.pl <vmatch_hits_file> <ebwt_hits_file>
#

use strict;
use warnings;
use Getopt::Std;

my %options=();
getopts("ovl",\%options);
my $ebwtOneHit = 0; $ebwtOneHit = $options{o} if defined $options{o};
my $verbose = 0; $verbose = $options{o} if defined $options{v};
my $lax = 0; $lax = $options{l} if defined $options{l};

print "ebwtOneHit: $ebwtOneHit\n";
print "lax: $lax\n";

my %vhitsFw = ();
my %vhitsRc = ();

defined($ARGV[0]) || die "Must specify vmatch hit file as first argument";
defined($ARGV[1]) || die "Must specify ebwt_search hit file as second argument";

my $vmatchHitsFile = $ARGV[0];
my $ebwtHitsFile   = $ARGV[1];
my $upto = -1;
$upto = $ARGV[2] if defined($ARGV[2]);
my $problems = 0;

my $numVmatchHits = 0;
my $numEbwtHits = 0;

if($vmatchHitsFile =~ /\.gz$/) {
    open(VMATCH, "zcat $vmatchHitsFile |") || die "Could not open gzipped vmatch hit file $vmatchHitsFile";
} else {
    open(VMATCH, $vmatchHitsFile) || die "Could not open vmatch hit file $vmatchHitsFile";
}
if($ebwtHitsFile =~ /\.gz$/) {
    open(EBWT, "zcat $ebwtHitsFile |") || die "Could not open gzipped ebwt hit file $ebwtHitsFile";
} else {
    open(EBWT, $ebwtHitsFile) || die "Could not open ebwt hit file $ebwtHitsFile";
}

# Read in vmatch hits and install in the global hashes
print "Reading Vmatch hits into hashes...\n";
while(<VMATCH>) {
	my @s = split;
	next if /^#/;
	my $vhits = ($s[3] eq "D") ? \%vhitsFw : \%vhitsRc;
	my $q = int($s[5]);
	unless(defined $vhits->{$q}) {
		$vhits->{$q} = {};
	}
	my $mm = ($s[7] eq "-1")? 1 : 0;
	my $key = (int($s[1]) + int($s[2]) + $mm);
	print "Vmatch query: $q, key: $key\n" if $verbose;
	$numVmatchHits++;
	$vhits->{$q}->{$key} = 1;
}

# Read in ebwt hits and compare against the global hashes
my $l = 0;
print "Checking Ebwt hits against Vmatch hits...\n";
while(<EBWT>) {
	next if /^Time/;
	chomp;
	my @s = split(/:/);
	my $q = $s[0];
	my $fw = (chop($q) eq '+') ? 1 : 0;
	my $vhits = ($fw == 1)? \%vhitsFw : \%vhitsRc;
	$q = int($q); # convert query index to int
	print "Ebwt query: $q\n" if $verbose;
	# Immediately check whether there are *any* hits for this same
	# query in the Vmatch output
	unless(defined($vhitsFw{$q}) || defined($vhitsRc{$q})) {
		print "No Vmatch hits for query: $q";
		if($fw) { print "+"; } else { print "-"; }
		print "\n  Ebwt output line: $_\n";
		$problems++;
		my $hstr = substr($s[1], 1);
		chop($hstr);
		my @hs = split(/>,</, $hstr);
		$numEbwtHits += ($#hs + 1);
		last if(++$l == $upto);
		next;
	}
	# Grab hit list and omit initial < and trailing >
	my $hstr = substr($s[1], 1);
	chop($hstr);
	# Split according to ">,<" ; if there's just one result, the split
	# will leave it alone
	my @hs = split(/>,</, $hstr);
	$numEbwtHits += ($#hs + 1);
	my @hInts = ();
	my @mms = ();
	# Turn the array of "index,offset" strings into integers equal to
	# int(index)+int(offset)+int(mms)
	for my $h (@hs) {
		my @tri = split(/,/, $h);
		print "  Try: @tri\n" if $verbose;
		my $r = int($tri[0]);
		my $o = int($tri[1]);
		my $mm = 0;
		$mm = int($tri[2]) if defined($tri[2]); # 0 or 1
		$mm == 0 || $mm == 1 || die "Bad mm: $mm";
		print "  r: $r, o: $o, mm: $mm\n" if $verbose;
		my $key = $r + $o + $mm;
		print "Ebwt query: $q, key: $key\n" if $verbose;
		push(@hInts, $key);
		push(@mms, $mm);
	}
	if($ebwtOneHit) {
		# One-Ebwt-against-all-Vmatch
		for(my $i = 0; $i <= $#hInts; $i++) {
			my $hInt = $hInts[$i];
			print "  Trying hInt: $hInt\n" if $verbose;
			my $hStr = $hs[$i];
			my $mm = $mms[$i];
			my $found = 0;
			foreach my $k (keys %{$vhits->{$q}}) {
			print "    comparing to vhits key: $k\n" if $verbose;
				if($k == $hInt) {
					print "    found!\n" if $verbose;
					$found = 1;
					# Delete the key since it's "covered" by an Ebwt
					# hit.  This same strategy does not work if we're
					# not in "oneHit" mode, since there may be another
					# set of Ebwt hits for this query later on (when
					# we process the transposed index).
					delete $vhitsFw{$q} if defined($vhitsFw{$q});
					delete $vhitsRc{$q} if defined($vhitsRc{$q});
					last;
				} else {
					print "    not found\n" if $verbose;
				}
			}
			if($lax != 0 && $found == 0) { # && $mm > 0) {
				# Count this as a match
				delete $vhitsFw{$q} if defined($vhitsFw{$q});
				delete $vhitsRc{$q} if defined($vhitsRc{$q});
			} elsif($found == 0) {
				print "No matching Vmatch hit for Ebwt hit: $q";
				if($fw) { print "+"; } else { print "-"; }
				print ":<$hStr>\n";
				$problems++;
			}
		}
	} else {
		# All-against-all
		my @vmatchHits = sort keys %{$vhits->{$q}};
		my @ebwtHits = sort @hInts;
		my $vi = 0;
		my $ei = 0;
		while($ei <= $#ebwtHits && $vi <= $#vmatchHits) {
			if($vmatchHits[$vi] == $ebwtHits[$ei]) {
				delete $vhits->{$q}->{$vmatchHits[$vi]};
				$vi++; $ei++;
			}
			# Unmatched Ebwt hits are a problem, but unmatched Vmatch hits
			# are not, since they may be matched later on when we process
			# the transposed index.
			elsif($vmatchHits[$vi] < $ebwtHits[$ei]) {
				$vi++; # This is OK
			} else {
				# Not OK; Ebwt has hit that vmatch doesn't have
				print "Unmatched Ebwt hit\n";
				$ei++; $problems++;
			}
		}
		# Unmatched Ebwt hits are a problem, but unmatched Vmatch hits
		# are not, since they may be matched later on when we process
		# the transposed index.
		for(; $ei <= $#ebwtHits; $ei++) {
			# Problem: unmatched Ebwt hit
			print "Unmatched Ebwt hit\n";
			$problems++;
		}
	}
	last if(++$l == $upto);
}
# Any remaining entries in vhitsFw and vhitsRc correspond to Vmatch
# hits that were not matched up with Ebwt hits
foreach my $q (keys %vhitsFw) {
	print "FW Vmatch query that wasn't covered by Ebwt: $q+:\n";
	$problems++;
}
foreach my $q (keys %vhitsRc) {
	print "RC Vmatch query that wasn't covered by Ebwt: $q-:\n";
	$problems++;
}
print "Found $problems discrepencies\n";
print "Total Ebwt hits: $numEbwtHits\n";
print "Total Vmatch hits: $numVmatchHits\n";
