#!/usr/bin/perl -w

# Input looks like this:
#
# Queries that hit: 6517315 (72.8%)
#    with at least 1 exact hit: 5524236 (61.7%)
#       1 exact hit: 4018084 (44.9%)
#         ...
#    with at least 1 1-mismatch hit (and no exact hits): 993079 (11.1%)
#       1 1-mismatch hit: 711113 (7.9%)
#         ...

my @trims = (25, 30, 35, 40);
my $reads = 8956597;
my $asPercent = 1;

open(EXACT,    ">exact.dat")    || die "Could not open >exact.dat";
open(EXACTU,   ">exactu.dat")   || die "Could not open >exactu.dat";
open(INEXACT,  ">inexact.dat")  || die "Could not open >inexact.dat";
open(INEXACTU, ">inexactu.dat") || die "Could not open >inexactu.dat";
open(MAQ1,     ">maq1.dat")     || die "Could not open >maq1.dat";
open(MAQ2,     ">maq2.dat")     || die "Could not open >maq2.dat";

for my $t (@trims) {
    open(STATS, "ebwt.$t.1tfra.arrows.stats") || die "Could not open stats ebwt.$t.1tfra.arrows.stats";
    print "$t bp: \n";
    print "  # reads: $reads\n";
    my $add = 0;
    while(<STATS>) {
	if(/  1 exact hit: ([0-9]+)/) {
	    my $foo = int($1);
	    $foo = $foo * 100.0 / $reads if $asPercent;
	    $foo /= 1000000 unless $asPercent;
	    print EXACTU "$t\t$foo\n";
	    print "  Unique exact hits: $1\n";
	}
	elsif(/with at least 1 exact hit: ([0-9]+)/) {
	    my $foo = int($1);
	    $add = $foo;
	    $foo = $foo * 100.0 / $reads if $asPercent;
	    $foo /= 1000000 unless $asPercent;
	    print EXACT "$t\t$foo\n";
	    print "  Exact hits: $1\n";
	}
	elsif(/1 1-mismatch hit: ([0-9]+)/) {
	    $add > 0 || 
		die "Expected to see \"with at least 1 exact hit\" before now";
	    my $foo = int($1) + $add;
	    $foo = $foo * 100.0 / $reads if $asPercent;
	    $foo /= 1000000 unless $asPercent;
	    print INEXACTU "$t\t$foo\n";
	    print "  Unique 1-mismatch hits: $1\n";
	}
	elsif(/at least 1 1-mismatch hit \(and no exact hits\): ([0-9]+)/) {
	    $add > 0 || 
		die "Expected to see \"with at least 1 exact hit\" before now";
	    my $foo = int($1) + $add;
	    $foo = $foo * 100.0 / $reads if $asPercent;
	    $foo /= 1000000 unless $asPercent;
	    print INEXACT "$t\t$foo\n";
	    print "  1-mismatch hits: $1\n";
	}
    } # while(<STATS>)
    {
	    my $hits = `wc -l ebwt.$t.maq.n1.tfr.hits | awk '{print \$1}'`;
		my $foo = $hits * 100.0 / $reads if $asPercent;
		$foo /= 1000000 unless $asPercent;
		print MAQ1 "$t\t$foo\n";
		print "  Other hits: $hits\n";
    }
    {
	    my $hits = `wc -l ebwt.$t.maq.tfr.hits | awk '{print \$1}'`;
		my $foo = $hits * 100.0 / $reads if $asPercent;
		$foo /= 1000000 unless $asPercent;
		print MAQ2 "$t\t$foo\n";
		print "  Other hits: $hits\n";
    }
}

close(EXACT);
close(EXACTU);
close(INEXACT);
close(INEXACTU);
close(MAQ1);
close(MAQ2);
