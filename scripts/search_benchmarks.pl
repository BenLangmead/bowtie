#/usr/bin/perl -w

my $do_ebwts = 1;
my $do_kg = 0;
my $do_sim = 1;

# Do several different numbers of reads
my $ebwt = "/fs/szasmg/langmead/ebwts/whole.t10o5l6";
my @ebwts = (
	"first100000000",
	"first250000000",
	"first500000000",
	"first1000000000",
	"first2000000000",
);

my %reads = ('kg'  => "/fs/szasmg/langmead/reads/SRR001113/s_7_0000_0255.fa",
             'sim' => "/fs/szasmg/langmead/reads/human_7mil_mms_sim.fa");
my @readks = ();
push(@readks, 'kg') if $do_kg;
push(@readks, 'sim') if $do_sim;

my @us = ( 1000000,
           2000000,
           4000000,
           7000000,
          10000000,
          14000000);
my @mms = (0, 1);

if($do_ebwts) {
	for my $eb (@ebwts) {
		my $ebf = "/fs/szasmg/langmead/ebwts/whole." . $eb;
		print "*** Doing Ebwt $eb\n";
		open(OUT, ">$eb.t10o5l6u$u.kg.1.out");
		my $u = 4000000;
		my $cmd = "./ebwt_search --stats -u $u -1 -frt $ebf /fs/szasmg/langmead/hs_ref_$eb.mfa.reads $eb.t10o5l6u$u.kg.1.hits";
		print "$cmd\n";
		open(CMD, "$cmd 2>&1 |");
		while(<CMD>) {
			if(/Overall time: ([0-9:]+)/) {
				print "Overall time: $1\n";
			}
			print OUT $_;
		}
		close(OUT);
		close(CMD);
		# Done; now let's print fraction of reads mapped
		my $exactMapped = `grep -c '0>' $eb.t10o5l6u$u.kg.1.hits`;
		my $inexactMapped = `grep -c '1>' $eb.t10o5l6u$u.kg.1.hits`;
		print "Fraction of reads exact-mapped: " . ($exactMapped * 1.0 / ($u/2)) . "\n";
		print "Fraction of reads inexact-mapped: " . ($inexactMapped * 1.0 / ($u/2)) . "\n";
	}
}

for my $rk (@readks) {
	print "*** Doing readset $rk\n";
	my $readfile = $reads{$rk};
	for my $mm (@mms) {
		print "  *** Doing $mm mismatches\n";
		for my $u (@us) {
			my $mmarg = "";
			print "    *** Doing $u patterns\n";
			$mmarg = "-1" if $mm == 1;
			open(OUT, ">whole.t10o5l6u$u.$rk.$mm.out");
			my $cmd = "./ebwt_search --stats -u $u $mmarg -frt $ebwt $readfile whole.t10o5l6u$u.$rk.$mm.hits";
			print "$cmd\n";
			open(CMD, "$cmd 2>&1 |");
			while(<CMD>) {
				if(/Overall time: ([0-9:]+)/) {
					print "Overall time: $1\n";
				}
				print OUT $_;
			}
			close(OUT);
			close(CMD);
			# Done; now let's print fraction of reads mapped
			my $exactMapped = `grep -c '0>' whole.t10o5l6u$u.$rk.$mm.hits`;
			my $inexactMapped = `grep -c '1>' whole.t10o5l6u$u.$rk.$mm.hits`;
			print "Fraction of reads exact-mapped: " . ($exactMapped * 1.0 / ($u/2)) . "\n";
			print "Fraction of reads inexact-mapped: " . ($inexactMapped * 1.0 / ($u/2)) . "\n";
		}
	}
}
