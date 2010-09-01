##
# inspect.pl
#
# Basic tests to ensure bowtie-inspect is behaving as expected.
#

use strict;
use warnings;

my $bowtie = "./bowtie";
my $bowtie_d = "${bowtie}-debug";
my $bowtie_build = "./bowtie-build";
my $bowtie_build_d = "${bowtie_build}-debug";
my $bowtie_inspect = "./bowtie-inspect";
my $bowtie_inspect_d = "${bowtie_inspect}-debug";

my @cases = (
	">\nN\n>\nATCTAG\n>\nN\n"
);

my $fn = ".inspect.pl.tmp.fa";
for my $c (@cases) {
	open(TMP, ">$fn") || die;
	print TMP $c;
	close(TMP);
	my $cmd = "$bowtie_build_d $fn $fn";
	system($cmd) == 0 || die "Exitlevel $? from command '$cmd'\n";
	$cmd = "$bowtie_inspect_d $fn $fn";
	system($cmd) == 0 || die "Exitlevel $? from command '$cmd'\n";
}
