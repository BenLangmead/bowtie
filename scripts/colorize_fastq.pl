#!/usr/bin/perl -w

##
# colorize_fastq.pl
#
# Convert nucleotide FASTQ input to colorspace CSFASTQ output.
# Colorspace versions of nucleotide sequnences are 1 character shorter.
# Names are unchanged.  No primer base is given.
#

##
# Given a string in nucleotide space, convert to colorspace.
#
sub colorize($$) {
	my ($s, $nucs) = @_;
	defined($s) || die;
	my %cmap = (
		"AA" => "0", "CC" => "0", "GG" => "0", "TT" => "0",
		"AC" => "1", "CA" => "1", "GT" => "1", "TG" => "1",
		"AG" => "2", "GA" => "2", "CT" => "2", "TC" => "2",
		"AT" => "3", "TA" => "3", "CG" => "3", "GC" => "3",
		"NA" => ".", "NC" => ".", "NG" => ".", "NT" => ".",
		"AN" => ".", "CN" => ".", "GN" => ".", "TN" => ".",
		"NN" => "."
	);
	my %nmap = ("0" => "A", "1" => "C", "2" => "G", "3" => "T", "." => "N");
	my $ret = "";
	for(my $i = 0; $i < length($s)-1; $i++) {
		my $di = uc substr($s, $i, 2);
		$di =~ tr/-NnMmRrWwSsYyKkVvHhDdBbXx/N/;
		defined($cmap{$di}) || die "Bad dinuc: $di\n";
		$ret .= ($nucs ? $nmap{$cmap{$di}} : $cmap{$di});
	}
	return $ret;
}

while(<>) {
	next if /^[;#]/;
	my $name  = $_;
	my $seq   = <>;
	my $name2 = <>;
	my $qual  = <>;
	chomp($seq);
	print $name.colorize($seq, 0)."\n$name2$qual";
}
