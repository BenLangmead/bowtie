#/usr/bin/perl -w

while(<>) {
	chomp;
	my $r = reverse($_);
	$r =~ tr/aAcCgGtT/tTgGcCaA/;
	print "$r\n";
}
