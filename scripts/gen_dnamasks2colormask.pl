#!/usr/bin/perl

print "uint8_t dnamasks2colormask[16][16] = {\n";
print "\t         /* 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 */\n";
my %color = (
	"0,0" => 0, "0,1" => 1, "0,2" => 2, "0,3" => 3,
	"1,0" => 1, "1,1" => 0, "1,2" => 3, "1,3" => 2,
	"2,0" => 2, "2,1" => 3, "2,2" => 0, "2,3" => 1,
	"3,0" => 3, "3,1" => 2, "3,2" => 1, "3,3" => 0
);
for(my $i = 0; $i < 16; $i++) {
	printf "\t/* %2d */ { ", $i;
	for(my $j = 0; $j < 16; $j++) {
		my $mask = 0;
		for(my $x = 0 ; $x < 4; $x++) {
			for(my $y = 0 ; $y < 4; $y++) {
				if(($i & (1<<$x)) != 0 && ($j & (1<<$y)) != 0) {
					$mask |= (1 << $color{"$x,$y"});
				}
			}
		}
		$mask >= 0 || die;
		printf "%2d", $mask;
		if($j == 15) {
			print " }";
			print "," if $i < 15;
			print "\n";
		} else {
			print ", ";
		}
	}
}
print "};\n";
