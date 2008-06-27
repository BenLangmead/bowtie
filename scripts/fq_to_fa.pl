#/usr/bin/perl -w

# ------------------------------------------
# Symbol       Meaning      Nucleic Acid
# ------------------------------------------
# A            A           Adenine
# C            C           Cytosine
# G            G           Guanine
# T            T           Thymine
# U            U           Uracil
# M          A or C
# R          A or G
# W          A or T
# S          C or G
# Y          C or T
# K          G or T
# V        A or C or G
# H        A or C or T
# D        A or G or T
# B        C or G or T
# X      G or A or T or C
# N      G or A or T or C

open(FQ, "egrep '^[ACGTUMRWSYKVHDBXN]{40}\$' $ARGV[0] |") || die "Could not open .fq";
my $i = 0;
while(<FQ>) {
    print ">fakelab$i\n";
    print $_;
    $i++;
}
