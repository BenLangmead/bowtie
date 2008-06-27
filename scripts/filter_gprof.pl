#!/usr/bin/perl -w

# Remove the excessive C++ cruft that makes it so difficult to read
# gprof results

use strict;
use warnings;

while(<>) {
	s/seqan:://g;
	s/std:://g;
	s/unsigned /u/g;
	s/, Alloc<void> //g;
	s/basic_string<char, char_traits<char>, allocator<char> >/string/g;
	s/, allocator<string > //g;
	s/SimpleType<uchar, _Dna>/Dna/g;
	s/basic_istream<char, char_traits<char> >/istream/g;
	s/Tag<TagFasta_>/Fasta/g;
	s/basic_ifstream<char, char_traits<char> >/ifstream/g;
	s/pair<uint, uint>/U32Pair/g;
	s/ >/>/g;
	s/, allocator<String<Dna>>//g;
	s/, allocator<Hit>//g;
	s/Fasta const/Fasta/g;
	s/Size<String<Dna>>::Type/size_t/g;
	s/Tag<TagExact_>/Exact/g;
	s/basic_ostream<char, char_traits<char>>/ostream/g;
	print $_;
}
