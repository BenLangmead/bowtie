/*
 *  inexact_extend.cpp
 *  Bowtie
 *
 *  Created by Cole Trapnell on 6/18/08.
 *  
 */
#include <vector>
#include "inexact_extend.h"
#include "LVKernel.h"
#include "ebwt.h"

using namespace std;

DnaWord pack_dna_string(String<Dna,Packed<> >& dna, bool pack_left)
{
	uint64_t w = 0;
	//cerr << "converting: " << dna << endl;
	if (length(dna) > 16)
	{
		w = host(dna)[0];
		w |= (uint64_t)(host(dna)[1]) << 32;
	}
	else
	{
		w = host(dna)[0];
	}
	
	//mask off extraneous bits 
	w <<= ((32 - length(dna)) << 1);
	w >>= ((32 - length(dna)) << 1);
	int l = length(dna);
	
	if (pack_left)
	{
		w <<= 64 - (length(dna) << 1);
	}
	
	return DnaWord(w,l);
	
}

bitset<max_read_bp> mismatching_bases(const DnaWord& w1, const DnaWord& w2, bool left_extend)
{
	bitset<max_read_bp> diffs = 0;
	int l = min(w2.len, w1.len);
	
	uint64_t w1_word = w1.word;
	uint64_t w2_word = w2.word;
	int w1_len = w1.len;
	int w2_len = w2.len;
	int shift = 0;
	while (shift < l)
	{
		int match_chars = 0;
		if (left_extend)
			match_chars = get_left_matching_chars(w1_word, w2_word);
		else
			match_chars = get_right_matching_chars(w1_word, w2_word);
		
		if (match_chars == -1)
		{
			match_chars = min(w2_len, w1_len);
			//diffs += abs(w2_len - w1_len);
			if (w2_len < w1_len)
			{
				for (int i = w2_len + 1; i < w1_len; ++i)
					 diffs.set(i);
			}
			else
			{
				 for (int i = w2_len + 1; i < w1_len; ++i)
					  diffs.set(i);
			}
			shift = l;
		}
		else
		{
			match_chars = min(w2_len, match_chars);
			match_chars = min(w1_len, match_chars);
			int shift_chars = (match_chars + 1);

			w2_len -= shift_chars;
			w1_len -= shift_chars;
			
			int shift_bits = shift_chars << 1;
			if (left_extend)
			{
				w1_word <<= (shift_bits);
				w2_word <<= (shift_bits);
			}
			else
			{
				w1_word >>= (shift_bits);
				w2_word >>= (shift_bits);
			}
			
			shift += shift_chars;
			diffs.set(shift - 1);
		}
	}
	return diffs;
}




