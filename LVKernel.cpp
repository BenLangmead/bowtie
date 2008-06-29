#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "LVKernel.h"

#define OS_X

using namespace std;

LVPyramid::LVPyramid(unsigned int Max_Diffs) : _max_diffs(Max_Diffs)
{
	// TODO: profile and check that this constructor isn't popping up
	// there are better ways to store the DP pyramid, but this is a
	// clear way.
	
	_cells.resize(_max_diffs + 2);
	for (unsigned int i = 0; i < _max_diffs + 2; ++i)
	{
		_cells[i] = vector<LVCell>(2 * (_max_diffs + 1) + 1);
	}
}

LVPyramid::~LVPyramid()
{
	
}

LVCell& LVPyramid::cell(int diffs, int offset)
{
	return _cells[diffs][offset];
	//return *(_cells + (_max_diffs * offset) + diffs);
}

/****************************************************************************/
/*                             Public Functions                             */
/****************************************************************************/

// Pilfered from googlecode:
// Find the first (least significant) set bit in a 64-bit integer.  The return
// value ranges from 0 (for no bit set), to 1 (for the least significant bit
// set), to 64 (for only the most significant bit set).
int find_64_lsm(uint64_t n)
{	
#if defined(OS_X) || defined(WINDOWS)
	n &= -n;
	int shift = (uint64_t) n <= 0xFFFFFFFFULL ? 0 : 32;
#endif
	
#if defined(LINUX)
	return ffsll(n);
#elif defined(OS_X)
	return ffs(n >> shift) + shift;
#elif defined(WINDOWS)
	return find_32(n >> shift) + shift;
#endif
}

// Find the first (most significant) set bit in a 64-bit integer.  The return
// value ranges from 0 (for no bit set), to 1 (for the least significant bit
// set), to 64 (for only the most significant bit set).
int find_64_msm(uint64_t n)
{	
	
	int pos = 0;
	uint64_t tmp;
	tmp = n >> 32;
	if (tmp != 0) { n = tmp; pos = pos + 32; }
	tmp = n >> 16;
	if (tmp != 0) { n = tmp; pos = pos + 16; }
	tmp = n >> 8;
	if (tmp != 0) { n = tmp; pos = pos + 8; }
	tmp = n >> 4;
	if (tmp != 0) { n = tmp; pos = pos + 4; }
	tmp = n >> 2;
	if (tmp != 0) { n = tmp; pos = pos + 2; }
	tmp = n >> 1;
	if (tmp != 0) { n = tmp; pos = pos + 1; }
	return pos + n - 1;
}


enum Dna { A, C, G, T };

DnaWord::DnaWord(const char* s, bool pack_left)
{
	len = strlen(s);
	assert (len <= 32);
	word = pack_dna_string(s, pack_left);
}

bool DnaWord::operator==(const DnaWord& rhs) const
{
	return (this->word == rhs.word) && (this->len == rhs.len);
}


// Packs a dna string into a 64 bit unsigned int.  The string dna must be 
// no more than 32 bp (excluding null terminator)
uint64_t pack_dna_string(const char* dna, bool pack_left)
{
	char c;
	uint64_t w = 0;
	unsigned int len = 0;
	while ((c = *dna++))
	{
		++len;
		c = toupper(c);
		switch (c)
		{
			case 'A': c = A; break;
			case 'C': c = C; break;
			case 'G': c = G; break;
			case 'T': c = T; break;
			default:  c = A; break;
		}
		w <<= 2;
		w |= c;
	}
	
	if (pack_left)
	{
		w <<= 64 - (len << 1);
	}
	return w;
}


int get_right_matching_chars(uint64_t w1, uint64_t w2)
{
	//find the least significant mismatching bit between w1 and w2
	int mismatch_bit = find_64_lsm(w1 ^ w2);
	
	if (!mismatch_bit)
		return -1;
	
	mismatch_bit -= 1;
	mismatch_bit -= ((mismatch_bit) & 1);
	mismatch_bit >>= 1;
	
	return mismatch_bit;
}

// Given two left-packed 64 bit packed dna strings,
// returns the number of bases that match between their
// left ends, or -1 if the strings are equal
int get_left_matching_chars(uint64_t w1, uint64_t w2)
{
	//find the most significant mismatching bit between w1 and w2
	int mismatch_bit = find_64_msm(w1 ^ w2);
	
	if (mismatch_bit < 0)
		return -1;
	
	int matching_bits = 64 - mismatch_bit;
	matching_bits -= ((matching_bits - 1) & 1);
	
	matching_bits >>= 1;
	
	return matching_bits;
}

// This routine assumes w1 != w2, and returns a pair (row, col) in the 
// pyramid where the alignment terminated.
void compute_pyramid(LVPyramid& py, 
					 const DnaWord& w1, 
					 const DnaWord& w2, 
					 bool left_extend,
					 int* row,
					 int* col,
					 int* w1_remaining)
{
	// diff_row tracks the number of edits so far in the alignment
	// since there are MAX_DIFFS + 1 rows, the number of actual edits is 
	// diff_row - 1 at any point in the alignment.
	unsigned int diff_row = 0;
	unsigned int max_diffs = py.max_diffs();
	// When shift_row changes, a gap is introduced. 
	unsigned int shift_col = max_diffs + 1;
	
	LVCell& start = py.cell(diff_row,shift_col);
	if (left_extend)
		start.match_chars = get_left_matching_chars(w1.word, w2.word);
	else
		start.match_chars = get_right_matching_chars(w1.word, w2.word);
	
	if (start.match_chars == -1)
		start.match_chars = min(w2.len, w1.len);
	else
	{
		start.match_chars = min(w2.len, start.match_chars);
		start.match_chars = min(w1.len, start.match_chars);
	}
	start.w1_shift = start.w2_shift = (start.match_chars << 1);
	unsigned int row_start = shift_col - 1;
	++diff_row;
	
	for (; diff_row <= max_diffs + 1; ++diff_row)
	{	
		shift_col = row_start--;
		for (; shift_col <= diff_row + max_diffs + 1; ++shift_col)
		{
			int total_matched_chars = py.cell(diff_row - 1,shift_col).match_chars;
			int w1_shift = py.cell(diff_row - 1,shift_col).w1_shift + 2;
			int w2_shift = py.cell(diff_row - 1,shift_col).w2_shift + 2;
			
			TracePtr back_ptr = UP;
			if (total_matched_chars < py.cell(diff_row - 1,shift_col - 1).match_chars)
			{
				// Gap in w2 (insertion in w1)
				total_matched_chars = py.cell(diff_row - 1,shift_col - 1).match_chars;
				w1_shift = py.cell(diff_row - 1,shift_col - 1).w1_shift + 2;
				w2_shift = py.cell(diff_row - 1,shift_col - 1).w2_shift;
				back_ptr = LEFT;
			}
			if (total_matched_chars < py.cell(diff_row - 1,shift_col + 1).match_chars)
			{
				// Gap in w1 (insertion in w2)
				total_matched_chars = py.cell(diff_row - 1,shift_col + 1).match_chars;
				w1_shift = py.cell(diff_row - 1,shift_col + 1).w1_shift;
				w2_shift = py.cell(diff_row - 1,shift_col + 1).w2_shift + 2;
				back_ptr = RIGHT;
			}
			
			uint64_t shifted_w1;
			uint64_t shifted_w2;
			
			if (!left_extend)
			{
				shifted_w1 = w1.word >> w1_shift;
				shifted_w2 = w2.word >> w2_shift;
			}
			else
			{
				shifted_w1 = w1.word << w1_shift;
				shifted_w2 = w2.word << w2_shift;	
			}
			int matching;
			int w2_chars_remaining = w2.len - (w2_shift >> 1);
			int w1_chars_remaining = w1.len - (w1_shift >> 1);
			
			if (w2_chars_remaining == 0 || w1_chars_remaining == 0)
			{
				matching = 0;
			}
			else if (shifted_w1 == shifted_w2 && 
					 w2_chars_remaining == w1_chars_remaining)
			{
				matching = w1.len - ((w1_shift >> (w1_shift & 1)) >> 1);
			}
			else
			{
				if (left_extend)
					matching = get_left_matching_chars(shifted_w1, shifted_w2);
				else
					matching = get_right_matching_chars(shifted_w1, shifted_w2);
				if (matching == -1)
					matching = min(w2_chars_remaining, w1_chars_remaining);
				else
				{
					matching = min(matching, w2_chars_remaining);
					matching = min(matching, w1_chars_remaining);
				}
			}
			
			total_matched_chars += matching;
			w1_shift += (matching << 1);
			w2_shift += (matching << 1);
			
			py.cell(diff_row,shift_col).match_chars = total_matched_chars;
			py.cell(diff_row,shift_col).w1_shift = w1_shift;
			py.cell(diff_row,shift_col).w2_shift = w2_shift;
			py.cell(diff_row,shift_col).back_ptr = back_ptr;
			
			w2_chars_remaining -= matching;
			w1_chars_remaining -= matching;
			
			if (w2_chars_remaining <= 0 /*&& w1_chars_remaining == 0*/ &&
				total_matched_chars + (int)(diff_row) >= min(w1.len, w2.len))
			{
				*row = diff_row;
				*col = shift_col;
				*w1_remaining = w1_chars_remaining;
				return;
			}
		}

	}
	
	*row = max_diffs + 1;
	*col = max_diffs + 1;
	*w1_remaining = -1;
}

void edit_distance(const DnaWord& w1,
				   const DnaWord& w2,
				   int max_diffs,
				   bool left_extend,
				   int* dist, 
				   int* w1_chars_remaining)
{
	if (w1 == w2)
	{
		*dist = 0;
		*w1_chars_remaining = 0;
		return;
	}
	
	// The edit distance between w1 and w2 is at least the difference in
	// their lengths
	if (abs(w1.len - w2.len) > max_diffs)
	{
		*dist = max_diffs + 1;
		*w1_chars_remaining = -1;
		return;
	}
	
	
	
	LVPyramid py(max_diffs);
	
	int col;
	*w1_chars_remaining = -1;
	
	compute_pyramid(py, w1, w2, left_extend, dist, &col, w1_chars_remaining);

	return;
}

void edit_distance(const char* s1, 
				  const char* s2, 
				  int max_diffs, 
				  bool left_extend,
				  int* dist, 
				  int* w1_chars_remaining)
{
	DnaWord w1(s1, left_extend);
	DnaWord w2(s2, left_extend);
	
	edit_distance(w1, w2, max_diffs, left_extend, dist, w1_chars_remaining);
}
