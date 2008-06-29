#ifndef LV_KERNEL_H
#define LV_KERNEL_H

#include <vector>
using std::vector;

uint64_t pack_dna_string(const char* s, bool pack_left);

/**
 * This structure stores a packed representation of a Dna sequence up to 63bp
 * long.  It may be left- or right-packed, allowing the use of shift-and ops
 * in higher-level alignment routines.
 */
struct DnaWord
{
	DnaWord(uint64_t w, int l) : word(w), len(l) {}
	DnaWord(const DnaWord& rhs) : word(rhs.word), len(rhs.len) {}
	DnaWord(const char* s, bool pack_left); 
	bool operator==(const DnaWord& rhs) const;
	uint64_t word;
	int len;
};

/** 
 * Computes the edit distance between w2 and the left (or right, if left_extend
 * is false) [0,len(w2)] characters of w1.  Returns the distance, and the 
 * number of characters leftover (the left or right overhang) from w1.
 */
void edit_distance(const DnaWord& w1, 
				   const DnaWord& w2, 
				   int max_diffs, 
				   bool left_extend,
				   int* dist, 
				   int* w2_chars_remaining);

// Convenience function for use with unpacked c-style dna strings
void edit_distance(const char* s1, 
				  const char* s2, 
				  int max_diffs, 
				  bool left_extend,
				  int* dist, 
				   int* w2_chars_remaining);

/*
 * A backpointer for use with the Landau-Vishkin dynamic programming table 
 */
enum TracePtr {UP, LEFT, RIGHT};

/*
 * A single cell from the Landau-Vishkin dynamic programming table
 */
struct LVCell
{
	LVCell() : match_chars(-1), w1_shift(0), w2_shift(0), back_ptr(UP) {}
	
	// How many total matching characters have been found, given the shifts
	// this cell in the DP table implies.
	int match_chars; 
	
	// The remains of w1 after the shifts this cell implies
	unsigned int w1_shift;
	
	// The remains of w1 after the shifts this cell implies
	unsigned int w2_shift;
	
	// Backpointer to the cell we came from in order to get to this cell
	TracePtr back_ptr;
};

/**
 * The dynamic programming table used by the Landau-Vishking algorithm for 
 * finding the best alignment with up to k differences between two strings
 * We can use this for seed-and-extend, if we want to allow indels in the 
 * extension.
 */
class LVPyramid
{
private:
	vector<vector<LVCell> > _cells;
	//LVCell* _cells;
	unsigned int _max_diffs;
public:
	LVPyramid(unsigned int Max_Diffs);
	~LVPyramid();
	LVCell& cell(int diffs, int offset);
	unsigned int max_diffs() { return _max_diffs; }
};

/**
 * Given two right-packed 64 bit packed dna strings,
 * returns the number of bases that match between their
 * right ends, or -1 if the strings are equal
 */
int get_right_matching_chars(uint64_t w1, uint64_t w2);

/**
 * Given two left-packed 64 bit packed dna strings,
 * returns the number of bases that match between their
 * left ends, or -1 if the strings are equal
 */
int get_left_matching_chars(uint64_t w1, uint64_t w2);

/**
 * Compute the results of the Landau-Vishkin k-difference algorithm
 */
void compute_pyramid(LVPyramid& py, 
					 const DnaWord& w1, 
					 const DnaWord& w2, 
					 bool left_extend,
					 int* row,
					 int* col,
					 int* w2_remaining);


#endif /*LV_KERNEL_H*/
