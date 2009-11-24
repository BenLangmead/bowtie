/*
 * color_dec.h
 *
 *  Created on: Oct 14, 2009
 *      Author: Ben Langmead
 */

#ifndef COLOR_DEC_H_
#define COLOR_DEC_H_

#include <stdint.h>
#include <string>
#include <utility>
#include "alphabet.h"

void decodeHit(
		const char *read, // ASCII colors, '0', '1', '2', '3', '.'
		const char *qual, // ASCII quals, Phred+33 encoded
		size_t readi, // offset of first character within 'read' to consider
		size_t readf, // offset of last char (exclusive) in 'read' to consider
		const char *ref, // reference sequence, as masks
		size_t refi, // offset of first character within 'ref' to consider
		size_t reff, // offset of last char (exclusive) in 'ref' to consider
		int snpPhred, // penalty incurred by a SNP
		char *ns,  // decoded nucleotides are appended here
		char *cmm, // where the color mismatches are in the string
		char *nmm, // where nucleotide mismatches are in the string
		int& cmms, // number of color mismatches
		int& nmms);// number of nucleotide mismatches

#endif /* COLOR_DEC_H_ */
