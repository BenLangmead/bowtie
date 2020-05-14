#ifndef FORMATS_H_
#define FORMATS_H_

#include <iostream>

/**
 * File-format constants and names
 */

enum file_format {
	FASTA = 1,
	FASTA_CONT,
	FASTQ,
	INTERLEAVED,
	TAB_MATE,
	RAW,
	CMDLINE,
	INPUT_CHAIN,
	RANDOM
};

static const std::string file_format_names[] = {
	"Invalid!",
	"FASTA",
	"FASTA sampling",
	"FASTQ",
	"Tabbed mated",
	"Raw",
	"Command line",
	"Chained",
	"Random"
};

#endif /*FORMATS_H_*/
