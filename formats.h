#ifndef FORMATS_H_
#define FORMATS_H_

/**
 * File-format constants and names
 */

typedef enum file_format {
	FASTA = 1,
	FASTQ,
	CMDLINE
};

static const std::string file_format_names[] = {
	"Invalid!",
	"FASTA",
	"FASTQ",
	"Command line"
};

#endif /*FORMATS_H_*/
