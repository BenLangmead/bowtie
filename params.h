#ifndef PARAMS_H_
#define PARAMS_H_

typedef enum algorithm_type {
	SKEW7 = 1,
	SKEW3,
	LARSSON_SADAKANE,
	
};

typedef enum sequence_type {
	DNA = 1,
	DNA5,
	RNA,
	RNA5,
	IUPAC,
	AMINO_ACID
};

static const std::string sequence_type_names[] = {
	"Invalid!",
	"Dna",
	"Dna5",
	"RNA",
	"Rna5",
	"Iupac",
	"Protein"
};

typedef enum file_format {
	FASTA = 1,
	FASTQ,
	PACKED,
	EMBL,
	GENBANK,
	BFQ,
	SOLEXA,
	RAW,
	CMDLINE
};

static const std::string file_format_names[] = {
	"Invalid!",
	"FASTA",
	"FASTQ",
	"Packed",
	"Embl",
	"Genbank",
	"Bfq/Maq",
	"Solexa",
	"Raw",
	"Command line"
};

#endif /*PARAMS_H_*/
