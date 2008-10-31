#ifndef MAQMAP_H_
#define MAQMAP_H_

#define MAX_NAMELEN 36

#define MAQMAP_FORMAT_OLD 0
#define MAQMAP_FORMAT_NEW -1


#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "const.h"

/*
  name: read name
  size: the length of the read
  seq: read sequence (see also below)
  seq[MAX_READLEN-1]: single end mapping quality (equals to map_qual if not paired)
  map_qual: the final mapping quality
  alt_qual: the lower quality of the two ends (equals to map_qual if not paired)
  flag: status of the pair
  dist: offset of the mate (zero if not paired)
  info1: mismatches in the 24bp (higher 4 bits) and mismatches (lower 4 bits)
  info2: sum of errors of the best hit
  c[2]: count of all 0- and 1-mismatch hits on the reference
 */
template<int MAXLEN>
struct aln_t {
	uint8_t seq[MAXLEN];
	uint8_t size, map_qual, info1, info2, c[2], flag, alt_qual;
	uint32_t seqid, pos;
	int32_t dist;
	char name[MAX_NAMELEN];
};

template<int MAXLEN>
struct header_t {
	int32_t n_ref, format;
	char **ref_name;
	uint64_t n_mapped_reads;
	aln_t<MAXLEN> *mapped_reads;
};
//#define maqmap_read1(fp, m1) gzread((fp), (m1), sizeof(maqmap1_t))

template<int MAXLEN>
header_t<MAXLEN>* maq_init_header()
{
	header_t<MAXLEN>* mm = (header_t<MAXLEN>*)calloc(1, sizeof(header_t<MAXLEN>));
	mm->format = MAQMAP_FORMAT_NEW;
	return mm;
}

template<int MAXLEN>
void maq_destroy_header(header_t<MAXLEN>* mm)
{
	int i;
	if (mm == 0) return;
	for (i = 0; i < mm->n_ref; ++i)
		free(mm->ref_name[i]);
	free(mm->ref_name);
	free(mm);
}

template<int MAXLEN>
void maq_write_header(gzFile fp, const header_t<MAXLEN>* mm)
{
	int i, len;
	gzwrite(fp, &mm->format, sizeof(int));
	gzwrite(fp, &mm->n_ref, sizeof(int));
	for (i = 0; i != mm->n_ref; ++i) {
		len = strlen(mm->ref_name[i]) + 1;
		gzwrite(fp, &len, sizeof(int));
		gzwrite(fp, mm->ref_name[i], len);
	}
	gzwrite(fp, &mm->n_mapped_reads, sizeof(bit64_t));
}

template<int MAXLEN>
header_t<MAXLEN>* maq_read_header(gzFile fp)
{
	header_t<MAXLEN> *mm;
	int k, len;
	mm = maq_init_header<MAXLEN>();
	gzread(fp, &mm->format, sizeof(int));
	if (mm->format != MAQMAP_FORMAT_NEW) {
		fprintf(stderr, "Corrupted format. Abort!\n");
		abort();
	}
	gzread(fp, &mm->n_ref, sizeof(int));
	mm->ref_name = (char**)calloc(mm->n_ref, sizeof(char*));
	for (k = 0; k != mm->n_ref; ++k) {
		gzread(fp, &len, 4);
		mm->ref_name[k] = (char*)malloc(len);
		gzread(fp, mm->ref_name[k], len);
	}
	gzread(fp, &mm->n_mapped_reads, 8);
	return mm;
}


#endif

