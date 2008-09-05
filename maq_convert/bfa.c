#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "bfa.h"

nst_bfa1_t *nst_new_bfa1()
{
	nst_bfa1_t *bfa1;
	bfa1 = (nst_bfa1_t*)malloc(sizeof(nst_bfa1_t));
	bfa1->name = 0;
	bfa1->seq = bfa1->mask = 0;
	bfa1->ori_len = bfa1->len = 0;
	return bfa1;
}
void nst_delete_bfa1(nst_bfa1_t *bfa1)
{
	if (bfa1 == 0) return;
	free(bfa1->name);
	free(bfa1->seq);
	free(bfa1->mask);
	free(bfa1);
}
nst_bfa1_t *nst_load_bfa1(FILE *fp)
{
	int len;
	nst_bfa1_t *bfa1;
	if (fread(&len, sizeof(int), 1, fp) == 0) return 0;
	bfa1 = nst_new_bfa1();
	bfa1->name = (char*)malloc(sizeof(char) * len);
	fread(bfa1->name, sizeof(char), len, fp);
	fread(&bfa1->ori_len, sizeof(int), 1, fp);
	fread(&bfa1->len, sizeof(int), 1, fp);
	bfa1->seq = (bit64_t*)malloc(sizeof(bit64_t) * bfa1->len);
	fread(bfa1->seq, sizeof(bit64_t), bfa1->len, fp);
	bfa1->mask = (bit64_t*)malloc(sizeof(bit64_t) * bfa1->len);
	fread(bfa1->mask, sizeof(bit64_t), bfa1->len, fp);
	return bfa1;
}
nst_bfa_t *nst_new_bfa()
{
	return (nst_bfa_t*)calloc(1, sizeof(nst_bfa_t));
}
void nst_delete_bfa(nst_bfa_t * bfa)
{
	int i;
	for (i = 0; i != bfa->n; ++i)
		nst_delete_bfa1(bfa->bfa1[i]);
	free(bfa->bfa1);
	free(bfa);
}
nst_bfa_t *nst_load_bfa(FILE *fp)
{
	nst_bfa_t *bfa = nst_new_bfa();
	nst_bfa1_t *bfa1;
	int n = 0;
	while ((bfa1 = nst_load_bfa1(fp)) != 0) {
		bfa->bfa1 = (nst_bfa1_t**)realloc(bfa->bfa1, sizeof(nst_bfa1_t*) * (n + 1));
		bfa->bfa1[n++] = bfa1;
	}
	bfa->n = n;
	return bfa;
}
