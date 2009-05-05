#include <iostream>
#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "bfa.h"

using namespace std;

nst_bfa1_t *nst_new_bfa1()
{
	nst_bfa1_t *bfa1;
	bfa1 = (nst_bfa1_t*)malloc(sizeof(nst_bfa1_t));
	if(bfa1 == NULL) {
		cerr << "Exhausted memory allocating space for the .bfa file" << endl;
		exit(1);
	}
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
static void bfa_read_error() {
	fprintf(stderr, "Error reading from .bfa file\n");
	exit(1);
}
nst_bfa1_t *nst_load_bfa1(FILE *fp)
{
	int len;
	nst_bfa1_t *bfa1;
	if (fread(&len, sizeof(int), 1, fp) == 0) return 0;
	bfa1 = nst_new_bfa1();
	bfa1->name = (char*)malloc(sizeof(char) * len);
	if(bfa1->name == NULL) {
		cerr << "Exhausted memory allocating space for the .bfa file name" << endl;
		exit(1);
	}
	/*
	 * BTL: I had to add in these return-value checks to keep gcc 4.3.2
	 * from complaining.
	 */
	if(fread(bfa1->name, sizeof(char), len, fp) != (size_t)len) {
		bfa_read_error();
	}
	if(fread(&bfa1->ori_len, sizeof(int), 1, fp) != 1) {
		bfa_read_error();
	}
	if(fread(&bfa1->len, sizeof(int), 1, fp) != 1) {
		bfa_read_error();
	}
	bfa1->seq = (bit64_t*)malloc(sizeof(bit64_t) * bfa1->len);
	if(bfa1->seq == NULL) {
		cerr << "Exhausted memory allocating space for the .bfa file sequence" << endl;
		exit(1);
	}
	if(fread(bfa1->seq, sizeof(bit64_t), bfa1->len, fp) != (size_t)bfa1->len) {
		bfa_read_error();
	}
	bfa1->mask = (bit64_t*)malloc(sizeof(bit64_t) * bfa1->len);
	if(bfa1->mask == NULL) {
		cerr << "Exhausted memory allocating space for the .bfa file mask" << endl;
		exit(1);
	}
	if(fread(bfa1->mask, sizeof(bit64_t), bfa1->len, fp) != (size_t)bfa1->len) {
		bfa_read_error();
	}
	return bfa1;
}
/*
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
*/
