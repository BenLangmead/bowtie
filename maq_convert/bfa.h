#ifndef NST_BFA_H
#define NST_BFA_H

#include <stdio.h>
#include "const.h"

typedef struct
{
	char *name;
	int ori_len, len;
	bit64_t *seq, *mask;
} nst_bfa1_t;

typedef struct
{
	int n;
	nst_bfa1_t **bfa1;
} nst_bfa_t;

#ifdef __cplusplus
extern "C" {
#endif
	nst_bfa1_t *nst_new_bfa1();
	void nst_delete_bfa1(nst_bfa1_t*);
	nst_bfa1_t *nst_load_bfa1(FILE *fp);
	nst_bfa_t *nst_new_bfa();
	void nst_delete_bfa(nst_bfa_t*);
	nst_bfa_t *nst_load_bfa(FILE *fp);
#ifdef __cplusplus
}
#endif

#endif
