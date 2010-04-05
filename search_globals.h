/*
 * search_globals.h
 *
 *  Created on: Dec 5, 2009
 *      Author: Ben Langmead
 */

#ifndef SEARCH_GLOBALS_H_
#define SEARCH_GLOBALS_H_

#include "threading.h"

// declared in ebwt_search.cpp
extern bool color;
extern bool colorExEnds;
extern bool colorSeq;
extern bool colorQual;
extern int  snpPhred;
extern bool showSeed;
extern bool quiet;

extern MUTEX_T gLock;

#endif /* SEARCH_GLOBALS_H_ */
