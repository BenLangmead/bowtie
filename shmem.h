/*
 * shmem.h
 *
 *  Created on: Aug 13, 2009
 *      Author: Ben Langmead
 */

#ifndef SHMEM_H_
#define SHMEM_H_

#ifdef BOWTIE_SHARED_MEM

#include <string>
#include <sys/shm.h>

extern bool allocSharedMem(
	std::string fname,
    size_t len,
    void ** dst,
    const char *memName,
    bool verbose);

extern void notifySharedMem(void *mem, size_t len);

extern void waitSharedMem(void *mem, size_t len);

#define ALLOC_SHARED allocSharedMem
#define FREE_SHARED shmdt
#define NOTIFY_SHARED notifySharedMem
#define WAIT_SHARED waitSharedMem

#else

#define ALLOC_SHARED(...)
#define FREE_SHARED(...)
#define NOTIFY_SHARED(...)
#define WAIT_SHARED(...)

#endif /*BOWTIE_SHARED_MEM*/

#endif /* SHMEM_H_ */
