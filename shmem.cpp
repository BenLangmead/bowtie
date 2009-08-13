/*
 * shmem.cpp
 *
 *  Created on: August 13, 2009
 *      Author: Ben Langmead
 */

#ifdef BOWTIE_SHARED_MEM

#include <iostream>
#include <string>
#include <unistd.h>
#include <sys/shm.h>
#include <errno.h>
#include "str_util.h"

using namespace std;

/**
 * Tries to allocate a shared-memory chunk for a given file of a given size.
 */
bool allocSharedMem(string fname,
                    size_t len,
                    void ** dst,
                    const char *memName,
                    bool verbose)
{
	int shmid = -1;
	// Calculate key given string
	key_t key = (key_t)hash_string(fname);
	shmid_ds ds;
	int ret;
	// Reserve 4 bytes at the end for silly synchronization
	size_t shmemLen = len + 4;
	if(verbose) {
		cout << "Reading " << len << "+4 bytes into shared memory for " << memName << endl;
	}
	void *ptr = NULL;
	while(true) {
		// Create the shrared-memory block
		if((shmid = shmget(key, shmemLen, IPC_CREAT | 0666)) < 0) {
			if(errno == ENOMEM) {
				cerr << "Out of memory allocating shared area " << memName << endl;
			} else if(errno == EACCES) {
				cerr << "EACCES" << endl;
			} else if(errno == EEXIST) {
				cerr << "EEXIST" << endl;
			} else if(errno == EINVAL) {
				cerr << "Warning: shared-memory chunk's segment size doesn't match expected size (" << (shmemLen) << ")" << endl
					 << "Deleteing old shared memory block and trying again." << endl;
				shmid = shmget(key, 0, 0);
				if((ret = shmctl(shmid, IPC_RMID, &ds)) < 0) {
					cerr << "shmctl returned " << ret
						 << " for IPC_RMID, errno is " << errno
						 << ", shmid is " << shmid << endl;
					exit(1);
				} else {
					cerr << "Deleted shared mem chunk with shmid " << shmid << endl;
				}
				continue;
			} else if(errno == ENOENT) {
				cerr << "ENOENT" << endl;
			} else if(errno == ENOSPC) {
				cerr << "ENOSPC" << endl;
			} else {
				cerr << "shmget returned " << shmid << " for and errno is " << errno << endl;
			}
			exit(1);
		}
		ptr = shmat(shmid, 0, 0);
		if(ptr == (void*)-1) {
			cerr << "Failed to attach " << memName << " to shared memory with shmat()." << endl;
			exit(1);
		}
		if(ptr == NULL) {
			cerr << memName << " pointer returned by shmat() was NULL." << endl;
			exit(1);
		}
		// Did I create it, or did I just attach to one created by
		// another process?
		if((ret = shmctl(shmid, IPC_STAT, &ds)) < 0) {
			cerr << "shmctl returned " << ret << " for IPC_STAT and errno is " << errno << endl;
			exit(1);
		}
		if(ds.shm_segsz != shmemLen) {
			cerr << "Warning: shared-memory chunk's segment size (" << ds.shm_segsz
				 << ") doesn't match expected size (" << shmemLen << ")" << endl
				 << "Deleteing old shared memory block and trying again." << endl;
			if((ret = shmctl(shmid, IPC_RMID, &ds)) < 0) {
				cerr << "shmctl returned " << ret << " for IPC_RMID and errno is " << errno << endl;
				exit(1);
			}
		} else {
			break;
		}
	} // while(true)
	*dst = ptr;
	if(ds.shm_cpid == getpid()) {
		if(verbose) {
			cout << "  I (pid = " << getpid() << ") created the "
			     << "shared memory for " << memName << endl;
		}
		// Set this value just off the end of the chunk to
		// indicate that the data hasn't been read yet.
		((volatile uint32_t*)((char*)ptr + len))[0] = 0xffffffff;
		return true;
	} else {
		if(verbose) {
			cout << "  I (pid = " << getpid()
			     << ") did not create the shared memory for "
			     << memName << ".  Pid " << ds.shm_cpid << " did." << endl;
		}
		return false;
	}
}

/**
 * Notify other users of a shared-memory chunk that the leader has
 * finished initializing it.
 */
void notifySharedMem(void *mem, size_t len) {
	((volatile uint32_t*)((char*)mem + len))[0] = 0x0f0f0f0f;
}

/**
 * Wait until the leader of a shared-memory chunk has finished
 * initializing it.
 */
void waitSharedMem(void *mem, size_t len) {
	while(((volatile uint32_t*)((char*)mem + len))[0] != 0x0f0f0f0f) {
		// spin
	}
}

#endif
