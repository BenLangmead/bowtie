#ifndef THREADING_H_
#define THREADING_H_

#include <iostream>
#include "tinythread.h"
#include "fast_mutex.h"

#ifdef NO_SPINLOCK
#   define MUTEX_T tthread::mutex
#else
#   define MUTEX_T tthread::fast_mutex
#endif /* NO_SPINLOCK */


/**
 * Wrap a lock; obtain lock upon construction, release upon destruction.
 */
class ThreadSafe {
public:
    ThreadSafe(MUTEX_T* ptr_mutex, bool locked = true) {
		if(locked) {
		    this->ptr_mutex = ptr_mutex;
		    ptr_mutex->lock();
		}
		else
		    this->ptr_mutex = NULL;
	}

	~ThreadSafe() {
	    if (ptr_mutex != NULL)
	        ptr_mutex->unlock();
	}
    
private:
	MUTEX_T *ptr_mutex;
};

#endif

