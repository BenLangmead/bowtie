#ifndef THREADING_H_
#define THREADING_H_

#include <iostream>
#ifdef WITH_TBB
# include <tbb/mutex.h>
# include <tbb/spin_mutex.h>
# include <tbb/task_group.h>
#else
# include "tinythread.h"
# include "fast_mutex.h"
#endif

#ifdef NO_SPINLOCK
# ifdef WITH_TBB
#   define MUTEX_T tbb::mutex
# else
#   define MUTEX_T tthread::mutex
# endif
#else
# ifdef WITH_TBB
#  	define MUTEX_T tbb::spin_mutex
# else
#  	define MUTEX_T tthread::fast_mutex
# endif
#endif /* NO_SPINLOCK */

#ifdef WITH_TBB
# define GUARD_LOCK(lock) tbb::spin_mutex::scoped_lock tbb_guard(lock)
#else
# define GUARD_LOCK(lock) tthread::lock_guard<MUTEX_T> tt_guard(lock)
#endif

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

