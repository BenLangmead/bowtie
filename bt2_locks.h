#ifndef __MCSLOCK_H__
#define __MCSLOCK_H__

#include <atomic>
#include <sched.h>

class cpu_backoff {
public:
	cpu_backoff(): count(1) {}
	void pause() {
		if (count <= LOOPS_BEFORE_YIELD) {
			for (int32_t i = 0; i < count; i++)
				__asm__ __volatile__("pause;");
			count *= 2;
		} else {
			sched_yield();
		}
	}
private:
	static const int32_t LOOPS_BEFORE_YIELD = 16;
	int32_t count;
};

class spin_lock {
	std::atomic_flag flag;
public:
	spin_lock() : flag(false) {}

	void lock();
	void unlock();
};

class mcs_lock {
public:
	mcs_lock(): q(nullptr) {}
        struct mcs_node {
		mcs_node *next;
		std::atomic_bool unlocked;
        };

	void lock();
	void unlock();
	typedef mcs_node* mcs_node_ptr;
private:
	void spin_while_eq(const volatile mcs_node_ptr& value, mcs_node *expected) {
		cpu_backoff backoff;
		while (value == expected)
			backoff.pause();
	}

	void spin_while_eq(const volatile std::atomic_bool& value, bool expected) {
	 	cpu_backoff backoff;
	 	while (value.load(std::memory_order_acquire) == expected)
	 		backoff.pause();
	}
        std::atomic<mcs_node *> q;
        static thread_local mcs_node node;
};
#endif
