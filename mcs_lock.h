#ifndef __MCSLOCK_H__
#define __MCSLOCK_H__

#include <atomic>
#include <iostream>

class MCSMutex {
public:
        struct MCSNode {
		MCSNode *next;
		std::atomic_bool unlocked;
        };

	void lock();
	void unlock();

private:
        std::atomic<MCSNode *> q{nullptr};
        static thread_local MCSNode node;
};

#endif
