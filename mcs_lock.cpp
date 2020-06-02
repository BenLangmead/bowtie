#include "mcs_lock.h"

thread_local MCSMutex::MCSNode MCSMutex::node;

void MCSMutex::lock() {
	node.next = nullptr;
	node.unlocked.store(false, std::memory_order_release);

	MCSNode *pred = q.exchange(&node, std::memory_order_acquire);
	if (pred != nullptr) {
		pred->next = &node;
		while (node.unlocked.load(std::memory_order_acquire) == false)
			;
	}
}

void MCSMutex::unlock() {
	if (node.next == nullptr) {
		MCSNode *node_ptr = &node;
		if (q.compare_exchange_strong(node_ptr,
					      nullptr,
					      std::memory_order_release,
					      std::memory_order_relaxed))
			return;
		while (node.next == nullptr)
			;
	}
	node.next->unlocked.store(true, std::memory_order_release);
}
