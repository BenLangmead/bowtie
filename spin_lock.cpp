#include <atomic>

class spin_lock {
	std::atomic_flag flag;
public:
	spin_lock() : flag(false) {}

	void lock() {
		while (flag.test_and_set(std::memory_order_acquire));
	}

	void unlock() {
		flag.clear(std::memory_order_release);
	}
};
