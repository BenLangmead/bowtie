#ifndef __MCSLOCK_H__
#define __MCSLOCK_H__

#if (__cplusplus >= 201103L)

#include <atomic>
#include <thread>

class MCSMutex {
public:
        struct MCSNode {
                MCSNode *next;
                uintptr_t unlocked;
        };

        void lock() {
                node.next = nullptr;
                node.unlocked = 0;

                MCSNode *pred = q.exchange(&node, std::memory_order_acquire);
                if (pred != nullptr) {
                        pred->next = &node;
                        while (node.unlocked == 0)
                                ;
                }
        }

        void unlock() {
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
                node.next->unlocked = 1;
        }
private:
        std::atomic<MCSNode *> q{nullptr};
        static thread_local MCSNode node;
};

thread_local auto MCSMutex::node = MCSMutex::MCSNode();

#endif

#endif
