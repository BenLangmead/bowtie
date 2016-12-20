#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include "threading.h"

class SyncLogger {
public:
	SyncLogger() : mutex_m() { }

	void msg(const char *s) {
		ThreadSafe ts(&mutex_m);
		std::cout << s << std::endl;
	}

	void msg(const std::string& s) {
		ThreadSafe ts(&mutex_m);
		std::cout << s << std::endl;
	}

private:
	MUTEX_T mutex_m;
};

extern SyncLogger glog;

#endif /*LOG_H_*/
