/*
 * auto_array.h
 *
 *  Created on: Oct 12, 2009
 *      Author: Ben Langmead
 */

#include <cstring>

#ifndef AUTO_ARRAY_H_
#define AUTO_ARRAY_H_

/**
 * A simple fixed-length array of type T, automatically freed in the
 * destructor.
 */
template<typename T>
class AutoArray {
public:
	AutoArray(size_t sz) {
		t_ = NULL;
		t_ = new T[sz];
		memset(t_, 0, sz*sizeof(T));
		sz_ = sz;
	}
	~AutoArray() { if(t_ != NULL) delete[] t_; }
	T& operator[](size_t sz) {
		return t_[sz];
	}
	const T& operator[](size_t sz) const {
		return t_[sz];
	}
private:
	T *t_;
	size_t sz_;
};

#endif /* AUTO_ARRAY_H_ */
