/*
 * ebwt_search.h
 *
 *  Created on: Apr 7, 2015
 *      Author: vantones
 */

#ifndef EBWT_SEARCH_H_
#define EBWT_SEARCH_H_
#ifdef WITH_TBB

#include <tbb/tbb.h>
#include <tbb/task_group.h>

class exactSearchWorkerStateful {
	int tid;

public:
	exactSearchWorkerStateful(const exactSearchWorkerStateful& W): tid(W.tid) {};
	exactSearchWorkerStateful(int id):tid(id) {};
	void operator()();
};

class exactSearchWorker {
	int tid;

public:
	exactSearchWorker(const exactSearchWorker& W): tid(W.tid) {};
	exactSearchWorker(int id):tid(id) {};
	void operator()();

};

class mismatchSearchWorkerFull {
	int tid;

public:
	mismatchSearchWorkerFull(const mismatchSearchWorkerFull& W): tid(W.tid) {};
	mismatchSearchWorkerFull(int id):tid(id) {};
	void operator()();

};

class mismatchSearchWorkerFullStateful {
	int tid;

public:
	mismatchSearchWorkerFullStateful(const mismatchSearchWorkerFullStateful& W): tid(W.tid) {};
	mismatchSearchWorkerFullStateful(int id):tid(id) {};
	void operator()();

};
class twoOrThreeMismatchSearchWorkerStateful {
	int tid;

public:
	twoOrThreeMismatchSearchWorkerStateful(const twoOrThreeMismatchSearchWorkerStateful& W): tid(W.tid) {};
	twoOrThreeMismatchSearchWorkerStateful(int id):tid(id) {};
	void operator()();

};
class twoOrThreeMismatchSearchWorkerFull {
	int tid;

public:
	twoOrThreeMismatchSearchWorkerFull(const twoOrThreeMismatchSearchWorkerFull& W): tid(W.tid) {};
	twoOrThreeMismatchSearchWorkerFull(int id):tid(id) {};
	void operator()();

};

class seededQualSearchWorkerFullStateful {
	int tid;

public:
	seededQualSearchWorkerFullStateful(const seededQualSearchWorkerFullStateful& W): tid(W.tid) {};
	seededQualSearchWorkerFullStateful(int id):tid(id) {};
	void operator()();

};
class seededQualSearchWorkerFull {
	int tid;

public:
	seededQualSearchWorkerFull(const seededQualSearchWorkerFull& W): tid(W.tid) {};
	seededQualSearchWorkerFull(int id):tid(id) {};
	void operator()();

};

#endif /* WITH_TBB */
#endif /* EBWT_SEARCH_H_ */
