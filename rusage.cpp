#include <iostream>
#include <stdexcept>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

/**
 * Get resource usage via a call to getrusage and output relevant
 * results.  Note that not all of these values are set properly on all
 * platforms, e.g., the Mac seems to always report 0 resident size.
 */
void printResourceUsage(ostream& out, bool verbose = false) {
	struct rusage ru;
	int ret = getrusage(0, &ru);
	if(ret == -1) {
		throw runtime_error("getrusage returned error");
	}
	time_t uhrs  = (ru.ru_utime.tv_sec / 60) / 60;
	time_t umins = (ru.ru_utime.tv_sec / 60) % 60;
	time_t usecs =  ru.ru_utime.tv_sec % 60;
	time_t shrs  = (ru.ru_stime.tv_sec / 60) / 60;
	time_t smins = (ru.ru_stime.tv_sec / 60) % 60;
	time_t ssecs =  ru.ru_stime.tv_sec % 60;
	out << "Usr time: " << uhrs << ":" << umins << ":" << usecs << endl;
	if(verbose) {
		out << "  The total amount of time spent executing in user mode." << endl; 
	}
	out << "Sys time: " << shrs << ":" << smins << ":" << ssecs << endl;
	if(verbose) {
		out << "  The total amount of time spent in the system executing on" << endl
            << "  behalf of the process(es)." << endl; 
	}
	out << "Max resident size: " << ru.ru_maxrss << " KB (" << (ru.ru_maxrss >> 10) << " MB)" << endl;
	if(verbose) {
		out << "  The maximum resident set size utilized." << endl; 
	}
	out << "Minor page faults: " << ru.ru_minflt << endl;
	if(verbose) {
		out << "  The number of page faults serviced without any I/O activity;" << endl
            << "  here I/O activity is avoided by reclaiming a page frame from" << endl
            << "  the list of pages awaiting reallocation.." << endl; 
	}
	out << "Major page faults: " << ru.ru_majflt << endl;
	if(verbose) {
		out << "  The number of page faults serviced that required I/O activity." << endl; 
	}
	out << "Process swaps: "     << ru.ru_nswap << endl;
	if(verbose) {
		out << "  The number of times a process was swapped out of main memory." << endl; 
	}
	out << "FS Inputs: "   << ru.ru_inblock << endl;
	if(verbose) {
		out << "  Real I/Os only; I/Os serviced from the cache don't count" << endl;
	}
	out << "FS Outputs: "  << ru.ru_oublock << endl;
	if(verbose) {
		out << "  Real I/Os only; I/Os serviced from the cache don't count" << endl;
	}
	out << "Voluntary context switches: "  << ru.ru_nvcsw << endl;
	if(verbose) {
		out << "  The number of times a context switch resulted due to a" << endl
            << "  process voluntarily giving up the processor before its time" << endl
            << "  slice was completed (usually to await availability of a" << endl
            << "  resource)." << endl;
	}
	out << "Involuntary context switches: "  << ru.ru_nivcsw << endl;
	if(verbose) {
		out << "  The number of times a context switch resulted due to a" << endl
            << "  higher priority process becoming runnable or because the" << endl
            << "  current process exceeded its time slice." << endl;
	}
}

#ifdef RUSAGE_MAIN

#include <sys/time.h>
#include <sys/resource.h>

int main(void) {
	struct rlimit r;
	int ret = getrlimit(RLIMIT_DATA, &r);
	if(ret == -1) {
		cerr << "Error getting RLIMIT_DATA limit" << endl;
	} else {
		cout << "RLIMIT_DATA: " << r.rlim_max << endl;
	}
	
	ret = getrlimit(RLIMIT_AS, &r);
	if(ret == -1) {
		cerr << "Error getting RLIMIT_AS limit" << endl;
	} else {
		cout << "RLIMIT_AS: " << r.rlim_max << endl;
	}
}

#endif

