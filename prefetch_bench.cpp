#include <iostream>
#include <pthread.h>
#include <getopt.h>
#include "random_source.h"
#include "timer.h"

/**
 * Example application showing how prefetching and "pseudo-threads" can
 * hide latency and increase throughput in a manner similar to hardware
 * multithreadedness (e.g. on Niagara, NPUs, GPGPUs), but in software.
 * 
 * Here is some example output showing 16 pseudo-threads ("under" a
 * single normal thread) giving substantially better performance than
 * one normal thread.
 * 
 * 1 normal thread (thread_main()):
 * 
 * sycamore $ ./prefetch_bench -v -r 150000000 -t 1 -w 2 -p 1
 * # threads: 1
 * Prefetch width: 1
 * Work units: 2
 * Time: 00:00:26
 * Final result: 456047716
 * 
 * 16 pseudo-threads (thread_main_prefetch()):
 * 
 * sycamore $ ./prefetch_bench -v -r 150000000 -t 1 -w 2 -p 16
 * # threads: 1
 * Prefetch width: 16
 * Work units: 2
 * Time: 00:00:10
 * Final result: 456047716
 */

using namespace std;

static int threads   = 2;
static int prefetchWidth = 1;
static int workUnits = 8;
static int rounds    = 1000000;
static bool verbose  = false;

static struct option long_options[] = {
	/* None for now */
	{0, 0, 0, 0}
};

static const char *short_options = "vt:p:w:r:";

/**
 * Print a friendly usage message.
 */
static void printUsage(ostream& out) {
	out << "Usage: prefetch_bench [-v] [-t <int>] [-r <int>] [-p <int>] [-w <int>]" << endl;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', then output the given error message and
 * exit with an error and a usage message.
 */
static int parseInt(int lower, const char *errmsg) {
	long l;
	char *endPtr= NULL;
	l = strtol(optarg, &endPtr, 10);
	if (endPtr != NULL) {
		if (l < lower) {
			cerr << errmsg << endl;
			printUsage(cerr);
			exit(1);
		}
		return (int32_t)l;
	}
	cerr << errmsg << endl;
	printUsage(cerr);
	exit(1);
	return -1;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, char **argv) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
		switch (next_option) {
			case 't':
				threads = parseInt(1, "-t arg must be at least 1");
	   			break;
			case 'r':
				rounds = parseInt(1, "-r arg must be at least 1");
	   			break;
			case 'p':
				prefetchWidth = parseInt(1, "-p arg must be at least 1");
				break;
	   		case 'w':
				workUnits = parseInt(1, "-w arg must be at least 1");
				break;
	   		case 'v':
	   			verbose = true;
	   			break;
			case -1: /* Done with options. */
				break;
			case 0:
				if (long_options[option_index].flag != 0)
					break;	
			default: 
				cerr << "Unknown option: " << (char)next_option << endl;
				printUsage(cerr);
				exit(1);
		}
	} while(next_option != -1);
}

#define ACC_ARR_SIZE (128 * 1024 * 1024)
#define NEXT_RAND(r) max(0, abs((int)r.nextU32()))

static char *arr;

/**
 * Not pseudo-threaded.
 */
void *thread_main(void *arg) {
	pair<int,int> *p = (pair<int,int>*)arg;
	int subresult = 0;
	for(int i = p->first; i < p->second; i++) {
		int cur = 0;
		RandomSource rand(i);
		if(cur == 0) cur = NEXT_RAND(rand);
		int off = cur % ACC_ARR_SIZE;
		cur = arr[off]; // expensive miss
		cur += NEXT_RAND(rand);
		// Do some busy-work
		for(int w = 0; w < workUnits; w++) {
			cur *= NEXT_RAND(rand);
			subresult ^= cur;
		}
	}
	return reinterpret_cast<void*>(subresult);
}

/**
 * Pseudo-threded.
 */
void *thread_main_prefetch(void *arg) {
	pair<int,int> *p = (pair<int,int>*)arg;
	int cur = 0;
	int subresult = 0;
	int *curs = new int[prefetchWidth];
	RandomSource *rands = new RandomSource[prefetchWidth];
	for(int j = p->first; j < p->second; j += prefetchWidth) {
		int lower = j;
		int upper = j + prefetchWidth;
		if(upper > p->second) upper = p->second;
		bzero(curs, prefetchWidth * sizeof(int));
		// Evaluate loop up until prefetch
		for(int i = lower; i < upper; i++) {
			int& cur = curs[i - lower];
			rands[i-lower].init(i);
			if(cur == 0) cur = NEXT_RAND(rands[i-lower]);
			// If we access arr[cur % ACC_ARR_SIZE] directly now, we'll
			// almost certainly miss in the L2; instead, issue a
			// prefetch and we'll come back to it in the next loop
			__builtin_prefetch(&arr[cur % ACC_ARR_SIZE], 0, 0);
		}
		// Evaluate rest of loop
		for(int i = lower; i < upper; i++) {
			int& cur = curs[i - lower];
			// Hopefully the prefetch has eliminated some or all of the
			// stall here
			cur = arr[cur % ACC_ARR_SIZE];
			cur += NEXT_RAND(rands[i-lower]);
			// Do some busy-work
			for(int w = 0; w < workUnits; w++) {
				cur *= NEXT_RAND(rands[i-lower]);
				subresult ^= cur;
			}
		}
	}
	delete[] curs;
	delete[] rands;
	return reinterpret_cast<void*>(subresult);
}

/**
 * Set up and run all worker threads, i.e., one or more normal threads
 * (pthreads threads), each running one or more pseudo-threads (using
 * thread_main_prefetch()).  If only one normal thread is needed, then
 * pthreads is bypassed entirely.  If only one pseudo-thread is needed
 * in a particular thread, then a simpler non-pseudo-threaded version
 * of the worker function is used (thread_main()).
 */
int main(int argc, char **argv) {
	parseOptions(argc, argv);
	if(verbose) {
		cout << "# threads: " << threads << endl;
		cout << "Prefetch width: " << prefetchWidth << endl;
		cout << "Work units: " << workUnits << endl;
	}
	pthread_t *thread_arr = new pthread_t[threads];
	pair<int,int> *pair_arr = new pair<int,int>[threads];
	arr = new char[ACC_ARR_SIZE];
	bzero(arr, ACC_ARR_SIZE);

	// Initialize pair_arr[]
	int off = 0;
	for(int t = 0; t < threads; t++) {
		int nextoff = (t == threads-1)? rounds : off + rounds/threads;
		pair_arr[t] = make_pair(off, nextoff);
		off = nextoff;
	}

	int result = 0;
	{
		Timer t(cout, "Time: ", true);
		if(threads > 1) {
			if(verbose) cout << "Creating pthreads threads\n";
			// Create and run threads
			for(int t = 0; t < threads; t++){
				int rc;
				if(prefetchWidth > 1) {
					rc = pthread_create(&thread_arr[t], NULL, thread_main_prefetch, (void *)&pair_arr[t]);
				} else {
					rc = pthread_create(&thread_arr[t], NULL, thread_main, (void *)&pair_arr[t]);
				}
				if(rc) {
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
			}
			// Reduce result
			for(int t = 0; t < threads; t++) {
				int i;
				pthread_join(thread_arr[t], (void**)&i);
				result ^= i;
			}
		} else {
			if(verbose) cout << "Not using pthreads; calling thread_main directly\n";
			if(prefetchWidth > 1) {
				// don't care about losing precision
				result = (int)reinterpret_cast<int64_t>(thread_main_prefetch(reinterpret_cast<void*>(&pair_arr[0])));
			} else {
				// don't care about losing precision
				result = (int)reinterpret_cast<int64_t>(thread_main(reinterpret_cast<void*>(&pair_arr[0])));
			}
		}
	}
	cout << "Final result: " << result << endl;
	pthread_exit(NULL);
	delete[] arr;
	return 0;
}
