/*
 * getusage.c, version 0.4
 *
 * Usage statistic checker for Linux 2.6:
 *   Prints out an approximation of a process's peak memory usage, CPU usage
 *   and run time.
 *
 * Copyright (C) 2007 Bernard Blackham (bernard NOSPAM at largestprime dot net)
 * Licenced under the GPLv2.
 *
 *
 * 
 * To compile:
 *   gcc getusage.c -shared -o getusage.so -ldl
 *
 * To use for a one-off program:
 *   LD_LIBRARY_PATH=path_of_getusage.so LD_PRELOAD=getusage.so ./program args
 *
 *   where path_of_getusage.so does *not* include getusage.so, eg:
 *   LD_LIBRARY_PATH=/tmp LD_PRELOAD=getusage.so ./program args
 *
 * To use for all processes spawned from your shell:
 *   $ export LD_LIBRARY_PATH=path_of_getusage.so
 *   $ export LD_PRELOAD=getusage.so
 *   $ ./program args
 *
 * To turn it off again:
 *   $ unset LD_PRELOAD LD_LIBRARY_PATH
 *
 *
 * NOTES:
 *
 *  Why is memory usage approximate?
 *   On 2.6 kernels that don't provide VmPeak in /proc/$pid/stat, we record the
 *   memory usage after every trapped invocation of libc's malloc. Memory could
 *   be allocated at other times, such as by mmap(), brk(), SysV IPC, or by
 *   programs which have their own malloc implementations or don't use the
 *   dynamic linker (in which case this hack is useless).
 */

#include <time.h>
#include <assert.h>
#include <ctype.h>
#include <malloc.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <unistd.h>
#include <sys/times.h>

/* Track maximum memory usage in kB */
static int max_mem_usage;

/* Store the process's short name (as appears in ps) */
static char process_name[1024];

/* Keep track of the old malloc hook function. */
static void *(*old_malloc_hook)(size_t, const void *);

/* Start and end times of this process. */
static clock_t start_time, end_time;

/* Recent 2.6 kernels offer the VmPeak statistic, which means we don't need to
 * interrogate the process after every malloc, and is *much* more accurate. */
static int have_vm_peak;

static void get_mem_usage() {
	FILE *f;
	char procfn[32];
	char s[1024];

	snprintf(procfn, sizeof(procfn), "/proc/%d/status", getpid());
	procfn[sizeof(procfn)-1] = '\0';

	if ((f = fopen(procfn, "r")) == NULL) {
		fprintf(stderr, "get_mem_usage: Could not open %s: %s\n",
			procfn, strerror(errno));
		return;
	}

	/* Search /proc/$pid/status for the lines we care about. */
	while (fgets(s, sizeof(s), f)) {
		/* Grab the process name, if we haven't already. */
		if (process_name[0] == '\0' && strncmp(s, "Name:", 5) == 0) {
			char *p;
			strncpy(process_name, s+6, sizeof(process_name));
			for (p = process_name; *p; p++)
				if (*p == '\r' || *p == '\n') {
					*p = '\0';
					break;
				}
		}
		else if (strncmp(s, "VmPeak:", 7) == 0 ||
			   strncmp(s, "VmSize:", 7) == 0)
		{
			int usage;
			char usage_str[32];
			char *p;

			/* s is either "VmPeak: ..." or "VmSize:..." */
			if (!have_vm_peak && s[3] == 'P')
				have_vm_peak = 1;
			
			/* Strip leading whitespace. */
			for (p = s+8; isspace(*p); p++)
				;
			strncpy(usage_str, p, sizeof(usage_str));

			/* Strip trailing whitespace */
			for (p = usage_str; *p && !isspace(*p); p++)
				;
			*p = '\0';

			/* We only deal with kB at the moment. */
			assert(p[1] == 'k' && p[2] == 'B');

			usage = atoi(usage_str);
			if (max_mem_usage < usage)
				max_mem_usage = usage;
			
			break;
		}
	}

	fclose(f);
}

static void *my_malloc(size_t size, const void *caller) {
	void *ret;

	__malloc_hook = old_malloc_hook;

	ret = malloc(size);
	if (ret != NULL)
		get_mem_usage();

	__malloc_hook = my_malloc;

	return ret;
}

static void __attribute__((constructor)) my_init() {
	start_time = times(NULL);

	/* Call this to find out if we have VmPeak or not. */
	get_mem_usage();

	/* If we have VmPeak, we don't need to bother hooking malloc. */
	if (!have_vm_peak) {
		old_malloc_hook = __malloc_hook;
		__malloc_hook = my_malloc;
	}
}

static void __attribute__((destructor)) print_usage() {
	FILE *f;
	long cps;
	struct tms tms;

	/* Get final memory usage. */
	get_mem_usage();

	/* Get user, system and total times. */
	end_time = times(&tms);
	cps = sysconf(_SC_CLK_TCK);

	/* Print to /dev/tty, as stdout/stderr maybe redirected to somewhere
	 * useless. */
	f = fopen("/dev/tty", "w");
	if (f == NULL)
		f = stderr;
	if (f != NULL)
		fprintf(f,
			"[*] Usage for %-10s (%5d): %6d kB, "
			"%2.2fs user, %2.2fs sys, %2.2f tot\n",
			process_name, getpid(), max_mem_usage,
			tms.tms_utime*1.0/cps, tms.tms_stime*1.0/cps,
			(end_time - start_time)*1.0/cps);
	if (f != stderr)
		fclose(f);
}
