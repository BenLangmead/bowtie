#ifndef MM_H_
#define MM_H_

/**
 * mm.h:
 *
 * Defines that make it easier to handle files in the two different MM
 * contexts: i.e. on Linux and Mac where MM is supported and POSIX I/O
 * functions work as expected, and on Windows where MM is not supported
 * and where there isn't POSIX I/O,
 */

#ifdef BOWTIE_MM
#define MM_FILE_CLOSE close
#define MM_READ_RET ssize_t
#define MM_READ read
#define MM_FILE int
#define MM_FILE_INIT -1
#else
#define MM_FILE_CLOSE fclose
#define MM_READ_RET size_t
#define MM_READ(file, dest, sz) fread(dest, 1, sz, file)
#define MM_FILE FILE*
#define MM_FILE_INIT NULL
#endif

#endif /* MM_H_ */
