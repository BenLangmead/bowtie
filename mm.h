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

#define MM_IS_IO_ERR(file_hd, ret, count) is_fread_err(file_hd, ret, count)
#define MM_READ(file, dest, sz) fread(dest, 1, sz, file)

#endif /* MM_H_ */
