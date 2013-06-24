#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([bowtie], [1.0.0], [])
#AM_INIT_AUTOMAKE([bowtie], [1.0.0])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([aligner.h])
#AC_CONFIG_HEADER([config.h])
AC_GNU_SOURCE
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CPP

CXXFLAGS=""
AC_ARG_ENABLE([debug],
  [AC_HELP_STRING([--enable-debug], [build with debug symbols])],
  [ case "${enableval}" in
    yes)
	CXXFLAGS+="-O0 -g3 ";
	AC_DEFINE([NDEBUG], [1], [])
	;;
    no) 
	CXXFLAGS+="-O3 ";
	;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-debug);;
	esac ],		
  [debug=yes])
AM_CONDITIONAL([ENABLE_DEBUG], [test x$enable_debug = xyes])

if test "$CXX" = "g++"; then
	gcc_version=`$CXX -dumpversion`
 	CXXFLAGS+="-fno-strict-aliasing -m64 -Wl,-hash-style=both " 	
	case $gcc_version in
		4.1.2)
			#AC_MSG_ERROR("CXX verions $gcc_version is")
			;;
	esac
fi

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h memory.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_CHECK_FUNCS([floor gettimeofday memmove memset sqrt strrchr strstr strtol])

AC_DEFINE([BOWTIE_MM], [1], [])
AC_DEFINE([BOWTIE_SHARED_MEM], [1], [])
AC_DEFINE_UNQUOTED([BOWTIE_VERSION], ["`cat VERSION`"], [])
AC_DEFINE_UNQUOTED([BUILD_HOST], ["`hostname`"], [])
AC_DEFINE_UNQUOTED([BUILD_TIME], ["`date "+%Y%m%d"`"], [])
AC_DEFINE_UNQUOTED([COMPILER_VERSION], ["$gcc_version"], [])
AC_DEFINE_UNQUOTED([COMPILER_OPTIONS], ["$CXXFLAGS"], [])
AC_DEFINE([PREFETCH_LOCALITY], [2], [])
#AC_DEFINE([_LARGEFILE_SOURCE], [], [])
AC_DEFINE([_FILE_OFFSET_BITS], [64], [])

AC_CANONICAL_HOST
#AC_USE_SYSTEM_EXTENSIONS
AC_FUNC_FSEEKO
AC_SYS_LARGEFILE

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
