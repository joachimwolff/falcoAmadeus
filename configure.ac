dnl  falco: quality control for sequencing read files
dnl
dnl  Copyright (C) 2019-2024 Guilherme De Sena Brandine and
dnl                          Andrew D. Smith
dnl  Authors: Guilherme De Sena Brandine, Andrew Smith
dnl
dnl  This program is free software: you can redistribute it and/or
dnl  modify it under the terms of the GNU General Public License as
dnl  published by the Free Software Foundation, either version 3 of the
dnl  License, or (at your option) any later version.
dnl
dnl  This program is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl  General Public License for more details.

AC_INIT([FalcoAmadeus], [1], [],
        [falcoAmadeus], [https://github.com/joachimwolff/falco])
dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX

# Set CXXFLAGS only if the user hasn't set them
AS_IF([test "x$CXXFLAGS_set" != "xset"], [
  CXXFLAGS="-O3 -DNDEBUG"
])
AC_MSG_NOTICE([CXXFLAGS is $CXXFLAGS])

AC_DEFINE_UNQUOTED([PROGRAM_PATH], ["$srcdir"],
  [The PROGRAM_PATH variable])

AX_CXX_COMPILE_STDCXX_17([noext], [mandatory])
dnl OpenMP powers the multi-threaded (-t) gzip parser; if it is unavailable the
dnl code falls back to a single worker (still correct, just serial).
AC_OPENMP

dnl check for the Zlib library
AC_CHECK_LIB([z], [zlibVersion], [],
  [AC_MSG_FAILURE([Zlib must be installed to build falco])], [])

dnl check for HTSLib if requested
hts_fail_msg="

Failed to locate HTSLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library
and headers can be found.
"
AC_ARG_ENABLE([hts],
  [AS_HELP_STRING([--enable-hts], [Enable HTSLib @<:@yes@:>@])],
  [enable_hts=yes], [enable_hts=no])
AS_IF([test "x$enable_hts" = "xyes"],
  [AC_CHECK_LIB([hts], [hts_version], [],
  [AC_MSG_FAILURE([$hts_fail_msg])])]
)
AM_CONDITIONAL([ENABLE_HTS], [test "x$enable_hts" = "xyes"])

dnl With HTSLib enabled, optionally use libdeflate directly for parallel BGZF
dnl block decompression (both are external libraries). Falls back gracefully.
have_libdeflate=no
AS_IF([test "x$enable_hts" = "xyes"],
  [AC_CHECK_HEADER([libdeflate.h],
    [AC_CHECK_LIB([deflate], [libdeflate_gzip_decompress],
      [have_libdeflate=yes
       LIBS="-ldeflate $LIBS"], [])], [])]
)
AM_CONDITIONAL([HAVE_LIBDEFLATE], [test "x$have_libdeflate" = "xyes"])

dnl Optionally use Intel ISA-L igzip for fast streaming gzip inflate. Independent
dnl of HTSLib -- it accelerates the plain (non-BGZF) gzip path, whose single
dnl DEFLATE stream cannot be parallel-decompressed. Falls back to zlib/htslib.
have_isal=no
AC_CHECK_HEADER([isa-l/igzip_lib.h],
  [AC_CHECK_LIB([isal], [isal_inflate],
    [have_isal=yes
     LIBS="-lisal $LIBS"], [])], [])
AM_CONDITIONAL([HAVE_ISAL], [test "x$have_isal" = "xyes"])

dnl Optional rapidgzip (header-only, MIT): parallel decompression of arbitrary
dnl single-stream gzip. Point --with-rapidgzip at the librapidarchive/src
dnl directory of a rapidgzip checkout. Requires zlib (already mandatory).
AC_ARG_WITH([rapidgzip],
  [AS_HELP_STRING([--with-rapidgzip=DIR],
     [rapidgzip source dir (librapidarchive/src) for parallel gzip inflate])],
  [rapidgzip_dir=$withval], [rapidgzip_dir=no])
have_rapidgzip=no
AS_IF([test "x$rapidgzip_dir" != "xno"],
  [AS_IF([test -f "$rapidgzip_dir/rapidgzip/ParallelGzipReader.hpp"],
     [have_rapidgzip=yes
      RAPIDGZIP_CPPFLAGS="-I$rapidgzip_dir -I$rapidgzip_dir/core -I$rapidgzip_dir/huffman -I$rapidgzip_dir/indexed_bzip2"
      AC_SUBST([RAPIDGZIP_CPPFLAGS])],
     [AC_MSG_FAILURE([rapidgzip/ParallelGzipReader.hpp not found under $rapidgzip_dir])])])
AM_CONDITIONAL([HAVE_RAPIDGZIP], [test "x$have_rapidgzip" = "xyes"])

AC_CONFIG_FILES([Makefile])
dnl make the test data files available in the build tree
AC_CONFIG_LINKS([
test_build/md5sum.txt:test/md5sum.txt
test_build/test_data.tgz:test/test_data.tgz
])

AC_OUTPUT
