# falco: quality control for sequencing read files
#
# Copyright (C) 2019 Guilherme De Sena Brandine and
#                    Andrew D. Smith
# Authors: Guilherme De Sena Brandine, Andrew Smith
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

EXTRA_DIST = README.md LICENSE AUTHORS Configuration example.fq test

ACLOCAL_AMFLAGS = -I m4

# ADS: fix the order of user and project set vars; shouldn't need to worry about
# setting CXXFLAGS in configure.ac
# -pthread: std::thread overlaps gzip decompression with parsing
# OPENMP_CXXFLAGS: OpenMP powers the multi-threaded (-t) parser (from AC_OPENMP)
AM_CXXFLAGS = -Wall -Wextra -Wpedantic -pthread $(OPENMP_CXXFLAGS)
AM_LDFLAGS = -pthread $(OPENMP_CXXFLAGS)
if ENABLE_HTS
CPPFLAGS += -DUSE_HTS
endif
if HAVE_LIBDEFLATE
CPPFLAGS += -DUSE_LIBDEFLATE
endif
if HAVE_ISAL
CPPFLAGS += -DUSE_ISAL
endif
if HAVE_RAPIDGZIP
CPPFLAGS += -DUSE_RAPIDGZIP @RAPIDGZIP_CPPFLAGS@
endif

bin_PROGRAMS = falcoAmadeus

falcoAmadeus_SOURCES = \
	src/falco.cpp \
	src/FastqStats.cpp \
	src/HtmlMaker.cpp \
	src/Module.cpp \
	src/StreamReader.cpp \
	src/FalcoConfig.cpp \
	src/OptionParser.cpp \
	src/smithlab_utils.cpp \
	src/PairDupAware.cpp \
	src/Complexity.cpp \
	src/InsertSize.cpp \
	src/Screen.cpp \
	src/Screen.hpp \
	src/InsertSize.hpp \
	src/Complexity.hpp \
	src/PairDupAware.hpp \
	src/Module.hpp \
	src/FastqStats.hpp \
	src/HtmlMaker.hpp \
	src/StreamReader.hpp \
	src/FalcoConfig.hpp \
	src/html_template.hpp \
	src/OptionParser.hpp \
	src/smithlab_utils.hpp \
	src/dup_hash_map.hpp \
	src/aux.hpp

TESTS = test/falco.test
TEST_EXTENSIONS = .test
CLEANFILES = tests_build/test_data.tgz

# EXTRA_DIST packs the whole test/ directory; keep stray in-tree test-suite
# artifacts out of release tarballs (a shipped falco.log breaks the VPATH
# test harness of make distcheck).
dist-hook:
	rm -f $(distdir)/test/*.log $(distdir)/test/*.trs
