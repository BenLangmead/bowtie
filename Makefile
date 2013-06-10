#
# Makefile for bowtie, bowtie-build, bowtie-inspect
#

SEQAN_DIR = SeqAn-1.1
SEQAN_INC = -I $(SEQAN_DIR)
INC = $(SEQAN_INC)
CPP = g++
CXX = $(CPP)
CC = gcc
HEADERS = $(wildcard *.h)
BOWTIE_MM = 1
BOWTIE_SHARED_MEM = 1
EXTRA_FLAGS =
EXTRA_CFLAGS =
EXTRA_CXXFLAGS =
CFLAGS += $(EXTRA_CFLAGS)
CXXFLAGS += $(EXTRA_CXXFLAGS)

# Detect Cygwin or MinGW
WINDOWS = 0
CYGWIN = 0
MINGW = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
    WINDOWS = 1
    CYGWIN = 1
    # POSIX memory-mapped files not currently supported on Windows
    BOWTIE_MM = 0
    BOWTIE_SHARED_MEM = 0
else
    ifneq (,$(findstring MINGW,$(shell uname)))
	WINDOWS = 1
	MINGW = 1
	# POSIX memory-mapped files not currently supported on Windows
	BOWTIE_MM = 0
	BOWTIE_SHARED_MEM = 0
    endif
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
    MACOS = 1
endif

LINUX = 0
ifneq (,$(findstring Linux,$(shell uname)))
    LINUX = 1
    EXTRA_FLAGS += -Wl,--hash-style=both
endif

MM_DEF = 
ifeq (1,$(BOWTIE_MM))
    MM_DEF = -DBOWTIE_MM
endif
SHMEM_DEF = 
ifeq (1,$(BOWTIE_SHARED_MEM))
    SHMEM_DEF = -DBOWTIE_SHARED_MEM
endif
PTHREAD_PKG =
PTHREAD_LIB =
PTHREAD_DEF =

ifeq (1,$(MINGW))
	PTHREAD_LIB = 
	EXTRA_FLAGS += -static-libgcc -static-libstdc++
else
    PTHREAD_LIB = -lpthread
endif

POPCNT_CAPABILITY ?= 0
ifeq (1, $(POPCNT_CAPABILITY))
    EXTRA_FLAGS += -DPOPCNT_CAPABILITY -msse4.2
    ifeq (1,$(MACOS))
        INC += -I third_party/macos
    endif
endif

PREFETCH_LOCALITY = 2
PREF_DEF = -DPREFETCH_LOCALITY=$(PREFETCH_LOCALITY)

LIBS = $(PTHREAD_LIB)
SEARCH_LIBS = 
BUILD_LIBS =
INSPECT_LIBS = 

ifeq (1,$(MINGW))
    BUILD_LIBS = 
    INSPECT_LIBS = 
endif

OTHER_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
             edit.cpp ebwt.cpp tinythread.cpp
SEARCH_CPPS = qual.cpp pat.cpp ebwt_search_util.cpp ref_aligner.cpp \
              log.cpp hit_set.cpp refmap.cpp annot.cpp sam.cpp \
              color.cpp color_dec.cpp hit.cpp
SEARCH_CPPS_MAIN = $(SEARCH_CPPS) bowtie_main.cpp

BUILD_CPPS =
BUILD_CPPS_MAIN = $(BUILD_CPPS) bowtie_build_main.cpp

SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
VERSION = $(shell cat VERSION)

BITS=32
# cygwin will stay 32 bit for now.
ifeq (1,$(MINGW))
    # msys will always be 32 bit so look at the cpu arch.
    ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
        BITS=64
    else
        ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITECTURE)))
            BITS=64
        endif
    endif
endif

ifeq (1,$(LINUX))
    ifeq (x86_64, $(shell uname -p))
        BITS=64
    endif
endif

# Convert BITS=?? to a -m flag
BITS_FLAG =
ifeq (32,$(BITS))
    BITS_FLAG = -m32
endif
ifeq (64,$(BITS))
    BITS_FLAG = -m64
endif

DEBUG_FLAGS = -O0 -g3 $(BITS_FLAG)
RELEASE_FLAGS = -O3 $(BITS_FLAG)
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

BIN_LIST = bowtie-build \
           bowtie \
           bowtie-inspect
BIN_LIST_AUX = bowtie-build-debug \
               bowtie-debug \
               bowtie-inspect-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
               $(wildcard scripts/*.pl) \
               $(wildcard indexes/e_coli*) \
               $(wildcard genomes/NC_008253.fna) \
               $(wildcard reads/e_coli_1000.*) \
               $(wildcard reads/e_coli_1000_*) \
               SeqAn-1.1 \
               doc/manual.html \
               doc/README \
               doc/style.css \
               reads/e_coli_10000snp.fa \
               reads/e_coli_10000snp.fq \
               $(PTHREAD_PKG) \
               AUTHORS \
               LICENSE \
               NEWS \
               MANUAL \
               MANUAL.markdown \
               TUTORIAL \
               VERSION

SRC_PKG_LIST = $(wildcard *.h) \
               $(wildcard *.hh) \
               $(wildcard *.c) \
               $(wildcard *.cpp) \
               doc/strip_markdown.pl \
               Makefile \
               $(GENERAL_LIST)

BIN_PKG_LIST = $(GENERAL_LIST)

all: $(BIN_LIST)

allall: $(BIN_LIST) $(BIN_LIST_AUX)

DEFS=-fno-strict-aliasing \
     -DBOWTIE_VERSION="\"`cat VERSION`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
     $(PTHREAD_DEF) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

ALL_FLAGS = $(EXTRA_FLAGS) $(CFLAGS) $(CXXFLAGS)
DEBUG_DEFS = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(ALL_FLAGS)\""
RELEASE_DEFS = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(ALL_FLAGS)\""

#
# bowtie-build targets
#

bowtie-build: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(ALL_FLAGS)  \
		$(DEFS) $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie-build_prof: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) -pg -p -g3 $(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie-build-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) -Wall \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

#
# bowtie targets
#

bowtie: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie_prof: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) -pg -p -g3 $(ALL_FLAGS) \
		$(DEFS) $(NOASSERT_FLAGS) -Wall \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie-debug: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) -Wall \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

#
# bowtie-inspect targets
#

bowtie-inspect: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(OTHER_CPPS) \
		$(LIBS)

bowtie-inspect-debug: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) -Wall \
		$(INC) -I . \
		-o $@ $< \
		$(OTHER_CPPS) \
		$(LIBS)

bowtie-src.zip: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/bowtie-$(VERSION)
	zip -r tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/bowtie-$(VERSION)
	cd .src.tmp/bowtie-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r $@ bowtie-$(VERSION)
	cp .src.tmp/$@ .
	rm -rf .src.tmp

bowtie-bin.zip: $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) 
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/bowtie-$(VERSION)
	if [ -f bowtie.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/bowtie-$(VERSION)
	cd .bin.tmp/bowtie-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r $@ bowtie-$(VERSION)
	cp .bin.tmp/$@ .
	rm -rf .bin.tmp

.PHONY: doc
doc: doc/manual.html MANUAL

doc/manual.html: MANUAL.markdown
	echo "<h1>Table of Contents</h1>" > .tmp.head
	pandoc -T "Bowtie Manual" -B .tmp.head \
	       --css style.css -o $@ \
	       --from markdown --to HTML \
	       --table-of-contents $^

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

.PHONY: clean
clean:
	rm -f $(BIN_LIST) $(BIN_LIST_AUX) \
	bowtie_prof \
	$(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX) bowtie_prof) \
	bowtie-src.zip bowtie-bin.zip
	rm -f core.*
