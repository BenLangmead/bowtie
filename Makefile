#
# Makefile for bowtie, bowtie-build, and bowtie-maqconvert
#

SEQAN_DIR = SeqAn-1.1
SEQAN_INC = -I $(SEQAN_DIR)
INC = $(SEQAN_INC)
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC = $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP = $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX = $(CPP)
HEADERS = $(wildcard *.h)

# Detect Cygwin or MinGW
WINDOWS = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
WINDOWS = 1
else
ifneq (,$(findstring MINGW,$(shell uname)))
WINDOWS = 1
endif
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
MACOS = 1
endif

BOWTIE_PTHREADS = 1
PTHREAD_PKG =
PTHREAD_LIB =
PTHREAD_DEF =
ifeq (1,$(BOWTIE_PTHREADS))
PTHREAD_DEF = -DBOWTIE_PTHREADS
ifeq (1,$(WINDOWS))
# pthreads for windows forces us to be specific about the library
PTHREAD_LIB = -lpthreadGC2
PTHREAD_PKG = pthreadGC2.dll
else
# There's also -pthread, but that only seems to work on Linux
PTHREAD_LIB = -lpthread
endif
endif

SHMEM_DEF = -DBOWTIE_SHARED_MEM
ifeq (1,$(WINDOWS))
# No shared-mem facility in Windows
SHMEM_DEF =
endif

CHUD=0
CHUD_DEF =
ifeq (1,$(CHUD))
ifeq (1,$(MACOS))
CHUD_DEF = -F/System/Library/PrivateFrameworks -weak_framework CHUD -DCHUD_PROFILING
endif
endif

PREFETCH_LOCALITY = 2
PREF_DEF = -DPREFETCH_LOCALITY=$(PREFETCH_LOCALITY)

LIBS = 
SEARCH_LIBS = $(PTHREAD_LIB)
BUILD_LIBS =

SEARCH_CPPS = qual.cpp pat.cpp ebwt_search_util.cpp ref_aligner.cpp
OTHER_CPPS = ccnt_lut.cpp hit.cpp ref_read.cpp alphabet.c
SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
MAQ_H   = $(wildcard maq_convert/*.h)
MAQ_CPP	= maq_convert/const.c \
		  maq_convert/bfa.c
# bowtie-maqconvert requires zlib because maq's format is compressed
MAQ_LIB = -lz
VERSION = $(shell cat VERSION)

# For a universal 32/64-bit x86 Mac binary compatible w/ Tiger:
# EXTRA_FLAGS="-arch i386 -arch x86_64 -mmacosx-version-min=10.4"
# For a 32-bit x86 binary on Linux (if not the default):
# EXTRA_FLAGS="-m32"
# For a 64-bit x86 binary on Linux (if not the default):
# EXTRA_FLAGS="-m64"
EXTRA_FLAGS =
DEBUG_FLAGS = -O0 -g3
DEBUG_DEFS = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS = -O3
RELEASE_DEFS = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
BIN_LIST = bowtie-build \
           bowtie \
           bowtie-maqconvert \
           bowtie-inspect \
           bowtie-maptool
BIN_LIST_AUX = bowtie-build-debug \
               bowtie-debug \
               bowtie-inspect-debug \
               bowtie-maptool-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
               $(wildcard scripts/*.pl) \
               $(wildcard indexes/e_coli*) \
               $(wildcard genomes/NC_008253.fna) \
               $(wildcard reads/e_coli_1000.*) \
               $(wildcard reads/e_coli_1000_*) \
               reads/e_coli_10000snp.fa \
               reads/e_coli_10000snp.fq \
               $(PTHREAD_PKG) \
               AUTHORS \
               COPYING \
               NEWS \
               MANUAL \
               TUTORIAL \
               VERSION

SRC_PKG_LIST = $(wildcard *.h) \
               $(wildcard *.hh) \
               $(wildcard *.c) \
               $(wildcard *.cpp) \
               $(wildcard maq_convert/*.h) \
               $(wildcard maq_convert/*.hh) \
               $(wildcard maq_convert/*.c) \
               $(wildcard maq_convert/*.cpp) \
               $(shell find SeqAn-1.1 -name '*.h') \
               $(shell find SeqAn-1.1 -name '*.txt') \
               Makefile \
               $(GENERAL_LIST)

BIN_PKG_LIST = $(GENERAL_LIST)

all: $(BIN_LIST)

allall: $(BIN_LIST) $(BIN_LIST_AUX)

DEFS=-DBOWTIE_VERSION="\"`cat VERSION`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(PTHREAD_DEF) \
     $(SHMEM_DEF) \
     $(PREF_DEF) \
     $(CHUD_DEF)

bowtie-build: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie-build_prof: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) -pg -p -g3 $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie-build-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) $(SEARCH_LIBS)

# 'bowtie' forced to be 32-bit using -m32
bowtie32: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) -m32 -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) $(SEARCH_LIBS)

# 'bowtie' forced to be 64-bit using -m64
bowtie64: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) -m64 -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) $(SEARCH_LIBS)

bowtie_prof: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) -pg -p -g3 $(EXTRA_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) $(SEARCH_LIBS)

bowtie_prof0: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) -fno-inline $(RELEASE_DEFS) -pg -p -g3 $(EXTRA_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) $(SEARCH_LIBS)

bowtie-debug: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) $(SEARCH_LIBS)

bowtie-asm: ebwt_asm.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_ASM_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

# 'bowtie-asm' forced to be 32-bit using -m32
bowtie-asm32: ebwt_asm.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) -m32 -DEBWT_ASM_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

# 'bowtie-asm' forced to be 64-bit using -m64
bowtie-asm64: ebwt_asm.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) -m64 -DEBWT_ASM_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-asm-debug: ebwt_asm.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_ASM_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-maptool: map_tool.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_MAPTOOL_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-maptool-debug: map_tool.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_MAPTOOL_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-maqconvert: maq_convert/bowtie_convert.cpp $(HEADERS) $(MAQ_H) $(MAQ_CPP)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) $(DEFS) -Wall $(INC) -I . -o $@ $< $(MAQ_CPP) $(LIBS) $(MAQ_LIB)

bowtie-maqconvert-debug: maq_convert/bowtie_convert.cpp $(HEADERS) $(MAQ_H) $(MAQ_CPP)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) $(DEFS) -Wall $(INC) -I . -o $@ $< $(MAQ_CPP) $(LIBS) $(MAQ_LIB)

bowtie-inspect: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_INSPECT_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -I . -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-inspect-debug: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS) 
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_INSPECT_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -I . -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-src.zip: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/bowtie-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
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

bowtie-all-bin.zip: $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) bowtie-asm bowtie-asm-debug
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/bowtie-$(VERSION)
	if [ -f bowtie.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX) bowtie-asm bowtie-asm-debug) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) bowtie-asm bowtie-asm-debug ; \
	fi
	mv tmp.zip .bin.tmp/bowtie-$(VERSION)
	cd .bin.tmp/bowtie-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r $@ bowtie-$(VERSION)
	cp .bin.tmp/$@ .
	rm -rf .bin.tmp

.PHONY: all-hadoop
all-hadoop: bowtie32 bowtie64 bowtie-asm32 bowtie-asm64 

.PHONY: clean
clean:
	rm -f $(BIN_LIST) $(BIN_LIST_AUX) \
	bowtie_prof bowtie-asm bowtie-asm-debug \
	$(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX) bowtie_prof) \
	bowtie-src.zip bowtie-bin.zip
	rm -f core.*
