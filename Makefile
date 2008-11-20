#
# Makefile for bowtie, bowtie-build, and bowtie-maqconvert
#

SEQAN_DIR = SeqAn-1.1
SEQAN_INC = -I $(SEQAN_DIR)
INC = $(SEQAN_INC)
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX = -4.2
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

LIBS = 
SEARCH_LIBS = $(PTHREAD_LIB)
BUILD_LIBS =

SEARCH_CPPS = qual.cpp pat.cpp
OTHER_CPPS = ccnt_lut.cpp hit.cpp ref_read.cpp alphabet.c
SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
MAQ_H   = $(wildcard maq_convert/*.h)
MAQ_CPP	= maq_convert/const.c \
		  maq_convert/bfa.c
# bowtie-maqconvert requires zlib because maq's format is compressed
MAQ_LIB = -lz
VERSION = $(shell cat VERSION)

EXTRA_FLAGS =
DEBUG_FLAGS = -O0 -g3
DEBUG_DEFS = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS = -O3
RELEASE_DEFS = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
BIN_LIST = bowtie-build \
           bowtie-build-packed \
           bowtie \
           bowtie-maqconvert \
		   bowtie-inspect
BIN_LIST_AUX = bowtie-build-debug \
               bowtie-build-packed-debug \
               bowtie-debug \
			   bowtie-inspect-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
               $(wildcard scripts/*.pl) \
               $(wildcard indexes/e_coli*) \
               $(wildcard genomes/NC_008253.fna) \
               $(wildcard reads/e_coli*) \
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
     $(PTHREAD_DEF)

bowtie-build: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie-build_prof: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) -pg -p -g3 $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie-build-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie-build-packed: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) -DPACKED_STRINGS $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie-build-packed-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) -DPACKED_STRINGS -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS) $(BUILD_LIBS)

bowtie: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) $(SEARCH_LIBS)

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

bowtie-inspect: maq_convert/bowtie_inspect.cpp $(HEADERS) $(MAQ_H) $(MAQ_CPP) $(OTHER_CPPS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_INSPECT_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -I . -o $@ $< $(MAQ_CPP) $(OTHER_CPPS) $(LIBS) $(MAQ_LIB)

bowtie-inspect-debug: maq_convert/bowtie_inspect.cpp $(HEADERS) $(MAQ_H) $(MAQ_CPP) $(OTHER_CPPS) 
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_INSPECT_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -I . -o $@ $< $(MAQ_CPP) $(OTHER_CPPS) $(LIBS) $(MAQ_LIB)


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

.PHONY: clean
clean:
	rm -f $(BIN_LIST) $(BIN_LIST_AUX) bowtie_prof \
	$(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX) bowtie_prof) \
	bowtie-src.zip bowtie-bin.zip
