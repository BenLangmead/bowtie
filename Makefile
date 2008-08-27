#
# Makefile for bowtie, bowtie-build, and bowtie-convert
#

SEQAN_DIR = SeqAn-1.1
SEQAN_INC = -I $(SEQAN_DIR)
INC = $(SEQAN_INC)
GCC_PREFIX = $(shell dirname `which gcc`)
CC = $(GCC_PREFIX)/gcc
CPP = $(GCC_PREFIX)/g++
CXX = $(CPP)
HEADERS = $(wildcard *.h)
LIBS =
OTHER_CPPS = ccnt_lut.cpp hit.cpp ref_read.cpp
MAQ_H   = $(wildcard maq_convert/*.h)
MAQ_CPP	= maq_convert/maqmap.c \
          maq_convert/const.c
# bowtie-convert requires zlib because maq's format is compressed
MAQ_LIB = -lz

EXTRA_FLAGS =
DEBUG_FLAGS = -O0 -g3
DEBUG_DEFS = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS = -O3
RELEASE_DEFS = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
BIN_LIST = bowtie-build \
           bowtie-build-packed \
           bowtie \
           bowtie-convert
BIN_LIST_AUX = bowtie-build-debug \
               bowtie-build-packed-debug \
               bowtie-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
               $(wildcard scripts/*.pl) \
               $(wildcard indexes/e_coli*) \
               $(wildcard genomes/NC_008253.fna) \
               $(wildcard reads/e_coli*) \
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
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\""

bowtie-build: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-build-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-build-packed: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) -DPACKED_STRINGS $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie-build-packed-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` $(DEFS) -DPACKED_STRINGS -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(LIBS)

bowtie: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) $(NOASSERT_FLAGS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS)

bowtie-debug: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) 
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` $(DEFS) -Wall $(INC) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS) $(LIBS) 

bowtie-convert: maq_convert/bowtie_convert.cpp $(HEADERS) $(MAQ_H) $(MAQ_CPP)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) -Wall $(INC) -I . -o $@ $< $(MAQ_CPP) $(LIBS) $(MAQ_LIB)

bowtie-src.zip: $(SRC_PKG_LIST)
	zip $@ $(SRC_PKG_LIST)

bowtie-bin.zip: $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) 
	if [ -f bowtie.exe ] ; then \
		zip $@ $(BIN_PKG_LIST) $(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX)) ; \
	else \
		zip $@ $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) ; \
	fi

.PHONY: clean
clean:
	rm -f $(BIN_LIST) $(BIN_LIST_AUX) \
	$(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX)) \
	bowtie-src.zip bowtie-bin.zip
