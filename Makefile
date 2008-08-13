#
# Makefile for bowtie, bowtie-build, and bowtie-convert.  Adjust the
# SEQAN_INC variable to point to your SeqAn installation.
#

SEQAN_INC = -I../SeqAn-1.0
INC = $(SEQAN_INC)
PREFIX = /usr/bin/
CC = ${PREFIX}gcc
CPP = ${PREFIX}g++
CXX = ${CPP}
HEADERS = $(wildcard *.h)
LIBS =
OTHER_CPPS = ccnt_lut.cpp \
             hit.cpp \
             ref_read.cpp

MAQ_HEADERS = $(wildcard maq_convert/*.h)
MAQ_CPPS	= maq_convert/maqmap.c \
              maq_convert/const.c
MAQ_LIBS    = -lz

DEBUG_FLAGS = -O0 -g3
RELEASE_FLAGS = -O3 
NOASSERT_FLAGS = -DNDEBUG
BIN_LIST = bowtie-build \
           bowtie-build-debug \
           bowtie-build-packed \
           bowtie-build-packed-debug \
           bowtie \
           bowtie-debug \
           bowtie-convert

all: $(BIN_LIST)

bowtie-build: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` -DEBWT_BUILD_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

bowtie-build-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` -DEBWT_BUILD_MAIN -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

bowtie-build-packed: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` -DEBWT_BUILD_MAIN -DPACKED_STRINGS $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

bowtie-build-packed-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) -DEBWT_BUILD_HASH=`cat .$@.cksum` -DEBWT_BUILD_MAIN -DPACKED_STRINGS -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

bowtie: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(RELEASE_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` -DEBWT_SEARCH_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS)

bowtie-debug: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) 
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	$(CXX) $(DEBUG_FLAGS) -DEBWT_SEARCH_HASH=`cat .$@.cksum` -DEBWT_SEARCH_MAIN -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS)

bowtie-convert: maq_convert/bowtie_convert.cpp $(HEADERS) $(MAQ_HEADERS) $(MAQ_CPPS)
	$(CXX) $(DEBUG_FLAGS) -Wall $(LIBS) $(MAQ_LIBS) $(INC) -I . -o $@ $< $(MAQ_CPPS)

.PHONY: clean
clean:
	rm -f $(BIN_LIST)
