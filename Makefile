INC = -I../SeqAn-1.0
GCC_VERSION = 3.4.6
PREFIX = /usr/bin/
#PREFIX = $HOME/software/gcc-4.2.4/targ/
CC = ${PREFIX}gcc
CPP = ${PREFIX}g++
CXX = ${CPP}
HEADERS = $(wildcard *.h)
OTHER_CPPS = endian.cpp \
             word_io.cpp \
             ccnt_lut.cpp \
             tokenize.cpp \
             blockwise_sa.cpp \
             rusage.cpp \
             diff_sample.cpp \
             hit.cpp \
             ref_read.cpp

MAQ_HEADERS = maq/maqmap.h \
              maq/const.h
			  
MAQ_CPPS	= maq/maqmap.c \
              maq/const.c
			 
SEARCH_CPPS = LVKernel.cpp \
              inexact_extend.cpp \
			  packed_io.cpp

DEBUG_FLAGS = -O0 -g3
RELEASE_FLAGS = -O3 
NOASSERT_FLAGS = -DNDEBUG
64BIT_FLAGS = -m64
BIN_LIST = ebwt_build \
           ebwt_build_packed \
           ebwt_build-debug \
           ebwt_build_packed-debug \
           ebwt_build-with-asserts \
           ebwt_build_packed-with-asserts \
           ebwt_search \
           ebwt_search_prof \
           ebwt_search-debug \
           ebwt_search-with-asserts \
           multikey_qsort \
           blockwise_sa \
           diff_sample \
           diff_sample-with-asserts \
           lcp \
           rusage \
           pack_fasta \
           simreads
LIBS = -lz

all: $(BIN_LIST)

ebwt_build: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -DEBWT_BUILD_HASH=$$m -DEBWT_BUILD_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build_prof0: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -DEBWT_BUILD_HASH=$$m -fno-inline -pg -g -DEBWT_BUILD_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build_prof1: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -DEBWT_BUILD_HASH=$$m -pg -g -DEBWT_BUILD_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build-with-asserts: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(DEBUG_FLAGS) -DEBWT_BUILD_HASH=$$m -DEBWT_BUILD_MAIN -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(DEBUG_FLAGS) -DEBWT_BUILD_HASH=$$m -DEBWT_BUILD_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build_packed: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -DEBWT_BUILD_HASH=$$m -DEBWT_BUILD_MAIN -DPACKED_STRINGS $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build_packed-with-asserts: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(DEBUG_FLAGS) -DEBWT_BUILD_HASH=$$m -DEBWT_BUILD_MAIN -DPACKED_STRINGS -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build_packed-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(DEBUG_FLAGS) -DEBWT_BUILD_HASH=$$m -DEBWT_BUILD_MAIN -DPACKED_STRINGS $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_build_packed_prof1: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -DEBWT_BUILD_HASH=$$m -pg -g -DEBWT_BUILD_MAIN -DPACKED_STRINGS $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

ebwt_search: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -DEBWT_SEARCH_HASH=$$m -DEBWT_SEARCH_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS)

ebwt_search_prof: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -DEBWT_SEARCH_HASH=$$m -DEBWT_SEARCH_MAIN $(NOASSERT_FLAGS) -g -pg -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS)

# Just like ebwt_search_prof but with inlines turned off so that it's
# easier to interpret the gprof results
ebwt_search_prof0: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(RELEASE_FLAGS) -fno-inline -DEBWT_SEARCH_HASH=$$m -DEBWT_SEARCH_MAIN $(NOASSERT_FLAGS) -g -pg -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS)

ebwt_search-with-asserts: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) 
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(DEBUG_FLAGS) -DEBWT_SEARCH_HASH=$$m -DEBWT_SEARCH_MAIN -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS)

ebwt_search-debug: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS)
	cat $^ | cksum | sed 's/[01-9][01-9] .*//' > .$@.cksum
	m=`cat .$@.cksum` && \
	$(CXX) $(DEBUG_FLAGS) -DEBWT_SEARCH_HASH=$$m -DEBWT_SEARCH_MAIN $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS) $(SEARCH_CPPS)

multikey_qsort: multikey_qsort.cpp diff_sample.cpp tokenize.cpp endian.cpp $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) -DMULTIKEY_QSORT_MAIN -Wall $(INC) $(LIBS) -o $@ $< diff_sample.cpp tokenize.cpp endian.cpp

blockwise_sa: blockwise_sa.cpp $(HEADERS) diff_sample.cpp
	$(CXX) $(DEBUG_FLAGS) -DBLOCKWISE_SA_MAIN -Wall $(INC) $(LIBS) -o $@ $< diff_sample.cpp

lcp: lcp.cpp $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) -DLCP_MAIN -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

diff_sample: diff_sample.cpp $(HEADERS) tokenize.cpp endian.cpp
	$(CXX) $(RELEASE_FLAGS) $(NOASSERT_FLAGS) -DDIFF_SAMPLE_MAIN -Wall $(INC) $(LIBS) -o $@ $< tokenize.cpp endian.cpp

diff_sample-with-asserts: diff_sample.cpp $(HEADERS) tokenize.cpp endian.cpp
	$(CXX) $(DEBUG_FLAGS) -DDIFF_SAMPLE_MAIN -Wall $(INC) $(LIBS) -o $@ $< tokenize.cpp endian.cpp

rusage: rusage.cpp rusage.h
	$(CXX) $(DEBUG_FLAGS) -DRUSAGE_MAIN -Wall $(INC) $(LIBS) -o $@ $<

bwt: bwt.cpp blockwise_sa.h
	$(CXX) $(RELEASE_FLAGS) $(NOASSERT_FLAGS) -Wall $(INC) $(LIBS) -o $@ $<

pack_fasta: pack_fasta.cpp packed_io.h endian.cpp tokenize.h tokenize.cpp packed_io.cpp
	$(CXX) $(RELEASE_FLAGS) -DPACK_FASTA_MAIN -Wall $(INC) $(LIBS) -o $@ $< endian.cpp tokenize.cpp packed_io.cpp

bowtie_convert: bowtie_convert.cpp tokenize.h tokenize.cpp pat.h hit.h params.h $(MAQ_HEADERS) $(MAQ_CPPS)
	$(CXX) $(DEBUG_FLAGS) -Wall $(LIBS) $(INC) -o $@ $< $(MAQ_CPPS) tokenize.cpp

simreads: simreads.cpp tokenize.h tokenize.cpp
	$(CXX) $(DEBUG_FLAGS) -DSIMREADS_MAIN -Wall $(INC) $(LIBS) -o $@ $< tokenize.cpp

txt_to_fastq: txt_to_fastq.cpp pat.h
	$(CXX) $(DEBUG_FLAGS) -DTXT_TO_FASTQ_MAIN -Wall $(INC) $(LIBS) -o $@ $< $(OTHER_CPPS)

cvg_islands: cvg_islands.cpp tokenize.cpp tokenize.h $(MAQ_HEADERS) $(MAQ_CPPS)
	$(CXX) $(DEBUG_FLAGS) -Wall $(LIBS) $(INC) -o $@ $< $(MAQ_CPPS) tokenize.cpp

.PHONY: clean
clean:
	rm -f $(BIN_LIST)
