AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -ISeqAn-1.1

CFLAGS=
LDFLAGS= -lpthread

OTHER_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
             edit.cpp ebwt.cpp tinythread.cpp
SEARCH_CPPS = qual.cpp pat.cpp ebwt_search_util.cpp ref_aligner.cpp \
              log.cpp hit_set.cpp refmap.cpp annot.cpp sam.cpp \
              color.cpp color_dec.cpp hit.cpp
SEARCH_CPPS_MAIN = $(SEARCH_CPPS) bowtie_main.cpp

if ENABLE_DEBUG
bin_PROGRAMS = bowtie-build-debug bowtie-debug bowtie-inspect-debug
else
bin_PROGRAMS = bowtie-build bowtie bowtie-inspect
endif

bowtie_build_SOURCES = ebwt_build.cpp $(OTHER_CPPS) bowtie_build_main.cpp
bowtie_build_debug_SOURCES = ebwt_build.cpp $(OTHER_CPPS) bowtie_build_main.cpp

bowtie_SOURCES = ebwt_search.cpp $(OTHER_CPPS) $(SEARCH_CPPS_MAIN) 
bowtie_debug_SOURCES = ebwt_search.cpp $(OTHER_CPPS) $(SEARCH_CPPS_MAIN) 

bowtie_inspect_SOURCES = $(OTHER_CPPS) bowtie_inspect.cpp
bowtie_inspect_debug_SOURCES = $(OTHER_CPPS) bowtie_inspect.cpp

