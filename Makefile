#
# Makefile for bowtie, bowtie-build, bowtie-inspect
#

prefix = /usr/local
bindir = $(prefix)/bin

SEQAN_DIR = ./SeqAn-1.1
# treat SeqAn as a sysdir to suppress warnings
SEQAN_INC = -isystem $(SEQAN_DIR)
INC = $(if $(RELEASE_BUILD),-I$(CURDIR)/.include) $(SEQAN_INC) -I third_party
CPP = g++
CXX = $(CPP)
CC = gcc
LIBS = $(LDFLAGS) $(if $(RELEASE_BUILD),-L$(CURDIR)/.lib) -lz
HEADERS = $(wildcard *.h)
BOWTIE_MM = 1
BOWTIE_SHARED_MEM = 1
EXTRA_FLAGS =
EXTRA_CFLAGS =
EXTRA_CXXFLAGS =
CFLAGS += $(EXTRA_CFLAGS)
CXXFLAGS += $(EXTRA_CXXFLAGS)
WARNING_FLAGS = -Wall -Wno-unused-parameter -Wno-reorder \
				-Wno-unused-local-typedefs

RELEASE_DEPENDENCIES = $(if $(RELEASE_BUILD),static-libs)

# Detect Cygwin or MinGW
WINDOWS =
CYGWIN =
MINGW =
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

MACOS =
ifneq (,$(findstring Darwin,$(shell uname)))
    MACOS = 1
	ifneq (,$(findstring 13,$(shell uname -r)))
		CPP = clang++
		CC = clang
		override EXTRA_FLAGS += -stdlib=libstdc++
	endif
	ifeq (1, $(RELEASE_BUILD))
		EXTRA_FLAGS += -mmacosx-version-min=10.9
	endif
endif

LINUX =
ifneq (,$(findstring Linux,$(shell uname)))
    LINUX = 1
    override EXTRA_FLAGS += -Wl,--hash-style=both
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

#if we're not using TBB, then we can't use queuing locks
ifeq (1,$(NO_TBB))
	NO_QUEUELOCK=1
endif

ifeq (1,$(MINGW))
	PTHREAD_LIB = 
	override EXTRA_FLAGS += -static-libgcc -static-libstdc++
else
    PTHREAD_LIB = -lpthread
endif

ifeq (1,$(NO_SPINLOCK))
	override EXTRA_FLAGS += -DNO_SPINLOCK
endif

ifneq (1,$(NO_TBB))
	LIBS += $(PTHREAD_LIB) -ltbb
	LIBS += -ltbbmalloc$(if $(RELEASE_BUILD),,_proxy)
	override EXTRA_FLAGS += -DWITH_TBB
else
	LIBS += $(PTHREAD_LIB)
endif

POPCNT_CAPABILITY ?= 1
ifeq (1, $(POPCNT_CAPABILITY))
    override EXTRA_FLAGS += -DPOPCNT_CAPABILITY
    INC += -I third_party
endif

PREFETCH_LOCALITY = 2
PREF_DEF = -DPREFETCH_LOCALITY=$(PREFETCH_LOCALITY)


SEARCH_LIBS = 
BUILD_LIBS =
INSPECT_LIBS = 

ifeq (1,$(MINGW))
    BUILD_LIBS = 
    INSPECT_LIBS = 
endif

ifeq (1,$(WITH_THREAD_PROFILING))
	override EXTRA_FLAGS += -DPER_THREAD_TIMING=1
endif

ifeq (1,$(WITH_AFFINITY))
	override EXTRA_FLAGS += -DWITH_AFFINITY=1
endif

ifeq (1,$(WITH_QUEUELOCK))
	override EXTRA_FLAGS += -DWITH_QUEUELOCK=1
endif

ifeq (1,$(WITH_FINE_TIMER))
	override EXTRA_FLAGS += -DUSE_FINE_TIMER=1
endif

OTHER_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
             edit.cpp ebwt.cpp

ifeq (1,$(WITH_COHORTLOCK))
	override EXTRA_FLAGS += -DWITH_COHORTLOCK=1
	OTHER_CPPS += cohort.cpp cpu_numa_info.cpp
endif

ifeq (1,$(NO_TBB))
	OTHER_CPPS += tinythread.cpp
endif

SEARCH_CPPS = qual.cpp pat.cpp ebwt_search_util.cpp ref_aligner.cpp \
              log.cpp hit_set.cpp sam.cpp \
              color.cpp color_dec.cpp hit.cpp
SEARCH_CPPS_MAIN = $(SEARCH_CPPS) bowtie_main.cpp

BUILD_CPPS =
BUILD_CPPS_MAIN = $(BUILD_CPPS) bowtie_build_main.cpp

SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
VERSION = $(shell cat VERSION)

BITS=32
ifeq (x86_64,$(shell uname -m))
	BITS=64
endif
# msys will always be 32 bit so look at the cpu arch instead.
ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
	ifeq (1,$(MINGW))
		BITS=64
	endif
endif

ifeq (1,$(LINUX))
    ifeq (x86_64, $(shell uname -p))
        BITS=64
    endif
endif

ifeq (32,$(BITS))
    $(error bowtie2 compilation requires a 64-bit platform )
endif

DEBUG_FLAGS = -O0 -g3 -m64
RELEASE_FLAGS = -O3 -m64
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

BIN_LIST = bowtie-build-s \
           bowtie-build-l \
           bowtie-align-s \
           bowtie-align-l \
           bowtie-inspect-s \
           bowtie-inspect-l
BIN_LIST_AUX = bowtie-build-s-debug \
               bowtie-build-l-debug \
               bowtie-align-s-debug \
               bowtie-align-l-debug \
               bowtie-inspect-s-debug \
               bowtie-inspect-l-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
               $(wildcard scripts/*.pl) \
               $(wildcard scripts/*.py) \
               $(wildcard indexes/e_coli*) \
               $(wildcard genomes/NC_008253.fna) \
               $(wildcard reads/e_coli_1000.*) \
               $(wildcard reads/e_coli_1000_*) \
               SeqAn-1.1 \
               bowtie \
               bowtie-build \
               bowtie-inspect \
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
               $(wildcard third_party/*.h) \
               $(wildcard third_party/*.hh) \
               $(wildcard third_party/*.c) \
               $(wildcard third_party/*.cpp) \
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

bowtie-build-s: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(ALL_FLAGS)  \
		$(DEFS) $(NOASSERT_FLAGS) $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie-build-l: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(ALL_FLAGS)  \
		$(DEFS) -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie-build_prof: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) -pg -p -g3 $(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(NOASSERT_FLAGS) $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie-build-s-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

bowtie-build-l-debug: ebwt_build.cpp $(OTHER_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) -DBOWTIE_64BIT_INDEX $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(BUILD_CPPS_MAIN) \
		$(LIBS) $(BUILD_LIBS)

#
# bowtie targets
#

bowtie-align-s: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(NOASSERT_FLAGS) $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie-align-l: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(NOASSERT_FLAGS) -DBOWTIE_64BIT_INDEX $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie_prof: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) -pg -p -g3 $(ALL_FLAGS) \
		$(DEFS) $(NOASSERT_FLAGS) $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie-align-s-debug: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

bowtie-align-l-debug: ebwt_search.cpp $(SEARCH_CPPS) $(OTHER_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) -DBOWTIE_64BIT_INDEX $(WARNING_FLAGS) \
		$(INC) \
		-o $@ $< \
		$(OTHER_CPPS) $(SEARCH_CPPS_MAIN) \
		$(LIBS) $(SEARCH_LIBS)

#
# bowtie-inspect targets
#

bowtie-inspect-s: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(WARNING_FLAGS) \
		$(INC) -I . \
		-o $@ $< \
		$(OTHER_CPPS) \
		$(LIBS)

bowtie-inspect-l: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
		$(RELEASE_DEFS) $(ALL_FLAGS) \
		$(DEFS) -DBOWTIE_64BIT_INDEX $(WARNING_FLAGS) \
		$(INC) -I . \
		-o $@ $< \
		$(OTHER_CPPS) \
		$(LIBS)

bowtie-inspect-s-debug: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) $(WARNING_FLAGS) \
		$(INC) -I . \
		-o $@ $< \
		$(OTHER_CPPS) \
		$(LIBS)

bowtie-inspect-l-debug: bowtie_inspect.cpp $(HEADERS) $(OTHER_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
		$(DEBUG_DEFS) $(ALL_FLAGS) \
		$(DEFS) -DBOWTIE_64BIT_INDEX $(WARNING_FLAGS) \
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

bowtie-bin.zip: $(RELEASE_DEPENDENCIES) $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) 
	$(eval PKG_DIR=bowtie-$(VERSION)-$(if $(MACOS),macos,$(if $(MINGW),mingw,linux))-x86_64)
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir -p .bin.tmp/$(PKG_DIR)
	if [ -f bowtie-align-s.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(BIN_LIST) $(BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/$(PKG_DIR)
	cd .bin.tmp/$(PKG_DIR) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r $(PKG_DIR).zip $(PKG_DIR)
	cp .bin.tmp/$(PKG_DIR).zip .
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

.PHONY: install
install: all
	mkdir -p $(DESTDIR)$(bindir)
	for file in $(BIN_LIST) bowtie-inspect bowtie-build bowtie ; do \
		cp -f $$file $(DESTDIR)$(bindir) ; \
	done

.PHONY: simple-test
simple-test: allall perl-deps
	eval `perl -I $(CURDIR)/.perllib.tmp/lib/perl5 -Mlocal::lib=$(CURDIR)/.perllib.tmp` ; \
	./scripts/test/simple_tests.pl --bowtie=./bowtie --bowtie-build=./bowtie-build

.PHONY: random-test
random-test: all perl-deps
	eval `perl -I $(CURDIR)/.perllib.tmp/lib/perl5 -Mlocal::lib=$(CURDIR)/.perllib.tmp` ; \
	./scripts/test/random_bowtie_tests.sh

.PHONY: perl-deps
perl-deps:
	if [ ! -e .perllib.tmp ]; then \
		DL=$$([ `which wget` ] && echo wget -O- || echo curl -L) ; \
		mkdir .perllib.tmp ; \
		$$DL http://cpanmin.us | perl - -l $(CURDIR)/.perllib.tmp App::cpanminus local::lib ; \
		eval `perl -I $(CURDIR)/.perllib.tmp/lib/perl5 -Mlocal::lib=$(CURDIR)/.perllib.tmp` ; \
		cpanm --force Math::Random Clone Test::Deep Sys::Info -n --quiet; \
	fi

static-libs:
	if [[ ! -d $(CURDIR)/.lib || ! -d $(CURDIR)/.inc ]]; then \
		mkdir $(CURDIR)/.lib $(CURDIR)/.include ; \
	fi ; \
	if [[ `uname` = "Darwin" ]]; then \
		export CFLAGS=-mmacosx-version-min=10.9 ; \
		export CXXFLAGS=-mmacosx-version-min=10.9 ; \
	fi ; \
	DL=$$([ `which wget` ] && echo "wget --no-check-certificate --content-disposition" || echo "curl -LJkO") ; \
	cd /tmp ; \
	$$DL https://zlib.net/zlib-1.2.11.tar.gz && tar xzf zlib-1.2.11.tar.gz && cd zlib-1.2.11 ; \
	$(if $(MINGW), mingw32-make -f win32/Makefile.gcc, ./configure --static && make) && cp libz.a $(CURDIR)/.lib && cp zconf.h zlib.h $(CURDIR)/.include ; \
	cd .. ; \
	$$DL https://github.com/01org/tbb/archive/2017_U8.tar.gz && tar xzf tbb-2017_U8.tar.gz && cd tbb-2017_U8; \
	$(if $(MINGW), mingw32-make compiler=gcc arch=ia64 runtime=mingw, make) extra_inc=big_iron.inc -j4 \
	&& cp -r include/tbb $(CURDIR)/.include && cp build/*_release/*.a $(CURDIR)/.lib


.PHONY: clean
clean:
	rm -f $(BIN_LIST) $(BIN_LIST_AUX) \
	bowtie_prof \
	$(addsuffix .exe,$(BIN_LIST) $(BIN_LIST_AUX) bowtie_prof) \
	bowtie-src.zip bowtie-bin.zip
	rm -f core.*
	rm -f bowtie-align-s-master* bowtie-align-s-no-io* 
	rm -rf .lib .include
