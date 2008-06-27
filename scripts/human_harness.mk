#!/bin/sh

OUTDIR=.
dummy := $(shell mkdir -p $(OUTDIR))

#
# Sanity tests for various bucket sizes.
#

all:         $(OUTDIR)/experiment.chr21_unpack_dc1024.out \
             $(OUTDIR)/experiment.chr11_unpack_dc1024.out \
             $(OUTDIR)/experiment.chr1_unpack_dc1024.out \
             $(OUTDIR)/experiment.chr21_unpack_nodc.out \
             $(OUTDIR)/experiment.chr11_unpack_nodc.out \
             $(OUTDIR)/experiment.chr1_unpack_nodc.out \
             $(OUTDIR)/experiment.chr21_unpack_dc4096.out \
             $(OUTDIR)/experiment.chr11_unpack_dc4096.out \
             $(OUTDIR)/experiment.chr1_unpack_dc4096.out \
             $(OUTDIR)/experiment.chr21_unpack_dc2048.out \
             $(OUTDIR)/experiment.chr11_unpack_dc2048.out \
             $(OUTDIR)/experiment.chr1_unpack_dc2048.out \
             $(OUTDIR)/experiment.chr21_unpack_dc512.out \
             $(OUTDIR)/experiment.chr11_unpack_dc512.out \
             $(OUTDIR)/experiment.chr1_unpack_dc512.out \
             $(OUTDIR)/experiment.chr21_unpack_dc256.out \
             $(OUTDIR)/experiment.chr11_unpack_dc256.out \
             $(OUTDIR)/experiment.chr1_unpack_dc256.out \
             $(OUTDIR)/experiment.chr21_unpack_dc128.out \
             $(OUTDIR)/experiment.chr11_unpack_dc128.out \
             $(OUTDIR)/experiment.chr1_unpack_dc128.out \
             $(OUTDIR)/experiment.chr21_pack_dc1024.out \
             $(OUTDIR)/experiment.chr11_pack_dc1024.out \
             $(OUTDIR)/experiment.chr1_pack_dc1024.out

all_asserts: $(OUTDIR)/experiment.chr21_unpack_dc1024_asserts.out \
             $(OUTDIR)/experiment.chr11_unpack_dc1024_asserts.out \
             $(OUTDIR)/experiment.chr1_unpack_dc1024_asserts.out

all_whole:   $(OUTDIR)/experiment.whole_unpack_dc1024.out \
             $(OUTDIR)/experiment.whole_unpack_dc512.out \
             $(OUTDIR)/experiment.whole_unpack_dc256.out \
             $(OUTDIR)/experiment.whole_pack_dc2048.out

# DCV=1024, asserts, chromosomes 1, 11, 21, variety of bucket sizes

$(OUTDIR)/experiment.chr21_unpack_dc1024_asserts.out:
	perl scripts/human.pl -c 21  -r "--dcv 1024" -a -b 64000000,32000000,16000000,8000000 chr21_unpack_dc1024_asserts $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_dc1024_asserts.out:
	perl scripts/human.pl -c 11  -r "--dcv 1024" -a -b 128000000,64000000,32000000,16000000 chr11_unpack_dc1024_asserts $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_dc1024_asserts.out:
	perl scripts/human.pl -c  1  -r "--dcv 1024" -a -b 256000000,128000000,64000000,32000000 chr1_unpack_dc1024_asserts $(OUTDIR) >$@ 2>&1

# No DC, no asserts, chromosomes 1, 11, 21, smallish bucket size

$(OUTDIR)/experiment.chr21_unpack_nodc.out:
	perl scripts/human.pl -c 21 -r "--noDc" -b 16000000 chr21_unpack_nodc $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_nodc.out:
	perl scripts/human.pl -c 11 -r "--noDc" -b 32000000 chr11_unpack_nodc $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_nodc.out:
	perl scripts/human.pl -c  1 -r "--noDc" -b 64000000 chr1_unpack_nodc $(OUTDIR) >$@ 2>&1


# DCV=128, no asserts, chromosomes 1, 11, 21, smallish bucket size

$(OUTDIR)/experiment.chr21_unpack_dc128.out:
	perl scripts/human.pl -c 21 -r "--dcv 128" -b 16000000 chr21_unpack_dc128 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_dc128.out:
	perl scripts/human.pl -c 11 -r "--dcv 128" -b 32000000 chr11_unpack_dc128 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_dc128.out:
	perl scripts/human.pl -c  1 -r "--dcv 128" -b 64000000 chr1_unpack_dc128 $(OUTDIR) >$@ 2>&1


# DCV=256, no asserts, chromosomes 1, 11, 21, smallish bucket size

$(OUTDIR)/experiment.chr21_unpack_dc256.out:
	perl scripts/human.pl -c 21 -r "--dcv 256" -b 16000000 chr21_unpack_dc256 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_dc256.out:
	perl scripts/human.pl -c 11 -r "--dcv 256" -b 32000000 chr11_unpack_dc256 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_dc256.out:
	perl scripts/human.pl -c  1 -r "--dcv 256" -b 64000000 chr1_unpack_dc256 $(OUTDIR) >$@ 2>&1


# DCV=512, no asserts, chromosomes 1, 11, 21, smallish bucket size

$(OUTDIR)/experiment.chr21_unpack_dc512.out:
	perl scripts/human.pl -c 21 -r "--dcv 512" -b 16000000 chr21_unpack_dc512 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_dc512.out:
	perl scripts/human.pl -c 11 -r "--dcv 512" -b 32000000 chr11_unpack_dc512 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_dc512.out:
	perl scripts/human.pl -c  1 -r "--dcv 512" -b 64000000 chr1_unpack_dc512 $(OUTDIR) >$@ 2>&1


# DC=1024, no asserts, chromosomes 1, 11, 21, variety of bucket sizes

$(OUTDIR)/experiment.chr21_unpack_dc1024.out:
	perl scripts/human.pl -c 21 -r "--dcv 1024" -b 32000000,16000000,8000000,4000000,2000000,1000000,500000 chr21_unpack_dc1024 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_dc1024.out:
	perl scripts/human.pl -c 11 -r "--dcv 1024" -b 128000000,64000000,32000000,16000000,8000000,4000000,2000000 chr11_unpack_dc1024 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_dc1024.out:
	perl scripts/human.pl -c  1 -r "--dcv 1024" -b 128000000,64000000,32000000,16000000,8000000,4000000,2000000 chr1_unpack_dc1024 $(OUTDIR) >$@ 2>&1


# DCV=2048, no asserts, chromosomes 1, 11, 21, smallish bucket size

$(OUTDIR)/experiment.chr21_unpack_dc2048.out:
	perl scripts/human.pl -c 21 -r "--dcv 2048" -b 16000000 chr21_unpack_dc2048 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_dc2048.out:
	perl scripts/human.pl -c 11 -r "--dcv 2048" -b 32000000 chr11_unpack_dc2048 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_dc2048.out:
	perl scripts/human.pl -c  1 -r "--dcv 2048" -b 64000000 chr1_unpack_dc2048 $(OUTDIR) >$@ 2>&1

# DCV=4096, no asserts, chromosomes 1, 11, 21, smallish bucket size

$(OUTDIR)/experiment.chr21_unpack_dc4096.out:
	perl scripts/human.pl -c 21 -r "--dcv 4096" -b 16000000 chr21_unpack_dc4096 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_unpack_dc4096.out:
	perl scripts/human.pl -c 11 -r "--dcv 4096" -b 32000000 chr11_unpack_dc4096 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_unpack_dc4096.out:
	perl scripts/human.pl -c  1 -r "--dcv 4096" -b 64000000 chr1_unpack_dc4096 $(OUTDIR) >$@ 2>&1


# DCV=1024, no asserts, chromosomes 1, 11, 21, smallish bucket size, packed

$(OUTDIR)/experiment.chr21_pack_dc1024.out:
	perl scripts/human.pl -p -c 21 -r "--dcv 1024" -b 16000000 chr21_pack_dc1024 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr11_pack_dc1024.out:
	perl scripts/human.pl -p -c 11 -r "--dcv 1024" -b 32000000 chr11_pack_dc1024 $(OUTDIR) >$@ 2>&1

$(OUTDIR)/experiment.chr1_pack_dc1024.out:
	perl scripts/human.pl -p -c  1 -r "--dcv 1024" -b 64000000 chr1_pack_dc1024 $(OUTDIR) >$@ 2>&1


#
# Whole-human-genome experiments
#

# Workstation (<2GB peak VM usage)
$(OUTDIR)/experiment.whole_pack_dc4096.out:
	perl scripts/human.pl -p -r "--dcv 4096" -b  128000000   whole_pack_dc4096 $(OUTDIR) >$@ 2>&1

# Server (~10GB peak VM usage)
$(OUTDIR)/experiment.whole_unpack_dc1024.out:
	perl scripts/human.pl    -r "--dcv 1024" -b 1280000000,960000000 whole_unpack_dc1024 $(OUTDIR) >$@ 2>&1

# Server
$(OUTDIR)/experiment.whole_unpack_dc512.out:
	perl scripts/human.pl    -r "--dcv 512" -b 1280000000,960000000 whole_unpack_dc512 $(OUTDIR) >$@ 2>&1

# Server
$(OUTDIR)/experiment.whole_unpack_dc256.out:
	perl scripts/human.pl    -r "--dcv 256" -b 1280000000,960000000 whole_unpack_dc256 $(OUTDIR) >$@ 2>&1

clean:
	rm -f $(OUTDIR)/*.chr21_*.out
	rm -f $(OUTDIR)/*.chr11_*.out
	rm -f $(OUTDIR)/*.chr1_*.out
	rm -f $(OUTDIR)/*.whole_*.out
