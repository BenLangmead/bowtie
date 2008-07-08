#!/bin/sh

BASE_ARGS="-tfr --concise"
INEXACT_ARGS="-1"

make -C ~/workspace/bowtie ebwt_search
cp ~/workspace/bowtie/ebwt_search .

for t in 0 5 10 ; do

TRIM_ARGS="-3 $t"

echo "Doing whole genome trimming $t bases..."
if [ ! -f whole.sim_reads.t$t.hits ] ; then
    echo "  Exact..."
    sh wrap.sh whole.sim_reads.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} \
	/fs/szasmg/langmead/ebwts/whole \
	whole.sim_reads.fa \
	whole.sim_reads.t$t.hits 2>&1 | tee whole.sim_reads.t$t.out
fi
if [ ! -f whole.sim_reads.1.t$t.hits ] ; then
    echo "  Inexact..."
    sh wrap.sh whole.sim_reads.1.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
	/fs/szasmg/langmead/ebwts/whole \
	whole.sim_reads.fa \
	whole.sim_reads.1.t$t.hits 2>&1 | tee whole.sim_reads.1.t$t.out
fi

echo "Doing half-genome trimming $t bases..."
if [ ! -f whole.1435609516.sim_reads.t$t.hits ] ; then
    echo "  Exact..."
    sh wrap.sh whole.1435609516.sim_reads.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.1435609516 \
	whole.1435609516.sim_reads.fa \
	whole.1435609516.sim_reads.t$t.hits 2>&1 | tee whole.1435609516.sim_reads.t$t.out
fi
if [ ! -f whole.1435609516.sim_reads.1.t$t.hits ] ; then
    echo "  Inexact..."
    sh wrap.sh whole.1435609516.sim_reads.1.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.1435609516 \
	whole.1435609516.sim_reads.fa \
	whole.1435609516.sim_reads.1.t$t.hits 2>&1 | tee whole.1435609516.sim_reads.1.t$t.out
fi

echo "Doing one-quarter-genome trimming $t bases..."
if [ ! -f whole.717804758.sim_reads.t$t.hits ] ; then
    echo "  Exact..."
    sh wrap.sh whole.717804758.sim_reads.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.717804758 \
	whole.717804758.sim_reads.fa \
	whole.717804758.sim_reads.t$t.hits 2>&1 | tee whole.717804758.sim_reads.t$t.out
fi
if [ ! -f whole.717804758.sim_reads.1.t$t.hits ] ; then
    echo "  Inexact..."
    sh wrap.sh whole.717804758.sim_reads.1.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.717804758 \
	whole.717804758.sim_reads.fa \
	whole.717804758.sim_reads.1.t$t.hits 2>&1 | tee whole.717804758.sim_reads.1.t$t.out
fi

echo "Doing one-eighth-genome trimming $t bases..."
if [ ! -f whole.358902379.sim_reads.t$t.hits ] ; then
    echo "  Exact..."
    sh wrap.sh whole.358902379.sim_reads.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.358902379 \
	whole.358902379.sim_reads.fa \
	whole.358902379.sim_reads.t$t.hits 2>&1 | tee whole.358902379.sim_reads.t$t.out
fi
if [ ! -f whole.358902379.sim_reads.1.t$t.hits ] ; then
    echo "  Inexact..."
    sh wrap.sh whole.358902379.sim_reads.1.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.358902379 \
	whole.358902379.sim_reads.fa \
	whole.358902379.sim_reads.1.t$t.hits 2>&1 | tee whole.358902379.sim_reads.1.t$t.out
fi

echo "Doing one-sixteenth-genome trimming $t bases..."
if [ ! -f whole.179451189.sim_reads.t$t.hits ] ; then
    echo "  Exact..."
    sh wrap.sh whole.179451189.sim_reads.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.179451189 \
	whole.179451189.sim_reads.fa \
	whole.179451189.sim_reads.t$t.hits 2>&1 | tee whole.179451189.sim_reads.t$t.out
fi
if [ ! -f whole.179451189.sim_reads.1.t$t.hits ] ; then
    echo "  Inexact..."
    sh wrap.sh whole.179451189.sim_reads.1.t$t
    ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
	/fs/szasmg/langmead/ebwts/whole.179451189 \
	whole.179451189.sim_reads.fa \
	whole.179451189.sim_reads.1.t$t.hits 2>&1 | tee whole.179451189.sim_reads.1.t$t.out
fi

done
