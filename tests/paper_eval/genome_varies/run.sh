#!/bin/sh

BASE_ARGS="-qfr --concise"
EXACT_ARGS="-0"
ONEMM_ARGS="-1"
MAQ_ARGS="--maq"

make -C ~/workspace/bowtie ebwt_search
cp ~/workspace/bowtie/ebwt_search .

# For each 3' trim length
for t in 0 5 10 ; do

	TRIM_ARGS="-3 $t"
	
	# Whole
	
	echo "Doing whole genome trimming $t bases..."
	if [ ! -f whole.sim_reads.t$t.hits ] ; then
	    echo "  Exact..."
	    sh wrap.sh whole.sim_reads.t$t \
	       ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${EXACT_ARGS} \
		   /fs/szasmg/langmead/ebwts/whole \
		   whole.sim_reads.fq \
		   whole.sim_reads.t$t.hits 2>&1 | tee whole.sim_reads.t$t.out
	fi
	if [ ! -f whole.sim_reads.1.t$t.hits ] ; then
	    echo "  Inexact..."
	    sh wrap.sh whole.sim_reads.1.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole \
		  whole.sim_reads.fq \
		  whole.sim_reads.1.t$t.hits 2>&1 | tee whole.sim_reads.1.t$t.out
	fi
	if [ ! -f whole.sim_reads.m.t$t.hits ] ; then
	    echo "  Maq-like..."
	    sh wrap.sh whole.sim_reads.m.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${MAQ_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole \
		  whole.sim_reads.fq \
		  whole.sim_reads.m.t$t.hits 2>&1 | tee whole.sim_reads.m.t$t.out
	fi
	
	# 1/2
	
	echo "Doing half-genome trimming $t bases..."
	if [ ! -f whole.1435609516.sim_reads.t$t.hits ] ; then
	    echo "  Exact..."
	    sh wrap.sh whole.1435609516.sim_reads.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${EXACT_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole.1435609516 \
		  whole.1435609516.sim_reads.fq \
		  whole.1435609516.sim_reads.t$t.hits 2>&1 | tee whole.1435609516.sim_reads.t$t.out
	fi
	if [ ! -f whole.1435609516.sim_reads.1.t$t.hits ] ; then
	    echo "  Inexact..."
	    sh wrap.sh whole.1435609516.sim_reads.1.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole.1435609516 \
		  whole.1435609516.sim_reads.fq \
		  whole.1435609516.sim_reads.1.t$t.hits 2>&1 | tee whole.1435609516.sim_reads.1.t$t.out
	fi
	if [ ! -f whole.1435609516.sim_reads.m.t$t.hits ] ; then
	    echo "  Maq-like..."
	    sh wrap.sh whole.1435609516.sim_reads.m.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${MAQ_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole.1435609516 \
		  whole.1435609516.sim_reads.fq \
		  whole.1435609516.sim_reads.m.t$t.hits 2>&1 | tee whole.1435609516.sim_reads.m.t$t.out
	fi
	
	# 1/4
	
	echo "Doing one-quarter-genome trimming $t bases..."
	if [ ! -f whole.717804758.sim_reads.t$t.hits ] ; then
	    echo "  Exact..."
	    sh wrap.sh whole.717804758.sim_reads.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${EXACT_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole.717804758 \
		  whole.717804758.sim_reads.fq \
		  whole.717804758.sim_reads.t$t.hits 2>&1 | tee whole.717804758.sim_reads.t$t.out
	fi
	if [ ! -f whole.717804758.sim_reads.1.t$t.hits ] ; then
	    echo "  Inexact..."
	    sh wrap.sh whole.717804758.sim_reads.1.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole.717804758 \
		  whole.717804758.sim_reads.fq \
		  whole.717804758.sim_reads.1.t$t.hits 2>&1 | tee whole.717804758.sim_reads.1.t$t.out
	fi
	if [ ! -f whole.717804758.sim_reads.m.t$t.hits ] ; then
	    echo "  Maq-like..."
	    sh wrap.sh whole.717804758.sim_reads.m.t$t \
	      ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${MAQ_ARGS} \
		  /fs/szasmg/langmead/ebwts/whole.717804758 \
		  whole.717804758.sim_reads.fq \
		  whole.717804758.sim_reads.m.t$t.hits 2>&1 | tee whole.717804758.sim_reads.m.t$t.out
	fi
	
	# 1/8
	
	echo "Doing one-eighth-genome trimming $t bases..."
	if [ ! -f whole.358902379.sim_reads.t$t.hits ] ; then
	    echo "  Exact..."
	    sh wrap.sh whole.358902379.sim_reads.t$t \
	       ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${EXACT_ARGS} \
		   /fs/szasmg/langmead/ebwts/whole.358902379 \
		   whole.358902379.sim_reads.fq \
		   whole.358902379.sim_reads.t$t.hits 2>&1 | tee whole.358902379.sim_reads.t$t.out
	fi
	if [ ! -f whole.358902379.sim_reads.1.t$t.hits ] ; then
	    echo "  Inexact..."
	    sh wrap.sh whole.358902379.sim_reads.1.t$t \
	       ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
		   /fs/szasmg/langmead/ebwts/whole.358902379 \
		   whole.358902379.sim_reads.fq \
		   whole.358902379.sim_reads.1.t$t.hits 2>&1 | tee whole.358902379.sim_reads.1.t$t.out
	fi
	if [ ! -f whole.358902379.sim_reads.m.t$t.hits ] ; then
	    echo "  Maq-like..."
	    sh wrap.sh whole.358902379.sim_reads.m.t$t \
	       ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${MAQ_ARGS} \
		   /fs/szasmg/langmead/ebwts/whole.358902379 \
		   whole.358902379.sim_reads.fq \
		   whole.358902379.sim_reads.m.t$t.hits 2>&1 | tee whole.358902379.sim_reads.m.t$t.out
	fi

	# 1/16
	
	echo "Doing one-sixteenth-genome trimming $t bases..."
	if [ ! -f whole.179451189.sim_reads.t$t.hits ] ; then
	    echo "  Exact..."
	    sh wrap.sh whole.179451189.sim_reads.t$t ${EXACT_ARGS} \
	       ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} \
		   /fs/szasmg/langmead/ebwts/whole.179451189 \
		   whole.179451189.sim_reads.fq \
		   whole.179451189.sim_reads.t$t.hits 2>&1 | tee whole.179451189.sim_reads.t$t.out
	fi
	if [ ! -f whole.179451189.sim_reads.1.t$t.hits ] ; then
	    echo "  Inexact..."
	    sh wrap.sh whole.179451189.sim_reads.1.t$t \
	       ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${INEXACT_ARGS} \
		   /fs/szasmg/langmead/ebwts/whole.179451189 \
		   whole.179451189.sim_reads.fq \
		   whole.179451189.sim_reads.1.t$t.hits 2>&1 | tee whole.179451189.sim_reads.1.t$t.out
	fi
	if [ ! -f whole.179451189.sim_reads.m.t$t.hits ] ; then
	    echo "  Inexact..."
	    sh wrap.sh whole.179451189.sim_reads.m.t$t \
	       ./ebwt_search ${BASE_ARGS} ${TRIM_ARGS} ${MAQ_ARGS} \
		   /fs/szasmg/langmead/ebwts/whole.179451189 \
		   whole.179451189.sim_reads.fq \
		   whole.179451189.sim_reads.m.t$t.hits 2>&1 | tee whole.179451189.sim_reads.m.t$t.out
	fi

done
