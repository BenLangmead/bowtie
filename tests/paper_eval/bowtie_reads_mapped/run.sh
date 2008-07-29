#!/bin/sh

BOWTIE_HOME=$HOME/workspace/bowtie
KG_READS=/fs/szasmg/langmead/reads/SRR001115/SRR001115/s_7_0000_0255
BOWTIE_ARGS="-1tfra --concise --arrows"
BOWTIE_MAQ_N1_ARGS="-n 1 -tfr --concise"
BOWTIE_MAQ_ARGS="-tfr --concise"
DO_CVS_UPDATE=0

# Optionally do a cvs update in the Bowtie home
if [ "$DO_CVS_UPDATE" = "1" ] ; then
	pushd ${BOWTIE_HOME}
	cvs update -d
	popd
fi

cp ${BOWTIE_HOME}/scripts/arrows_to_stats.pl .

make -C ${BOWTIE_HOME} ebwt_search
if [ ! -f ${BOWTIE_HOME}/ebwt_search ] ; then
   echo "Failed to build ebwt_search in ${BOWTIE_HOME}; aborting..."
   exit 1
fi
cp ${BOWTIE_HOME}/ebwt_search .

# Run Bowtie in 1-mismatch end-to-end mode with arrow reporting; this
# gives us the information needed to draw the bottom 4 lines 
if [ ! -f ebwt.40.1tfra.arrows.hits ] ; then
	./ebwt_search -3  0 ${BOWTIE_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.40.1tfra.arrows.hits
fi
if [ ! -f ebwt.35.1tfra.arrows.hits ] ; then
	./ebwt_search -3  5 ${BOWTIE_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.35.1tfra.arrows.hits
fi
if [ ! -f ebwt.30.1tfra.arrows.hits ] ; then
	./ebwt_search -3 10 ${BOWTIE_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.30.1tfra.arrows.hits
fi
if [ ! -f ebwt.25.1tfra.arrows.hits ] ; then
	./ebwt_search -3 15 ${BOWTIE_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.25.1tfra.arrows.hits
fi

# Process the arrows down to some statistics
if [ ! -f ebwt.40.1tfra.arrows.stats ] ; then
	cat ebwt.40.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
	  > ebwt.40.1tfra.arrows.stats
fi
if [ ! -f ebwt.35.1tfra.arrows.stats ] ; then
	cat ebwt.35.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
	  > ebwt.35.1tfra.arrows.stats
fi
if [ ! -f ebwt.30.1tfra.arrows.stats ] ; then
	cat ebwt.30.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
	  > ebwt.30.1tfra.arrows.stats
fi
if [ ! -f ebwt.25.1tfra.arrows.stats ] ; then
	cat ebwt.25.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
	  > ebwt.25.1tfra.arrows.stats
fi

# Run Bowtie in Maq mode (both -n 1 and the default of -n 2); this
# gives us the information needed to draw the 5th and 6th lines
if [ ! -f ebwt.40.maq.n1.tfr.hits ] ; then
	./ebwt_search -3  0 ${BOWTIE_MAQ_N1_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.40.maq.n1.tfr.hits
fi
if [ ! -f ebwt.35.maq.n1.tfr.hits ] ; then
	./ebwt_search -3  5 ${BOWTIE_MAQ_N1_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.35.maq.n1.tfr.hits
fi
if [ ! -f ebwt.30.maq.n1.tfr.hits ] ; then
	./ebwt_search -3 10 ${BOWTIE_MAQ_N1_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.30.maq.n1.tfr.hits
fi
if [ ! -f ebwt.25.maq.n1.tfr.hits ] ; then
	./ebwt_search -3 15 ${BOWTIE_MAQ_N1_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.25.maq.n1.tfr.hits
fi

if [ ! -f ebwt.40.maq.tfr.hits ] ; then
	./ebwt_search -3  0 ${BOWTIE_MAQ_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.40.maq.tfr.hits
fi
if [ ! -f ebwt.35.maq.tfr.hits ] ; then
	./ebwt_search -3  5 ${BOWTIE_MAQ_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.35.maq.tfr.hits
fi
if [ ! -f ebwt.30.maq.tfr.hits ] ; then
	./ebwt_search -3 10 ${BOWTIE_MAQ_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.30.maq.tfr.hits
fi
if [ ! -f ebwt.25.maq.tfr.hits ] ; then
	./ebwt_search -3 15 ${BOWTIE_MAQ_ARGS} \
	    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.25.maq.tfr.hits
fi
    