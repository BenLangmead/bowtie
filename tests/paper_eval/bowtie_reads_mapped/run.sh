#!/bin/sh

BOWTIE_HOME=$HOME/workspace/bowtie
KG_READS=/fs/szasmg/langmead/reads/SRR001115/s_7_0000_0255
BOWTIE_ARGS="-1tfra --concise --arrows"
BOWTIE_MAQ_ARGS="--maq -tfr --concise"
cp ${BOWTIE_HOME}/scripts/arrows_to_stats.pl .

make -C ${BOWTIE_HOME} ebwt_search
cp ${BOWTIE_HOME}/ebwt_search .

# Run Bowtie in 1-mismatch end-to-end mode with arrow reporting; this
# gives us the information needed to draw the bottom 4 lines 
./ebwt_search -3  0 ${BOWTIE_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.40.1tfra.arrows.hits
./ebwt_search -3  5 ${BOWTIE_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.35.1tfra.arrows.hits
./ebwt_search -3 10 ${BOWTIE_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.30.1tfra.arrows.hits
./ebwt_search -3 15 ${BOWTIE_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.25.1tfra.arrows.hits

# Process the arrows down to some statistics
cat ebwt.40.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
  > ebwt.40.1tfra.arrows.stats
cat ebwt.35.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
  > ebwt.35.1tfra.arrows.stats
cat ebwt.30.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
  > ebwt.30.1tfra.arrows.stats
cat ebwt.25.1tfra.arrows.hits | perl arrows_to_stats.pl -s -r 8956597 \
  > ebwt.25.1tfra.arrows.stats

# Run Bowtie in Maq mode; this gives us the information needed to draw
# the 5th line
./ebwt_search -3  0 ${BOWTIE_MAQ_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.40.maq.tfr.hits
./ebwt_search -3  5 ${BOWTIE_MAQ_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.35.maq.tfr.hits
./ebwt_search -3 10 ${BOWTIE_MAQ_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.30.maq.tfr.hits
./ebwt_search -3 15 ${BOWTIE_MAQ_ARGS} \
    ../../whole.t10o5l6 ${KG_READS}.fa ebwt.25.maq.tfr.hits
