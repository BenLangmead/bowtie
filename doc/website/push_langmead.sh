#!/bin/sh

# If you're me (sourceforge: ben_langmead), use this to push the
# contents of this directory to the appropriate place on sourceforge.
# Changes will be instantaneous on http://bowtie-bio.sf.net.
scp *.ssi *.shtml *.html *.txt ben_langmead,bowtie-bio@web.sourceforge.net:/home/groups/b/bo/bowtie-bio/htdocs
scp -r css images ben_langmead,bowtie-bio@web.sourceforge.net:/home/groups/b/bo/bowtie-bio/htdocs
