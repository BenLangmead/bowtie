#!/bin/sh

echo "Short test emphasizing searching..."
perl scripts/random_bowtie_tests.pl 77 1000 200 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; fi
echo "Short test emphasizing searching with short patterns..."
perl scripts/random_bowtie_tests.pl 77 1000 200 300 200 6 6
if [ "$?" != "0" ] ; then echo "Error!" ; fi
echo "Short test emphasizing searching with long patterns..."
perl scripts/random_bowtie_tests.pl 77 1000 200 300 200 30 20
if [ "$?" != "0" ] ; then echo "Error!" ; fi
echo "Short test emphasizing building..."
perl scripts/random_bowtie_tests.pl 77 1200 10 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; fi

echo "Long test emphasizing searching..."
perl scripts/random_bowtie_tests.pl 77 5000 200 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; fi
echo "Long test emphasizing building..."
perl scripts/random_bowtie_tests.pl 77 6000 10 300 200 8 30
if [ "$?" != "0" ] ; then echo "Error!" ; fi
