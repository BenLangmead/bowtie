#!/bin/sh

cvs -d :ext:langmead@junkfood.cs.umd.edu:/fs/pugh/langmead/cvs/research co CSAMapper
cvs -d :ext:langmead@junkfood.cs.umd.edu:/fs/pugh/langmead/cvs/research co SeqAn-1.0

rm -rf SeqAn-1.0/doc
cp CSAMapper/README .
cp CSAMapper/papers/CMSC858P_Report.pdf .

tar jcvf langmead_proj2.tar.bz --exclude '*CVS*' *
