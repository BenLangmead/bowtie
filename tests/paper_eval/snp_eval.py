#!/usr/bin/env python
# encoding: utf-8
"""
snp_eval.py

Accepts a transcript from polymorph and a maq .snp file and then calculates
the sensitivity and specificity of the SNP calls.

Created by Cole Trapnell on 2008-07-04.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt


help_message = '''
Accepts a transcript from polymorph and a maq .snp file and then calculates
the sensitivity and specificity of the SNP calls.

Usage:
    snp_eval.py <polymorph_transcript> <cns.snp>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
        
        if len(args) != 2:
            raise(Usage(help_message))
    
        true_snps = open(args[0])
        pred_snps = open(args[1])
        
        snps = set([])
        
        fp = []
        tp = []
        
        for line in true_snps.readlines():
            cols = line.strip().split()
            seq = cols[1]
            pos = int(cols[2])
            snps.add((seq,pos))
        
        for line in pred_snps.readlines():
            cols = line.strip().split()
            seq = cols[0]
            pos = int(cols[1])
            #print (seq, pos)
            if (seq, pos) in snps:
                tp.append((seq,pos))
            else:
                fp.append((seq,pos))
        
        tp = set(tp)
        fp = set(fp)
        
        print "SNP calls in", args[1]
        print len(tp), "true positive predictions"
        print len(fp), "false positive predictions"
        print len(snps), "Possible true positives"
        if len(snps):
            print "sensitivity:", len(tp) / float(len(snps))
        else:
            print "No true snps in donor!"
            
        if len(tp) + len(fp):    
            print "specificity:", len(tp) / float(len(tp) + len(fp))
        else:
            print "No predictions made!"
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
