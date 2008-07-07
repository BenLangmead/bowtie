#!/usr/bin/env python
# encoding: utf-8
"""
polymorph

Created by Cole Trapnell on 2007-04-26.
Copyright (c) 2007 __MyCompanyName__. All rights reserved.
"""

import sys
import getopt
import random
import math

help_message = """
polymorph reads a multi-fasta file and for each sequence, prints to standard 
output a mutated version of it. Prints a transcript of polymorphims in 1-based
coords

Usage:
    polymorph [options] <multi-fasta reference> 
    
Options:
    [TODO]
"""


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hm:s:vm:d:i:", 
                ["help", "seed=", "mismatch-rate", "deletion-rate", "insertion-rate"])
        except getopt.error, msg:
            raise Usage(msg)
    
        mismatch_rate = 0.0001
        insert_rate = 0.0001
        delete_rate = 0.0001
        
        seed_val = 0
        sorted = False
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-s", "--seed"):
                seed_val = long(value)
            if option in ("-m", "--mismatch-rate"):
                mismatch_rate = float(value)
            if option in ("-d", "--deletion-rate"):
                delete_rate = float(value)
            if option in ("-i", "--insertion-rate"):
                insert_rate = float(value)
            
        random.seed(seed_val)
    
        fasta_input = args[0]
        
        f = open(fasta_input, "r")
        
        seqs = {}
        curr_seq = []
        curr_seq_name = None
        for line in f:
            if line[0] == ">":
                if curr_seq_name != None:
                    seqs[curr_seq_name] = curr_seq
                    curr_seq = []
                    curr_seq_name = None
                curr_seq_name = line.strip()
            else:
                curr_seq.extend(list(line.strip()))
                
        if curr_seq_name != None:
            seqs[curr_seq_name] = curr_seq
            curr_seq = []
            curr_seq_name = None
        
        bases = ['A','C','T','G']
        base_mismatcher = {
            "A" : ['C','T','G'],
            "C" : ['A','T','G'],
            "G" : ['A','C','T'],
            "T" : ['A','C','G']
        }
        
        for (name, seq) in seqs.iteritems():
            i = 0

            num_ins = math.floor(insert_rate * len(seq))
            num_dels = math.floor(delete_rate * len(seq))
            num_mismats = math.floor(mismatch_rate * len(seq))
            
            name = name[1:]
            
            # while num_ins > 0:
            #   ins = random.randint(0, len(seq))
            #   seq.insert(ins, random.choice(bases))
            #   num_ins -= 1
            #   
            # while num_dels > 0:
            #   deletion = random.randint(0, len(seq))
            #   del seq[deletion]
            #   num_dels -= 1
            while num_mismats > 0:
                mis = random.randint(0, len(seq))
                print >> sys.stderr, "MISMATCH\t%s\t%d" % (name, mis + 1)
                seq[mis] = random.choice(base_mismatcher[seq[mis]])
                num_mismats -= 1
            
            
            print name, "(polymorphed)"
            seq = "".join(seq)
            while i < len(seq):
                print seq[i:i + 60]
                i += 60

            
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2


if __name__ == "__main__":
    sys.exit(main())
