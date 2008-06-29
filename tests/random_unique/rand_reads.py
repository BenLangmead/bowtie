#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Cole Trapnell on 2008-06-19.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt
import random

help_message = '''
The help message goes here.
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
        num_reads = 0
        read_len = 0
        if len(args) >= 2:
            num_reads = int(args[0])
            read_len = int(args[1])
        else:
            raise Usage(help_message)
        bases = ['A','C','G','T']
        
        for i in range(0, num_reads):
            read = []
            for j in range(0, read_len):
                read.append(random.choice(bases))
            print ">%d" % i
            seq = "".join(read)
            j = 0
            while j + 60< len(seq):
                print seq[j:j + 60]
                j += 60
            print seq[j:]
                
            
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
