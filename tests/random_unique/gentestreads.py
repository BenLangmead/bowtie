#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Cole Trapnell on 2007-04-26.
Copyright (c) 2007 __MyCompanyName__. All rights reserved.
"""

import sys
import getopt
import random


help_message = """
"""


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def get_reference(ref_file_name):
	f = open(ref_file_name, "r")
	lines = f.readlines()
	
	if lines[0].find(">") == -1:
		raise Usage("File is not FASTA format")

	seq = "".join([line.strip() for line in lines[1:]])
	return seq
	
	
def combine(*seqin):
    '''returns a list of all combinations of argument sequences.
for example: combine((1,2),(3,4)) returns
[[1, 3], [1, 4], [2, 3], [2, 4]]'''
    def rloop(seqin,listout,comb):
        '''recursive looping function'''
        if seqin:                       # any more sequences to process?
            for item in seqin[0]:
                newcomb=comb+[item]     # add next item to current comb
                # call rloop w/ rem seqs, newcomb
                rloop(seqin[1:],listout,newcomb)
        else:                           # processing last sequence
            listout.append(comb)        # comb finished, add to list
    listout=[]                      # listout initialization
    rloop(seqin,listout,[])         # start recursive process
    return listout

def xcombine(*seqin):
    '''returns a generator which returns combinations of argument sequences
for example xcombine((1,2),(3,4)) returns a generator; calling the next()
method on the generator will return [1,3], [1,4], [2,3], [2,4] and
StopIteration exception.  This will not create the whole list of 
combinations in memory at once.'''
    def rloop(seqin,comb):
        '''recursive looping function'''
        if seqin:                   # any more sequences to process?
            for item in seqin[0]:
                newcomb=comb+[item]     # add next item to current combination
                # call rloop w/ remaining seqs, newcomb
                for item in rloop(seqin[1:],newcomb):   
                    yield item          # seqs and newcomb
        else:                           # processing last sequence
            yield comb                  # comb finished, add to list
    return rloop(seqin,[])
	
def apply_diff(orig, diff):
	(pos, diff_type) = diff
	new_seq = [orig[:pos]]
	
	base_mismatcher = {
		"A" : "C",
		"C"	: "G",
		"G" : "T",
		"T" : "A"
	}
		
	if diff_type == 'READ_MISMATCH':
		new_seq.append(base_mismatcher[orig[pos]])
		new_seq.append(orig[pos+1:])
	elif diff_type == 'READ_INSERT':
		new_seq.append('A')
		new_seq.append(orig[pos:])
	else:
		new_seq.append(orig[pos+1:])
	
	return ''.join(new_seq)	
	
def main(argv=None):

	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hs:v", ["help", "seed="])
		except getopt.error, msg:
			raise Usage(msg)
	
		seed_val = 0
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-s", "--seed"):
				seed_val = long(value)
			
		random.seed(seed_val)
		
		if len(args) != 2:
			raise Usage(help_message)
	
		fasta_input = args[0]
		read_length = int(args[1])
		
		seq = get_reference(fasta_input)
		
		L = len(seq)
		
		rid = 0

		reads = []
		diff_types = ['READ_MISMATCH','READ_INSERT','READ_DELETE']
		
		input_file = open(fasta_input + "_in.fna", "w")
		output_file = open(fasta_input + "_out.fna", "w")
					
		mods = combine(range(0,read_length), range(0,read_length))
		for (p1,p2) in mods:
			p1_mods = combine([p1],diff_types)
			p2_mods = combine([p2],diff_types)
			p1_p2_mods = combine(p1_mods, p2_mods)
			
			for ((p1pos,p1diff), (p2pos,p2diff)) in p1_p2_mods:
				if p1pos == p2pos:
					continue
				start = random.randint(0, L - read_length)
				end = start + read_length
				
				orig = seq[start:end]
				new = apply_diff(orig, (p1pos, p1diff))
				
				if p1diff == 'READ_DELETE':
					if p2pos > p1pos:
						p2pos -= 1
				elif p1diff == 'READ_INSERT':
					if p2pos > p1pos:
						p2pos += 1
				
				new = apply_diff(new, (p2pos, p2diff))
				
				if len(new) < read_length:
					new += seq[end:end+(read_length - len(new))]
				
				new = new[:read_length]
				
				rid +=1
				#print ">" + "rid" + str(rid), "@", start, "len = ", len(new), ((p1pos,p1diff), (p2pos,p2diff))
				#print new
				
				header = ">" + " rid" + str(rid) + " "+ str(((p1pos,p1diff), (p2pos,p2diff)))
				print >> input_file, header
				print >> input_file, new
				
				print >> output_file, header
				print >> output_file, "@ %d: 2 diffs" % (start + 1)
				
				assert len(new) >= 30 and len(new) <= 34
		
		input_file.close()
		output_file.close()
		# for i in range(0, num_reads):
		# 			start = random.randint(0, L - read_length)
		# 			end = start + read_length
		# 
		# 			rid +=1
		# 			#print ">" + "rid" + str(rid)
		# 			#print seq[start:end]
			
			
			
			
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		return 2


if __name__ == "__main__":
	sys.exit(main())
