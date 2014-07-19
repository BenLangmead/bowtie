#!/usr/bin/env python

import os
import logging
import subprocess


class BowtieSuite(object):
    
    def __init__(self,bt_path=''):
        curr_path           = os.path.realpath(bt_path)
        self.bowtie_bin     = os.path.join(curr_path,'bowtie')
        self.bowtie_build   = os.path.join(curr_path,'bowtie-build')
        self.bowtie_inspect = os.path.join(curr_path,'bowtie-inspect')


    def run(self, *args):
        cmd = self.bowtie_bin + " " + " ".join([i for i in args])
        return(subprocess.call(cmd,shell=True))        
        
        
    def silent_run(self, *args):
        cmd = self.bowtie_bin + " " + " ".join([i for i in args])
        return(subprocess.call(cmd,shell=True,stderr=open(os.devnull, 'w')))        
        
        
    def build(self, *args):
        cmd = self.bowtie_build + " " + " ".join([i for i in args])
        curr_dir = os.getcwd()
        os.chdir(os.path.dirname(cmd.split()[-2]))
        ret = subprocess.check_call(cmd,shell=True)
        os.chdir(curr_dir)
        return(ret)
        
    

