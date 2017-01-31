#!/usr/bin/env python 

import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import multiprocessing
import subprocess
#sys.path.insert(0,'/usr/local/bin/')
import functions
import argparse

__version__ = 'v01'
# created on 17/08/2016

if __name__ == '__main__':

    """ This script creates a file which needs to be filled with information required for a methylSeq project.
    - It takes one argument, the 'outfile', which is the name of the output file. The default is 'analysis_info.txt'""" 

    parser=argparse.ArgumentParser(prog='remove_bases_file_info.py', description='Creates remove_bases_file.txt')
    parser.add_argument('-v','--version',action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--outfile', help='Name of the output file.')
    args=parser.parse_args()

    outfile_name=args.outfile
    lines=['sample\t5R1\t3R1\t5R2\t3R2']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()
    
