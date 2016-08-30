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

    parser=argparse.ArgumentParser(prog='analysis_info.py', description='Creates analysis_info.txt')
    parser.add_argument('-v','--version',action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--outfile', help='Name of the output file. Default=analysis_info.txt', default='analysis_info.txt')
    args=parser.parse_args()

    outfile_name=args.outfile
    lines=['working_directory = ', 'run_folder = ', 'run_samplesheet = ', 'bcl2fastq_output = fastq/ ', 'readType = pairedEnd', 'reference_genome = ', 'bismark_params = --bowtie2; --bam; -N 0; -L 20; -D 15; -R 2; --score_min L,0,-0.2; ', 'methyl_extract_params= --bedGraph; --gzip; --merge_non_CpG;', 'target_regions_bed = ', 'ncores = 8']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()
    
