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

def mbias_make_tables(i, ai,allFiles):
    sample= in_dir + [allFiles[y] for y,x in enumerate(allFiles) if re.findall(i,x)][0]
    os.system("python bin/split_Mbias_data.py " + sample + " all " + out_dir + '/' + i )


#########################

__version__='v01'
# created 24/08/2016

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='trimmingReads.py',description = 'QC and adapter trimming using Trim Galore')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing trimmed fastq files. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=alignedReads/QC/', default='alignedReads/QC/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names.txt', default='sample_names.txt')
    parser.add_argument('--pattern', help='String to match with the name of the file. Default=bismark_bt2_pe.deduplicated.M-bias.txt', default='bismark_bt2_pe.deduplicated.M-bias.txt')
    parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    args=parser.parse_args()

    # Get path of working directory
    ai=functions.read_analysis_info_file(args.analysis_info_file)
    path=os.getcwd()

    #Ncores
    ncores=int(args.ncores)

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories if not 'rawReads/'
    in_dir=path + '/' + args.in_dir
    out_dir=path + '/' + args.out_dir

    # Run table 
    functions.make_sure_path_exists(out_dir)
    allFiles=os.listdir(in_dir+'/')
    allFiles=[allFiles[y] for y,x in enumerate(allFiles) if re.findall(args.pattern,x)]
    Parallel(n_jobs=ncores)(delayed(mbias_make_tables)(i,ai,allFiles) for i in sampleNames)

