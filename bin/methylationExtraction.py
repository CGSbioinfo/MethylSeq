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


def methylationExtraction(i, ai):
    params=ai['methyl_extract_params'].replace(';','')
    alignedReads = os.listdir(in_dir + "/")
    alignedReads = [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall(i, x)]
    if args.dedup == True:
        bamFile = in_dir + [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall("bismark_bt2_.*deduplicated.bam$", x)][0]
    else:
        bamFile = in_dir + [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall("val_1_bismark_bt2_pe.bam$", x)][0]
    os.system("bismark_methylation_extractor  " + params + " --output " + out_dir + " " + bamFile + " &>" + out_dir + "/methylExtract_log_"+i+".txt")


#########################

__version__='v01'
# created 18/08/2016

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='methylationExtraction.py',description = 'Methylation extraction from bamfiles.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing bam files. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names.txt', default='sample_names.txt')
    parser.add_argument('--dedup', action='store_true')
    parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    args=parser.parse_args()

    # Set path of working directory
    ai=functions.read_analysis_info_file(args.analysis_info_file)
    path=os.getcwd()

    #Ncores
    ncores=int(args.ncores)

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories if not 'rawReads/'
    in_dir= args.in_dir
    out_dir= args.out_dir

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Run bismark alignment
    functions.make_sure_path_exists(out_dir)
    Parallel(n_jobs=ncores)(delayed(methylationExtraction)(i,ai) for i in sampleNames)

