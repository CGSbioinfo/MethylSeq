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


def alignment(i, ai):
    # Bismark params
    bismark_params=ai['bismark_params'].replace(';','')
    trimmedReads = os.listdir(in_dir + "/")
    trimmedReads = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall(i, x)]
    r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2_.*val_2.fq", x)]
    if r2:
        r1 = in_dir + [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1_.*val_1.fq", x)][0]
        r2 = in_dir + [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2_.*val_2.fq", x)][0]
        os.system("bismark " + bismark_params + " --output_dir " + out_dir + " " + ai['reference_genome']    +" -1 " + r1 + " -2 " + r2 + " &>" + out_dir + "/bismark_log_"+i+".txt ")
    else:
        r1 = in_dir + [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1_.*val_1.fq", x)][0]
        os.system("bismark " + bismark_params + " --output_dir " + out_dir + " " + ai['reference_genome'] +" -1 " + r1 + " -2 " + r2 + " &>" + out_dir + "/bismark_log_"+i+".txt")


def deduplicate(i):
    alignedReads = os.listdir(out_dir + "/")
    alignedReads = [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall(i, x)]
    bamFile = out_dir + [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall( "bismark_bt2_.*.bam$", x)][0]
    os.system("deduplicate_bismark --bam " + bamFile + " &>" + out_dir + "/dedup_log_"+i+".txt")

#########################

__version__='v01'
# created 18/08/2016

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='trimmingReads.py',description = 'QC and adapter trimming using Trim Galore')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing trimmed fastq files. Default=trimmedReads/', default='trimmedReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names.txt', default='sample_names.txt')
    parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    parser.add_argument('--run', help='Choose a section of the script to run. Possible options: bismark_alignment; deduplicate_bismark; all. Default=all', default='all')
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
    in_dir=path + '/' + args.in_dir
    out_dir=path + '/' + args.out_dir

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Run bismark alignment
    if args.run == 'all' or args.run =='bismark_alignment': 
        functions.make_sure_path_exists(out_dir)
        Parallel(n_jobs=ncores)(delayed(alignment)(i, ai) for i in sampleNames)

    # Run deduplication
    if args.run == 'all' or args.run == 'deduplicate_bismark':
        Parallel(n_jobs=ncores)(delayed(deduplicate)(i) for i in sampleNames)
