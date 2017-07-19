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
import functions
import argparse

def trimming(i):
    '''
    Trim sequencing reads using trim_galore
    '''
    allFiles = os.listdir(in_dir + "/" + i )
    pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_R2", x)]
    if pairedReads_temp:
        functions.runAndCheck("srun -c "+ncores+" trim_galore --gzip --paired --fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + 
            in_dir + i + "/" + i + "*_R1*.fastq" + gz + " " + 
            in_dir + i + "/" + i + "*_R2*.fastq" + gz, "Error in trim_galore")
    else:
        functions.runAndCheck("srun -c "+ncores+" trim_galore --gzip --fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + in_dir + i + "/" + i + "*_R1*.fastq" + gz, "Error in trim_galore")


#########################

__version__='v02'
# created 18/08/2016

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='trimmingReads.py',
        description = 'QC and adapter trimming using Trim Galore')
    parser.add_argument('-v', '--version', 
        action='version', 
        version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', 
        help='Text file with details of the analysis. Default=analysis_info.txt', 
        default='analysis_info.txt')
    parser.add_argument('--in_dir', 
        help='Path to folder containing fastq files. Default=rawReads/', 
        default='rawReads/')
    parser.add_argument('--out_dir', 
        help='Path to out put folder. Default=trimmedReads/', 
        default='trimmedReads/')
    parser.add_argument('--out_dir_report', 
        help='Path to out put folder. Default=Report/figure/data/', 
        default='Report/figure/data/')
    parser.add_argument('--sample_names_file', 
        help='Text file with sample names. Default=sample_names_info.txt', 
        default='sample_names.txt')

    args=parser.parse_args()


    path=os.getcwd()
    os.chdir(path)
    ai=functions.read_analysis_info_file(args.analysis_info_file)    

    #Ncores
    ncores=int(ai['ncores'])
    ninstances=int(ai['ninstances'])

    # Read sample names text file
    sampleNames = functions.read_sample_names(args.sample_names_file)

    # Set input and output directories if not 'rawReads/'
    in_dir        =path + '/' + args.in_dir
    out_dir       =path + '/' + args.out_dir
    out_dir_report=path + '/' + args.out_dir_report

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Run trimm galore
    functions.make_sure_path_exists(out_dir)
    Parallel(n_jobs=ninstances)(delayed(trimming)(i) for i in sampleNames)
    functions.make_sure_path_exists(out_dir_report)
    
    # Generate summary stats
    cmd_str = ("srun /usr/bin/Rscript bin/trimming_summary.R %s %s %s " %
        (in_dir, out_dir, out_dir_report))
    functions.runAndCheck(cmd_str, "Error making trimming summary")

