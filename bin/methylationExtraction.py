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


def methylationExtraction(i, ai, remove_bases_dict):
    '''
    Run the bismark methylation extraction on the given sample.in_dir

    If a non-empty dictionary is supplied in remove_bases_dict, the corresponding
    bases will be trimmed
    '''
    params=ai['methyl_extract_params'].replace(';','')
    alignedReads = os.listdir(in_dir + "/")
    alignedReads = [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall(i, x)]
    if not alignedReads:
        print "No aligned reads found for %s." % i
        return

    # Filter to filenames containing the sample name
    dedupregex = re.compile("bismark_bt2_.*deduplicated.bam$")
    bamregex   = re.compile("val_1.fq.gz_bismark_bt2_pe.bam$")

    dedup = filter(dedupregex.search, alignedReads)
    bams = filter(bamregex.search, alignedReads)

    if args.dedup == True:
        if not dedup:
            print "No deduplicated bam files found for %s" % i
            return;
        else: 
            bamFile = dedup[0]
    else:
        if not bams:
            print "No bam files found for %s" % i
            return;
        else:
            bamFile = bams[0]
    
    if not bamFile:
        print "No bam file to process for %s." % i
        return

    bamFile = in_dir + bamFile

    # Create the ignore string if needed
    if(remove_bases_dict):
        basesInfo=remove_bases_dict[i]
        logFile = ("%s/methylExtract_REMOVEDBASES_log_%s.txt"
            % (out_dir, i))
        ignoreString = (" --ignore %s --ignore_3prime %s --ignore_r2 %s --ignore_3prime_r2 %s "
            % (basesInfo[0], basesInfo[1], basesInfo[2], basesInfo[3]))
    else:
        logFile = ("%s/methylExtract_log_%s.txt"
            % (out_dir, i))
        ignoreString = ""

    
    cmdStr = ('srun -c %i bismark_methylation_extractor %s --output %s %s %s &>%s'
            % (ncores, params, out_dir, ignoreString, bamFile, logFile ))
    print("Extracting methylation from "+i)
    functions.runAndCheck(cmdStr, "Error extracting methylation")


def createRemoveBaseDict(remove_bases_file):
    '''
    Read the given file and make a dictionary 
    of bases to remove.

    If the remove bases file does not exist, an
    empty dictionary will be returned. 
    '''
    remove_bases_dict = {}
    if os.path.exists(remove_bases_file):
        with open(remove_bases_file, 'r') as f:
            next(f)
            for line in f:
                (key, ignore, ignore_3prime, ignore_r2, ignore_3prime_r2)=line.split()
                remove_bases_dict[key]=(ignore, ignore_3prime, ignore_r2, ignore_3prime_r2)
    return remove_bases_dict

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
    # parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    parser.add_argument('--remove_bases_file', help='Text file with bases to remove from each sample. Default=mbias_remove_bases.txt', default='mbias_remove_bases.txt')
    args=parser.parse_args()

    # Set path of working directory
    ai=functions.read_analysis_info_file(args.analysis_info_file)

    path = ai['working_directory']
    # path=os.getcwd()

    #Ncores
    ncores=int(ai['ncores'])
    ninstances=int(ai['ninstances'])

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories
    in_dir = path + args.in_dir
    out_dir= path + args.out_dir

    functions.make_sure_path_exists(out_dir)

    if not os.path.exists(in_dir):
        print "Input directory was not found: " + in_dir
        sys.exit(1)

    # Read remove bases file
    remove_bases_dict = createRemoveBaseDict(remove_bases_file)

    Parallel(n_jobs=ninstances)(delayed(methylationExtraction)(i,ai, remove_bases_dict) for i in sampleNames)

