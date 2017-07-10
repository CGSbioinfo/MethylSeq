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

def alignment(sampleName, ai):
    '''Perform alignment of reads using Bismark params

    This function check that read 1 and read 2 files exist
    and runs bismark using the parameters detailed in the 
    analysis_info file.

    sampleName -- the name of the sample to process
    ai -- the settings
    '''
    bismark_params=ai['bismark_params'].replace(';','')
    
    # Get all files in the input dir as a list
    trimmedReads = os.listdir(in_dir + "/")

    # Filter to filenames containing the sample name
    regex = re.compile(sampleName)
    trimmedReads = filter(regex.search, trimmedReads)

    r1regex = re.compile("_R1_.*val_1.fq"+gz)
    r2regex = re.compile("_R2_.*val_2.fq"+gz)

    r1 = filter(r1regex.search, trimmedReads)
    r2 = filter(r2regex.search, trimmedReads)

    # Ensure read 1 is present
    if not r1:
        print "No read 1 file found. Ignoring %s." % sampleName
        return

    r1 = in_dir + r1[0]

    # Check read 2 is present when needed
    if(ai['readType']=='pairedEnd'):
        if not r2:
            print "No read 2 file found when readType is paired end. Ignoring %s." % sampleName
            return

    if r2: 
        r2 = in_dir + r2[0]
        cmdStr = ('srun -c %i bismark %s --output_dir %s %s -1 %s -2 %s &>%s/bismark_log_%s.txt' 
            % (ncores, bismark_params, out_dir, ai['reference_genome'], r1, r2, out_dir, sampleName ))
        print "Invoking: %s" % cmdStr
        functions.runAndCheck(cmdStr, "Error in bismark read mapping")
    else:
        print "No read 2 for %s." % sampleName
        print "Single reads not currently implemented"
        # r1 = in_dir + [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1_.*val_1.fq", x)][0]
        # print "Found "+str(len(r1))+" r1 files"
        # cmdStr = str("srun -c %(ncores)i -mem=20G bismark " + bismark_params + " --output_dir " + out_dir + " " + ai['reference_genome'] +" -1 " + r1 + " -2 " + r2 + " &>" + out_dir + "/bismark_log_"+sampleName+".txt ")
        # print "Invoking: "+cmdStr
        # functions.runAndCheck(cmdStr, "Error in bismark read mapping")

def deduplicate(sampleName):
    '''Perform deduplication of reads using Bismark params.

    The function checks if the required bam file exists before running.
    sampleName -- the name of the sample to process
    '''
    alignedReads = os.listdir(out_dir + "/")

     # Filter to filenames containing the sample name
    regex = re.compile(sampleName)
    alignedReads = filter(regex.search, alignedReads)

    if not alignedReads:
        print "No aligned reads found for %s." % sampleName
        return

    regex = re.compile("bismark_bt2_.*.bam$")

    bamFiles = filter(regex.search, alignedReads)
    if not bamFiles:
        print "No bam file found for %s." % sampleName
        return

    bamFile = out_dir + bamFiles[0]

    cmdStr = ('srun -c %i deduplicate_bismark --bam %s &>%s/dedup_log_%s.txt' 
            % (ncores, bamFile, out_dir, sampleName ))

    print "Invoking: "+cmdStr
    alignedReads = [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall(sampleName, x)]
    # bamFile = out_dir + [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall( "bismark_bt2_.*.bam$", x)][0]
    functions.runAndCheck(cmdStr, "Error deduplicating")

#########################

__version__='v01'
# created 18/08/2016

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='mappingReads.py',description = 'Read mapping using Bismark')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing trimmed fastq files. Default=trimmedReads/', default='trimmedReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names.txt', default='sample_names.txt')
    # parser.add_argument('--ncores', help='Number of cores to use. Default=3', default='3')
    parser.add_argument('--run', help='Choose a section of the script to run. Possible options: bismark_alignment; deduplicate_bismark; all. Default=all', default='all')
    args=parser.parse_args()

    # Set path of working directory
    ai=functions.read_analysis_info_file(args.analysis_info_file)
    path=os.getcwd()

    #Ncores
    ncores = int(ai['ncores'])
    ninstances=int(ai['ninstances'])

    # Read sample names text file
    sampleNames = functions.read_sample_names(args.sample_names_file)

    # Set input and output directories if not 'rawReads/'
    in_dir =path + '/' + args.in_dir
    out_dir=path + '/' + args.out_dir

    # Detect if files are gz
    gz = functions.check_gz(in_dir, "fq")
    print ("Appending gz string '%s'" % gz)

    # Display config options
    print ("Setting to run in batches of %i instances with %i cores per instance" % (ninstances, ncores))

    # Run bismark alignment
    if args.run == 'all' or args.run =='bismark_alignment': 
        print "Running bismark alignment"
        functions.make_sure_path_exists(out_dir)

        # Use Parallel to batch the task into groups of ninstances
        # alignments so the cluster is not overloaded.
        Parallel(n_jobs=ninstances)(delayed(alignment)(i, ai) for i in sampleNames)
        # Parallel(n_jobs=ncores)(delayed(alignment)(i, ai) for i in sampleNames)

    # Run deduplication
    if args.run == 'all' or args.run == 'deduplicate_bismark':
        print "Running deduplication"
        Parallel(n_jobs=ninstances)(delayed(deduplicate)(i) for i in sampleNames)
        # Parallel(n_jobs=ncores)(delayed(deduplicate)(i) for i in sampleNames)
