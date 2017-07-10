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

__version__='v02'
# created 17/08/2016
# modified 2017-06-13 BS - will move files to existing empty raw reads folder

def moveReads(sampleNames, fastq):
    ''' Move the sample reads

    Move the fastq files for the given sample names
    into the rawReads folder.

    '''
    functions.make_sure_path_exists('rawReads')
    sampleDir = []
    for sample in sampleNames:
        reads = [fastq[i] for i,x in enumerate(fastq) if re.findall(sample,x)]
        if sample not in sampleDir:
            functions.make_sure_path_exists('rawReads/'+sample)
        for r in reads:
            os.system('mv ' + '"' + r + '"' + ' rawReads/' + sample)
        sampleDir.append(sample)

if __name__ == '__main__':

    # Parser
    parser = argparse.ArgumentParser(prog='organizeWorkingDirectory.py',description = 'Organize working directory of the analysis')
    parser.add_argument('-v','--version', action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    parser.add_argument('--in_dir', help='directory with fastq files. Default= corresponding to bcl2fastq_output', default='bcl2fastq_output')
    args=parser.parse_args()

    # Read analysis info file
    ai=functions.read_analysis_info_file(args.analysis_info_file)

    # Read sample names
    sample_names_file = args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)
        
    # Check if rawReads exists
    folders = os.listdir('.')
    readsFiles = [folders[i] for i, x in enumerate(folders) if re.findall('rawReads',x)]
    
     # Collect fastq files analysis_info_file
    if args.in_dir == 'bcl2fastq_output':
        allFiles=functions.get_filepaths(ai[args.in_dir])
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[fastq[y] for y,x in enumerate(fastq) if not re.findall('Undetermined', x)]
    elif args.in_dir != 'bcl2fastq_output':
        allFiles=os.listdir(args.in_dir)
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[args.in_dir + x for x in fastq]

    


    # Move reads
    if not readsFiles:
        # The rawReads directory does not exist.
        # Create the directory and move the files.
        print "rawReads/ folder does not exist. Creating folder and moving files."
        moveReads(sampleNames, fastq)
    else:
        # If the rawReads directory exists and is empty
        # - for example if step 1 failed before bcl2fastq 
        # completed, then we need to continue moving files. 
        # Otherwise fail with meaningful error.
        if os.listdir('rawReads') == []: 
            print "rawReads/ folder exists and is empty. Moving files."
            moveReads(sampleNames, fastq)
        else:
            print "rawReads/ folder exists and contains files. Not moving files."

