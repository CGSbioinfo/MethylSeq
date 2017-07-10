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

__version__ = 'v02'
# created on 17/08/2016

if __name__ == '__main__':

    '''This script runs bcl2fastq.

    If uses the parameters in the analysis_info file.analysis_info_file.
    If bcl2fastq fails, the script will quit with exit code 1, otherwise
    it will quit with exit code 0
    '''

    parser=argparse.ArgumentParser(prog='analysis_info.py', description='Creates analysis_info.txt')
    parser.add_argument('-v','--version',action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    args=parser.parse_args()

    # Collect info from analysis_info_file
    ai=functions.read_analysis_info_file(args.analysis_info_file)

    functions.runAndCheck("bcl2fastq -R " + ai['run_folder'] 
        + " -o " + ai['bcl2fastq_output'] 
        + " --no-lane-splitting -l NONE --sample-sheet " + ai['run_samplesheet'])
    # functions.runAndCheck("echo 'Moo'", "Moose")
