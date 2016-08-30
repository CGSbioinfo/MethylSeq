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
# created on 30/08/2016

if __name__ == '__main__' :
	# Script information
	parser = argparse.ArgumentParser(prog='runMethylationAnalysis.py',description = '')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)
    
    # Analysis information
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--sample_names_file', help='File with sample names. Default=sample_names.txt', default='sample_names.txt')

    # Choose a section of the pipeline to run
    parser.add_argument('--run', help='Choose a section of the pipeline to run. Possible options: ', default='')

    # Bcl2fastq
    

