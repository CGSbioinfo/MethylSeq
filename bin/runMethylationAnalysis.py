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

if __name__ == '__main__':
	# Script information
	parser = argparse.ArgumentParser(prog='runMethylationAnalysis.py',description = '')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)

	# Analysis information
	parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
	parser.add_argument('--sample_names_file', help='File with sample names. Default=sample_names.txt', default='sample_names.txt')

	# Choose a section of the pipeline to run
	parser.add_argument('--run', help='Choose a section of the pipeline to run. Possible options: step1_prepare_analysis, step2_qc_and_trimming', default='')

	# Parse arguments
	args=parser.parse_args()

	# Step 1, prepare analysis
	if args.run == 'step1_prepare_analysis':

		print "Started bcl2fastq"
		os.system("python bin/run_bcl2fastq.py ")
		print "Finished bcl2fastq. Output is in the fastq/ folder"

		print "Started creating sample names file"
		os.system("python bin/create_sampleNames.py ")
		print "Created sample_names_file"

		print "Started organizing working directory"
		os.system("python bin/organizeWorkingDirectory.py ")
		print "Finished organizing working directory"

	if args.run == 'step2_qc_and_trimming':
		print "Running fastqc"
		os.system("python bin/qcReads.py ")
		print "Finished fastqc"

		print "Creating plots and tables raw data"
		os.system("python bin/fastqc_tables_and_plots.py --suffix_name _raw --out_dir_report Report/figure/rawQC --plot_device pdf")
		print "Finished plots and tables raw data"

		print "Running trimgalore"
		os.system("python bin/trimmingReads.py --in_dir rawReads/ --out_dir trimmedReads/")
		print "Finished trimgalore"

		print "Creating plots and tables trimmed data"
		os.system("python bin/fastqc_tables_and_plots.py --in_dir trimmedReads/ --out_dir trimmedReads/ --suffix_name _trimmed --out_dir_report Report/figure/trimmedQC --plot_device pdf")
		print "Finished plots and tables trimmed data"




