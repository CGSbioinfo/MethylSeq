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
# created on 30/08/2016

def step_1(analysis_info_file):
    '''Run the first part of the analysis pipeline.

    Checks the info files, runs bcl2fastq, and creates the sample names file
    for future scripts.

    --analysis_info_file - the file of analysis options
    '''

    cmd_str = "python bin/validateAnalysisInfo.py --analysis_info_file %s" % analysis_info_file
    functions.runAndCheck(cmd_str, "Error in analysis info")

    cmd_str = "python bin/validateSampleSheet.py --analysis_info_file %s" % analysis_info_file
    functions.runAndCheck(cmd_str, "Error in in sample sheet")

    cmd_str = "python bin/run_bcl2fastq.py --analysis_info_file %s" % analysis_info_file
    functions.runAndCheck(cmd_str, "Error in bcl2fastq")

    print "Finished bcl2fastq. Output is in the fastq/ folder"

    functions.runAndCheck("python bin/create_sampleNames.py ", "Error creating samples names")

    print "Created sample_names_file"

    functions.runAndCheck("python bin/organizeWorkingDirectory.py ", "Error organizing working directory")

    print "Finished running step 1"

def step_2():
    '''Run the second part of the analysis pipeline.

    Runs fastqc, trims reads, and generates the QC plots for 
    raw and trimmed data.
    '''
    print "Running fastqc"
    functions.runAndCheck("python bin/qcReads.py ", "Error in fastqc")
    print "Finished fastqc"

    print "Creating plots and tables raw data"
    functions.runAndCheck("python bin/fastqc_tables_and_plots.py --suffix_name _raw --out_dir_report Report/figure/rawQC --plot_device pdf", "Error making plots")
    print "Finished plots and tables raw data"

    print "Running trimgalore"
    functions.runAndCheck("python bin/trimmingReads.py --in_dir rawReads/ --out_dir trimmedReads/", "Error in trimgalore")
    print "Finished trimgalore"

    print "Creating plots and tables for trimmed data"
    functions.runAndCheck("python bin/fastqc_tables_and_plots.py --in_dir trimmedReads/ --out_dir trimmedReads/ --suffix_name _trimmed --out_dir_report Report/figure/trimmedQC --plot_device pdf", "Error making plots")
    print "Finished plots and tables trimmed data"
    print "Finished running step 2"

def step_3(analysis_info_file):
    '''Run the third part of the analysis pipeline.

    Maps and deduplicates reads with bismark, and creates mapping 
    QC plots.

    --analysis_info_file - the file of analysis options
    '''
    functions.runAndCheck("python bin/validateAnalysisInfo.py --analysis_info_file "+analysis_info_file, "Error in analysis info")

    functions.runAndCheck("python bin/mappingReads.py --run bismark_alignment " + "--analysis_info_file " + analysis_info_file + " --sample_names_file " + args.sample_names_file, "Error mapping reads")

    print "Finished mapping"

    print "Running deduplication"
    functions.runAndCheck("python bin/mappingReads.py --run deduplicate_bismark " + "--analysis_info_file " + analysis_info_file + " --sample_names_file " + args.sample_names_file, "Error deduplicating")
    print "Finished deduplication"

    print "Making plots"

    ai=functions.read_analysis_info_file(analysis_info_file)

    path        = ai['working_directory']
    read_dir    = path + "alignedReads/"
    sample_file = path + "sample_names.txt"
    out_dir     = path + "Report/figure/mappingQC/"
    functions.make_sure_path_exists(out_dir)

    cmd_str = ('srun /usr/bin/Rscript bin/mappingQC.R %s %s  _bismark_bt2_PE_report.txt %s' 
        % (read_dir, path+"sample_names.txt", out_dir))

    functions.runAndCheck(cmd_str, "Error making mapping QC plots")
    print "Finished running step 3"


if __name__ == '__main__':
    # Script information
    parser = argparse.ArgumentParser(prog='runMethylationAnalysis.py',description = '')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)

    # Analysis information
    parser.add_argument('--analysis_info_file', 
        help='Text file with details of the analysis. Default=analysis_info.txt', 
        default='analysis_info.txt')

    # # Note - this argument is never used at the moment.
    # parser.add_argument('--sample_names_file', 
    #     help='File with sample names. Default=sample_names.txt', 
    #     default='sample_names.txt')

    # Choose a section of the pipeline to run
    parser.add_argument('--run', 
        help='Choose a section of the pipeline to run. Possible options: step1_prepare_analysis, step2_qc_and_trimming, step3_mapping_and_deduplication', 
        default='')

    # Parse arguments
    args=parser.parse_args()

    if args.run == '':
        print "No run argument was supplied!"
        print "Use --help to see valid parameters."
        sys.exit(1)

    if args.run == 'step1_prepare_analysis':
        step_1(args.analysis_info_file)

    if args.run == 'step2_qc_and_trimming':
        step_2()

    if args.run == 'step3_mapping_and_deduplication':
        step_3(args.analysis_info_file)
        
    sys.exit(0)
