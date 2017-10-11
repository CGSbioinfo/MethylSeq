#!/usr/bin/env python 

import os
import os.path
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
import validation
import invokations
import plotting
import mapping

__version__ = 'v02'

def step_1():
    '''Run the first part of the analysis pipeline.

    Checks the info files, runs bcl2fastq, and creates the sample names file
    for future scripts.
    '''
    if not validation.validateSampleSheet(args.analysis_info_file):
        logger.error("Sample sheet failed validation")
        sys.exit(1)

    invokations.invoke_bcl2fastq(args.analysis_info_file)

    logger.info("Finished bcl2fastq. Output is in the fastq/ folder")

    invokations.create_sample_names(args.analysis_info_file)

    invokations.organize_working_directory()

    logger.info("Finished running step 1")

def step_2():
    '''Run the second part of the analysis pipeline.

    Runs fastqc, trims reads, and generates the QC plots for 
    raw and trimmed data.
    '''
    logger.info("Running fastqc")
    invokations.qcReads()

    logger.info("Creating plots and tables for raw data")
    plotting.fastqc_tables_and_plots(suffix_name='_raw', out_dir_report='Report/figure/rawQC', plot_device='pdf')

    logger.info("Trimming reads")
    invokations.trim_reads()

    logger.info("Creating plots and tables for trimmed data")
    plotting.fastqc_tables_and_plots(in_dir='trimmedReads/', out_dir='trimmedReads/', suffix_name='_trimmed', out_dir_report='Report/figure/trimmedQC', plot_device='pdf')

    logger.info("Finished running step 2")

def step_3():
    '''Run the third part of the analysis pipeline.

    Maps and deduplicates reads with bismark, and creates mapping 
    QC plots.
    '''
    logger.info("Running alignment")
    mapping.map_reads(run='bismark_alignment')

    logger.info("Running deduplication")
    mapping.map_reads(run='deduplicate_bismark')

    logger.info("Making plots")

    path        = functions.getWorkingDir(args.analysis_info_file)
    read_dir    = path + "alignedReads/"
    sample_file = path + "sample_names.txt"
    out_dir     = path + "Report/figure/mappingQC/"
    functions.make_sure_path_exists(out_dir)

    cmd_str = ('/usr/bin/Rscript bin/mappingQC.R %s %s  _bismark_bt2_PE_report.txt %s' 
        % (read_dir, sample_file, out_dir))

    functions.srun(cmd_str)
    logger.info("Finished running step 3")

def step_4():
    '''Extract methylation data from the mapped reads
    '''
    logger.info("Extracting methlyation data")
    invokations.extract_methylation(args.analysis_info_file)
    invokations.create_remove_bases_file_info("remove_bases.txt")
    logger.info("Finished running step 4")

def step_5():
    '''Extract methylation data from the mapped reads using a remove bases file
    '''
    logger.info("Extracting methlyation data with read trimming")
    invokations.extract_methylation(analysis_info_file=args.analysis_info_file)
    invokations.calculate_coverage(args.analysis_info_file)
    logger.info("Finished running step 5")

def step_6():
    '''Run differential methylation analysis
    '''
    logger.info("Running differential methylation analysis")
    invokations.analyse_methylation(args.analysis_info_file)
    logger.info("Finished running step 6")


def setupLogger():
    '''Configure the logger. 

    Use a log file in the analysis working directory to store
    all messages of level DEBUG or higher. Send all messages of
    level INFO or higher to sdout.
    '''
    logger = logging.getLogger("runMethylationAnalysis")
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(functions.getWorkingDir(args.analysis_info_file)+"analysis.log")
 
    formatter = logging.Formatter('%(asctime)s\t%(name)s\t%(funcName)s\t%(levelname)s\t%(message)s')
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)
    return(logger)

def check_info_file():
    ''' Check that the analysis info file is valid.
        Quit if file is not valid
    '''
    if not validation.validateAnalysisInfo(args.analysis_info_file):
        logger.error("Analysis info file failed validation")
        sys.exit(1)

def runTests():
    ''' Run tests
    '''
    logger.info("Testing function")
    functions.testLogging() 

def noArgs():
    ''' Run when no --run argument was input
    '''
    print "No run argument was supplied!"
    print "Use --help to see valid parameters."
    sys.exit(1)  

def makeRemoveFile():
    invokations.create_remove_bases_file_info("remove_bases.txt")

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='runMethylationAnalysis.py',description = '')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s-'+__version__)

    # Analysis information
    parser.add_argument('--analysis_info_file', 
        help='Text file with details of the analysis. Default=analysis_info.txt', 
        default='analysis_info.txt')

    # Choose a section of the pipeline to run
    parser.add_argument('--run', 
        help='Choose a section of the pipeline to run. Possible options: step0_create_info, step1_prepare_analysis, step2_qc_and_trimming, step3_mapping_and_deduplication, step4_extract_methylation, step5_extract_methylation, step6_analyse_methylation, create_remove_file', 
        default='')

    args=parser.parse_args()

    if(args.run=='step0_create_info'):
        invokations.create_analysis_info_file(args.analysis_info_file)
        print "Analysis file created"
        print "Fill in the file before running step 1"
        sys.exit(0)


    if(not os.path.isfile(args.analysis_info_file)):
        print "Analysis file '%s' not found" % args.analysis_info_file
        print "Run step 0 and try again"
        sys.exit(1)

    logger = setupLogger()
    check_info_file()

    # Run the command input
    argDict = { 'step1_prepare_analysis': step_1,
                'step2_qc_and_trimming': step_2,
                'step3_mapping_and_deduplication': step_3,
                'step4_extract_methylation': step_4,
                'step5_extract_methylation': step_5,
                'step6_analyse_methylation': step_6,
                'test': runTests,
                '': noArgs,
                'create_remove_file': makeRemoveFile }

    for i in argDict:
        if(args.run==i):
            logger.info("Running command '"+i+"'")
            argDict[i]()

    sys.exit(0)
