#!/usr/bin/env python 

import argparse
from datetime import datetime
import logging
import os.path
import sys
import functions
import invokations
import mapping
import pipeline 
import plotting
import validation
import socket


__version__ = 'v02'

def step_1():
    '''Run the first part of the analysis pipeline.

    Checks the info files, runs bcl2fastq, and creates the sample names file
    for future scripts.
    '''
    _logger = logging.getLogger(__name__)
    if not validation.validateSampleSheet(args.analysis_info_file):
        _logger.error("Sample sheet failed validation")
        sys.exit(1)

    invokations.invoke_bcl2fastq(args.analysis_info_file)
    
    _logger.info("Finished bcl2fastq. Output is in the fastq/ folder")

    invokations.create_sample_names(args.analysis_info_file)

    invokations.organize_working_directory()

    _logger.info("Finished running step 1")

def step_2():
    '''Run the second part of the analysis pipeline.

    Runs fastqc, trims reads, and generates the QC plots for 
    raw and trimmed data.
    '''
    _logger = logging.getLogger(__name__)
    _logger.info("Running fastqc")
    invokations.qcReads()

    _logger.info("Creating plots and tables for raw data")
    plotting.fastqc_tables_and_plots(suffix_name='_raw', out_dir_report='Report/figure/rawQC', plot_device='pdf')

    _logger.info("Trimming reads")
    invokations.trim_reads()

    _logger.info("Creating plots and tables for trimmed data")
    plotting.fastqc_tables_and_plots(in_dir='trimmedReads/', out_dir='trimmedReads/', suffix_name='_trimmed', out_dir_report='Report/figure/trimmedQC', plot_device='pdf')

    _logger.info("Finished running step 2")

def step_3():
    '''Run the third part of the analysis pipeline.

    Maps and deduplicates reads with bismark, and creates mapping 
    QC plots.
    '''
    _logger = logging.getLogger(__name__)
    _logger.info("Running alignment")
    mapping.map_reads(run='bismark_alignment')

    _logger.info("Running deduplication")
    mapping.map_reads(run='deduplicate_bismark')

    _logger.info("Making plots")

    path        = functions.getWorkingDir(args.analysis_info_file)
    read_dir    = path + "alignedReads/"
    sample_file = path + "sample_names.txt"
    out_dir     = path + "Report/figure/mappingQC/"
    functions.make_sure_path_exists(out_dir)

    cmd_str = ('/usr/bin/Rscript src/mappingQC.R %s %s  _bismark_bt2_PE_report.txt %s' 
        % (read_dir, sample_file, out_dir))

    pipeline.srun(cmd_str, ncores=1, mem=20)
    _logger.info("Finished running step 3")

def step_4():
    '''Extract methylation data from the mapped reads
    '''
    _logger = logging.getLogger(__name__)
    _logger.info("Extracting methlyation data")
    invokations.extract_methylation(args.analysis_info_file)
    invokations.create_remove_bases_file_info("remove_bases.txt")
    _logger.info("Finished running step 4")

def step_5():
    '''Extract methylation data from the mapped reads using a remove bases file
    '''
    _logger = logging.getLogger(__name__)
    _logger.info("Extracting methlyation data with read trimming")
    invokations.extract_methylation(analysis_info_file=args.analysis_info_file)
    invokations.calculate_coverage(args.analysis_info_file)
    _logger.info("Finished running step 5")

def step_6():
    '''Run differential methylation analysis via BiSeq
    '''
    _logger = logging.getLogger(__name__)
    _logger.info("Running differential methylation analysis using BiSeq")
    invokations.analyse_methylation(args.analysis_info_file)
    _logger.info("Finished running step 6")

def step_7():
    '''Run differential methylation analysis via SVA
    '''
    _logger = logging.getLogger(__name__)
    _logger.info("Running differential methylation analysis using SVA")
    invokations.run_sva(args.analysis_info_file)
    _logger.info("Finished running step 7")


def setup_logger():
    '''Configure the logger. 
    Use a log file to store all messages of level DEBUG 
    or higher. Send all messages of level INFO or higher
    to sdout.
    '''
    _logger = logging.getLogger()
    pipeline.make_sure_path_exists("logs/")
    _logger.setLevel(logging.DEBUG) 
    logfile = "logs/"+datetime.now().strftime('analysis.%Y-%m-%d_%H-%M.log')
    formatter = logging.Formatter('%(asctime)s\t%(name)s\t%(funcName)s\t%(levelname)s\t%(message)s')

    fh = logging.FileHandler(logfile)
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG) # Log ALL THE THINGS to file
    _logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO) # Log important things to console
    _logger.addHandler(sh)       

def check_info_file():
    ''' Check that the analysis info file is valid.
        Quit if file is not valid
    '''
    if not validation.validateAnalysisInfo(args.analysis_info_file):
        _logger.error("Analysis info file failed validation")
        sys.exit(1)
    _logger.debug("Analysis info file passed validation")

def runTests():
    ''' Run tests
    '''
    _logger = logging.getLogger(__name__)
    _logger.info("Testing pipeline module")
    pipeline.testLogging() 
    plotting.test()

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
        help='Choose a section of the pipeline to run. Possible options: step0_create_info, step1_prepare_analysis, step2_qc_and_trimming, step3_mapping_and_deduplication, step4_extract_methylation, step5_extract_methylation, step6_analyse_methylation, step7_SVA_analysis, create_remove_file', 
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

    setup_logger()
    _logger = logging.getLogger(__name__)
    _logger.debug(("Started script on %s" % socket.gethostname()))
    check_info_file()

    # Run the command input
    argDict = { 'step1_prepare_analysis': step_1,
                'step2_qc_and_trimming': step_2,
                'step3_mapping_and_deduplication': step_3,
                'step4_extract_methylation': step_4,
                'step5_extract_methylation': step_5,
                'step6_analyse_methylation': step_6,
                'step7_SVA_analysis': step_7,
                'test': runTests,
                '': noArgs,
                'create_remove_file': makeRemoveFile }

    if(argDict[args.run]):
        _logger.info(("Running command '%s' " % args.run))
        argDict[args.run]()
    sys.exit(0)
