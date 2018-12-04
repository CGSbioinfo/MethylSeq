''' Read mapping.

    This module contains methods used to invoke Bismark to map
    and deduplicate bisulphite converted sequence reads 
'''

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
import pipeline

def _alignment(sampleName, ai, in_dir, out_dir, gz):
    '''Perform alignment of reads using Bismark params

    This function check that read 1 and read 2 files exist
    and runs bismark using the parameters detailed in the 
    analysis_info file.

    --sampleName - the name of the sample to process
    --ai - the settings
    --in_dir - the input directory
    --out_dir - the output directory
    --gz - the gzip suffix to use when searching for files (empty string if files are not compressed)
    '''
    logger = logging.getLogger(__name__)
    bismark_params=ai['bismark_params'].replace(';','')
    
    # Get all files in the input dir as a list
    trimmedReads = os.listdir(in_dir)

    # Filter to filenames containing the sample name
    regex = re.compile(sampleName)
    trimmedReads = filter(regex.search, trimmedReads)

    r1regex = re.compile("_R1_.*val_1.fq"+gz)
    r2regex = re.compile("_R2_.*val_2.fq"+gz)

    r1 = filter(r1regex.search, trimmedReads)
    r2 = filter(r2regex.search, trimmedReads)

    # Ensure read 1 is present
    if not r1:
        logger.info("No read 1 file found. Ignoring %s." % sampleName)
        return

    r1 = in_dir + r1[0]

    # Check read 2 is present when needed
    if(ai['readType']=='pairedEnd'):
        if not r2:
            logger.warning("No read 2 file found when readType is paired end. Ignoring %s." % sampleName)
            return

    if r2: 
        r2 = in_dir + r2[0]
        cmdStr = ('bismark %s --output_dir %s %s -1 %s -2 %s &>%sbismark_log_%s.txt' 
            % ( bismark_params, out_dir, ai['reference_genome'], r1, r2, out_dir, sampleName ))
        logger.info("Aligning %s" % sampleName)
        pipeline.srun(cmd=cmdStr, ncores=ai['ncores'], mem=40)
    else:
        logger.info("No read 2 for %s" % sampleName)
        logger.warning("Single reads not currently implemented")
    logger.info("Alignment complete for %s" % sampleName)

def _deduplicate(sampleName, ai, in_dir, out_dir):
    '''Perform deduplication of reads using Bismark params.

        The function checks if the required bam file exists before running.
        --sampleName - the name of the sample to process
        --ai - the analysis options
        --in_dir - the input directory
        --out_dir - the output directory
    '''
    logger = logging.getLogger(__name__)
    alignedReads = os.listdir(out_dir)

     # Filter to filenames containing the sample name
    regex = re.compile(sampleName)
    alignedReads = filter(regex.search, alignedReads)

    if not alignedReads:
        logger.warning("No aligned reads found for %s." % sampleName)
        return

    regex = re.compile("bismark_bt2_.*.bam$")

    bamFiles = filter(regex.search, alignedReads)
    if not bamFiles:
        logger.warning("No bam file found for %s." % sampleName)
        return

    bamFile = out_dir + bamFiles[0]

    cmdStr = ('deduplicate_bismark --bam %s &>%sdedup_log_%s.txt' 
            % (bamFile, out_dir, sampleName ))

    alignedReads = [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall(sampleName, x)]
    
    logger.info("Deduplicating %s" % sampleName)
    pipeline.srun(cmd=cmdStr, ncores=ai['ncores'], mem=40)
    logger.info("Deduplication complete for %s" % sampleName)

def map_reads(analysis_info_file='analysis_info.txt', 
    in_dir='trimmedReads/', 
    out_dir='alignedReads/', 
    out_dir_report='Report/figure/data/', 
    sample_names_file='sample_names.txt',
    run='all'):
    '''Run mapping of reads using Bismark

        --analysis_info_file - the info data
        --in_dir - the input directory (default trimmedReads/)
        --out_dir - the output directory (default alignedReads/)
        --out_dir_report - the directory for reports (default Report/figure/data/)
        --sample_names_file - the sample names file (default sample_names.txt)
        --run - the command to run (options: all, bismark_alignment, deduplicate_bismark; default: all)
    '''

    logger = logging.getLogger(__name__)

    # Set path of working directory
    ai   = functions.read_analysis_info_file(analysis_info_file)
    path = functions.getWorkingDir(analysis_info_file)

    #Ncores
    ncores = int(ai['ncores'])
    ninstances=int(ai['ninstances'])

    # Read sample names text file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories
    path    = functions.slash_terminate(path)
    in_dir  = functions.slash_terminate(in_dir)
    out_dir = functions.slash_terminate(out_dir)

    in_dir = path + in_dir
    out_dir= path + out_dir

    # Detect if files are gz
    gz = functions.check_gz(in_dir, "fq")
    logger.debug("Appending gz string '%s'" % gz)

    # Run bismark alignment
    if run == 'all' or run =='bismark_alignment': 
        pool = multiprocessing.Pool(processes=ninstances)
        logger.info("Running bismark alignment")
        functions.make_sure_path_exists(out_dir)

        logger.info("Invoking apply_async with %d processes" % ninstances)
        for i in sampleNames:
            pool.apply_async(_alignment, [i, ai, in_dir, out_dir, gz])
        pool.close()
        pool.join()
        # Use Parallel to batch the task into groups of ninstances
        # alignments so the cluster is not overloaded.
        # Parallel(n_jobs=ninstances)(delayed(_alignment)(i, ai, in_dir, out_dir, gz) for i in sampleNames)
    
    # Run deduplication
    if run == 'all' or run == 'deduplicate_bismark':
        pool = multiprocessing.Pool(processes=ninstances)
        logger.info("Running deduplication")
        logger.info("Invoking apply_async with %d processes" % ninstances)
        for i in sampleNames:
            pool.apply_async(_deduplicate, [i, ai, in_dir, out_dir])
        pool.close()
        pool.join()   
