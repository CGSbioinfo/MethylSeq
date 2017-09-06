''' Plotting methods for methylSeq analysis.

    This module contains methods used in sequence
    analysis pipelines for making QC plots
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
import argparse

def _tables(i):
    logger = logging.getLogger("runMethylationAnalysis.plotting")
    logger.info("Making tables for "+i)
    outdir = re.sub('fastqc_data.txt', '', i)

    fastq_data = i
    plot = 'all' # all, per_base_sequence_content, per_base_sequence_quality, per_sequence_gc_content, per_sequence_quality_scores, seq_dup_levels, kmer_content

    dictionary = {'seq_length':'>>Sequence Length Distribution','per_base_sequence_content':'>>Per base sequence content', 'per_base_sequence_quality':'>>Per base sequence quality', 'kmer_content':'>>Kmer Content', 'seq_dup_levels':'>>Sequence Duplication Levels', 'per_sequence_quality_scores':'>>Per sequence quality scores', 'per_sequence_gc_content': '>>Per sequence GC content'}

    if plot == 'all':
        for d in dictionary.keys():
            f = open(fastq_data,'r')
            i=0
            for line in f:
                line = line.strip()
                if line.startswith(dictionary[d]):
                    name = d + '.txt'
                    line_number = i
                i += 1
            f.close()
            output = open(outdir + name,'w')
            f = open(fastq_data,'r')
            i=0
            for line in f:
                line = line.strip()
                if i > line_number:
                    if line != '>>END_MODULE':
                        output.write(line + '\n')
                    else:
                        break
                i += 1
            f.close()
    else:
        f = open(fastq_data,'r')
        i=0
        for line in f:
            line = line.strip()
            if line.startswith(dictionary[plot]):
                name = plot + '.txt'
                line_number = i
            i += 1
        f.close()
        output = open(outdir + name,'w')
        f = open(fastq_data,'r')
        i=0
        for line in f:
            line = line.strip()
            if i > line_number:
                if line != '>>END_MODULE':
                    output.write(line + '\n')
                else:
                    break
            i += 1
    # functions.runAndCheck("srun python bin/create_fastqcTables.py " + i + " all " + outdir, "Error making fastqc tables")

def _plots(i, in_dir, readType, out_dir_report, suffix_name, plot_device):
    logger = logging.getLogger("runMethylationAnalysis.plotting")
    logger.info("Making plots for "+i)
    functions.srun('/usr/bin/Rscript bin/create_fastqcPlots_perSample.R ' + in_dir + ' ' + i + ' ' + readType + ' ' + out_dir_report + ' ' + suffix_name + ' ' + plot_device)

def fastqc_tables_and_plots(analysis_info_file='analysis_info.txt', 
    in_dir='rawReads/', 
    out_dir='rawReads/', 
    out_dir_report='Report/figure', 
    suffix_name= '_plot',
    sample_names_file='sample_names.txt',
    plot_device = 'png'):
    
    logger = logging.getLogger("runMethylationAnalysis.plotting")

    # Get wd
    path=os.getcwd()

    # Read sample names text file
    ai=functions.read_analysis_info_file(analysis_info_file)
        #Ncores
    ncores=int(ai['ncores'])

    sampleNames = functions.read_sample_names(sample_names_file)

    ninstances=int(ai['ninstances'])

    # Set input and output directories
    path    = functions.slash_terminate(path)
    in_dir  = functions.slash_terminate(in_dir)
    out_dir_report = functions.slash_terminate(out_dir_report)

    in_dir = path + in_dir
    out_dir_report= path + out_dir_report
    readType=ai['readType']

    # Create tables
    files=functions.get_filepaths(in_dir)
    files = [files[y] for y, x in enumerate(files) if re.findall("fastqc_data.txt", x)] 
    Parallel(n_jobs=ninstances)(delayed(_tables)(i) for i in files)


    # Create plots
    functions.make_sure_path_exists(out_dir_report)
    Parallel(n_jobs=ninstances)(delayed(_plots)(i, in_dir, readType, out_dir_report, suffix_name, plot_device) for i in sampleNames)

    logger.info("Making plots for all samples")
    functions.srun('/usr/bin/Rscript bin/create_fastqcPlots_allSamples.R ' + in_dir + ' ' + sample_names_file + ' ' + readType + ' ' + out_dir_report + ' ' + suffix_name + ' ' + plot_device)

