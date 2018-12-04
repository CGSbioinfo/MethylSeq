''' System commands to invoke for methylSeq analysis.

    This module contains methods used to invoke tools for
    QC and processing of bisulphite converted sequence reads 
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


def create_analysis_info_file(analysis_info_file="analysis_info.txt"):
    lines=['working_directory = ', 'run_folder = ', 'run_samplesheet = ', 'bcl2fastq_output = fastq/ ', 
    'readType = pairedEnd', 'reference_genome = ', 
    'bismark_params = --bowtie2; --bam; -N 0; -L 20; -D 15; -R 2; --score_min L,0,-0.2; ', 
    'methyl_extract_params= --bedGraph; --gzip; --merge_non_CpG;', 'target_regions_bed = ', 'ncores = 10', 
    'ninstances = 3', 'clean_files = False', 'sill = 1', 'dmr_bandwidth = 50', 'gtf_file = ', 'cx_extract = False']
    outfile=open(analysis_info_file,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()

def create_remove_bases_file_info(outfile_name, sample_names_file='./sample_names.txt'):
    """ Creates a file of bases to remove.

        Creates a file for base trimming in a methylSeq project. This is tab delimited
        and provides for a number of bases to remove at 3' and 5' ends of each read.

        --outfile_name - the name of the output file.
    """ 
    lines=['sample\t5R1\t3R1\t5R2\t3R2\n']
    sampleNames = functions.read_sample_names(sample_names_file)
    for name in sampleNames:
        lines.append(name+'\t\n')
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l)
    outfile.close()

def invoke_bcl2fastq(analysis_info_file):

    '''This function runs bcl2fastq.

    It uses the parameters in the analysis_info file.
    If bcl2fastq fails, the script will quit with exit code 1, otherwise
    it will quit with exit code 0
    '''

    ai=functions.read_analysis_info_file(analysis_info_file)

    cmdStr = "bcl2fastq -R " + ai['run_folder'] + " -o " + ai['bcl2fastq_output'] + " --no-lane-splitting -l NONE --sample-sheet " + ai['run_samplesheet']

    pipeline.srun(cmdStr, ncores=ai['ncores'], mem=20)


def create_sample_names(analysis_info_file, in_dir='bcl2fastq_output', out_file='./sample_names.txt'):
    """ 
    This function creates a file which needs to be filled with information 
    required for a methylSeq project.

    It takes one argument, the 'outfile', which is the name of the output file. 
    The default is 'analysis_info.txt'
    """ 
    if in_dir == 'bcl2fastq_output': 
        ai=functions.read_analysis_info_file(analysis_info_file)
        allFiles=functions.get_filepaths(ai[in_dir])
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[fastq[y] for y,x in enumerate(fastq) if not re.findall('Undetermined', x)]
    elif in_dir != 'bcl2fastq_output':
        allFiles=os.listdir(in_dir)
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
    
    fastq=[re.sub('.*/','', f) for f in fastq]
    fastq=[re.sub('_R[0-9_].*','',f) for f in fastq]
    fastq=list(set(fastq))

    outfile=open(out_file, 'w')
    for f in fastq:
        outfile.write(f.split('/')[-1] + '\n')
    outfile.close()       

def _moveReads(sampleNames, fastq):
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

def organize_working_directory(analysis_info_file='analysis_info.txt', sample_names_file='sample_names.txt', in_dir='bcl2fastq_output'):

    # Read analysis info file
    ai=functions.read_analysis_info_file(analysis_info_file)

    # Read sample names
    sampleNames = functions.read_sample_names(sample_names_file)
        
    # Check if rawReads exists
    folders = os.listdir('.')
    readsFiles = [folders[i] for i, x in enumerate(folders) if re.findall('rawReads',x)]
    
     # Collect fastq files analysis_info_file
    if in_dir == 'bcl2fastq_output':
        allFiles=functions.get_filepaths(ai[in_dir])
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[fastq[y] for y,x in enumerate(fastq) if not re.findall('Undetermined', x)]
    elif in_dir != 'bcl2fastq_output':
        allFiles=os.listdir(in_dir)
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[in_dir + x for x in fastq]

    # Move reads
    if not readsFiles:
        # The rawReads directory does not exist.
        # Create the directory and move the files.
        print "rawReads/ folder does not exist. Creating folder and moving files."
        _moveReads(sampleNames, fastq)
    else:
        # If the rawReads directory exists and is empty
        # - for example if step 1 failed before bcl2fastq 
        # completed, then we need to continue moving files. 
        # Otherwise fail with meaningful error.
        if os.listdir('rawReads') == []: 
            print "rawReads/ folder exists and is empty. Moving files."
            _moveReads(sampleNames, fastq)
        else:
            print "rawReads/ folder exists and contains files. Not moving files."


def _qc_check(i, in_dir, out_dir, gz):
    allFiles = os.listdir(in_dir + i )
    pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_R2", x)]
    functions.make_sure_path_exists(out_dir+i)

    cmd_str = "fastqc "  + in_dir + i + "/" + i + "*_R1*.fastq" + gz + " --outdir=" + out_dir + i + " --nogroup --extract "
    pipeline.srun(cmd_str)
    if pairedReads_temp:
        cmd_str = "fastqc " + in_dir + i + "/" + i + "*_R2*.fastq" + gz + " --outdir=" + out_dir + i + " --nogroup --extract"
        pipeline.srun(cmd_str)

def qcReads(analysis_info_file='analysis_info.txt', 
    in_dir='rawReads/', 
    out_dir='rawReads/', 
    out_dir_report='Report/figure/data/', 
    sample_names_file='sample_names.txt' ):
    """ This method will run fastqc on all samples in parallel. 
     - You can specify the number of cores used
     - The input and output default is rawReads/samplesubfolder/
     - It also returns a table and a plot with the number of reads for each sample, the output by default is Report/figure/data 
     """
    _logger = logging.getLogger(__name__)
    
    # Set path of working directory
    ai=functions.read_analysis_info_file(analysis_info_file)
    path=os.getcwd()

    #Ncores
    ncores=int(ai['ncores'])
    ninstances=int(ai['ninstances'])

    # Read sample names text file
    sampleNames = functions.read_sample_names(sample_names_file)


    # Set input and output directories
    path    = functions.slash_terminate(path)
    in_dir  = functions.slash_terminate(in_dir)
    out_dir = functions.slash_terminate(out_dir)
    out_dir_report = functions.slash_terminate(out_dir_report)

    in_dir = path + in_dir
    out_dir= path + out_dir
    out_dir_report = path + out_dir_report

    # Create out_dir_report
    functions.make_sure_path_exists(out_dir_report)

    # Detect if files are gz 
    gz = functions.check_gz(in_dir)

    # Run fastqc
    # print "Running fastqc"
    _logger.info("Running fastqc")

    pool = multiprocessing.Pool(processes=ninstances)

    _logger.info("Invoking apply_async with %d processes" % ninstances)

    for i in sampleNames:
        pool.apply_async(_qc_check, [i, in_dir, out_dir, gz])

    pool.close()
    pool.join()
    
    # Parallel(n_jobs=ninstances)(delayed(_qc_check)(i, in_dir, out_dir, gz) for i in sampleNames)

    # Number of reads per sample
    _logger.info("Running index QC")
    pipeline.runAndCheck("/usr/bin/Rscript src/indexQC.R " + in_dir + " " + out_dir_report, "Error in index QC") 


def _trim(i, in_dir, out_dir, gz, ncores):
    '''
    Trim sequencing reads using trim_galore
    '''
    allFiles = os.listdir(in_dir + i )
    pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_R2", x)]

    cmdStr ="trim_galore --gzip"
    if pairedReads_temp:
        cmdStr = cmdStr + " --paired --fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + in_dir + i + "/" + i + "*_R1*.fastq" + gz + " " + in_dir + i + "/" + i + "*_R2*.fastq" + gz
        pipeline.srun(cmd=cmdStr, ncores=1, mem=20)
    else:
        pipeline.srun(cmdStr + "--fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + in_dir + i + "/" + i + "*_R1*.fastq" + gz, ncores=1, mem=20)


def trim_reads(analysis_info_file='analysis_info.txt', 
    in_dir='rawReads/', 
    out_dir='trimmedReads/', 
    out_dir_report='Report/figure/data/', 
    sample_names_file='sample_names.txt'):
    
    _logger = logging.getLogger(__name__)
    
    path=os.getcwd()
    os.chdir(path)
    ai=functions.read_analysis_info_file(analysis_info_file)    

    #Ncores
    ncores=int(ai['ncores'])
    ninstances=int(ai['ninstances'])

    # Read sample names text file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories'
    path    = functions.slash_terminate(path)
    in_dir  = functions.slash_terminate(in_dir)
    out_dir = functions.slash_terminate(out_dir)
    out_dir_report = functions.slash_terminate(out_dir_report)

    in_dir = path + in_dir
    out_dir = path + out_dir
    out_dir_report=path + out_dir_report

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Run trimm galore
    functions.make_sure_path_exists(out_dir)
    _logger.info("Trimming reads")

    pool = multiprocessing.Pool(processes=ninstances)

    _logger.info("Invoking apply_async with %d processes" % ninstances)

    for i in sampleNames:
        pool.apply_async(_trim, [i, in_dir, out_dir, gz, ncores])

    pool.close()
    pool.join()

    functions.make_sure_path_exists(out_dir_report)
    
    # Generate summary stats
    cmd_str = ("/usr/bin/Rscript src/trimming_summary.R %s %s %s " %
        (in_dir, out_dir, out_dir_report))
    pipeline.srun(cmd_str)


def _methylationExtraction(i, ai, in_dir, out_dir, remove_bases_dict, isDedup, cx):
    '''
    Run the bismark methylation extraction on the given sample.in_dir

    If a non-empty dictionary is supplied in remove_bases_dict, the corresponding
    bases will be trimmed
    '''
    _logger = logging.getLogger(__name__)

    params=ai['methyl_extract_params'].replace(';','')
    alignedReads = os.listdir(in_dir)
    alignedReads = [alignedReads[y] for y, x in enumerate(alignedReads) if re.findall(i, x)]
    if not alignedReads:
        _logger.warning("No aligned reads found for %s." % i)
        return

    # Filter to filenames containing the sample name
    dedupregex = re.compile("bismark_bt2_.*deduplicated.bam$")
    bamregex   = re.compile("val_1.fq.gz_bismark_bt2_pe.bam$")

    dedup_files = filter(dedupregex.search, alignedReads)
    bam_files   = filter(bamregex.search, alignedReads)

    # Sanity check: if we requested deduplication, make sure the files exist
    if isDedup == True:
        _logger.debug("Attempting to extract from deduplicated file")
        if not dedup_files:
            _logger.warning("No deduplicated bam files found for %s" % i)
            return;
        else: 
            bamFile = dedup_files[0]
    else:
        _logger.debug("Not attempting to extract from deduplicated file")
        if not bam_files:
            _logger.warning("No bam files found for %s" % i)
            return;
        else:
            bamFile = bam_files[0]
    
    if not bamFile:
        _logger.warning("No bam file to process for %s." % i)
        return

    bamFile = in_dir + bamFile

    # Create the ignore string if needed
    if(remove_bases_dict):
        _logger.info("Clipping bases for %s" % i)
        basesInfo=remove_bases_dict[i]
        logFile = ("%smethylExtract_REMOVEDBASES_log_%s.txt"
            % (out_dir, i))
        ignoreString = (" --ignore %s --ignore_3prime %s --ignore_r2 %s --ignore_3prime_r2 %s "
            % (basesInfo[0], basesInfo[1], basesInfo[2], basesInfo[3]))
    else:
        logFile = ("%smethylExtract_log_%s.txt"
            % (out_dir, i))
        ignoreString = ""


    cx_string = "--cytosine_report --CX_context --genome_folder "+ai['reference_genome']+" " if cx else ""

    param_string = ('--multicore %s %s %s --output %s %s' %
        ( ai['ncores'], params, cx_string, out_dir, ignoreString ))

    cmdStr = ('bismark_methylation_extractor %s %s &>%s' % (param_string, bamFile, logFile ))
    _logger.info("Extracting methylation from "+i)
    pipeline.srun(cmd=cmdStr, ncores=ai['ncores'], mem=20)
    _logger.info("Extracted methylation from "+i)


def _createRemoveBaseDict(remove_bases_file):
    '''
    Read the given file and make a dictionary 
    of bases to remove.

    If the remove bases file does not exist, an
    empty dictionary will be returned. 
    '''
    remove_bases_dict = {}
    if os.path.exists(remove_bases_file):
        with open(remove_bases_file, 'r') as f:
            next(f)
            for line in f:
                (key, ignore, ignore_3prime, ignore_r2, ignore_3prime_r2)=line.split()
                remove_bases_dict[key]=(ignore, ignore_3prime, ignore_r2, ignore_3prime_r2)
    return remove_bases_dict


def _makeMbiasPlots(remove_bases_dict, out_dir, sample_names_file):
    '''Call the R script to generate methylation bias plots.
    '''
    _logger = logging.getLogger(__name__)
    clipped = "_clipped" if remove_bases_dict else ""
    out_file = "Report/figure/methExtractQC/Mbias_plot%s.pdf" % clipped

    cmdStr = ('/usr/bin/Rscript src/methylExtractQC_mbias_plot.R %s %s .M-bias.txt %s'
            % (out_dir, sample_names_file, out_file))
    _logger.info("Making M-bias plots")
    pipeline.srun(cmd=cmdStr, mem=5)
       
def extract_methylation(analysis_info_file='analysis_info.txt',
    in_dir='alignedReads/', 
    out_dir='methylExtraction/', 
    sample_names_file='sample_names.txt',
    dedup=True,
    remove_bases_file='remove_bases.txt'):
    '''
    Run the bismark methylation extraction
    '''
    _logger = logging.getLogger(__name__)
    ai=functions.read_analysis_info_file(analysis_info_file)

    cx = ai['cx_extract']=='True'

    path = functions.getWorkingDir(analysis_info_file)

    ninstances=int(ai['ninstances'])


    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories
    path    = functions.slash_terminate(path)
    in_dir  = functions.slash_terminate(in_dir)
    out_dir = functions.slash_terminate(out_dir)
    
    in_dir = path + in_dir

    if not os.path.exists(in_dir):
        _logger.error("Input directory was not found: " + in_dir)
        sys.exit(1)

    out_dir = path + out_dir

    functions.make_sure_path_exists(out_dir)

    # Read remove bases file if present or create an empty dict otherwise
    remove_bases_dict = _createRemoveBaseDict(remove_bases_file)
    pool = multiprocessing.Pool(processes=ninstances)

    _logger.info("Invoking apply_async with %d processes" % ninstances)

    for i in sampleNames:
        pool.apply_async(_methylationExtraction, [i, ai, in_dir, out_dir,remove_bases_dict, dedup, cx])

    pool.close()
    pool.join()

    # Create the QC plots
    _makeMbiasPlots(remove_bases_dict, out_dir, sample_names_file)
    return


def calculate_coverage(analysis_info_file='analysis_info.txt',
    cov_dir='methylExtraction/',
    sample_names_file='sample_groups.txt'):
    '''
    Calculate the sequence coverage of the genome and target region
    '''
    _logger = logging.getLogger(__name__)

    if not os.path.exists(cov_dir):
        _logger.error("Input directory was not found: " + cov_dir)
        sys.exit(1)

    ai=functions.read_analysis_info_file(analysis_info_file)

    wg_output = 'Report/figure/methExtractQC/coverage_rawData_wg.pdf'
    tr_output = 'Report/figure/methExtractQC/coverage_rawData_tr.pdf'

    cmdStr = ('/usr/bin/Rscript src/coverage_methExtractedData.r %s %s %s %s %s '
        %(sample_names_file, cov_dir, wg_output, tr_output, ai['target_regions_bed']))

    logging.info("Calculating coverage in genome and target regions")
    pipeline.srun(cmdStr)



def analyse_methylation(analysis_info_file='analysis_info.txt',
    cov_dir='methylExtraction/', 
    sample_group_file='sample_groups.txt',
    temp_folder_name='biseqAnalysis/',
    chunk_size=20000):
    '''
    Run differential methylation analysis in BiSeq.
    '''
    _logger = logging.getLogger(__name__)
    if not os.path.exists(cov_dir):
        logging.error("Input directory was not found: " + cov_dir)
        sys.exit(1)
        
    ai=functions.read_analysis_info_file(analysis_info_file)

    temp_folder_name = functions.slash_terminate(temp_folder_name)

    cmdStr = ('/usr/bin/Rscript src/analyseMethylationPatterns.r %s %s %s %s %s %s %s %s %s %s' 
        % (sample_group_file, cov_dir, ai['target_regions_bed'], temp_folder_name, ai['ncores'],  chunk_size, ai['gtf_file'], ai['dmr_bandwidth'], ai['sill'], ai['tissue']))

    _logger.info("Running methylation analysis")
    pipeline.srun(cmdStr, ncores=ai['ncores'], mem=120)






