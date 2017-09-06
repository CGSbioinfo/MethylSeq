''' Generic methods for methylSeq analysis.

    This module contains general methods used in sequence
    analysis pipelines
'''
import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
import subprocess
from joblib import Parallel, delayed
import multiprocessing
# version v01

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raisec

def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths

def read_analysis_info_file(info_file):
    
    """ Reads the analysis_info_file and creates variables for each element """

    analysis_info_dict = {}
    with open(info_file, 'r') as f:
        for line in f:
            entry = line.strip().split("=")
            analysis_info_dict[entry[0].strip()]= entry[1].strip()
    return(analysis_info_dict)

def create_rawReads_folder(sampleNames, reads_dir):
    """ Creates a rawReads/ folder and subfolders for each sample. Moves the fastq reads from reads_dir/ to their corresponding sample subfolder"""

    # Check if rawReads exists
    folders = os.listdir('.')
    readsFiles = [folders[i] for i, x in enumerate(folders) if re.findall('rawReads',x)]
    
    # Get reads from reads directory
    allFiles = os.listdir(reads_dir)
    allFiles= [allFiles[i] for i, x in enumerate(allFiles) if re.findall("_R\d.*.fastq", x)]

    # Move reads
    if not readsFiles:
        make_sure_path_exists('rawReads')
        sampleDir = []
        for sample in sampleNames:
            reads = [allFiles[i] for i,x in enumerate(allFiles) if re.findall(sample,x)]
            if sample not in sampleDir:
                make_sure_path_exists('rawReads/'+sample)
            for r in reads:
               r = reads_dir + '/' + r
               os.system('mv ' + '"' + r + '"' + ' rawReads/' + sample)
            sampleDir.append(sample)

def read_sample_names(sample_names_file):
    """ Read the sample_names_file and returns a list of sample names"""
    sampleNames = []
    sample_names_file = open(sample_names_file,'r')
    for line in sample_names_file:
        sampleNames.append(line.strip())
    return(sampleNames)

def check_gz(dir, suffix="fastq"):
    """ Determines if files with the given extension are gzipped and returns a gz string which 
    will be added as a suffix when file names are called in further steps
    """
    readFiles = []
    for root, dir, files in os.walk(dir):
        readFiles.extend(files)
    indicesgzFiles = [i for i, x in enumerate(readFiles) if re.findall("."+suffix+".gz", x)]
    if indicesgzFiles:
        gz =".gz"
    else:
        gz =""
    return(gz)

def get_strand(strand):
    """ Reads the strand information and returns the specific terms used in picard and htseq """

    if strand == 'NONE' or strand == 'FIRST_READ_TRANSCRIPTION_STRAND' or strand == 'SECOND_READ_TRANSCRIPTION_STRAND':
        strand_picard = strand
        if strand_picard == 'NONE':
            strand_htseq = 'no'
        elif strand_picard == 'FIRST_READ_TRANSCRIPTION_STRAND':
            strand_htseq = 'yes'
        elif strand_picard == 'SECOND_READ_TRANSCRIPTION_STRAND':
            strand_htseq = 'reverse'
    elif strand == 'no' or strand == 'yes' or strand == 'reverse':
        strand_htseq = strand
        if strand_htseq == 'no':
            strand_picard = 'NONE'
        elif strand_htseq == 'yes':
            strand_picard = 'FIRST_READ_TRANSCRIPTION_STRAND'
        elif strand_htseq == 'reverse':
            strand_picard = 'SECOND_READ_TRANSCRIPTION_STRAND'
    return([strand_picard, strand_htseq])

def slash_terminate(s):
    '''Add a trailing '/' to the string if not present
    '''
    return(s if s.endsWith('/') else s + '/')

def runAndCheck(cmd, msg):
    '''Run the given command and check the return value.

        The script will print the given error message and quit 
        with exit code 1 if the command did not return
        exit code 0

        --cmd - the command to be run
        --msg - the message to print if the command returns an error
    '''
    logger = logging.getLogger("runMethylationAnalysis.functions")
    logger.debug("Invoking: " + cmd)
    try:
        retcode = subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
        if retcode != 0:
            logger.debug("System call returned "+str(retcode))
            logger.error(msg+". Quitting.")
            sys.exit(1)
    except OSError as e:
        logger.exception("Error invoking command")
        print >>sys.stderr, "Execution exception:", e
        sys.exit(1)

def srun(cmd, ncores="1", mem=""):
    '''Run the given command through srun with the given number of
    cores and memory in Gb. If srun is not installed, invokes the 
    command directly.
    '''
    logger = logging.getLogger("runMethylationAnalysis.functions")
    if(srun_is_installed()):
        s = "srun -c "+str(ncores)+" "
        s = s + " --mem="+str(mem)+"G " if mem else s
        cmd = s+cmd
    else:
        logger.debug("srun not found")

    runAndCheck(cmd, "Error in command") 

def getWorkingDir(analysis_info_file):
    '''Get the working directory specified in the analysis info file
    '''
    ai=read_analysis_info_file(analysis_info_file)
    return(ai['working_directory'])

def getNCores(analysis_info_file):
    '''Get the number of cores specified in the analysis info file
    '''
    ai=read_analysis_info_file(analysis_info_file)
    return(int(ai['ncores']))

def testLogging():
    '''Testing for the logger
    ''' 
    logger = logging.getLogger("runMethylationAnalysis.functions")
    logger.info("Testing logging info level")
    logger.debug("Testing logging debug level")
    logger.error("Testing logging error level")
    try:
        logger.info("Testing stack trace")
        raise RuntimeError
    except Exception, err:
        logger.exception("Expected exception")

def srun_is_installed():
    '''Tests if srun is installed. 

        Uses 'command -v' for POSIX compliance; works in sh and bash.
        When srun is not found, the result will be 1, else null.
    '''
    cmd = 'command -v srun >/dev/null 2>&1 || { echo "1" >&2; }'
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    output = p.stdout.read()
    return(output != "1\n")