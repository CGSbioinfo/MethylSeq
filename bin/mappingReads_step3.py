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

def alignment(i):
    trimmedReads = os.listdir("./trimmedReads")
    trimmedReads = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall(i, x)]
    r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2_val_2.fq", x)]
    if r2:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1_val_1.fq", x)]
        r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2_val_2.fq", x)]
        os.system("bismark -n 1 --bowtie2 --bam --chunkmbs 1000 "+"-o "+path+"alignedReads/1pass_ann/ "+refGenome +" -1 "+path+"/trimmedReads/" + r1[0] + " -2 "+path+"/trimmedReads/" + r2[0]+" &>log_"+i+".txt")
        logging.info('Sample ' + i + ' running, pairedEnd mode')
    else:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_all_lane_R1_trimmed.fq", x)]
        os.system("bismark -n 1 --bowtie2  --bam --chunkmbs 1000 "+"-o "+path+"alignedReads/1pass_ann/ "+ refGenome +" "+path+"/trimmedReads/" + r1[0] + ' >'+path+'alignedReads/1pass_ann/'+" &>log_"+i+".txt" )
        logging.info('Sample ' + i + ' running, singleEnd mode')


def deduplicate(i):
    os.system("deduplicate_bismark --paired --bam "+path+"/alignedReads/1pass_ann/" + i+"_bismark_bt2_pe.bam ")
    logging.info("..duplicate removal finish for sample: "+i)

def methylCall(i):
    os.system("bismark_methylation_extractor -p --no_overlap --comprehensive -o "+path+"/MethylCall "+i+"_bismark_bt2_pe.deduplicated.bam")
    logging.info("Methylation call done for sample: "+i)



###################
#  * Rscripts used:
# mapping_summary.R
# mapping_distribution.R
###################

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
print(logFilename)
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log.txt",
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')

logging.info(" ")
logging.info(" ")
logging.info("***************************************")
logging.info("*************MAPPING READS*************")
logging.info(" ")
logging.info("User command: " + str(sys.argv))


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... READING PARAMETERS FILE")
params_file = sys.argv[1]
args = {}
with open(params_file, 'r') as f:
    for line in f:
        entry = line.strip().split("=")
        if entry[0]:
            args[entry[0].strip(" ")] = entry[1].strip()
logging.info('-- Input parameters:')
path = args['Working directory'].replace("\\", "")
logging.info('Working Directory = ' + path)
refGenome = args['Reference Genome']
logging.info('Reference Genome = ' + refGenome)
nsamples = int(args['Number of samples'])
#logging.info("\n -- Input parameters: \n Working Directory = " + path + "\n GTF file = " + gtfFile + "\n Reference Genome = " + refGenome + "\n Number of samples = " + str(nsamples) + "\n BedFile = " + bedFile + "\n BedFile_10k = " + bedFile_10k + "\n refFlat = " + refFlat + "\n rRNA_interval_list = " + rRNA_interval_list + "\n strand = " + strand )
os.chdir(path)

logging.info(" ")
logging.info(" ")
logging.info("#################################")
# Read in sampleNames
logging.info('... Sample names:')
sampleNames = []
sample_names_file = open('sample_names.txt','r')
for line in sample_names_file:
    sampleNames.append(line.strip())
    logging.info(line.strip())


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Mapping of reads")
logging.info("Samples to be mapped: ")
make_sure_path_exists("./alignedReads/1pass_ann")
alignmentOut = os.listdir(path+"/alignedReads/1pass_ann") # Check if aligned reads files exist
bamFiles = get_filepaths(path)
bamFiles= [bamFiles[i] for i, x in enumerate(bamFiles) if re.findall(".bam", x)]

sampleToRun=list()
for sample in sampleNames:
    found=False
    for files in bamFiles:
        #print(files)
        file=files.split("_all_lane")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToRun.append(sample)

logging.info("Sample listing done.... Starting mapping...")
if len(sampleToRun)>0:
#    for i in sampleToRun:
#       print(i)
#
#    sys.exit("done")
    Parallel(n_jobs=4)(delayed(alignment)(i) for i in sampleToRun)
logging.info("... Finished Mapping")

# Removing duplicate
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... removing duplicate reads")
logging.info("Samples to be processed: ")
dedupOut = os.listdir(path+"/alignedReads/1pass_ann")
dedupFiles = get_filepaths(path)
dedupFiles= [dedupFiles[i] for i, x in enumerate(dedupFiles) if re.findall(".deduplicated", x)]
sampleToDedup=[]
for sample in sampleNames:
    found=False
    for files in dedupFiles:
        #print(files)
        file=files.split("_all_lane")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToDedup.append(sample)

logging.info("Sample listing done.... Starting duplicates removal: ")
if len(sampleToDedup)>0:
    Parallel(n_jobs=8)(delayed(deduplicate)(i) for i in sampleToDedup)

logging.info("Duplicates removal completed.... Moving the files ")
dedupFiles = get_filepaths(path)
dedupFiles= [dedupFiles[i] for i, x in enumerate(dedupFiles) if re.findall(".deduplicated", x)]
make_sure_path_exists("./alignedReads/Deduplicated")
for dedup in dedupFiles:
    os.system("mv "+dedup+" ./alignedReads/Deduplicated")

logging.info("... Finished moving")
logging.info(" ")

logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Methylation Calls")
logging.info("Sample to call: ")
make_sure_path_exists("./MethylCall")
dedupFiles = get_filepaths(path)
dedupFiles= [dedupFiles[i] for i, x in enumerate(dedupFiles) if re.findall(".deduplicated", x)]
sampleToMethCall=[]
for sample in sampleNames:
    found=False
    for files in dedupFiles:
        #print(files)
        file=files.split("_all_lane")[0]
        file2=file.split("/")[-1]
        if sample == file2:
            found=True
    if not found:
        logging.info(sample)
        sampleToMethCall.append(sample)

if len(sampleToMethCall)>0:
    Parallel(n_jobs=8)(delayed(methylCall)(i) for i in sampleToMethCall)

logging.info("Methylation call done")



##extract methylation call
# logging.info(" ")
# logging.info(" ")
# logging.info("#################################")
# logging.info("... remove duplicate")
# dedupOut = os.listdir(path+"/alignedReads/1pass_ann")
# indicesdedupFiles = [i for i, x in enumerate(dedupOut) if re.findall(".dedup*", x)] # Check if aligned reads files exist
# if not indicesdedupFiles:
#     Parallel(n_jobs=8)(delayed(deduplicate)(i) for i in sampleNames)
# logging.info("... Finished removing duplicates")
# logging.info('Mapping summary done')
# logging.info(" ")


# Index Bam files
# logging.info(" ")
# logging.info(" ")
# logging.info("#################################")
# logging.info("... SAM 2 BAM files")
# alignmentOut = os.listdir(path+"/alignedReads/1pass_ann") # Check if aligned reads files exist
# indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".bam", x)] # Check if aligned reads files exist
# if not indicesAlignmentFiles:
#     Parallel(n_jobs=8)(delayed(sam2bam)(i) for i in sampleNames)
# logging.info("... Finished SAM 2 BAM")
#
#
#
#
# # Mapping QC
# logging.info(" ")
# logging.info(" ")
# logging.info("#################################")
# logging.info("... Mapping QC")
# qcOut = os.listdir(path+"/alignedReads/1pass_ann")
# indicesQCFiles = [i for i, x in enumerate(qcOut) if re.findall("_bam_stats.txt", x)] # Check if aligned reads files exist
# if not indicesQCFiles:
#     Parallel(n_jobs=8)(delayed(bamstats)(i) for i in sampleNames)
# logging.info("... Finished generating mapping QC")
# logging.info('Mapping summary done')
# logging.info(" ")
#
# #sort by coordinate
# logging.info(" ")
# logging.info(" ")
# logging.info("#################################")
# logging.info("... Sorting BAM files by coord ")
# # Sort bam files by name
# alignmentOut = os.listdir(path+"/alignedReads/1pass_ann") # Check if aligned reads files exist
# indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".sortedByCoord.out.bam", x)] # Check if aligned reads files exist
# if not indicesAlignmentFiles:
#     Parallel(n_jobs=8)(delayed(sortByCoord)(i) for i in sampleNames)
# logging.info("... Finished sorting BAM files by name")






# # Index Bam files
# logging.info(" ")
# logging.info(" ")
# logging.info("#################################")
# logging.info("... Indexing BAM files")
# alignmentOut = os.listdir(path+"/alignedReads/1pass_ann") # Check if aligned reads files exist
# indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".bai", x)] # Check if aligned reads files exist
#
# if not indicesAlignmentFiles:
#     Parallel(n_jobs=8)(delayed(indexing)(i) for i in sampleNames)
# logging.info("... Finished Indexing")


# logging.info(" ")
# logging.info(" ")
# logging.info("#################################")
# logging.info("... Convert bam to bed")
# bedOut = os.listdir(path+"/alignedReads/1pass_ann")
# indicesbedFiles = [i for i, x in enumerate(bedOut) if re.findall(".bed", x)] # Check if aligned reads files exist
# if not indicesbedFiles:
#     Parallel(n_jobs=8)(delayed(bam2bed)(i) for i in sampleNames)
# logging.info("... Finished converting bam to bed")
# logging.info('Mapping summary done')
# logging.info(" ")



logging.info(" ")
logging.info(" ")
logging.info(" ")
logging.info("##################################################################")
logging.info("------------------------------------------------------------------")
logging.info("##################################################################")
logging.info(" ")
logging.info(" ")
logging.info(" ")
# GeneBody Coverage
# make_sure_path_exists("./alignedReads/1pass_ann/QC")
# make_sure_path_exists("./Report/figure/")
# alignmentQCOut = os.listdir(path+"/alignedReads/1pass_ann/QC") # Check if aligned reads files exist
# indicesAlignmentQCFiles = [i for i, x in enumerate(alignmentQCOut) if re.findall(".geneBodyCoverage.", x)] # Check if aligned reads files exist
# if not indicesAlignmentQCFiles:
#     os.system("ls ./alignedReads/1pass_ann/*Aligned.sortedByCoord.out.bam > tempbamfiles.txt")
#     os.system("python ~/bin/geneBody_coverage.py -r " + bedFile_10k + " -i tempbamfiles.txt -o alignedReads/1pass_ann/QC/10KGenes")
#     os.system("rm tempbamfiles.txt")
# os.system('cp alignedReads/1pass_ann/QC/10KGenes.geneBodyCoverage.curves.pdf Report/figure/10KGenes_geneBodyCoverage_curves.pdf')
# logging.info('GeneBody Coverage done')
# logging.info(" ")
#
# # Junctions and junction saturation
# qcOut = os.listdir(path+"/alignedReads/1pass_ann/QC")
# qcFiles = [i for i, x in enumerate(qcOut) if re.findall(".junction.", x)] # Check if aligned reads files exist
# if not qcFiles:
#      Parallel(n_jobs=8)(delayed(junctions)(i) for i in sampleNames)
# os.system('grep "y=c(" alignedReads/1pass_ann/QC/*junctionSaturation*  | sed \'s/:y=c(/,/g\' | sed \'s/.junctionSaturation_plot.r//g\' | sed \'s/)//g\' | sed \"s/.*\///g\"  > alignedReads/1pass_ann/QC/junctionSat_all.csv')
# os.system('~/bin/junctionPlotAll.R .')
# logging.info(" ")
#
# # Collect Metrics
# qcOut = os.listdir(path+"/alignedReads/1pass_ann/QC") # Check if aligned reads files exist
# qcFiles = [i for i, x in enumerate(qcOut) if re.findall(".metrics.txt", x)] # Check if aligned reads files exist
# logging.info("Calculating RNASeq Metrics")
# if not qcFiles:
#     os.system("find alignedReads/1pass_ann/*.sortedByCoord.out.bam | sed \"s/.*\///g\" | sed \"s/.sortedByCoord.out.bam//g\" | " # get file names and format them
#           "parallel -j 3 --no-notice "
#           "\"java -jar ~/tools/picard-tools-1.127/picard.jar CollectRnaSeqMetrics "
#           "REF_FLAT=" + refFlat + " "
#           "RIBOSOMAL_INTERVALS=" + rRNA_interval_list + " "
#           "STRAND_SPECIFICITY=" + strand + " "
#           "INPUT=alignedReads/1pass_ann/{}.sortedByCoord.out.bam "
#           "OUTPUT=alignedReads/1pass_ann/QC/{}_metrics.txt \"")
#     logging.info("Calculated RNASeq Metrics")
# def pct(i):
#     os.system('mapping_distribution.R ./alignedReads/1pass_ann/QC/ ' + i)
# Parallel(n_jobs=8)(delayed(pct)(i) for i in sampleNames)
# logging.info("... Finished mapping QC")
#



