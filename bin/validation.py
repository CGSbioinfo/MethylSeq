'''Validation functions for methylation analysis
'''
import os
import sys
import functions
import argparse
import csv

__version__='v01'


def validateAnalysisInfo(analysis_info_file):
    '''Validate an analysis_info file for a methylation pipeline.

    This script will read an analysis_info.txt file and 
    check that key parameters have valid inputs.
    This includes checking directories exist. If any
    of the tests fail, the script will exit with status 1, 
    otherwise the script will exit with status 0.

    '''

    # Check that the directories exist:
    # working_directory
    # run_folder
    # run_samplesheet
    # reference_genome

    ai=functions.read_analysis_info_file(analysis_info_file)

    ok = True;

    print "Checking analysis info file"
    
    if not os.path.exists(ai['working_directory']):
        ok = False
        print "Working directory was not found: " + ai['working_directory']

    if not os.path.exists(ai['run_folder']):
        ok = False
        print "Run folder was not found: " + ai['run_folder']

    if not os.path.isfile(ai['run_samplesheet']):
        ok = False
        print "Run samplesheet was not found: " + ai['run_samplesheet']

    if not os.path.exists(ai['reference_genome']):
        ok = False
        print "Reference genome folder was not found: " + ai['reference_genome']

    # Check that the read type is a valid entry
    if not (ai['readType']=='pairedEnd' or ai['readType']=='singleEnd'):
        ok = False
        print "Read type is not valid: " + ai['readType']

    # Check that the number of cores per instance is a valid number
    try:
        ncores = int(ai['ncores'])
        if(ncores==0):
            ok = False
            print "Number of cores must be greater than 0"
    except Exception as e:
        ok = False
        print "Number of cores is not a valid integer: " + ai['ncores']

    # Check that the number of bismark instances is a valid number
    try:
        ninstances = int(ai['ninstances'])
        if(ninstances==0):
            ok = False
            print "Number of instances must be greater than 0"
    except Exception as e:
        ok = False
        print "Number of instances is not a valid integer: " + ai['ninstances']

    # Check the cleanup option
    clean = (ai['clean_files'])
    if(clean!='True' and clean!='False'):
        ok = False
        print "The clean_files option must be either True or False: "+ clean


    # Output result to system


    if ok:
        print "Analysis info file passed validation"
        # sys.exit(0)
    else:    
        print "Analysis info file failed validation"
        # sys.exit(1)
    return(ok)


    # Example analysis_info.txt:
    # working_directory = /mnt/research/bms41/methylSeq/aggressive/2017-06/
    # run_folder = /mnt/research/bms41/methylSeq/151018_NS500125_0141_AHL5VMBGXX/
    # run_samplesheet = /mnt/research/bms41/methylSeq/151018_NS500125_0141_AHL5VMBGXX/SampleSheet.csv
    # bcl2fastq_output = fastq/ 
    # readType = pairedEnd
    # reference_genome = /mnt/cgs-fs3/Sequencing/Genome/Pig/
    # bismark_params = --bowtie2; --bam; -N 0; -L 20; -D 15; -R 2; --score_min L,0,-0.2; 
    # methyl_extract_params= --bedGraph; --gzip; --merge_non_CpG;
    # target_regions_bed = 
    # ncores = 8
    # ninstances = 2
    # clean_files = False

    # Note: currently, target_regions_bed is only called in calculateCoverage.py, which is not in the main pipeline

def checkExperimentName(name):
    '''Check the validity of the given experiment name.

    Replace with appropriate regex at some point.

        name -- the experiment name to test
    '''
    illegalChars = [ "." ]
    ok = True
    for char in illegalChars:
        if char in name:
            ok = False
            print "Illegal character '"+char+"' in name: "+row[1]
    return ok

def validateSampleSheet(analysis_info_file):
    '''Validate an Illumina sample sheet for a methylation pipeline.

    This script will read a sample sheet and check that
    it contains values valid for bcl2fastq. If any
    of the tests fail, the script will exit with status 1, 
    otherwise the script will exit with status 0.
    '''

    ai=functions.read_analysis_info_file(analysis_info_file)

    ok = True;
    
    sampleSheet = ai['run_samplesheet']

    nameKey = "Experiment Name"
    
    with open(sampleSheet, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            if(len(row)>1):
                if(nameKey==row[0]): 
                    ok = checkExperimentName(row[1])

    # Output result to system
    if ok:
        print "Sample sheet passed validation"
    else:    
        print "Sample sheet file failed validation"

    return(ok)

    # Example sample sheet.
    # Known constraints: the experiment name cannot contain periods.

    # [Header],
    # IEMFileVersion,4
    # Experiment Name,NGS-C_Quilter-40078-C
    # Date,12/06/2017
    # Workflow,GenerateFASTQ
    # Application,NextSeq FASTQ Only
    # Assay,SureSelectXT
    # Description,
    # Chemistry,Default

    # [Reads],
    # 76,
    # 76,

    # [Settings],
    # Adapter,GTAGTCCGGCTGACTGACT
    # AdapterRead2,GCTGGCACACAATTACCATA,,,,,,

    # [Data],,,,,,,
    # Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
    # CQ_2169_cortex_aggr,CQ_2169_cortex_aggr,,,A01,ATGCCTAA,NGS-C_Quilter-40078-C,
    # CQ_4557_cortex_aggr,CQ_4557_cortex_aggr,,,B01,GAATCTGA,NGS-C_Quilter-40078-C,
    # CQ_5221_cortex_aggr,CQ_5221_cortex_aggr,,,C01,AACGTGAT,NGS-C_Quilter-40078-C,
    # CQ_1313_cortex_non-aggr,CQ_1313_cortex_non-aggr,,,D01,CACTTCGA,NGS-C_Quilter-40078-C,
    # CQ_1299_cortex_non-aggr,CQ_1299_cortex_non-aggr,,,E01,GCCAAGAC,NGS-C_Quilter-40078-C,
    # CQ_365_cortex_non-aggr,CQ_365_cortex_non-aggr,,,F01,GACTAGTA,NGS-C_Quilter-40078-C,