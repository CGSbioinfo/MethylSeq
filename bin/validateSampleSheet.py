#!/usr/bin/env python 

'''
Validate an Illumina sample sheet for a methylation pipeline.

This script will read a sample sheet and check that
it contains values valid for bcl2fastq. If any
of the tests fail, the script will exit with status 1, 
otherwise the script will exit with status 0.
'''

import os
import csv
import sys
import functions
import argparse

__version__='v01'

def checkExperimentName(name):
    '''Check the validity of the given experiment name.

    Replace with appropriate regex at some point.

        name -- the experiment name to test
    '''
    ok = True
    for char in illegalChars:
        if char in name:
            ok = False
            print "Illegal character '"+char+"' in name: "+row[1]
    return ok


if __name__ == '__main__':

    # Parser
    parser = argparse.ArgumentParser(prog='validateSampleSheet.py',
        description = 'Check the sample sheet was filled in correctly')
    parser.add_argument('-v','--version', 
        action='version',
        version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', 
        help='Text file with details of the analysis. Default=analysis_info.txt', 
        default='analysis_info.txt')
    args=parser.parse_args()

    # Read the sample sheet
    ai=functions.read_analysis_info_file(args.analysis_info_file)

    ok = True;
    
    sampleSheet = ai['run_samplesheet']

    nameKey = "Experiment Name"


    illegalChars = [ "." ]
        
    with open(sampleSheet, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            if(len(row)>1):
                if(nameKey==row[0]): 
                    ok = checkExperimentName(row[1])

    # Output result to system
    if ok:
        print "Sample sheet passed validation"
        sys.exit(0)
    else:    
        print "Sample sheet file failed validation"
        sys.exit(1)



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