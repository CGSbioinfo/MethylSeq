#!/usr/bin/env python
import sys

# USAGE: python fastq_data.txt per_base_sequence_content.py
fastq_data = sys.argv[1].replace('\\','')
plot = sys.argv[2] # all, per_base_sequence_content, per_base_sequence_quality, per_sequence_gc_content, per_sequence_quality_scores, seq_dup_levels, kmer_content
outdir = sys.argv[3] + '_'

dictionary = {'cpg_r1':'CpG context (R1)','chg_r1':'CHG context (R1)', 'chh_r1':'CHH context (R1)', 'cpg_r2':'CpG context (R2)', 'chg_r2':'CHG context (R2)', 'chh_r2':'CHH context (R2)'}

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
            if i > line_number+1:
                if line != '':
                    output.write(line + '\n')
                else:
                    break
            i += 1
        f.close()
else:
    f = open(fastq_data,'r')
    i=0
    next(f)
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
        if i > line_number+1:
            if line != '':
                output.write(line + '\n')
            else:
                break
        i += 1

