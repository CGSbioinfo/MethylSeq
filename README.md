# MethylSeq pipeline

----
#### Software and scripts required

Software: 
- bcl2fastq
- FastQC
- Trim Galore
- Bismark

R packages:
- ggplot2
- reshape
- grid
- xtable


---- Before running the pipeline
## Run the pipeline

Important note: Most python scripts have a help page, so if unsure on how to use it or require information about default and additional arguments, check. For example:
```bash
 $ analysis_info.py -h
```

#### Starting an analysis and downloading the scripts.
1) Go to your home folder. For example:
```bash
 $ cd /mnt/research/mjg225/
```

2) Create a directory where you will do the analysis and go to that folder. For example:
```bash
 $ mkdir methylSeq/pipelineTest/heroG/aug2016
 $ cd methylSeq/pipelineTest/heroG/aug2016
```

3) Create a bin/ folder. This folder will have the scripts you need to run the analysis:
```bash
 $ mkdir bin/ 
```

4) Go to https://github.com/CGSbioinfo/MethylSeq. Click on the green top right button "Clone or download". This will download the scripts to the Downloads/ folder in your home directory (for example: /home/mjg225/Downloads/). Make sure it exists by typing the following (Note that /home/mjg225/ should be replaced with your own home directory): 
```bash
 $ ls /home/mjg225/Downloads/MethylSeq-master
```
This should print on your terminal screen the MethylSeq-master folder

5) From your currect directory, unzip the file (Note that /home/mjg225/ should be replaced with your own home directory):
```bash
 $ unzip /home/mjg225/Downloads/MethylSeq-master
```
A folder named "MethylSeq-master" should appear. Check by just typing ls.

6) Move the contents of MethylSeq-master/bin to bin/
```bash
 $ mv MethylSeq-master/bin/* bin/
```

At this point you should have a directory where you will do the analysis and a bin/ folder in such directory with the analysis scripts copied from the github download.


--------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------

#### Providing information about the analysis, creating fastq files and setting up a working directory.

1\. Run the analysis_info.py script
```bash  
$ python bin/analysis_info.py
```

A file named 'analysis_info.txt' will be created in the folder. Open it in a text editor or vi and fill it. For example: 
>Working directory = /mnt/cgs-fs2/Bioinfo_pipeline/MethylSeq/test/aug2016/heroG/   
>run_folder = /mnt/cgs-fs3/Sequencing/NextSeq_Output/160711_NS500125_0298_AHFW35BGXY/  
>run_samplesheet = /mnt/cgs-fs3/Sequencing/NextSeq_Output/160711_NS500125_0298_AHFW35BGXY/SampleSheet.csv  
>bcl2fastq_output = /mnt/cgs-fs2/Bioinfo_pipeline/MethylSeq/test/aug2016/heroG/fastq/  
>readType = pairedEnd  
>reference_genome =  
>bismark_params = --bowtie2; --bam; --directional; -N 0; -L 20; --no-mixed; --no-discordant; -D 15; -R 2; --score_min L,0,-0.2;  
>methyl_extract_params= --bedGraph; --gzip; --merge_non_CpG;  
>target_regions_bed =  
>ncores = 8  

Explanation of 'analysis_info.txt':
>Working directory = *\<path to directory of the analysis\>*  
>run_folder = *\<path to the run folder\>*  
>run_samplesheet = *\<sample sheet used to generate fastq files. This is created using the Illumina Expert Manager\>*  
>bcl2fastq_output = *\<path to the desired output of bcl2fastq. The defaults is fastq/ and the folder will be created automatically\>*  
>readType = *\<either pairedEnd or singleEnd\>*
>reference_genome = *\<path to the bismark reference genome that will be used at the mapping step\>*  
>bismark_params = *\<parameters passed to the bismark script. Leave the default, and you can add additional parameters separated by ";"\>*  
>methyl_extract_params= *\<parameters passed to the methyl extractor script. Leave the default, and you can add additional parameters separated by ";"\>*  
>target_regions_bed = *\<path to bedfile\>*  
>ncores = *\<Number of cores to use to pararellize analysis\>*

#### Running the analysis
All the analysis scripts are wrapped in the main python script **runMethylationAnalysis.py**. This main script takes an argument *--run*, which is used to indicate which section of the main script to run. The following commands are used to run the analysis with the main script:  

##### Step 1
1\. Using the command *--run step1_prepare_analysis*, the main script will read the analysis_info_file and run bcl2fastq, create a sample names file, and organize the working directory.
```bash
$ python bin/runMethylationAnalysis.py --run step1_prepare_analysis
```
Once it finishes, there will be a folder named rawReads with fastq files sorted according to sample names. There will also be a sample_names.txt file with a list of sample names, one per line.  

##### Step 2 
2\. Using the command *--run step2_qc_and_trimming*, the main script will read the analysis_info_file, the sample_names file and run fastqc, create a folder Report/figure/rawQC with plots, create a folder Report/figure/data with tables, run trim galore, run fastqc on the trimmed reads, and create a folder Report/figure/trimmedQC with plots:
```bash
$ python bin/runMethylationAnalysis.py --run step2_qc_and_trimming
```

##### Step 3
3\. Using the command *--run step3_mapping_and_deduplication*, the main script will read the analysis_info_file, the sample_names file and run bismark and deduplication scripts. The output includes bam file (original and deduplicated), and log files of each sample, and will be saved in a folder alignedReads/ (created automatically).
```bash
$ python bin/runMethylationAnalysis.py --run step3_mapping_and_deduplication.
```



