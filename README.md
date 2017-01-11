# MethylSeq pipeline

----
## Software and scripts required

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
  
  
---- 
## Before running the pipeline

Important note: Most python scripts have a help page, so if unsure on how to use it or require information about default and additional arguments, check. For example:
```bash
 $ analysis_info.py -h
```

### Starting an analysis and downloading the scripts.
**1**\. Go to your home folder. For example:
```bash
 $ cd /mnt/research/mjg225/
```
  
  
**2**\. Create a directory where you will do the analysis and go to that folder. For example:
```bash
 $ mkdir methylSeq/pipelineTest/heroG/aug2016
 $ cd methylSeq/pipelineTest/heroG/aug2016
```
  
  
**3**\. Create a bin/ folder. This folder will have the scripts you need to run the analysis:
```bash
 $ mkdir bin/ 
```
  
  
**4**\. Go to https://github.com/CGSbioinfo/MethylSeq. Click on the green top right button "Clone or download". This will download the scripts to the Downloads/ folder in your home directory (for example: /home/mjg225/Downloads/). Make sure it exists by typing the following (Note that /home/mjg225/ should be replaced with your own home directory): 
```bash
 $ ls /home/mjg225/Downloads/MethylSeq-master
```
This should print on your terminal screen the MethylSeq-master folder
  
  
**5**\. From your currect directory, unzip the file (Note that /home/mjg225/ should be replaced with your own home directory):
```bash
 $ unzip /home/mjg225/Downloads/MethylSeq-master
```
A folder named "MethylSeq-master" should appear. Check by just typing ls.
  
  
**6**\. Move the contents of MethylSeq-master/bin to bin/
```bash
 $ mv MethylSeq-master/bin/* bin/
```
At this point you should have a directory where you will do the analysis and a bin/ folder in such directory with the analysis scripts copied from the github download.  
  
  

**7**\. Run the analysis_info.py script
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
  
  
  

----
## Running the analysis
All the analysis scripts are wrapped in the main python script **runMethylationAnalysis.py**. This main script takes an argument *--run*, which is used to indicate which section of the main script to run. The following commands are used to run the analysis with the main script:  
  
  

### Step 1
1\. Using the command *--run step1_prepare_analysis*, the main script will read the analysis_info_file and run bcl2fastq, create a sample names file, and organize the working directory.
```bash
$ python bin/runMethylationAnalysis.py --run step1_prepare_analysis
```
Once it finishes, there will be a folder named rawReads with fastq files sorted according to sample names. There will also be a sample_names.txt file with a list of sample names, one per line.  
  
  
### Step 2 
2\. Using the command *--run step2_qc_and_trimming*, the main script will read the analysis_info_file, the sample_names file and run fastqc, create a folder Report/figure/rawQC with plots, create a folder Report/figure/data with tables, run trim galore, run fastqc on the trimmed reads, and create a folder Report/figure/trimmedQC with plots:
```bash
$ python bin/runMethylationAnalysis.py --run step2_qc_and_trimming
```
  
  
### Step 3
3\. Using the command *--run step3_mapping_and_deduplication*, the main script will read the analysis_info_file, the sample_names file and run bismark and deduplication scripts. The output includes bam file (original and deduplicated), and log files of each sample, and will be saved in a folder alignedReads/ (created automatically).
```bash
$ python bin/runMethylationAnalysis.py --run step3_mapping_and_deduplication.
```

### Step 4
4\. Bismark outputs a bam file with the mapped reads and a report about the alignment.    
The following command uses an executable Rscript which summarises mapping QC metrics. 

Arguments:   
>Rscript bin/mappingQC.R *\<input folder containing bam files*\> *\<sample names file*\> *\<suffix pattern of report files output of bismark*\> *\<outdir*\>   

Example:
```bash
$ Rscript bin/mappingQC.R /mnt/research/jb393/MethylSeq_Pilot/Aligned_data/Raw_bam/ sample_names_all.txt _bismark_bt2_PE_report.txt Report/figure/mappingQC/
```


### Step 5
5\. Run the methylation extraction  
```bash 
 $ python bin/methylationExtraction.py
```
This creates 3 output files per sample: bedGraph.gz, bismark.cov.gz, and M-bias.txt.   
There is also a log file per sample: methylExtract_log_sampleName.txt.   
The M-bias.txt sample will be used in the next step to detect any bias in the %Methylation across the reads"'" positions.     


### Step 6
6\. Run the mbias plot   
Arguments:
>Rscript bin/methylExtractQC_mbias_plot.R *\<input folder containing .M-bias.txt files*\> *\<sample names file*\> *\<suffix pattern of M-bias.txt output of bismark*\> *\<outdir*\>

```bash
 $ Rscript bin/methylExtractQC_mbias_plot.R alignedReads/ sample_names.txt .M-bias.txt Report/figure/methExtractQC/
```
This creates a plot in the specified outdir with the %Methylation across reads"'" positions.   
Based on this plot, we need to decide whether or not to trim bases from 5p and 3p for each sample. 



### Step 7
7\. Create and fill a file *mbias_remove_bases.txt* with information about which bases to clip from reads.   
Arguments:
>python bin/remove_bases_file_info.py --outfile *\<name of output txt file*\>   

```bash
 $ python bin/remove_bases_file_info.py --outfile remove_bases.txt
```
A file with the specified name will be creates. Open the file and fill it with the following information (one sample per line):
>sample: *\<sample name*\>   
>5R1: *\<number of bases to clip from the 5 prime end from read 1 (forward read)*\>  
>3R1: *\<number of bases to clip from the 3 prime end from read 1 (reverse read)*\>   
>5R2: *\<number of bases to clip from the 5 prime end from read 2 (forward read)*\>   
>3R2: *\<number of bases to clip from the 3 prime end from read 2 (reverse read)*\>   


### Step 8
8\. Run the methyl extraction again, removing biased bases from reads with information in mbias_remove_bases.txt 
 $ python bin/methylationExtraction_removeBases.py --ncores 1 --sample_names_file sample_names_test.txt --out_dir alignedReads/clean --remove_bases_file mbias_remove_bases.txt
