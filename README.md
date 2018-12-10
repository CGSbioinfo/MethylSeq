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
- dplyr
- purrr
- tidyr
- BiSeq (optional)
- sva (optional)
- parallel

---- 
## Before running the pipeline

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
 $ mkdir src/ 
```
  
  
**4**\. Go to https://github.com/CGSbioinfo/MethylSeq. Click on the green top right button "Clone or download". This will download the scripts to the Downloads/ folder in your home directory (for example: /home/mjg225/Downloads/). Make sure it exists by typing the following (Note that /home/mjg225/ should be replaced with your own home directory): 
```bash
 $ ls /home/mjg225/Downloads/MethylSeq-master
```
This should print on your terminal screen the MethylSeq-master folder
  
  
**5**\. From your currect directory, unzip the file (Note that ```/home/mjg225/``` should be replaced with your own home directory):
```bash
 $ unzip /home/mjg225/Downloads/MethylSeq-master
```
A folder named ```MethylSeq-master``` should appear. Check by just typing ls.
  
  
**6**\. Move the contents of ```MethylSeq-master/src``` to ```src/```
```bash
 $ mv MethylSeq-master/src/* src/
```
At this point you should have a directory where you will do the analysis and a src/ folder in such directory with the analysis scripts copied from the github download.  
  
----
## Running the analysis

All the analysis scripts are wrapped in the main python script ```runMethylationAnalysis.py```. This main script takes an argument ```--run```, which is used to indicate which section of the main script to run. A log file is created in the working directory called ```analysis.log```. 

The analysis is intended to be run on a slurm cluster. The main script should be invoked without srun; this will then call ```ninstances``` of slurm, each with ```ncores``` cores. The reason for this is that if the main script is invoked through srun, there will only be ```ncores``` cores available for subprocesses, which will be shared between the ```ninstances```. If slurm is not installed, the script will still work, but analyses will be correspondingly slower.

The following commands are used to run the analysis pipeline:  


### Step 0

Run the setup to create an info file that contains parameters for the subsequent programs to be run.

```bash  
$ python src/runMethylationAnalysis.py --run step0_create_info
```
A file named ```analysis_info.txt``` will be created in the folder. Open it in a text editor or vi and fill it. For example: 
>working_directory = /mnt/cgs-fs2/Bioinfo_pipeline/MethylSeq/test/aug2016/heroG/   
>run_folder = /mnt/cgs-fs3/Sequencing/NextSeq_Output/160711_NS500125_0298_AHFW35BGXY/  
>run_samplesheet = /mnt/cgs-fs3/Sequencing/NextSeq_Output/160711_NS500125_0298_AHFW35BGXY/SampleSheet.csv  
>bcl2fastq_output = /mnt/cgs-fs2/Bioinfo_pipeline/MethylSeq/test/aug2016/heroG/fastq/  
>readType = pairedEnd  
>reference_genome =  /mnt/cgs-fs3/Sequencing/Genome/Sscrofa11.1/  
>bismark_params = --bowtie2; --bam; --directional; -N 0; -L 20; --no-mixed; --no-discordant; -D 15; -R 2; --score_min L,0,-0.2;  
>methyl_extract_params= --bedGraph; --gzip; --merge_non_CpG;  
>target_regions_bed =  /mnt/research/bms41/methylSeq/SureSelect/Sscrofa11.1.target.bed  
>ncores = 8  
>ninstances = 3  
>clean_files = False  
>sill = 1  
>dmr_bandwidth = 50  
>gtf_file = /mnt/cgs-fs3/Sequencing/Genome/Sscrofa11.1/annotation/Sus_scrofa.Sscrofa11.1.90.gtf  
>extract_cx = False  
tissue = Cortex

Field | Description | Values
------|------|------
working_directory | path to directory of the analysis | 
run_folder | path to the run folder |
run_samplesheet | Sample sheet used to generate fastq files. This is created using the Illumina Experiment Manager |
bcl2fastq_output | path to the desired output for bcl2fastq. The default is fastq/ and the folder will be created automatically | 
readType | The type of sequencing read | ```pairedEnd``` or ```singleEnd```
reference_genome | path to the bismark reference genome that will be used at the mapping step |
bismark_params | parameters passed to the bismark script. Leave the default, and you can add additional parameters separated by ";" |
methyl_extract_params | parameters passed to the methyl extractor script. Leave the default, and you can add additional parameters separated by ";" |
target_regions_bed | path to BED format file of targeted genomic regions |
ninstances | Number of slurm instances to run in parallel | Integer >= 1
ncores | Number of cores to use in each slurm instance to pararellize analysis | Integer >= 1
clean_files | Should files be cleaned | ```True``` or ```False```
sill | The sill value when smoothing variograms |  0-1
dmr_bandwidth | The DMR smoothing width | Integer >=10
gtf_file | The gene annotations for the reference genome in GTF format | 
extract_cx | Should methylation calls be extracted in CX context | ```True``` or ```False```  
tissue | the source of the sample | any sensible string

  
### Step 1
Using the command ```--run step1_prepare_analysis```, the main script will read the analysis_info_file and run bcl2fastq, create a sample names file, and organize the working directory.

```bash
$ python src/runMethylationAnalysis.py --run step1_prepare_analysis
```
Once it finishes, there will be a folder ```rawReads/``` with fastq files sorted according to sample names. There will also be a ```sample_names.txt``` file with a list of sample names, one per line.  

From the sample name file, create a sample group file. This is a tab separated file that contains extra information about the samples. The minimum required columns are ```Sample.Name```, ```Group``` and ```Tissue```. Other columns can be added as needed if modifying the downstream scripts. Example:

```text
Sample.Name Group   Tissue
CQ_1299_anteriorcortex_non-aggr_S6  Non-aggr    Cortex
CQ_1301_hypothalamus_aggrs_S3   Aggr    Hypothalamus
```
  

### Step 2 
Using the command ```--run step2_qc_and_trimming```, the main script will read the analysis_info_file, the sample_names file and run fastqc, create a folder ```Report/figure/rawQC``` with plots, create a folder ```Report/figure/data``` with tables, run trim galore, run fastqc on the trimmed reads, and create a folder ```Report/figure/trimmedQC``` with plots:

```bash
$ python src/runMethylationAnalysis.py --run step2_qc_and_trimming
```
  
  
### Step 3
Using the command ```--run step3_mapping_and_deduplication```, the main script will read the analysis_info_file, the sample_names file and run bismark and deduplication scripts. The output includes bam file (original and deduplicated), and log files of each sample, and will be saved in a folder ```alignedReads/``` (created automatically).
```bash
$ python src/runMethylationAnalysis.py --run step3_mapping_and_deduplication
```

### Step 4
After checking the QC metrics in ```Report/figure/mappingQC/```, extract the methylation data using the command ```--run step4_extract_methylation```.
```bash
$ python src/runMethylationAnalysis.py --run step4_extract_methylation
```

This creates 3 output files per sample: ```bedGraph.gz```, ```bismark.cov.gz```, and ```M-bias.txt```.   
There is also a log file per sample: ```methylExtract_log_sampleName.txt```.   
The M-bias.txt sample is used to detect any bias in the %Methylation across the reads positions.  The data is plotted in the files ```Report/figure/methExtractQC/Mbias_plot.pdf```. Based on this plot, we need to decide whether or not to trim bases from 5' and 3' for each sample. A file will be created in the working directory called ```remove_bases.txt```

Fill in the file ```remove_bases.txt``` with information about which bases to clip from reads.   

Open the file and fill it with the following columns (one sample per line):

>sample: *\<sample name*\>   
>5R1: *\<number of bases to clip from the 5 prime end from read 1 (forward read)*\>  
>3R1: *\<number of bases to clip from the 3 prime end from read 1 (reverse read)*\>   
>5R2: *\<number of bases to clip from the 5 prime end from read 2 (forward read)*\>   
>3R2: *\<number of bases to clip from the 3 prime end from read 2 (reverse read)*\>   

Example:

sample  | 5R1 | 3R1 | 5R2 | 3R2
------|-----|-----|-----|-----
CQ_4557_cortex_aggr_S2 | 9 | 6 | 10 | 7  
CQ_365_cortex_non-aggr_S6 |  9   |   6   |  10  |  7
CQ_5221_cortex_aggr_S3  |  9  |   6  |   10 |   7
CQ_2169_cortex_aggr_S1  |  9  |   6  |   10 |   7
CQ_1299_cortex_non-aggr_S5 |   9  |   6   |  10  |  7
CQ_1313_cortex_non-aggr_S4  |  9  |   6   |  10  |  7

### Step 5
Run the methyl extraction again. This will remove biased bases from reads using the information in  ```remove_bases.txt```.    

```bash
$ python src/runMethylationAnalysis.py --run step5_extract_methylation
```
The results will be output to ```Report/figure/methExtractQC/Mbias_plot_clipped.pdf```
Confirm that the clipping of bases worked.

The sequences coverage in the targeted regions and the whole genome will be output as pdfs to ```Report/figure/methExtractQC/```. Check that the coverage is enriched in target regions and sufficient for differential analysis.

### Step 6 - Optional

Run differential methylation analysis using BiSeq. Familiarity with BiSeq analysis is expected.

```bash
$ python bin/runMethylationAnalysis.py --run step6_analyse_methylation
```

Logging for the differential analysis is output to console and to the file ```biseqAnalysis/log.txt```.

The BiSeq analysis is quite long, so progress is saved as ```.Rdata``` files at various points. If the main script is cancelled at any step, it will detect the last save point when invoked again and resume from there.

There are instances where further action is needed within the script.

#### 1) Beta regression

The beta regression step is time-consuming, and so the script has been designed to allow greater parallelisation if desired. Once the regression has begun the script will report ```Parallel script can now be invoked```.

At this point multi-node parallelisation can be started using a secondary R script, ```parallelBetaRegression.R```. The main script saves out chunks of data that can be worked on independently; the parallel script works through any unfinished data chunks.

The parallel script takes some command line parameters; the required parameters have been generated by the main script and logged to the console for you to copy and paste. For reference, these paramters are as follows:

Field | Description | Values
------|------|------
temp_folder | the temporary folder for the regression | Use ```biseqAnalysis/``` unless you known why you should change this
ncores | number of cores to use in this invokation | Integer >=1
is_null_regression | is this is the null regression step | ```T``` if this is the null regression step, or ```F``` if this is the first beta regression
log_file | the file to write logging information to | a valid file path

Example:

```bash
 $ srun -c 20 --mem=20G /usr/bin/Rscript src/parallelBetaRegression.r biseqAnalysis 20 F ./biseqAnalysis/log.txt

```

Multiple invokations of the parallel script can be made; the main script will wait until all parallel instances have completed.

#### 2) Setting the sill value for variogram smoothing

Following variogram generation, the main script will exit to allow a sill value to be chosen. The variogram is saved to ```biseqAnalysis/Null_variogram.png```. Enter the desired sill in the analysis info file and rerun analysis step 6. The analysis will resume at the correct stage, and the smoothed variogram will be saved to ```biseqAnalysis/Smooth_variogram.png```

```bash
$ python src/runMethylationAnalysis.py --run step6_analyse_methylation
```

The detected DMRs will be exported to ```Report/figure/DMR_results/```. Tables with DMRs and DMRs overlapping genes are exported to ```Report/figure/DMR_results/DMRs_overlapping_target_regions.tsv``` and ```Report/figure/DMR_results/DMRs_overlapping_annotated_genes.tsv``` respectively.

### Step 7 - Optional

Run differential methylation analysis using SVA to correct for unknown cell type composition. Some familiarity with surrogate variable analysis is expected in order to interpret the results and adjust the script if necessary.

```bash
$ python src/runMethylationAnalysis.py --run step7_SVA_analysis
```

Logging for the differential analysis is output to console and to the file ```SVA/{date}.SVA.log.txt```.

Since SVA is performing a CpG by CpG comparison of groups using limma, it needs a more stringent read threshold than BiSeq. The default parameters are to filter out loci with fewer than 5 reads in any sample. It also removes the top 1% most frequent reads from each sample, to avoid bias from PCR amplifiction. The total read counts and filtered read counts are exported to images.

As with the BiSeq analysis, intermediate analysis steps are stored so you can modify and rerun sections if needed without repeating the whole time-consuming process.

The output of the SVA is a list of significant CpGs at the p=0.05 threshold after FDR correction, with a median methylation difference of at least 5% between groups. CpGs overlapping genes are also exported, plus a list of the unique genes.

If a BiSeq analysis was run, the DMRs are compared to the SVA loci, and any overlaps are also exported.