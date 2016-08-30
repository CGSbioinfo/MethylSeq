# MethylSeq pipeline


-----------------------------
#### Software and scripts required
-----------------------------




## Run the pipeline

Important note: All the main python scripts have a help page, so if unsure on how to use it or require information about default and additional arguments, check. For example:
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

#### Providing information about the analysis and setting up a working directory.

1) Run the analysis_info.py script  
 **$ python bin/analysis_info.py**  

A file named 'analysis_info.txt' will be created in the folder. Open it in a text editor or vi and fill it. For example: 
>Working directory = /mnt/cgs-fs2/Bioinfo_pipeline/MethylSeq/test/aug2016/heroG/   
>run_folder = /mnt/cgs-fs3/Sequencing/NextSeq_Output/160711_NS500125_0298_AHFW35BGXY/  
>run_samplesheet = /mnt/cgs-fs3/Sequencing/NextSeq_Output/160711_NS500125_0298_AHFW35BGXY/SampleSheet.csv  
>bcl2fastq_output = /mnt/cgs-fs2/Bioinfo_pipeline/MethylSeq/test/aug2016/heroG/fastq/  


2) Next run bcl2fastq.py. This script will run bcl2fastq with the information given in analysis_info.txt. By default, the script finds analysis_info.txt in the current directory. However, if the file has a different name you will need to specify it adding '--analysis_info_file whatever_name.txt'. Use the '-h' argument for details.  
 $ python bin/run_bcl2fastq.py 

9. Create sample_names.txt file. The following command will create a sample_names.txt file with the sample names of the project based on the fastq files that exist in the output from the run_bcl2fastq.py (previous step). You can also specify a different '--in_dir ' argument.
```bash
 $ python bin/create_sampleNames.py
```

10. Run organizeWorkingDirectory.py. This creates a folder called rawReads/ with subfolders corresponding to each sample, and fastq files sorted according to sample name. This structure is useful for further steps in the analysis
```bash
 $ python bin/organizeWorkingDirectory.py 
```

We are ready now to run the main steps of the pipeline
#### Quality control and trimming

11. Our fastq files are organized now in the rawReads/ folder. The next command will run fastqc in all our samples and create a csv file with the number of reads in each sample. The output of the fastqc is created in rawReads/, in the corresponding sample folder. The csv file is created in a Report/figure/dataQC/ folder 
```bash
 $ python bin/qcReads.py
```

12. The next line takes the output of fastqc and makes some tables and plots to include in the Report. The tables are saved in rawReads/sampleX/sampleX_fastqc/ , and the plots are saved by default in Report/figure/rawQC. 
```bash
 $ python bin/fastqc_tables_and_plots.py --in_dir rawReads/ --out_dir_report Report/figure/rawQC --suffix_name _raw --sample_names_file sample_names.txt --plot_device pdf
```

13. Next run the trimming script:
```bash
 $ python bin/trimmingReads.py --in_dir rawReads/ --out_dir trimmedReads --ncores 6 
```

14. The next line takes the output of trimmedReads fastqc and makes some tables and plots to include in the Report. The tables are saved in trimmedReads/, and the output report folder should be indicated.
```bash
 $ python bin/fastqc_tables_and_plots.py --in_dir trimmedReads/ --out_dir_report Report/figure/trimmedQC --suffix_name _trimmed --sample_names_file sample_names.txt --plot_device pdf
```

#### Mapping and de-duplication
15. First make sure that the reference genome is correctly indicated in the analysis_info.txt. Also see the bismark default parameters. You can change or add parameters.
```bash
 $ python bin/mappingReads.py --in_dir trimmedReads/ --out_dir alignedReads --ncores 2 
```









