# MethylSeq pipeline


-----------------------------
#### Software and scripts required
-----------------------------




## Run the pipeline

Important note: All the main python scripts have a help page, so if unsure on how to use it or require information about default and additional arguments, check. For example:
```bash
 $ analysis_info.py -h
```

#### Setting up the analysis
1. Go to your home folder. For example:
```bash
 $ cd /mnt/research/mjg225/
```

2. Create a directory where you will do the analysis and go to that folder. For example:
```bash
 $ mkdir methylSeq/pipelineTest/heroG/aug2016
 $ cd methylSeq/pipelineTest/heroG/aug2016
```

3. Create a bin/ folder. This folder will have the scripts you need to run the analysis:
```bash
 $ mkdir bin/
 $ 
```



1.	Go to the main folder of the project and run
```bash
 $ analysis_info.py
```
This will create a file named analysis_info.txt, which needs to be filled in a text editor. For more details about how to fill this file, see analysis_info_guide.txt

2.	Create a sample_names.txt file with the list of the sample names. If your fastq files are saved in the reads/ subfolder, you can run
```bash
 $ ls reads/*fastq.gz | sed 's/.*\///g' | sed 's/_R.*//g' | uniq > sample_names.txt
```
The names of the samples normally consist of sampleX_S1, sampleY_S2, etc.









