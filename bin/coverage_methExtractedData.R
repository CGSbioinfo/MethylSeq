#!/usr/bin/R

cat("Loading R packages ...\n")

#' This function loads the external functions file.
#' @title Load external functions
#' @export
loadFunctionsFile = function(){
  cat("Loading functions ...\n")
  file.arg.name = "--file="
  script.name   = sub(file.arg.name, "", commandArgs()[grep(file.arg.name, commandArgs())])
  script.folder = dirname(script.name)
  script.to.load = paste(sep="/", script.folder, "functions.r")
  source(script.to.load)
}

loadFunctionsFile()

install.missing(packages=c("dplyr", "ggplot2", "reshape"), 
  biopackages=c("BiSeq", "GenomicRanges", "rtracklayer"))

sampleGroupsFile           =commandArgs(TRUE)[1]
bedGraphFolder           =commandArgs(TRUE)[2]
wg.out.file              =commandArgs(TRUE)[3]
tr.out.file              =commandArgs(TRUE)[4]
annotation_target_regions=commandArgs(TRUE)[5]

##################
#
# Check validity of 
# file paths
#
##################

# functions::checkPath() will quit if file not found
quit.if.not.exists(sampleGroupsFile, "Sample group file") 
quit.if.not.exists(bedGraphFolder, "Coverage file folder")
quit.if.not.exists(annotation_target_regions, "Genome annotation")

##################
#
# Read Sample information
#
##################


cat("Reading sample names file ...\n")
sampleGroups = read.csv(sampleGroupsFile, sep="\t", header=T, stringsAsFactors=F)
rownames(sampleGroups) = as.character(sampleGroups$Sample.Name)

write.samples = function(s) cat(s, "\n")
invisible(lapply(sampleGroups$Sample.Name, write.samples))

# Locate coverage files
cat("Locating cov files ...\n")
cov.file.list = list.files(path=bedGraphFolder,pattern=".cov",full.names = T)
cov.file.list = grep(paste0(sampleGroups$Sample.Name,collapse='|'),
    cov.file.list, value=TRUE)

invisible(lapply(cov.file.list, write.samples))

##################
#
# QC raw data targeted region
#
##################

# Load annotation
cat("QC raw data target regions:\n")
cat("\tLoading annotation file ... \n")
annotation_capture = import.bed(annotation_target_regions)
# Only use full pig chromosomes
annotation_capture = annotation_capture[seqnames(annotation_capture) %in% c(1:18, "X", "Y")]

# Read bismark files
cat("Reading cov files ...\n")
methylDataRaw = readBismark(cov.file.list, sampleGroups)

cat("\tSubset target regions from methyl data ... \n")
methylDataRaw.rk = subsetByOverlaps(methylDataRaw, annotation_capture)

cat("\tCalculating total coverage ... \n")
totalCoverage = colSums(totalReads(methylDataRaw.rk))

cat("\tCalculating number of sites covered ...\n")
number_of_sites_covered = apply(totalReads(methylDataRaw.rk),2,function(x){table(x>0)['TRUE']})
melt.data = melt(totalReads(methylDataRaw))

# cat("\tCalculating distribution of coverage per site ...\n")
coverage_per_site_distribution=function(x){
  zero=table(x==0)['TRUE']
  one_ten=table(x>=1 & x<=10)['TRUE']
  ten_onehund=table(x>=11 & x<=100)['TRUE']
  onehund_onethous=table(x>=101 & x<=1000)['TRUE']
  onethous_max=table(x>=1001 & x<=max(x))['TRUE']
  max=max(x)
  vec=c(zero,one_ten,ten_onehund,onehund_onethous,onethous_max)
        names(vec)=c('zero','one_ten','ten_onehund','onehund_onethous','onethous_max')
  vec
}

cat("\tCalculating distribution of coverage per site ... \n")
cov_table = apply(totalReads(methylDataRaw.rk), 2, coverage_per_site_distribution)
cov_table[is.na(cov_table)] = 0
cov_table.melt = melt(cov_table)
cov_table.melt$X1 = factor( cov_table.melt$X1, levels=unique(cov_table.melt$X1))
cov_table.melt$X2 = as.character(cov_table.melt$X2)
colnames(cov_table.melt)[1:2] = c('Coverage','Sample')
cov_table.melt$percent = cov_table.melt$value/number_of_sites_covered
cov_table.melt$Region = 'Raw data targeted regions'
cov_table.melt=cov_table.melt %>% filter(Coverage!='zero')

cat("\tPlotting data ...\n")


# Plot the coverage as a line chart showing percentage
#
# data - the coverage table
plotCoverage = function(data){
    p1 = ggplot(data, aes(x=Coverage,y=value/number_of_sites_covered), group=Sample, color=Sample) + 
        geom_line(aes(group=Sample, color=Sample)) + 
        scale_y_continuous(labels=scales::percent) + 
        theme_bw() + 
        scale_x_discrete(labels=c('1-10X','10-100X','100-1000X','>1000X')) + 
        ylab('') + 
        xlab('Coverage per site')
    ggsave(filename=tr.out.file, plot=p1, device="pdf", height=7, width=7)
}

  plotCoverage(cov_table.melt)

cat("Done\n")
quit(save="no", status=0)
