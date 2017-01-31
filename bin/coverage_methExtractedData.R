#!/usr/bin/R

print("Loading R packages ...")
suppressMessages(library(BiSeq))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))

sampleInfoFile=commandArgs(TRUE)[1]
bedGraphFolder=commandArgs(TRUE)[2]
out_wholeGenome=commandArgs(TRUE)[3]
out_targetRegion=commandArgs(TRUE)[4]
annotation_target_regions=commandArgs(TRUE)[5]
# /mnt/research/jb393/MethylSeq_Pilot/Kit_annotation/S03770311_Covered_GR38.bed

# Read Sample information file
print("Reading Sample Info ...")
sampleInfo<- read.table(sampleInfoFile, sep="\t", header=T)
rownames(sampleInfo)=as.character(sampleInfo$Sample.Name)

# Locate coverage files
print("Locating cov files ... ")
listMethylFile<-list.files(path=bedGraphFolder,pattern=".cov",full.names = T)
listMethylFile=grep(paste0(c(as.character(sampleInfo$Sample.Name)),collapse='|'),listMethylFile, value=TRUE)

# Read bismark files
print("Reading cov files ...")
methylDataRaw<-readBismark(listMethylFile, sampleInfo)

# QC raw data whole genome
print("QC raw data whole genome:")
print("Calculating total coverage ...")
totalCoverage=colSums(totalReads(methylDataRaw))
print("Calculating number of sites covered ...")
number_of_sites_covered=apply(totalReads(methylDataRaw),2,function(x){table(x>0)['TRUE']})
melt.data=melt(totalReads(methylDataRaw))

print("Calculating distribution of coverage per site ...")
coverage_per_site_distribution=function(x){
	zero=table(x==0)['TRUE']
	one_ten=table(x>=1 & x<=10)['TRUE']
	ten_onehund=table(x>=11 & x<=100)['TRUE']
	onehund_onethous=table(x>=101 & x<=1000)['TRUE']
	onethous_max=table(x>=1001 & x<=max(x))['TRUE']
	max=max(x)
	vec=c(zero,one_ten,ten_onehund,onehund_onethous,onethous_max)
        names(vec)=c('zero','one_ten','ten_onehund','onehund_onethous','onethous_max')
	return(vec)
}
coverage_per_site_distribution_table=apply(totalReads(methylDataRaw), 2, coverage_per_site_distribution)
coverage_per_site_distribution_table[is.na(coverage_per_site_distribution_table)]=0
coverage_per_site_distribution_table.melt=melt(coverage_per_site_distribution_table)
coverage_per_site_distribution_table.melt$X1=factor( coverage_per_site_distribution_table.melt$X1, levels=unique(coverage_per_site_distribution_table.melt$X1))
coverage_per_site_distribution_table.melt$X2=as.character(coverage_per_site_distribution_table.melt$X2)
colnames(coverage_per_site_distribution_table.melt)[1:2]=c('Coverage','Sample')
coverage_per_site_distribution_table.melt$percent=coverage_per_site_distribution_table.melt$value/number_of_sites_covered
coverage_per_site_distribution_table.melt$Region='Raw data whole genome'
coverage_per_site_distribution_table.melt=coverage_per_site_distribution_table.melt %>% filter(Coverage!='zero')
print("Plotting data ...")
p1=ggplot(coverage_per_site_distribution_table.melt, aes(x=Coverage,y=value/number_of_sites_covered), group=Sample, color=Sample) + geom_line(aes(group=Sample, color=Sample)) + scale_y_continuous(labels=scales::percent) + theme_bw() + scale_x_discrete(labels=c('1-10X','10-100X','100-1000X','>1000X')) + ylab('') + xlab('Coverage per site')
ggsave(filename=out_wholeGenome,plot=p1,device='pdf', height=7,width=7)

# QC raw data targeted region
# Load annotation
print("QC raw data target regions:")
print("Loading annotation file ... ")
annotation_capture<-import.bed(annotation_target_regions)
print("Subset target regions from methyl data ... ")
methylDataRaw.rk<-subsetByOverlaps(methylDataRaw, annotation_capture)
print("Calculating total coverage ... ")
totalCoverage=colSums(totalReads(methylDataRaw.rk))
print("Calculating number of sites covered ...")
number_of_sites_covered=apply(totalReads(methylDataRaw.rk),2,function(x){table(x>0)['TRUE']})
melt.data=melt(totalReads(methylDataRaw))
print("Calculating distribution of coverage per site ... ")
coverage_per_site_distribution_table=apply(totalReads(methylDataRaw.rk), 2, coverage_per_site_distribution)
coverage_per_site_distribution_table[is.na(coverage_per_site_distribution_table)]=0
coverage_per_site_distribution_table.melt=melt(coverage_per_site_distribution_table)
coverage_per_site_distribution_table.melt$X1=factor( coverage_per_site_distribution_table.melt$X1, levels=unique(coverage_per_site_distribution_table.melt$X1))
coverage_per_site_distribution_table.melt$X2=as.character(coverage_per_site_distribution_table.melt$X2)
colnames(coverage_per_site_distribution_table.melt)[1:2]=c('Coverage','Sample')
coverage_per_site_distribution_table.melt$percent=coverage_per_site_distribution_table.melt$value/number_of_sites_covered
coverage_per_site_distribution_table.melt$Region='Raw data targeted regions'
coverage_per_site_distribution_table.melt=coverage_per_site_distribution_table.melt %>% filter(Coverage!='zero')
print("Plotting data ...")
p1=ggplot(coverage_per_site_distribution_table.melt, aes(x=Coverage,y=value/number_of_sites_covered), group=Sample, color=Sample) + geom_line(aes(group=Sample, color=Sample)) + scale_y_continuous(labels=scales::percent) + theme_bw() + scale_x_discrete(labels=c('1-10X','10-100X','100-1000X','>1000X')) + ylab('') + xlab('Coverage per site')
ggsave(filename=out_targetRegion, plot=p1, device="pdf", height=7, width=7)

q()

# QC raw data, methylation proportion and coverage at target regions. Not finished, Postponed for a different script
methReadsRaw.rk=methReads(methylDataRaw.rk)
totalReadsRaw.rk=totalReads(methylDataRaw.rk)

methProportion=melt(methReadsRaw.rk/totalReadsRaw.rk)
totalReadsRaw.rk=melt(totalReadsRaw.rk)

table=full_join(totalReadsRaw.rk,methProportion,by=c('X2','X1'))
colnames(table)=c('site','Sample','Coverage','methProportion')
table$Sample=as.factor(table$Sample)

###############################################################################################################


png('differentialMethylation/female_coord/covBoxplots.png')
covBoxplots(methylDataRaw.rk)
dev.off()

#######################
#     Define clusters #
#######################
cluster.rk<-clusterSites(methylDataRaw.rk, perc.samples=.5, min.sites = 20, max.dist = 100, mc.cores=30)

##Filtering highly covered region to avoid bias due to high weights
ind.cov<-totalReads(cluster.rk)>0
ninety_quant<-quantile(totalReads(cluster.rk)[ind.cov],0.9)
cluster.rk.f<-limitCov(cluster.rk, maxCov = ninety_quant)

###Smoothing methylation (will take some time!)
predictedMeth<-predictMeth(cluster.rk.f, mc.cores=30)

placebo=predictedMeth[,sampleInfo$Treated=='Placebo']
treated=predictedMeth[,sampleInfo$Treated=='Treated']

mean.placebo=rowMeans(methLevel(placebo))
mean.treated=rowMeans(methLevel(treated))

png('differentialMethylation/female_coord/meth_placebo_vs_treated.png')
plot(mean.placebo,mean.treated)
dev.off()

group=sampleInfo$Treated
betaResults=betaRegression(formula=~group, link='probit', object=predictedMeth, type='BR', mc.cores=30)



