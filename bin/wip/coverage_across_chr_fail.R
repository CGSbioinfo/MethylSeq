###################################################
###################################################
#      Coverage plots
###################################################
###################################################


library(GenomicRanges)
library(rtracklayer)
library(BiSeq)
library(ggplot2)
library(reshape)
library(dplyr)
library(caTools)

########################
# Arguments            #
########################

option_list = list(make_option(c("-s","--sampleInfoFile"), type="character", default=NULL, help="Sample information file", metavar="character"),
                   make_option(c("-c", "--comparisonFile"), type="character", default=NULL, help="Comparisons file", metavar = "character"),
                   make_option(c("-r", "--bedGraphDataFolder"), type="character", default=NULL, help="Folder containing the raw data (cels file)", metavar="character")
                   make_option(c("-b", "--bedFile"), type="character", default=NULL, help="Pah to annotation file", metavar="character")
)

opt_parser=OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (length(opt)<4){
  print_help(opt_parser)
  stop("All arguments must be supplied", call.=FALSE)
}

sampleInfoFile = opt$sampleInfoFile
groupInfoFile = opt$comparisonFile
bedGraphFolder = opt$bedGraphDataFolder
annotationBedPath = opt$bedFile

########################
#      Data loading    #
########################

sampleInfo<- read.table(sampleInfoFile, sep="\t", header=T)
listMethylFileName<-list.files(path=bedGraphFolder,pattern=".cov.gz$")
methylSampleName<-strsplit(listMethylFileName, split="_R1.*")
first=T
for(i in (1:length(listMethylFileName))){
  if(first){
    first=F
    sampleNameOrder<-methylSampleName[[i]][1]
  }else{
    sampleNameOrder<-c(sampleNameOrder,methylSampleName[[i]][1])
  }
}
sampleInfo<-sampleInfo[match(sampleNameOrder, sampleInfo$Sample.Name),]
listMethylFile<-list.files(path=bedGraphFolder,pattern=".cov.gz$",full.names = T)
listMethylFileName<-list.files(path=bedGraphFolder, pattern=".cov.gz$")
sample_names_table<-strsplit(listMethylFileName, split="_R1.*")
first=T
for(i in 1:length(listMethylFileName)){
  if(first){
    first=F
    sample_name<-sample_names_table[[i]][1]
  }else{
    sample_name<-c(sample_name, sample_names_table[[i]][1])
  }
}
sampleInfo<-sampleInfo[match(sample_name, sampleInfo$Sample.Name),]
methylDataRaw<-readBismark(listMethylFile, as.character(sampleInfo$Sample.Name))
genomic_ranges_methyldata<-as.data.frame(rowRanges(methylDataRaw))
genomic_ranges_coverage<-cbind(genomic_ranges_methyldata, totalReads(methylDataRaw))

##########################
#      Genome loading    #
##########################
x=read.table('QC/genomeCoverage/temp_header_CQ-hero-G-0277Wcordblood-group1_S4')
granges_genome=data.frame(seqnames=x$V1[1], position=1:x$V2[1])
granges_genome=GRanges(seqnames=granges_genome$seqnames,ranges=IRanges(start=granges_genome$position,end=granges_genome$position))

# sample granges
granges_samplecov=GRanges(seqnames=genomic_ranges_coverage$seqnames, ranges=IRanges(start=genomic_ranges_coverage$start, end=genomic_ranges_coverage$end), cov=genomic_ranges_coverage$`CQ-hero-G-0277Wcordblood-group1_S4`)

test=findOverlaps(granges_samplecov,granges_genome)

# sample and genome overlap
granges_genome=as.data.frame(granges_genome)
granges_genome$cov=0
granges_genome$cov[subjectHits(test)]=granges_samplecov[queryHits(test)]$cov

# sliding window mean
slideFunct=function(data,window,step){
	total=length(data)
	spots=seq(from=1, to=(total-window), by=step)
	result=sapply(spots,function(x){mean(data[x:x+window])})
	return(result)
}
timestamp()
test1_sw=slideFunct(granges_genome$cov, 100000, 1000)
timestamp()


sub=test1_sw_density$y
#plot(1:length(sub), sub)
sub.melt=melt(sub)
sub.melt$position=1:nrow(sub.melt)
ggplot(sub.melt,aes(x=position,y=value)) + geom_line()







slideFunct=function(data,window,step){
	total=length(data)
	spots=seq(from=1, to=(total-window), by=step)
	result=vector(length=length(spots))
	for (i in 1:length(spots)){
		result[i]=mean(data[spots[i]:(spots[i]+window)])
	}
	return(result)
}


test1_sw=slideFunct(granges_genome$cov, 10000, 10000)
spots=seq(from=1, to=(nrow(granges_genome)), by=100)
spots=rep(c(spots),each=100)

test_x=sapply(spots,function(x){mean(data[x:x+window])})

spots=rep(c(1:(248956422/10000000)),each=24)
spots=c(spots,rep((248956422/24)+1,248956422-length(spots)))
granges_genome$spots=spots
granges_genome_windows_meancov=granges_genome %>% group_by(spots) %>% summarise(mean_window=mean(cov))

	
	

