#!/usr/local/bin/Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(require(grid))
suppressMessages(library(xtable))

in_dir       = commandArgs(TRUE)[1]
sample.file  = commandArgs(TRUE)[2]
readType     = commandArgs(TRUE)[3]
outdir       = commandArgs(TRUE)[4]
suffix       = commandArgs(TRUE)[5]
plot_device  = commandArgs(TRUE)[6]
if (is.na(suffix)){
  suffix=""
}

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

pathExistsOrQuit(in_dir, "Sample input folder")
pathExistsOrQuit(sample.file, "Sample names file")
pathExistsOrQuit(outdir, "Output directory")


files        = list.files(path = in_dir, full.names = TRUE, recursive = TRUE)
sample.names = read.table(sample.file)[,1]


# Per sequence quality scores
#----------------------------
qual_scores = files[grep('per_sequence_quality_scores.txt',files)]


#' Read the quality scores for a sample
#' @title Read quality scores
#' @param sample the sample name to read
#' @return a table of scores
#' @export
getScores = function(sample){
  x = qual_scores[grep(sample,qual_scores)]
  x = x[grep('_R1_',x)]
  x = read.table(x, stringsAsFactors = FALSE)
  colnames(x) = c('x','y')
  x$Sample = sample
  # mr1 = rbind(mr1,x)
  x
}

score.list = lapply(sample.names, getScores)
mr1 = do.call("rbind", score.list)

# mr1=matrix(ncol=3)
# colnames(mr1)=c('x','y','Sample')
# for (i in 1:length(sample.names)){
#   x=qual_scores[grep(sample.names[i],qual_scores)]
#   x=x[grep('_R1_',x)]
#   x=read.table(x, stringsAsFactors = FALSE)
#   colnames(x)=c('x','y')
#   x$Sample=sample.names[i]
#   mr1=rbind(mr1,x)
# }

#' Create a quality score plot for paired end reads
#' @title Create a paired end read quality plot
#' @return a \code{\link{ggplot}} object
#' @export
makePairedEndPlot = function(mr1){
  mr2=matrix(ncol=3)
  colnames(mr2)=c('x','y','Sample')
  
  for (i in 1:length(sample.names)){
    y= qual_scores[grep(sample.names[i],qual_scores)]
    y= y[grep('_R2_',y)]
    y= read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('x','y')
    y$Sample=sample.names[i]
    mr2=rbind(mr2,y)
  }
  # mr2=mr2[-1,]
  mr2=cbind(mr2,Read='Read 2')
  
  d=rbind(mr1,mr2) 
  ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + 
    geom_line(size=0.3) + 
    geom_point(size=1) + 
    facet_wrap(~Read) +
    theme( axis.title.x =element_text(size=9), axis.title.y =element_text(size=6), 
           axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), legend.text=element_text(size=8),  
           legend.key.height=unit(.8,"line") ) + 
    # axis.title.y=element_blank()
    ylab("") + 
    xlab('Mean Quality score') + 
    scale_shape_manual(values=1:length(unique(d$Sample))) 
}



# mr1=mr1[-1,]
mr1=cbind(mr1,Read='Read 1')

if (readType=='pairedEnd') {
  p = makePairedEndPlot(mr1)
  # mr2=matrix(ncol=3)
  # colnames(mr2)=c('x','y','Sample')
  # for (i in 1:length(sample.names)){
  #   y=qual_scores[grep(sample.names[i],qual_scores)]
  #   y=y[grep('_R2_',y)]
  #   y=read.table(y, stringsAsFactors = FALSE)
  #   colnames(y)=c('x','y')
  #   y$Sample=sample.names[i]
  #   mr2=rbind(mr2,y)
  # }
  # mr2=mr2[-1,]
  # mr2=cbind(mr2,Read='Read 2')
  
  # d=rbind(mr1,mr2)
  # p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=0.3) + geom_point(size=1) + facet_wrap(~Read) +
  #   theme( axis.title.x =element_text(size=9), axis.title.y =element_text(size=6), 
  #          axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), legend.text=element_text(size=8),  
  #          legend.key.height=unit(.8,"line"), axis.title.y=element_blank()) + ylab("") + xlab('Mean Quality score') + scale_shape_manual(values=1:length(unique(d$Sample))) 
} else {
  d=mr1
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=0.3) + geom_point(size=1) + facet_wrap(~Read) +
    theme(axis.title.y =element_text(size=16), 
           axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=7),  
           legend.key.height=unit(.8,"line"), axis.title.y=element_blank()) + ylab("") + xlab('Mean Quality score') + scale_shape_manual(values=1:length(unique(d$Sample)))
}

ggsave(filename=paste0(outdir,'/per_sequence_quality_scores', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)


# Per sequence gc content
#----------------------------
gc=files[grep('per_sequence_gc_content.txt',files)]
mr1=matrix(ncol=3)
colnames(mr1)=c('x','y','Sample')
for (i in 1:length(sample.names)){
  x=gc[grep(sample.names[i],gc)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('x','y')
  x$Sample=sample.names[i]
  mr1=rbind(mr1,x)
}
mr1=cbind(mr1,Read='Read 1')
mr1=mr1[-1,]


if (readType=='pairedEnd') {
  mr2=matrix(ncol=3)
  colnames(mr2)=c('x','y','Sample')
  for (i in 1:length(sample.names)){
    y=gc[grep(sample.names[i],gc)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('x','y')
    y$Sample=sample.names[i]
    mr2=rbind(mr2,y)
  }
  mr2=cbind(mr2,Read='Read 2')
  mr2=mr2[-1,]
  
  d=rbind(mr1,mr2)
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) +  xlab('%GC') + ylab('') + theme(legend.key.height=unit(.8,"line")) + scale_shape_manual(values=1:length(unique(d$Sample)))

} else {
  d=mr1
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) + 
    xlab('%GC') + ylab('') + theme(legend.key.height=unit(.8,"line")) +  scale_shape_manual(values=1:length(unique(d$Sample)))

}
ggsave(filename=paste0(outdir,'/per_sequence_gc_content', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)


# Per sequence length distribution
#---------------------------------
length_dist=files[grep('seq_length.txt',files)]
mr1=matrix(ncol=4)
colnames(mr1)=c('Position','Frequency','Sample', 'Read')
for (i in 1:length(sample.names)){
  x=length_dist[grep(sample.names[i],length_dist)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Position','Frequency')
  x$Position <- factor(x$Position, as.character(x$Position))
  x$Read='Read 1'
  x$Sample=sample.names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
if (readType=='pairedEnd') {
  mr2=matrix(ncol=4)
  colnames(mr2)=c('Position','Frequency','Sample', 'Read')
  for (i in 1:length(sample.names)){
    y=length_dist[grep(sample.names[i],length_dist)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('Position','Frequency')
    y$Position <- factor(y$Position, as.character(y$Position))
    y$Read='Read 2'
    y$Sample=sample.names[i]
    mr2=rbind(mr2,y)
  }
  mr2=mr2[-1,]
  d=rbind(mr1,mr2)
  p=ggplot(d, aes(x = Position, y = Frequency, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) +
    theme( axis.title.x =element_text(size=12), axis.title.y =element_text(size=6), 
           axis.text.x=element_text(size=7,angle=90), axis.text.y=element_text(size=7))  + ylab("") + xlab('Length') + scale_shape_manual(values=1:length(unique(d$Sample)))
} else {
  d=mr1
  p=ggplot(d, aes(x = Position, y = Frequency, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) +
    scale_shape_manual(values=1:length(unique(d$Sample))) + 
    theme( axis.title.x =element_text(size=12), axis.title.y =element_text(size=12), legend.key.height=unit(.8,"line"), 
           axis.text.x=element_text(size=7,angle=90), axis.text.y=element_text(size=12))  + ylab("") + xlab('length')
}
ggsave(filename=paste0(outdir,'/sequence_length_distribution', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)


# Per duplication levels
#-----------------------
dup_levels=files[grep('seq_dup_levels.txt',files)]
mr1=matrix(ncol=4)
colnames(mr1)=c('Duplication_Level','Percentage', 'Sample','Read')
for (i in 1:length(sample.names)){
  x=dup_levels[grep(sample.names[i],dup_levels)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Duplication_Level','Percentage of Deduplicated','Percentage')
  x$Duplication_Level <- factor(x$Duplication_Level, as.character(x$Duplication_Level))
  x=x[,-which(colnames(x)=="Percentage of Deduplicated")]
  x$Read='Read 1'
  x$Sample=sample.names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
if (readType=='pairedEnd') {
  mr2=matrix(ncol=4)
  colnames(mr2)=c('Duplication_Level','Percentage', 'Sample','Read')
  for (i in 1:length(sample.names)){
    y=dup_levels[grep(sample.names[i],dup_levels)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('Duplication_Level','Percentage of Deduplicated','Percentage')
    y$Duplication_Level <- factor(y$Duplication_Level, as.character(y$Duplication_Level))
    y=y[,-which(colnames(y)=="Percentage of Deduplicated")]
    y$Read='Read 2'
    y$Sample=sample.names[i]
    mr2=rbind(mr2,y)
  }
  mr2=mr2[-1,]
  d=rbind(mr1,mr2)
  d$Duplication_Level=factor(d$Duplication_Level,levels=unique(d$Duplication_Level))
  p=ggplot(d, aes(x=Duplication_Level, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) + scale_shape_manual(values=1:length(unique(d$Sample))) +
    theme(axis.text.x = element_text(size = 9, angle=90), legend.key.height=unit(.8,"line"),
                        axis.title.y=element_blank())  + xlab("Number of copies per read")
} else {
  d=mr1
  d$Duplication_Level=factor(d$Duplication_Level,levels=unique(d$Duplication_Level))
  p=ggplot(d, aes(x=Duplication_Level, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.3) + geom_point(size=1) + facet_wrap(~Read) + 
    theme(axis.text.x = element_text(size = 9, angle=90), legend.key.height=unit(.8,"line"),
          axis.title.y=element_blank())  + xlab("Number of copies per read") + scale_shape_manual(values=1:length(unique(d$Sample)))
}
ggsave(filename=paste0(outdir,'/sequence_dup_levels', suffix, '.', plot_device), width=10, height=3.5, units='in', plot=p)



# Base content
#-----------------------
base_content=files[grep('per_base_sequence_content.txt',files)]
mr1=matrix(ncol=5)
colnames(mr1)=c('Base', 'Percentage', 'Nucleotide','Sample','Read')
for (i in 1:length(sample.names)){
  x=base_content[grep(sample.names[i],base_content)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Base','G','A', 'T', 'C')
  all=matrix(ncol=3)
  colnames(all)=c('Base','Percentage','Nucleotide')
  nt=c('G','A','T','C')
  for (j in 1:length(nt)){
    all=rbind(all,cbind(x[,1], 'Percentage'=x[,which(colnames(x)==nt[j])], 'Nucleotide'=nt[j]))
  }
  all=as.data.frame(all[-1,])
  all$Percentage=as.numeric(as.character(all$Percentage))
  all$Read='Read 1'
  all$Sample=sample.names[i]
  mr1=rbind(mr1,all)
}
mr1=mr1[-1,]
mr1$Base=factor(mr1$Base, levels=unique(mr1$Base))
p=ggplot(mr1, aes(x=Base, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.2) + geom_point(size=0.8) + facet_wrap(~Nucleotide) + 
    theme(axis.text.x = element_text(size = 7, angle=90, v=0.5), axis.text.y = element_text(size = 7), legend.text=element_text(size=6), legend.title=element_text(size=7),
          axis.title=element_text(size=8), legend.key.height=unit(.8,"line"), strip.text=element_text(size=7, face='bold'), panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.2, linetype=2), panel.grid.minor=element_line(colour='grey',size=.2,linetype=2)) + 
    xlab("Position") + scale_shape_manual(values=1:length(unique(mr1$Sample))) + scale_x_discrete(breaks = c('1','25','50','75','100'))
ggsave(filename=paste0(outdir,'/per_base_sequence_content_r1', suffix, '.', plot_device), width=8, height=3.5, units='in', plot=p)

if (readType=='pairedEnd') {
    mr2=matrix(ncol=5)
    colnames(mr2)=c('Base', 'Percentage', 'Nucleotide','Sample','Read')
    for (i in 1:length(sample.names)){
      x=base_content[grep(sample.names[i],base_content)]
      x=x[grep('_R2_',x)]
      x=read.table(x, stringsAsFactors = FALSE)
      colnames(x)=c('Base','G','A', 'T', 'C')
      all=matrix(ncol=3)
      colnames(all)=c('Base','Percentage','Nucleotide')
      nt=c('G','A','T','C')
      for (j in 1:length(nt)){
        all=rbind(all,cbind(x[,1], 'Percentage'=x[,which(colnames(x)==nt[j])], 'Nucleotide'=nt[j]))
      }
      all=as.data.frame(all[-1,])
      all$Percentage=as.numeric(as.character(all$Percentage))
      all$Read='Read 2'
      all$Sample=sample.names[i]
      mr2=rbind(mr2,all)
    }
mr2=mr2[-1,]
mr2$Base=factor(mr2$Base, levels=unique(mr2$Base))
p=ggplot(mr2, aes(x=Base, y=Percentage, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=.2) + geom_point(size=0.8) + facet_wrap(~Nucleotide) + 
    theme(axis.text.x = element_text(size = 7, angle=90, v=0.5), axis.text.y = element_text(size = 7), legend.text=element_text(size=6), legend.title=element_text(size=7),
          axis.title=element_text(size=8), legend.key.height=unit(.8,"line"), strip.text=element_text(size=7, face='bold'), panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.2, linetype=2), panel.grid.minor=element_line(colour='grey',size=.2,linetype=2)) + 
    xlab("Position") + scale_shape_manual(values=1:length(unique(mr2$Sample))) + scale_x_discrete(breaks = c('1','25','50','75','100'))

}
ggsave(filename=paste0(outdir,'/per_base_sequence_content_r2', suffix, '.', plot_device), width=8, height=3.5, units='in', plot=p)

quit(save="no", status=0)