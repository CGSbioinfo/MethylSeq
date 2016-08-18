#!/usr/local/bin/Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(require(grid))
suppressMessages(library(xtable))

in_dir=commandArgs(TRUE)[1]
sample_names=commandArgs(TRUE)[2]
readType=commandArgs(TRUE)[3]
outdir=commandArgs(TRUE)[4]
suffix=commandArgs(TRUE)[5]
plot_device=commandArgs(TRUE)[6]
if (is.na(suffix)){
  suffix=''
}

files=list.files(path = in_dir, full.names = TRUE, recursive = TRUE)
sample_names=read.table(sample_names)[,1]


# Per sequence quality scores
#----------------------------
qual_scores=files[grep('per_sequence_quality_scores.txt',files)]
mr1=matrix(ncol=3)
colnames(mr1)=c('x','y','Sample')
for (i in 1:length(sample_names)){
  x=qual_scores[grep(sample_names[i],qual_scores)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('x','y')
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
mr1=cbind(mr1,Read='Read 1')

if (readType=='pairedEnd') {
  mr2=matrix(ncol=3)
  colnames(mr2)=c('x','y','Sample')
  for (i in 1:length(sample_names)){
    y=qual_scores[grep(sample_names[i],qual_scores)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('x','y')
    y$Sample=sample_names[i]
    mr2=rbind(mr2,y)
  }
  mr2=mr2[-1,]
  mr2=cbind(mr2,Read='Read 2')
  
  d=rbind(mr1,mr2)
  p=ggplot(d, aes(x = x, y = y, group=Sample, colour=Sample, shape=Sample)) + geom_line(size=0.3) + geom_point(size=1) + facet_wrap(~Read) +
    theme( axis.title.x =element_text(size=9), axis.title.y =element_text(size=6), 
           axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), legend.text=element_text(size=8),  
           legend.key.height=unit(.8,"line"), axis.title.y=element_blank()) + ylab("") + xlab('Mean Quality score') + scale_shape_manual(values=1:length(unique(d$Sample))) 
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
for (i in 1:length(sample_names)){
  x=gc[grep(sample_names[i],gc)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('x','y')
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=cbind(mr1,Read='Read 1')
mr1=mr1[-1,]


if (readType=='pairedEnd') {
  mr2=matrix(ncol=3)
  colnames(mr2)=c('x','y','Sample')
  for (i in 1:length(sample_names)){
    y=gc[grep(sample_names[i],gc)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('x','y')
    y$Sample=sample_names[i]
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
for (i in 1:length(sample_names)){
  x=length_dist[grep(sample_names[i],length_dist)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Position','Frequency')
  x$Position <- factor(x$Position, as.character(x$Position))
  x$Read='Read 1'
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
if (readType=='pairedEnd') {
  mr2=matrix(ncol=4)
  colnames(mr2)=c('Position','Frequency','Sample', 'Read')
  for (i in 1:length(sample_names)){
    y=length_dist[grep(sample_names[i],length_dist)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('Position','Frequency')
    y$Position <- factor(y$Position, as.character(y$Position))
    y$Read='Read 2'
    y$Sample=sample_names[i]
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
for (i in 1:length(sample_names)){
  x=dup_levels[grep(sample_names[i],dup_levels)]
  x=x[grep('_R1_',x)]
  x=read.table(x, stringsAsFactors = FALSE)
  colnames(x)=c('Duplication_Level','Percentage of Deduplicated','Percentage')
  x$Duplication_Level <- factor(x$Duplication_Level, as.character(x$Duplication_Level))
  x=x[,-which(colnames(x)=="Percentage of Deduplicated")]
  x$Read='Read 1'
  x$Sample=sample_names[i]
  mr1=rbind(mr1,x)
}
mr1=mr1[-1,]
if (readType=='pairedEnd') {
  mr2=matrix(ncol=4)
  colnames(mr2)=c('Duplication_Level','Percentage', 'Sample','Read')
  for (i in 1:length(sample_names)){
    y=dup_levels[grep(sample_names[i],dup_levels)]
    y=y[grep('_R2_',y)]
    y=read.table(y, stringsAsFactors = FALSE)
    colnames(y)=c('Duplication_Level','Percentage of Deduplicated','Percentage')
    y$Duplication_Level <- factor(y$Duplication_Level, as.character(y$Duplication_Level))
    y=y[,-which(colnames(y)=="Percentage of Deduplicated")]
    y$Read='Read 2'
    y$Sample=sample_names[i]
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

