library(ggplot2)
library(reshape)
library(dplyr)

in_dir=commandArgs(TRUE)[1]
sample_names=commandArgs(TRUE)[2]
pattern=commandArgs(TRUE)[3]
outdir=commandArs(TRUE)[4]

allFiles=list.files(in_dir, recursive=TRUE, pattern=pattern)
sample_names=as.character(read.csv(sample_names, header=FALSE)[,1])

# Start matrix
all_data=matrix(ncol=8)

# Get colnames
x=readLines(grep(sample_names[1], allFiles, value=TRUE))
title=grep("^=",x)-1
end=which(x=='')
temp=x[(title[1]+2):(end[1]-1)]
coln=gsub(' ','_',unlist(strsplit(temp,'\t')[1]))
coln=gsub('%','pct',coln)
colnames(all_data)=c(coln,'sample', 'class', 'read')

for (i in 1:length(sample_names)) { 
	x=readLines(grep(sample_names[i], allFiles, value=TRUE))
	title=grep("^=",x)-1
	end=which(x=='')
	for (j in 1:length(title)){
		name=gsub(" ","_",x[title[j]])
		name=gsub(paste0(c('\\(','\\)'), collapse='|'),"",name)
		read=gsub(".*_","",name)
		class=gsub("_context_R\\d$","",name)
		temp=x[(title[j]+2):(end[j]-1)]
		coln=gsub(' ','_',unlist(strsplit(temp,'\t')[1]))
		coln=gsub('%','pct',coln)
		df=as.data.frame(strsplit(temp,'\t')[-1])
		df=as.data.frame(apply(df,1,as.numeric))
		colnames(df)=coln
		df$sample=sample_names[i]
		df$class=class
		df$read=read
		all_data=rbind(all_data,df)
		}
	}
all_data=all_data[-1,]
all_data$position=factor(all_data$position, levels=unique(all_data$position))

all_data %>% filter(read=='R1')  %>% ggplot(aes(x=position,y=pct_methylation, colour=sample, group=sample, shape=sample)) + geom_line() + geom_point() + facet_wrap(~class) + theme(axis.text.x=element_text(angle=90)) + scale_shape_manual(values=1:length(sample_names)) + ylab('% Methylation') + scale_x_discrete(breaks=as.character(c(1,seq(from=0,to=76,by=5)[2:length(seq(from=0,to=76,by=5))])))

