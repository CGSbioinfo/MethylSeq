library(ggplot2)
library(reshape)
library(dplyr)

# Read the input arguments
in_dir      =commandArgs(TRUE)[1]
sample_names=commandArgs(TRUE)[2]
pattern     =commandArgs(TRUE)[3]
outfile     =commandArgs(TRUE)[4]

dir.create(file.path("Report/figure/", "methExtractQC/"))

allFiles=list.files(in_dir, pattern=pattern)
sample_names=as.character(read.csv(sample_names, header=FALSE)[,1])

# Start matrix
all_data=matrix(ncol=8)

# Get colnames from the first sample
in_dir = gsub("/$","", in_dir)
toRead = paste0(in_dir,'/',grep(sample_names[2], allFiles, value=TRUE))

# Error occurs when deduplicated and non-deduplicated samples are present in the same folder.
# For now, just pick the first of the samples.
infile = toRead[1]
if (!file.exists(infile)){
	cat("File does not exist: ", infile, "\n")
	stop("Error, stopping")
} else {
	cat("Input file exists:\n'", infile, "'\n")
}

cat("Reading file:\n'", infile, "'\n")
x = readLines(infile, n = -1)
# x=readLines(con=infile, n=-1L, ok=TRUE, warn=TRUE, encoding = "unknown", skipNul = FALSE )
title=grep("^=",x)-1
end=which(x=='')
temp=x[(title[1]+2):(end[1]-1)]
coln=gsub(' ','_',unlist(strsplit(temp,'\t')[1]))
coln=gsub('%','pct',coln)
colnames(all_data)=c(coln,'sample', 'class', 'read')
cat("Successfully got column names\n")

# Get the value for each sample
for (i in 1:length(sample_names)) { 
	cat("Reading sample ",sample_names[i], "\n")
	infile = paste0(in_dir,'/',grep(sample_names[i], allFiles, value=TRUE))[1]
	cat("Reading file:\n'", infile, "'\n")
	if (!file.exists(infile)){
		cat("File does not exist: ", infile, "\n")
		stop("Error, stopping")
	}
	x=readLines(infile)
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

p1=all_data %>% 
	filter(class=='CpG') %>% 
	ggplot(aes(x=position,y=pct_methylation, colour=sample, group=sample, shape=sample)) + 
	geom_line(size=0.1) + geom_point(size=0.4) + 
	facet_wrap(~read) + 
	theme(axis.text.x=element_text(angle=90)) + 
	scale_shape_manual(values=1:length(sample_names)) + 
	ylab('% Methylation') + 
	scale_x_discrete(breaks=as.character(c(1,seq(from=0,to=76,by=5)[2:length(seq(from=0,to=76,by=5))]))) + 
	geom_vline(xintercept=c(1:10,65:76), colour='grey', linetype='dotted') + 
	ylim(0,70)
ggsave(p1,filename=outfile, width=14)
