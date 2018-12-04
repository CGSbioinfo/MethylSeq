suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(parallel))
suppressMessages(library(data.table))

## Some functions
read.data=function(sample,opt){
	allfiles=list.files(opt$in_dir)
	allfiles=grep(opt$pattern,allfiles,value=TRUE)
	file=grep(paste0(c(sample),collapse='|'), allfiles, value=TRUE)
	out=as.data.frame(fread(file))
	out$feature=paste(out$V1,out$V2,out$V3,sep='_')
	out=out[,-c(1:3)]
	out$sample=gsub(opt$pattern,'',gsub('.*\\/','',sample))
	return(out)
	}

# Parse arguments
option_list=list(
	make_option(c("--one_sample_file"), type="character", action='store_true', default=FALSE),
	make_option(c("--sample_names_file"), type="character", action='store_true', default=FALSE),
	make_option(c("--in_dir"), type="character", default='./'),
	make_option(c("--out_dir"), type="character", default='./'),
	make_option(c("--pattern"), type="character", default='_coverage.txt'))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Read data
if (opt$one_sample != FALSE){
	file=opt$one_sample_file
	file=read.table(file)
	file$feature=paste(file$V1,file$V2,file$V3,sep='_')
	file=file[,-c(1:3)]
	file$sample=gsub(opt$pattern,'',gsub('.*\\/','',opt$one_sample_file))
} else if (opt$sample_names_file != FALSE) {
	sample_names=read.table(opt$sample_names_file, stringsAsFactors=FALSE)[,1]
	#no_cores=detectCores()-1
	print(Sys.time())
	cl=makeCluster(6)
	clusterExport(cl,list('read.data','opt'))
	clusterEvalQ(cl, library(data.table))
	file=parLapply(cl,sample_names,function(x){read.data(sample=x,opt=opt)})
	stopCluster(cl)
	file=do.call(rbind,file)
	print(Sys.time())
	}
}
colnames(file)[1:4]=c('depth','number_of_positions_with_depth','feature_size','fraction')

# Calculate feature depth
file$depth_total=file$depth*file$number_of_positions_with_depth
feature_depth=file %>% group_by(sample,feature) %>% summarise(total=sum(depth_total))

# Normalize by size
feature_size=file%>% group_by(sample,feature) %>% summarise(size=unique(feature_size))
feature_depth=left_join(feature_depth,feature_size,by=c('sample','feature'))
feature_depth$normalised_cov=feature_depth$total/feature_depth$size

feature_depth_plot=feature_depth %>% ggplot(aes(x=log2(normalised_cov+.1), fill=sample, color=sample, group=sample)) + 
	geom_density(alpha=0.3) + 
	scale_x_continuous(breaks=log2(c(0.1,mean(feature_depth$normalised_cov),max(feature_depth$normalised_cov))), labels=paste(round(c(log2(c(0.1,round(mean(feature_depth$normalised_cov),digits=2),round(max(feature_depth$normalised_cov), digits=2)))),digits=2),c(0,round(mean(feature_depth$normalised_cov),digits=2),round(max(feature_depth$normalised_cov), digits=2)), sep='\n')) + 
	theme(panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.3,linetype=2), panel.grid.minor=element_line(colour='grey',size=.3,linetype=2), axis.text.x=element_text(angle=90,vjust=0.5)) + xlab('Log2 - Feature-length-Normalised coverage')
ggsave(feature_depth_plot,filename=paste0(opt$out_dir,'/', 'Normalised_feature_coverage.pdf'), device='pdf', width=15 )


# Calculate zscore
feature_depth=feature_depth%>% group_by(sample) %>% mutate(mean=mean(normalised_cov),sd=sd(normalised_cov), zscore=scale(normalised_cov))

feature_depth_plot_zscore=feature_depth %>% ggplot(aes(x=zscore,fill=sample, color=sample, group=sample)) + 
	geom_histogram(bins=500,alpha=0.3) + 
	scale_x_continuous(breaks=c(0,5,7,10,13,15,17,20,max(feature_depth$zscore)), labels=c(0,5,7,10,13,15,17,20,round(max(feature_depth$zscore),digits=2)))
ggsave(feature_depth_plot_zscore,filename=paste0(opt$out_dir,'/', 'Zscore_feature_coverage.pdf'), device='pdf', width=15 )

vals=c(3,4,5,6,7)
vals_matrix=matrix(ncol=length(vals),nrow=length(sample_names))
rownames(vals_matrix)=sample_names
colnames(vals_matrix)=vals

for (i in 1:length(vals)){
	temp=feature_depth %>% filter(zscore>=vals[i]) %>% group_by(sample) %>% summarise(nfeatures=length(feature))
	vals_matrix[,as.character(vals[i])]=temp[which(temp$sample %in% rownames(vals_matrix)),]$nfeatures
}
write.csv(vals_matrix,'Number_of_features_removed_zscore_threshold.csv')

# Output list

