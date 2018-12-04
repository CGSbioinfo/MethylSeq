suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(parallel))
suppressMessages(library(data.table))

# Check if we are provided with one sample or a text file with sample names
option_list=list(
	make_option(c("--one_sample_file"), type="character", action='store_true', default=FALSE),
	make_option(c("--sample_names_file"), type="character", action='store_true', default=FALSE),
	make_option(c("--in_dir"), type="character", default='./'),
	make_option(c("--pattern"), type="character", default='_coverage_per_position.txt'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# # Read data
# if (opt$one_sample != FALSE){
# 	file=opt$one_sample_file
# 	file=read.table(file)
# 	file$feature=paste(file$V1,file$V2,file$V3,sep='_')
# 	file=file[,-c(1:3)]
# 	file$sample=gsub(opt$pattern,'',gsub('.*\\/','',opt$one_sample_file))
# } else if (opt$sample_names_file != FALSE) {
# 	sample_names=read.table(opt$sample_names_file, stringsAsFactors=FALSE)[,1]
# 	allfiles=list.files(opt$in_dir)
# 	allfiles=grep(opt$pattern,allfiles,value=TRUE)
# 	allfiles=grep(paste0(c(sample_names),collapse='|'), allfiles, value=TRUE)
# 	file=read.table(paste0(opt$in_dir,'/',allfiles[1]))
# 	file$feature=paste(file$V1,file$V2,file$V3,sep='_')
# 	file=file[,-c(1:3)]
# 	file$sample=gsub(opt$pattern,'',gsub('.*\\/','',allfiles[1]))
# 	for (i in 2:length(sample_names)){
# 		temp=read.table(paste0(opt$in_dir,'/',allfiles[i]))
# 		temp$feature=paste(temp$V1,temp$V2,temp$V3,sep='_')
# 		temp=temp[,-c(1:3)]
# 		temp$sample=gsub(opt$pattern,'',gsub('.*\\/','',allfiles[i]))
# 		file=rbind(file,temp)
# 	}
# }



opt$sample_names='../../sample_names.txt'
sample_names=read.table(opt$sample_names_file, stringsAsFactors=FALSE)[,1]
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



no_cores=detectCores()-1
print(Sys.time())
cl=makeCluster(6)
clusterExport(cl,list('read.data','opt'))
clusterEvalQ(cl, library(data.table))
file=parLapply(cl,sample_names,function(x){read.data(sample=x,opt=opt)})
stopCluster(cl)
file=do.call(rbind,file)
colnames(file)[1:2]=c('Position','coverage')
print(Sys.time())
#res=lapply(res,function(x){x$feature=paste(x$V1,x$V2,x$V3,sep='_')})
#file=rbind(res[[1]],res[[2]])
#res=res[-c(1,2)]
#for (i in 1:length(res)){
#	file=rbind(file,res[[1]])
#	res=res[-1]
#}



# Coverage total
feature_depth=file %>% group_by(sample,feature) %>% summarise(total=sum(coverage))

# Normalize by size
feature_size=file%>% group_by(sample,feature) %>% summarise(size=max(Position))
feature_depth=left_join(feature_depth,feature_size,by=c('sample','feature'))
feature_depth$normalised_cov=feature_depth$total/feature_depth$size

feature_depth_plot=feature_depth %>% ggplot(aes(x=log2(normalised_cov+.1), fill=sample, color=sample, group=sample)) + 
	geom_density(alpha=0.3) + 
	scale_x_continuous(breaks=log2(c(0.1,mean(feature_depth$normalised_cov),max(feature_depth$normalised_cov))), labels=paste(round(c(log2(c(0.1,round(mean(feature_depth$normalised_cov),digits=2),round(max(feature_depth$normalised_cov), digits=2)))),digits=2),c(0,round(mean(feature_depth$normalised_cov),digits=2),round(max(feature_depth$normalised_cov), digits=2)), sep='\n')) + 
	theme(panel.background=element_rect(fill='white'), panel.grid.major=element_line(colour='grey',size=.3,linetype=2), panel.grid.minor=element_line(colour='grey',size=.3,linetype=2), axis.text.x=element_text(angle=90,vjust=0.5)) + xlab('Log2 - Feature-length-Normalised coverage')
ggsave(feature_depth_plot,filename='test_all.pdf', device='pdf', width=15 )


# zscore
feature_depth=feature_depth%>% group_by(sample) %>% mutate(mean=mean(normalised_cov),sd=sd(normalised_cov), zscore=scale(normalised_cov))

feature_depth_plot_zscore=feature_depth %>% ggplot(aes(x=zscore,fill=sample, color=sample, group=sample)) + 
	geom_histogram(bins=500,alpha=0.3) + 
	scale_x_continuous(breaks=c(0,5,7,10,13,15,17,20,max(feature_depth$zscore)), labels=c(0,5,7,10,13,15,17,20,round(max(feature_depth$zscore),digits=2)))
ggsave(feature_depth_plot_zscore,filename='test_all_zscore3.pdf', device='pdf', width=15 )

vals=c(3,4,5,6,7)
vals_matrix=matrix(ncol=length(vals),nrow=length(sample_names))
rownames(vals_matrix)=sample_names
colnames(vals_matrix)=vals

for (i in 1:length(vals)){
	temp=feature_depth %>% filter(zscore>=vals[i]) %>% group_by(sample) %>% summarise(nfeatures=length(feature))
	vals_matrix[,as.character(vals[i])]=temp[which(temp$sample %in% rownames(vals_matrix)),]$nfeatures
}
write.csv(vals_matrix,'Number_of_features_removed_zscore_threshold_2.csv')




# Feature bode coverage bias

# Calculate sd of feature coverage per feature per sample
print(Sys.time())
per_feature_sd=function(file=file, sample_name=sample_name){
	print(sample_name)
	test= file[file$sample==sample_name,]
	print(head(test))
	sd = test %>% group_by(feature) %>% summarise(per_feature_sd=sd(coverage))
	sd$sample=sample_name
	sd=as.data.frame(sd)
	print(head(sd))
	return(sd)
}
print(Sys.time())

#no_cores=detectCores()-1
#cl=makeCluster(6)
#clusterExport(cl,list('file','per_feature_sd'))
#clusterEvalQ(cl, library(dplyr))
#clusterEvalQ(cl, library(reshape))
#res=parLapply(cl,sample_names,function(x){per_feature_sd(file=file,sample_name=x)})
##stopCluster(cl)
print(Sys.time())
res=lapply(sample_names, function(x){per_feature_sd(file=file,sample_name=x)})
res=do.call(rbind,res)
res=res[!is.na(res$per_feature_sd),]
print(Sys.time())

per_feature_sd_plot=res %>% ggplot(aes(x=log2(per_feature_sd+.1),fill=sample, color=sample, group=sample)) + geom_density(alpha=0.3) + 
    scale_x_continuous(breaks=c(0,4,5,6,7,8,9,max(res$per_feature_sd)), 
    	labels=paste0(c(0,4,5,6,7,8,9,max(res$per_feature_sd)),'\n',2^(c(0,4,5,6,7,8,9,max(res$per_feature_sd)))))
ggsave(per_feature_sd_plot,filename='test_feature_sd2.pdf', device='pdf', width=15 )

res$per_feature_sd>7
outlier=res[res$per_feature_sd>128,]
outlier=unique(outlier$feature)
test=file %>% filter(feature%in%outlier)
not_outlier=unique(file$feature)[!unique(file$feature)%in%outlier]

per_feature_body_plot=test %>% ggplot(aes(x=Position,y=coverage,color=feature)) + geom_line() + facet_wrap(~sample) + theme(legend.position='none')
ggsave(per_feature_body_plot,filename='test_feature_bodycov.pdf', device='pdf', width=15 )


# Convert position to %position

calculate_coverage_per_chunk=function(x,df){
	val=mean(df[which(df$Position %in% x),]$coverage)
	return(val)
}

calculate_feature_values=function(table){
	pos=table$Position
	split_positions=split(pos, ceiling(seq_along(pos)/(length(pos)/100)))
	feature_values=round(unlist(lapply(split_positions,function(y){calculate_coverage_per_chunk(y,df=table)})),digits=2)
	feature_values_df=data.frame(Position=as.character(names(feature_values)), coverage=as.numeric(feature_values))
	return(feature_values_df)
}


res1=test  %>% group_by(sample,feature) %>% do(calculate_feature_values(.))
res1$Position=factor(res1$Position,levels=unique(res1$Position))
p= res1 %>% ggplot(aes(x=Position,y=coverage,color=feature,group=feature)) + geom_line() + facet_wrap(~sample) + theme(legend.position='none', axis.text.x=element_text(angle=90)) + ylim(0,1000)
ggsave(p,filename='test_feature_bodycov.pdf', device='pdf', width=15 )


univ=unique(res1$feature)
unc_sample=sample(univ,10,replace=FALSE)
try=res1%>%filter(feature%in%unc_sample)

p= try %>% ggplot(aes(x=Position,y=coverage,color=feature,group=feature)) + geom_line() + facet_wrap(~sample) + theme(legend.position='none', axis.text.x=element_text(angle=90)) + ylim(0,1000)
ggsave(p,filename='test_feature_bodycov.pdf', device='pdf', width=15 )


































