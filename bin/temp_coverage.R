library(ggplot2)
library(reshape)
library(dplyr)

file='temp'
file=read.table(file)
file$feature=paste(file$V1,file$V2,file$V3,sep='_')
file=file[,-c(1:3)]

colnames(file)[1:4]=c('depth','number_of_positions_with_depth','feature_size','fraction')
file$depth_total=file$depth*file$number_of_positions_with_depth

feature_depth=file%>% group_by(feature) %>% summarise(total=sum(depth_total))

feature_depth_plot=feature_depth %>% ggplot(aes(x=log2(total+.1))) + geom_histogram(bins=300) + scale_x_continuous(breaks=c(-3.32,0,5,10,15,20), labels=paste(c(-3.32,0,5,10,15,20),format(round(2^c(-3.32,0,5,10,15,20),digits=2),big.mark=','), sep='\n'))




expanded_file=untable(file,num=file$number_of_positions_with_depth)

mean = expanded_file %>% group_by(feature) %>% summarise(summary=mean(depth))
median = expanded_file %>% group_by(feature) %>% summarise(summary=median(depth))
min = expanded_file %>% group_by(feature) %>% summarise(summary=min(depth))
sd = expanded_file %>% group_by(feature) %>% summarise(summary=sd(depth))
pc1 = expanded_file %>% group_by(feature) %>% summarise(summary=quantile(depth, .25))
pc2 = expanded_file %>% group_by(feature) %>% summarise(summary=quantile(depth, .50))
pc3 = expanded_file %>% group_by(feature) %>% summarise(summary=quantile(depth, .75))
pc4 = expanded_file %>% group_by(feature) %>% summarise(summary=quantile(depth, .95))
pc5 = expanded_file %>% group_by(feature) %>% summarise(summary=quantile(depth, .99))
max = expanded_file %>% group_by(feature) %>% summarise(summary=max(depth))