library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(reshape)
library(ggplot2)

load('predictedMeth.RData')

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

# Coverage in targeted region
coverage_per_site_distribution_table=apply(totalReads(methylDataRaw.rk), 2, coverage_per_site_distribution)
coverage_per_site_distribution_table[is.na(coverage_per_site_distribution_table)]=0
coverage_per_site_distribution_table.melt=melt(coverage_per_site_distribution_table)
coverage_per_site_distribution_table.melt$X1=factor( coverage_per_site_distribution_table.melt$X1, levels=unique(coverage_per_site_distribution_table.melt$X1))
coverage_per_site_distribution_table.melt$X2=as.character(coverage_per_site_distribution_table.melt$X2)
colnames(coverage_per_site_distribution_table.melt)[1:2]=c('Coverage','Sample')
coverage_per_site_distribution_table.melt$percent=coverage_per_site_distribution_table.melt$value/nrow(totalReads(methylDataRaw.rk))
coverage_per_site_distribution_table.melt$Region='Sites in targeted regions'
p1=ggplot(coverage_per_site_distribution_table.melt, aes(x=Coverage,y=value/nrow(totalReads(methylDataRaw.rk))), group=Sample, color=Sample) + geom_line(aes(group=Sample, color=Sample)) + scale_y_continuous(labels=scales::percent) + theme_bw() + scale_x_discrete(labels=c('0X','1-10X','10-100X','100-1000X','>1000X')) + ylab('') + xlab('Coverage per site')
ggsave(p1, filename='coverage_per_site_kit.png')

# Coverage not in targeted regions
all_df=as.data.frame(rowRanges(methylDataRaw))
kit_df=as.data.frame(rowRanges(methylDataRaw.rk))
all_df=paste(all_df$seqnames, all_df$start, all_df$end, sep='_')
kit_df=paste(kit_df$seqnames, kit_df$start, kit_df$end, sep='_')
diff=setdiff(all_df,kit_df)
not_target_regions=GRanges(seqnames=gsub('_.*','',diff), ranges=IRanges(start=as.numeric(gsub('.*_','',gsub('_\\d*$','',diff))),end=as.numeric(gsub('.*_','', diff))))
methylDataRaw.not_target_regions<-subsetByOverlaps(methylDataRaw, not_target_regions)
coverage_per_site_distribution_table_not_targeted=apply(totalReads(methylDataRaw.not_target_regions), 2, coverage_per_site_distribution)
coverage_per_site_distribution_table_not_targeted[is.na(coverage_per_site_distribution_table_not_targeted)]=0
coverage_per_site_distribution_table_not_targeted.melt=melt(coverage_per_site_distribution_table_not_targeted)
coverage_per_site_distribution_table_not_targeted.melt$X1=factor( coverage_per_site_distribution_table_not_targeted.melt$X1, levels=unique(coverage_per_site_distribution_table_not_targeted.melt$X1))
coverage_per_site_distribution_table_not_targeted.melt$X2=as.character(coverage_per_site_distribution_table_not_targeted.melt$X2)
colnames(coverage_per_site_distribution_table_not_targeted.melt)[1:2]=c('Coverage','Sample')
coverage_per_site_distribution_table_not_targeted.melt$percent=coverage_per_site_distribution_table_not_targeted.melt$value/nrow(totalReads(methylDataRaw.not_target_regions))
coverage_per_site_distribution_table_not_targeted.melt$Region='Sites not in targeted regions'

# Rbind in targeted and not in targeted regions
coverage_per_site_distribution_table.melt=rbind(coverage_per_site_distribution_table.melt, coverage_per_site_distribution_table_not_targeted.melt)

# Plot 
p1=ggplot(coverage_per_site_distribution_table.melt, aes(x=Coverage,y=percent), group=Sample, color=Sample) + geom_line(aes(group=Sample, color=Sample)) + scale_y_continuous(labels=scales::percent) + theme_bw() + scale_x_discrete(labels=c('0X','1-10X','10-100X','100-1000X','>1000X')) + ylab('') + xlab('Coverage per site') + facet_wrap(~Region) + theme(panel.grid.major=element_line(colour='#333333', linetype='dotted'), axis.text.x=element_text(size=28), axis.text.y=element_text(size=28), strip.text=element_text(size=28), legend.text=element_text(size=28))
ggsave(p1,filename='test.png', width=15)


# Table f coverage at 1-10X and >10X
coverage_1_10X=rbind(in_target=coverage_per_site_distribution_table['one_ten',], not_in_target=coverage_per_site_distribution_table_not_targeted['one_ten',])
coverage_1_10X=t(t(coverage_1_10X)/colSums(coverage_1_10X))*100
write.csv(coverage_1_10X,'less_than_10X_coverage.csv')

coverage_more_than_10X=rbind(in_target=colSums(coverage_per_site_distribution_table[3:5,]), not_in_target=colSums(coverage_per_site_distribution_table_not_targeted[3:5,]))
coverage_more_than_10X=t(t(coverage_more_than_10X)/colSums(coverage_more_than_10X))*100
write.csv(coverage_more_than_10X,'more_than_10X_coverage.csv')








