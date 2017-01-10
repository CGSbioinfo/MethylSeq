library(ggplot2)
library(reshape)

in_dir=commandArgs(TRUE)[1]
sample_names=commandArgs(TRUE)[2]
pattern=commandArgs(TRUE)[3]
outdir=commandArgs(TRUE)[4]

allFiles=list.files(in_dir, recursive=TRUE, pattern=pattern)
sample_names=as.character(read.csv(sample_names, header=FALSE)[,1])

data_alignment=data.frame(matrix(ncol=6))
colnames(data_alignment)=c("Sequence pairs analysed in total:","Number of paired-end alignments with a unique best hit:","Mapping efficiency:" ,"Sequence pairs with no alignments under any condition:","Sequence pairs did not map uniquely:","Sequence pairs which were discarded because genomic sequence could not be extracted:")

data_aligned_pairs=data.frame(matrix(ncol=4))
colnames(data_aligned_pairs)=c("converted top strand","complementary to converted top strand","complementary to converted bottom strand","converted bottom strand")

data_methylation=data.frame(matrix(ncol=12))
colnames(data_methylation)=c("Total methylated C's in CpG context:", "Total unmethylated C's in CpG context:", "C methylated in CpG context:", "Total methylated C's in CHG context:" , "Total unmethylated C's in CHG context:", "C methylated in CHG context:", "Total methylated C's in CHH context:", "Total unmethylated C's in CHH context:", "C methylated in CHH context:", "Total methylated C's in Unknown context:", "Total unmethylated C's in Unknown context:", "C methylated in unknown context (CN or CHN):")

for (i in 1:length(sample_names)) { 
	x=readLines(paste0(in_dir,'/',grep(sample_names[i], allFiles, value=TRUE)))

	start_line=grep('Final Alignment report', x)
	x=x[start_line:length(x)]

# Alignment
	df1=t(as.data.frame(strsplit(x[3:(which(x=='')[1]-1)],'\t')))
	rownames(df1)=1:nrow(df1)
	df1=t(df1)
	colnames(df1)=df1[1,]
	df1=df1[-1,]

	data_alignment=rbind(data_alignment,df1)

	df2=t(as.data.frame(strsplit(x[(which(x=='')[1]+2):(which(x=='')[2]-1)],'\t')))
        rownames(df2)=1:nrow(df2)
        df2=t(df2[,2:3])
        colnames(df2)=df2[2,]
        df2=df2[-2,]
        names(df2)=gsub("\\(",'',names(df2))
	names(df2)=gsub("\\)",'',names(df2))

	data_aligned_pairs=rbind(data_aligned_pairs,df2)	

	# Final Cytosine Methylation Report
	
	start_line=grep('Final Cytosine Methylation Report', x)
	x=x[start_line:length(x)]

	# CpG context 
	cpg=grep('CpG',x, value=TRUE)
	cpg=t(as.data.frame(strsplit(cpg, '\t')))
	rownames(cpg)=1:nrow(cpg)
	cpg=t(cpg)
	colnames(cpg)=cpg[1,]
	cpg=cpg[-1,]
	
	# CHG
	chg=grep('CHG',x, value=TRUE)
	chg=t(as.data.frame(strsplit(chg, '\t')))
	rownames(chg)=1:nrow(chg)
	chg=t(chg)
	colnames(chg)=chg[1,]
	chg=chg[-1,]

	# CHG
	chh=grep('CHH',x, value=TRUE)
	chh=t(as.data.frame(strsplit(chh, '\t')))
	rownames(chh)=1:nrow(chh)
	chh=t(chh)
	colnames(chh)=chh[1,]
	chh=chh[-1,]

	# Unknown
	unkn=grep('nknown',x, value=TRUE)
	unkn=t(as.data.frame(strsplit(unkn, '\t')))
	rownames(unkn)=1:nrow(unkn)
	unkn=t(unkn)
	colnames(unkn)=unkn[1,]
	unkn=unkn[-1,]

	vec2=c(cpg,chg,chh,unkn)

	data_methylation=rbind(data_methylation,vec2)
}

data_alignment=data_alignment[-1,]
rownames(data_alignment)=sample_names

data_aligned_pairs=data_aligned_pairs[-1,]
rownames(data_aligned_pairs)=sample_names

data_methylation=data_methylation[-1,]
rownames(data_methylation)=sample_names


## Alignment plot
# aligned pairs
total=rowSums(apply(as.matrix(data_aligned_pairs),2,as.numeric))
data_aligned_pairs_pct=apply(as.matrix(data_aligned_pairs),2,as.numeric)/total*100
rownames(data_aligned_pairs_pct)=rownames(data_aligned_pairs)

plot_data=cbind(data_aligned_pairs_pct,Mapping_effiency=as.numeric(gsub('%','',data_alignment[,which(colnames(data_alignment)=='Mapping efficiency:')])))
colnames(plot_data)=c('converted\ntop strand','complementary\nto converted\ntop strand','complementary\nto converted\nbottom strand','converted\nbottom strand','Mapping_effiency')
melt_plot_data=melt(plot_data)
melt_plot_data$layer='Sequence pairs with unique best alignment'
melt_plot_data$layer[melt_plot_data$X2=='Mapping_effiency']='Mapping effiency'
melt_plot_data[melt_plot_data$layer=='Sequence pairs with unique best alignment',]$layer='Directionality'
pdf(paste0(outdir,'/mappingQC_efficiency_and_strand.pdf'),width=12)
ggplot(melt_plot_data,aes(x=X1,y=value,col=X2)) + geom_point(size=3) + facet_wrap(~layer) + theme_bw()  + theme(panel.grid.major=element_line(colour='#000000',linetype='dashed'), legend.title=element_blank(),axis.text.x=element_blank(), axis.text=element_text(size=20), axis.title=element_text(size=20), legend.text=element_text(size=20), strip.text=element_text(size=20), legend.position='bottom') + xlab('Sample') + ylab('Percentage')
dev.off()


# Methylated C's plot 1 
plot_data=data_methylation[,-grep('%',data_methylation)]
plot_data=data.frame(suppressWarnings(apply(plot_data,2,as.numeric)), row.names=rownames(plot_data))
all_cs=data.frame(matrix(0,ncol=4,nrow=nrow(plot_data)))
context=c('CpG.context','CHG.context','CHH.context','Unknown.context')
colnames(all_cs)=context
for (i in 1:length(context)){
	all_cs[,grep(context[i],colnames(all_cs))]=rowSums(plot_data[,grep(context[i],colnames(plot_data))])
}
rownames(all_cs)=rownames(plot_data)
all_cs_pct=(all_cs/rowSums(all_cs))*100
all_cs_pct.melt=melt(data.frame(Sample=rownames(all_cs_pct),all_cs_pct))
p=ggplot(all_cs_pct.melt,aes(x=Sample,y=value,color=variable, fill=variable)) + geom_bar(stat='identity') + theme(axis.text.x=element_text(angle=90,vjust=.5), legend.title=element_blank()) + ylab('Percentage')
ggsave(p,filename=paste0(outdir,'/mappingQC_methC_context.pdf'), width=12)

context_pct_range=round(apply(all_cs_pct,2,range))
context_pct_range=apply(context_pct_range,2,paste0,collapse='-')

# Methylated C's plot 2
plot_data=data_methylation[,grep('%',data_methylation)]
plot_data=apply(plot_data,2,gsub, patter='%', replacement='')
plot_data=apply(plot_data,2,as.numeric)
rownames(plot_data)=rownames(data_methylation)
melt_plot_data=melt(plot_data)
melt_plot_data$X2=as.character(melt_plot_data$X2)
melt_plot_data$X2[melt_plot_data$X2=='C methylated in CpG context:']=paste0('CpG context  (', context_pct_range[1],'%)')
melt_plot_data$X2[melt_plot_data$X2=='C methylated in CHG context:']=paste0('CHG context  (', context_pct_range[2],'%)')
melt_plot_data$X2[melt_plot_data$X2=='C methylated in CHH context:']=paste0('CHH context  (', context_pct_range[3],'%)')
melt_plot_data$X2[melt_plot_data$X2=='C methylated in unknown context (CN or CHN):']=paste0('Unknown context  (', context_pct_range[4],'%)')

p=ggplot(melt_plot_data,aes(x=X1,y=value,col=X2)) + geom_point() + theme(axis.text.x=element_text(angle=90), legend.position='none') + xlab('Sample') + ylab('Percentage of methylated C\'s') + facet_wrap(~X2)
ggsave(p,filename=paste0(outdir,'/mappingQC_methC.pdf'), width=12)
melt_plot_data=melt_plot_data[melt_plot_data$X2=='CpG context  (9-11%)',]
p=ggplot(melt_plot_data,aes(x=X1,y=value)) + geom_point(size=3) + theme_bw() + theme(panel.grid.major=element_line(linetype='dashed',colour='grey'),axis.text.y=element_text(size=22), axis.title=element_text(size=22),legend.position='none', axis.text.x=element_blank()) + xlab('Sample') + ylab('Percentage of methylated C\'s') + ylim(0,100)



#Tables
write.csv(data_alignment, paste0(outdir,'/alignment_report1.csv'),quote=F)
write.csv(data_aligned_pairs, paste0(outdir,'/alignment_report2.csv'),quote=F)
write.csv(data_methylation, paste0(outdir,'/cytosine_methylation_report.csv'),quote=F)
