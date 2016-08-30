library(ggplot2)
library(reshape)

in_dir=commandArgs(TRUE)[1]
sample_names=commandArgs(TRUE)[2]
pattern=commandArgs(TRUE)[3]
outdir=commandArs(TRUE)[4]

allFiles=list.files(in_dir, recursive=TRUE, pattern=pattern)
sample_names=as.character(read.csv(sample_names, header=FALSE)[,1])

data_alignment=data.frame(matrix(ncol=6))
colnames(data_alignment)=c("Sequence pairs analysed in total:","Number of paired-end alignments with a unique best hit:","Mapping efficiency:" ,"Sequence pairs with no alignments under any condition:","Sequence pairs did not map uniquely:","Sequence pairs which were discarded because genomic sequence could not be extracted:")

data_aligned_pairs=data.frame(matrix(ncol=4))
colnames(data_aligned_pairs)=c("converted top strand","complementary to converted top strand","complementary to converted bottom strand","converted bottom strand")

data_methylation=data.frame(matrix(ncol=12))
colnames(data_methylation)=c("Total methylated C's in CpG context:", "Total unmethylated C's in CpG context:", "C methylated in CpG context:", "Total methylated C's in CHG context:" , "Total unmethylated C's in CHG context:", "C methylated in CHG context:", "Total methylated C's in CHH context:", "Total unmethylated C's in CHH context:", "C methylated in CHH context:", "Total methylated C's in Unknown context:", "Total unmethylated C's in Unknown context:", "C methylated in unknown context (CN or CHN):")

for (i in 1:length(sample_names)) { 
	x=readLines(grep(sample_names[i], allFiles, value=TRUE))

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
melt_plot_data$layer='Sequence pairs with unique best alignment'
melt_plot_data$layer[melt_plot_data$X2=='Mapping_effiency']='Mapping effiency'
ggplot(melt_plot_data,aes(x=X1,y=value,col=X2)) + geom_point() + facet_wrap(~layer) + theme(axis.text=element_text(angle=90))


# Methylated C's plot
plot_data=data_methylation[,grep('%',data_methylation)]
plot_data=apply(plot_data,2,gsub, patter='%', replacement='')
plot_data=apply(plot_data,2,as.numeric)
rownames(plot_data)=rownames(data_methylation)
melt_plot_data=melt(plot_data)
p=ggplot(melt_plot_data,aes(x=X1,y=value,col=X2)) + geom_point() + theme(axis.text=element_text(angle=90))
ggsave(p,filename='plot.pdf')


#Tables
write.csv(data_alignment, 'alignment_report1.csv',quote=F)
write.csv(data_aligned_pairs, 'alignment_report2.csv',quote=F)
write.csv(data_methylation, 'cytosine_methylation_report.csv',quote=F)