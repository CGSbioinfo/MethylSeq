# Perform differential methylation analysis.
# Extended comments are copied from the BiSeq guide:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/BiSeq/inst/doc/BiSeq.pdf

cat("Loading R packages ...\n")
suppressMessages(library(BiSeq))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(parallel))

sampleInfoFile           =commandArgs(TRUE)[1]
sampleGroupFile          =commandArgs(TRUE)[2]
bedGraphFolder           =commandArgs(TRUE)[3]
out_wholeGenome          =commandArgs(TRUE)[4]
out_targetRegion         =commandArgs(TRUE)[5]
annotation_target_regions=commandArgs(TRUE)[6]

ncores = 20
dmr.bandwidth = 50

gtf.file = '/mnt/cgs-fs3/Sequencing/Genome/Pig/ensembl/gtf/Sscrofa10.2/release-85/Sus_scrofa.Sscrofa10.2.85.gtf'
bed.file = annotation_target_regions

##################
#
# Load external 
# functions
#
##################

cat("Loading functions ...\n")

file.arg.name   = "--file="
script.name     = sub(file.arg.name, "", commandArgs()[grep(file.arg.name, commandArgs())])
script.basename = dirname(script.name)
other.name      = paste(sep="/", script.basename, "functions.r")

source(other.name)

##################
#
# Check validity of 
# file paths
#
##################

# functions::checkPath() will quit if file not found
checkPath(sampleInfoFile, "Sample name file") 
checkPath(sampleGroupFile, "Sample group file") 
checkPath(bedGraphFolder, "BedGraph folder")
checkPath(annotation_target_regions, "Genome annotation")

checkPath(dirname(out_wholeGenome), "Whole genome output folder")
checkPath(dirname(out_targetRegion), "Target region output folder")

tempImagePath = paste0(dirname(sampleInfoFile),"/tmpImages/")
if(!dir.exists(tempImagePath)){
  dir.create(tempImagePath)
}

printTime = function(ptc){
  days  = floor(ptc / 86400)
  ptc   = ptc - (days*86400)
  hours = floor(ptc / 3600)
  ptc   = ptc - (hours*3600)
  mins  = floor(ptc / 60)
  secs  = ptc - (mins*60)
  cat("Completed in", days[3], "days", hours[3], "hours", mins[3],"minutes and", secs[3], "seconds\n" ) # elapsed time
}

# Check if the next analysis step output exists. If so, skip 
# this analysis step. If not, check if this step's output exists
# and if not, run the analysis.
runOrSkip = function(nextStepTempFile, tempFile, runFunction){
  if(!file.exists(paste0(tempImagePath, nextStepTempFile))){
    tempFile = paste0(tempImagePath, tempFile)
    if(file.exists(tempFile)){
      cat("Loading temporary analysis file", tempFile,"...\n")
      load(tempFile, envir = globalenv())
      cat("Loaded temporary analysis file", tempFile,"...\n")
    } else{
      ptm = proc.time()
      runFunction()
      printTime(proc.time() - ptm)
      cat("Saving temporary analysis file", tempFile, "...\n")
      save.image(tempFile) # save workspace
      
    }
  }
}

##################
#
# Read Sample information
#
##################

readSamples = function(){
  cat("Reading sample names file ...\n")
  sampleInfo = read.table(sampleInfoFile, sep="\t", header=F)
  rownames(sampleInfo) = as.character(sampleInfo$V1)
  colnames(sampleInfo) = c("Sample.Name")
  assign( "sampleInfo", sampleInfo, envir=globalenv())

  cat("Reading sample groups file ...\n")
  sampleGroups = read.table(sampleGroupFile, sep="\t", header=T)
  rownames(sampleGroups) = as.character(sampleGroups$Sample.Name)
  assign( "sampleGroups", sampleGroups, envir=globalenv())

  # Locate coverage files
  cat("Locating cov files ...\n")
  listMethylFile<-list.files(path=bedGraphFolder,pattern=".cov",full.names = T)
  listMethylFile=grep(paste0(c(as.character(sampleGroups$Sample.Name)),collapse='|'),
      listMethylFile, value=TRUE)
  assign( "listMethylFile", listMethylFile, envir=globalenv())

  print(listMethylFile, collapse="\n\t")
  methylDataRaw = readBismark(listMethylFile, sampleGroups)
  assign( "methylDataRaw", methylDataRaw, envir=globalenv())
}

runOrSkip("Filtered_methylation.RData", "Read_files.RData", readSamples)

##################
#
# Read the annotation 
# target regions
#
##################

readTargetRegions = function(){
  cat('Reading target regions\n')

  capture.region = import.bed(annotation_target_regions)
  # This creates a GRanges object with the chr, start, end of the agilent target regions:annotation_capture to look at it
  # Merge any overlapping regions, convert to data frame 
  capture.region<-reduce(capture.region) %>% # region distribution: distance between regions and region size
    as.data.frame %>%
    arrange(seqnames, start) %>% # sort the regions by chromosome (seqnames)and start position
    group_by(seqnames) %>%  # calculate distance between regions
    mutate(dist=c(start[-1], end[length(end)])-end)  %>%  # calculate distance from one region to the next one
    group_by(seqnames) %>% 
    mutate(size=end-start)  # Calculate size of each region

  capture.region = GRanges(seqnames=capture.region$seqnames,
    IRanges(start=capture.region$start, end=capture.region$end),
                              strand=capture.region$strand, dist = capture.region$dist, size= capture.region$size)

  capture.region$region = paste0('region_', 1:length(capture.region))

  assign( "capture.region", capture.region, envir=globalenv())

  #if you want to change the minimum coverage, change the minCov setting in the line below
  methylDataRaw.filter10 = filterByCov(methylDataRaw, minCov=10, global=F)
  assign( "methylDataRaw.filter10", methylDataRaw.filter10, envir=globalenv())
  # subset the methylDataRaw object that overlaps the target regions

  methylDataRaw.filter10.rk = subsetByOverlaps(methylDataRaw.filter10, capture.region)
  assign( "methylDataRaw.filter10.rk", methylDataRaw.filter10.rk, envir=globalenv())
}

runOrSkip("Predicted_methylation.RData", "Filtered_methylation.RData", readTargetRegions)

##################
#
# Cluster data
#
##################

# The next line is clustering the data, 
# you can set up the stringency of coverage using perc.sample, of clustering using max.dist and min.site,
# If this step crashes you might consider reducing mc.cores

defineCpGClusters = function(){
  cat('Defining Cpg clusters\n')

  data.clustered = clusterSites(methylDataRaw.filter10.rk,
    groups=colData(methylDataRaw.filter10.rk)$Group,
    perc.samples=0.25, 
    min.site=5, 
    max.dist=100, 
    mc.cores=ncores)

  # In the smoothing step CpG sites with high coverages get high weights. To
  # reduce bias due to unusually high coverages we limit the coverage.
  # Default is the 90% quantile
  ind.cov = totalReads(data.clustered) > 0
  quant   = quantile(totalReads(data.clustered)[ind.cov],0.9)
  methylDataRaw.filter10.rk.clustered.lim = limitCov(data.clustered, maxCov=quant)
  assign("methylDataRaw.filter10.rk.clustered.lim", methylDataRaw.filter10.rk.clustered.lim, envir=globalenv())

  #   We then smooth the methylation values of CpG sites within the clusters with
  # the given bandwidth (h)
  cat("Smoothing methylation ...\n")
  predictedMeth = predictMeth(methylDataRaw.filter10.rk.clustered.lim, h=dmr.bandwidth, mc.cores=ncores)

  assign( "predictedMeth", predictedMeth, envir=globalenv())

  # makeSmoothPlot = function(sample.name){
  #   png(file=paste0(tempImagePath,"Smoothed methylation", sample.name,".png"),width = 480, height=480)
  #   plotMeth(object.raw = methylDataRaw[,rownames(methylDataRaw)==sample.name], object.rel = predictedMeth[,6], region = region,
  #     lwd.lines = 2,
  #     col.points = "blue",
  #     cex = 1.5)
  #   dev.off()
  # }

}

runOrSkip("Beta_regression.RData", "Predicted_methylation.RData", defineCpGClusters)

##################
#
# Run regression
#
##################

runBetaRegression = function(){
  cat("Running beta regression ...\n")

  # To detect the CpG sites where the DNA methylation differs between case
  # and control samples we model the methylation within a beta regression with
  # the group as explanatory variable and use the Wald test to test if there is a
  # group effect
  # By setting type = "BR" the maximum likelihood with bias reduction is
  # called.  This is especially useful, when the sample size is small.
  betaResults = betaRegression(formula=~Group, link='probit', type='BR', object=predictedMeth, mc.cores=ncores)
  assign( "betaResults", betaResults, envir=globalenv())


  # The aim is to detect CpG clusters containing at least one differentially methylated location.
  # To do so the P values p from the Wald tests are transformed to Z scores which are normally 
  # distributed under Null hypothesis (no group effect). As cluster test statistic a standardized
  # Z score average is used. To estimate the standard deviation of the Z scores we have to estimate 
  # the correlation and hence the variogram of methylation between two CpG sites within a cluster.
  # The estimation of the standard deviation requires that the distribution of the Z scores 
  # follows a standard normal distribution.  However, if methylation in both groups differs for 
  # many CpG sites the density distribution of P values shows a peak near 0.  To ensure that the 
  # P values are roughly uniformly distributed to get a variance of the Z scores that is Gaussian
  # with variance 1 we recommend to estimate the variogram (and hence the correlation of Z scores)
  # under the null hypothesis. To do so we model the beta regression again for resampled data.

  cat("Creating resampled model for null hypothesis ...\n")

  # Since there are only a small number of samples, use an even division.
  # The group null labels should cover the number of samples in your predictmeth object
  # Do not use real group assignments. I.e. have an even split of case and control between
  # samples.

  # Since there are 3 case and 3 control, we need to select 2 and 2.
  # TODO - make dynamic based on group column
  predictedMethNull = predictedMeth[,c(1, 3, 2, 5)]

  colData(predictedMethNull)$group.null=rep(c(1,2), nrow(colData(predictedMethNull))/2)

  cat("Selected samples for modelling:\n")
  cat(colData(predictedMethNull))

  betaResultsNull = betaRegression(formula=~group.null, 
    link='probit', 
    object=predictedMethNull, 
    type='BR',
    mc.cores=ncores)

  assign( "betaResultsNull", betaResultsNull, envir=globalenv())

  # save.image(paste0(tempImagePath, "Beta_regression_null.RData")) # save workspace

  vario = makeVariogram(betaResultsNull)

  assign( "vario", vario, envir=globalenv())
  # save.image(paste0(tempImagePath, "Null_variogram.RData")) # save workspace

  # Plot the variogram
  png(file=paste0(tempImagePath,"Null_variogram.png"),width = 480, height=480)
  plot(vario$variogram$v)
  dev.off()
}

runOrSkip("DMRs_found.RData", "Beta_regression_null.RData", runBetaRegression)

##################
#
# Smooth
#
##################

smoothValues = function(){
  cat('Smoothing\n')
  ## Break, and call new script with desired sill value, or continue with default sill of 1?
  ## sill_value= readline("What is the value of sill to use (0 to 1)? please check graph ")
  sill_value= 1

  vario.sm = smoothVariogram(vario, sill=as.numeric(sill_value))
  assign( "vario.sm", vario.sm, envir=globalenv())

  # Plot the smoothed variogram
  png(file=paste0(tempImagePath,"Smooth_variogram.png"),width = 480, height=480)
  plot(vario.sm$variogram[,c("h", "v.sm")], type="l", col = "red", lwd = 1.5)
  dev.off()

  vario.aux = makeVariogram(betaResults, make.variogram=F)
  vario.sm$pValsList  = vario.aux$pValsList
  locCor = estLocCor(vario.sm)

  clusters.rej     = testClusters(locCor, FDR.cluster=0.1)
  clusters.trimmed = trimClusters(clusters.rej, FDR.loc=0.05)

  DMRs    = findDMRs(clusters.trimmed, max.dist=100, diff.dir=T)
  assign( "DMRs", DMRs, envir=globalenv())

  DMRs.df = data.frame(DMRs)
  assign( "DMRs.df", DMRs.df, envir=globalenv())
}

runOrSkip("Overlaps_with_GTF.RData", "DMRs_found.RData", smoothValues)

##################
#
# Read GTF files
#
##################

findOverlapsWithGTF = function(){
  cat('Finding overlaps with GTF genes\n')
  ##Plot DMRs
  merge_hits=function(x){
    data.frame(seqnames=x$seqnames[1], start=x$start[1], end=x$end[1], width=x$width[1],
                         median.p=x$median.p[1],
                         median.meth.group1=x$median.meth.group1[1],
                         median.meth.group2=x$median.meth.group2[1],
                         median.meth.diff=x$median.meth.diff[1],
                         annotation=paste(x$annotation, collapse=' | '))
  }

  # loading annotation, if using different kit/annotation resources, please change to the correct path
  GTF     = import.gff(gtf.file, format="gtf",feature.type="gene")
  bed.kit = import.bed(bed.file)

  assign( "GTF", GTF, envir=globalenv())
  assign( "bed.kit", bed.kit, envir=globalenv())

  tempoverlap1     = findOverlaps(DMRs, GTF)
  tempoverlap1.kit = findOverlaps(DMRs, bed.kit)
  tempoverlap2     = DMRs[queryHits(tempoverlap1)]
  tempoverlap2.kit = DMRs[queryHits((tempoverlap1.kit))]

  tempoverlap2$annotation=paste(GTF[subjectHits(tempoverlap1)]$gene_id,
                                GTF[subjectHits(tempoverlap1)]$type,
                                GTF[subjectHits(tempoverlap1)]$gene_biotype,
                                GTF[subjectHits(tempoverlap1)]$gene_name, sep='; ')

  tempoverlap2.kit$annotation = bed.kit[subjectHits(tempoverlap1.kit)]$name
  tempoverlap2     = as.data.frame(tempoverlap2)
  tempoverlap2.kit = data.frame(tempoverlap2.kit)
  tempoverlap2$seqnames2     = paste0(tempoverlap2$seqnames,'_', tempoverlap2$start)
  tempoverlap2.kit$seqnames2 = paste0(tempoverlap2.kit$seqnames,'_', tempoverlap2.kit$start)

  DMRs.annotated     = tempoverlap2 %>% group_by(seqnames2) %>% do(merge_hits(.)) %>% as.data.frame
  DMRs.annotated.kit = tempoverlap2.kit %>% group_by(seqnames2) %>% do(merge_hits(.)) %>% as.data.frame

  assign( "DMRs.annotated", DMRs.annotated, envir=globalenv())
  assign( "DMRs.annotated.kit", DMRs.annotated.kit, envir=globalenv())

  write.table(file=paste0(tempImagePath, "DMRs_annotated_kit.tsv"),   sep="\t",row.names = F, DMRs.annotated.kit)
  write.table(file=paste0(tempImagePath, "DMRs_annotated_genes.tsv"), sep="\t",row.names = F, DMRs.annotated)
}

runOrSkip("DMRs.RData", "Overlaps_with_GTF.RData", findOverlapsWithGTF)

##################
#
# Plot DMRs
#
##################

plotDMRs = function(){

  outputDir = paste0(dirname(sampleInfoFile),"/DMR_images/")
  if(!dir.exists(outputDir)){
    dir.create(outputDir)
  }

  for (i in 1:length(DMRs)){
    temp = subsetByOverlaps(predictedMeth, DMRs[i])
    t2   = rowRanges(temp) %>% as.data.frame
    t2   =c bind(index=as.numeric(rownames(t2)), t2)
    t2$seqnames2 = paste0(t2$seqnames, '_', t2$start)
    df = melt(methLevel(temp))

    annot = rowRanges(temp)
    annot$seqnames2 = paste0(seqnames(annot), '_', start(annot))
    annot = left_join(as.data.frame(annot), DMRs.annotated,by='seqnames2')
    annot = paste0(annot[complete.cases(annot),]$annotation, collapse= ' \n ')

    colnames(df) = c('index','Sample.Name','methLevel')
    df = left_join(df,as.data.frame(colData(temp)), by='Sample.Name')
    df = left_join(df,t2,by='index')

    positions = rowRanges(subsetByOverlaps(methylDataRaw.filter10.rk.clustered,DMRs[i])) %>% as.data.frame
    coverages = totalReads(subsetByOverlaps(methylDataRaw.filter10.rk.clustered,DMRs[i])) %>% as.data.frame
    positions = cbind(index=rownames(positions), positions)
    coverages = cbind(index=rownames(coverages), coverages)
    coverages = melt(coverages)
    cov = left_join(coverages,positions,by='index')
    cov$seqnames2 = paste0(cov$seqnames, '_', cov$start)
    colnames(cov)[2] = 'Sample.Name'
    cov = cov %>% select(Sample.Name,value,seqnames2)

    df = left_join(df,cov,by=c('Sample.Name','seqnames2'))
    df$seqnames2 = factor(df$seqnames2,levels=unique(df$seqnames2))

    p = ggplot(df,aes(x=seqnames2, y=methLevel, col=value)) +
     geom_point() + 
     facet_wrap(~Group) + 
     ggtitle(annot) + 
     scale_fill_gradient() + 
     theme(axis.text.x=element_text(angle=90,vjust=0.5))

    #change those 2 lines to change the location where the graphs are saved.
    ggsave(plot=p, filename=paste0(outputDir, 'DMR_',df$seqnames2[1],'.png'))
  }

}

runOrSkip("end.RData", "DMRs.RData", plotDMRs)
