# Perform differential methylation analysis using BiSeq.
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

meta.env = new.env() # parameters for the analysis
data.env = new.env() # data for the ongoing analysis
func.env = new.env() # functions defined in this file

meta.env$sample.info.file   = commandArgs(TRUE)[1]
meta.env$sample.group.file  = commandArgs(TRUE)[2]
meta.env$cov.file.folder    = commandArgs(TRUE)[3]
meta.env$target.region.file = commandArgs(TRUE)[4]
meta.env$temp.image.folder  = commandArgs(TRUE)[5]
meta.env$ncores             = as.numeric(commandArgs(TRUE)[6])
meta.env$chunk.size         = as.numeric(commandArgs(TRUE)[7])

meta.env$dmr.bandwidth = 50

meta.env$gtf.file = '/mnt/cgs-fs3/Sequencing/Genome/Pig/ensembl/gtf/Sscrofa10.2/release-85/Sus_scrofa.Sscrofa10.2.85.gtf'

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
checkPath(meta.env$sample.info.file, "Sample name file") 
checkPath(meta.env$sample.group.file, "Sample group file") 
checkPath(meta.env$cov.file.folder, "BedGraph folder")
checkPath(meta.env$target.region.file, "Genome annotation")

meta.env$tempImagePath = paste0(dirname(meta.env$sample.info.file),"/", meta.env$temp.image.folder, "/")
if(!dir.exists(meta.env$tempImagePath)){
  dir.create(meta.env$tempImagePath)
}

##################
#
# Define the functions
#
##################

func.env$runOrSkip = function(nextStepTempFile, tempFile, runFunction){
  # Check if the next analysis step output exists. If so, skip 
  # this analysis step. If not, check if this step's output exists
  # and if not, run the analysis.
  #
  # nextStepTempFile - the name of the temp file after this step
  # tempFile - the name of the temp file for this step
  # runFunction - the function to be run in this step
  if(!file.exists(paste0(meta.env$tempImagePath, nextStepTempFile))){
    tempFile = paste0(meta.env$tempImagePath, tempFile)
    if(file.exists(tempFile)){
      cat("Loading temporary analysis file", tempFile,"...\n")
      RunAndTime( f=function(){load(tempFile, envir = globalenv())} )
      cat("Loaded temporary analysis file", tempFile,"...\n")
    } else{
      RunAndTime(runFunction)
      cat("Saving temporary analysis file", tempFile, "...\n")

      # Only save the data values - ensures updated functions will not be overwritten 
      RunAndTime(f=function(){save(data.env, file = tempFile)})
      # RunAndTime(f=function(){save(vars.only, file=tempFile)})
      # RunAndTime(f=function(){save.image(tempFile)})
    }
  }
}

func.env$readSamples = function(){
  cat("Reading sample names file ...\n")
  sampleInfo = read.table(meta.env$sample.info.file, sep="\t", header=F)
  rownames(sampleInfo) = as.character(sampleInfo$V1)
  colnames(sampleInfo) = c("Sample.Name")
  assign( "sampleInfo", sampleInfo, envir=data.env)

  cat("Reading sample groups file ...\n")
  sampleGroups = read.table(meta.env$sample.group.file, sep="\t", header=T)
  rownames(sampleGroups) = as.character(sampleGroups$Sample.Name)
  assign( "sampleGroups", sampleGroups, envir=data.env)

  # Locate coverage files
  cat("Locating cov files ...\n")
  listMethylFile = list.files(path=meta.env$cov.file.folder,pattern=".cov",full.names = T)
  listMethylFile = grep(paste0(c(as.character(sampleGroups$Sample.Name)),collapse='|'),
      listMethylFile, value=TRUE)

  print(listMethylFile, collapse="\n\t")
  methylDataRaw = readBismark(listMethylFile, sampleGroups)
  assign( "methylDataRaw", methylDataRaw, envir=data.env)
}

func.env$readTargetRegions = function(){
  # Read target region annotations
  cat('Reading target regions\n')

  capture.region = import.bed(meta.env$target.region.file)

  # Merge overlapping regions
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

  #if you want to change the minimum coverage, change the minCov setting in the line below
  methylDataRaw.filter10 = filterByCov(data.env$methylDataRaw, minCov=10, global=F)

  # subset the methylDataRaw object that overlaps the target regions
  methylDataRaw.filter10.rk = subsetByOverlaps(methylDataRaw.filter10, capture.region)
  assign( "methylDataRaw.filter10.rk", methylDataRaw.filter10.rk, envir=data.env)
  rm(methylDataRaw, envir=data.env) # no longer needed
  rm(methylDataRaw.filter10, envir=data.env) # no longer needed
}

func.env$defineCpGClusters = function(){
  cat('Defining Cpg clusters\n')

  # Within a BSraw object clusterSites searches for agglomerations of CpG sites
  # across all samples. In a first step the data is reduced to CpG sites covered
  # in at least perc.samples samples. In a second step regions are detected where
  # not less than min.sites frequently covered CpG sites are sufficiantly close to each other
  # (max.dist). Note, that the frequently covered CpG sites are considered to define
  # the boundaries of the CpG clusters only. For the subsequent analysis the 
  # methylation data of all CpG sites within these clusters are used.
  data.clustered = BiSeq::clusterSites(data.env$methylDataRaw.filter10.rk,
    groups=colData(data.env$methylDataRaw.filter10.rk)$Group,
    perc.samples=0.25, 
    min.site=5, 
    max.dist=100, 
    mc.cores=meta.env$ncores)

  # In the smoothing step CpG sites with high coverages get high weights. To
  # reduce bias due to unusually high coverages we limit the coverage.
  # Default is the 90% quantile
  ind.cov = totalReads(data.clustered) > 0
  quant   = quantile(totalReads(data.clustered)[ind.cov],0.9)
  data.limited = limitCov(data.clustered, maxCov=quant)

  #   We then smooth the methylation values of CpG sites within the clusters with
  # the given bandwidth (h)
  cat("Smoothing methylation ...\n")

  # Using multiple cores through slurm does not work here - a single core is spread
  # across each instance.
  predictedMeth = predictMeth(data.limited, h=meta.env$dmr.bandwidth, mc.cores=meta.env$ncores)

  assign( "predictedMeth", predictedMeth, envir=data.env)

  rm(methylDataRaw.filter10.rk, envir=data.env)  #no longer needed

  # makeSmoothPlot = function(sample.name){
  #   png(file=paste0(tempImagePath,"Smoothed methylation", sample.name,".png"),width = 480, height=480)
  #   plotMeth(object.raw = methylDataRaw[,rownames(methylDataRaw)==sample.name], object.rel = predictedMeth[,6], region = region,
  #     lwd.lines = 2,
  #     col.points = "blue",
  #     cex = 1.5)
  #   dev.off()
  # }
}

func.env$cleanEnvironment = function(){
  # Remove any unneeded global variables left behind from past runs
  cat("Cleaning environment...\n")
  rm(methylDataRaw, envir=data.env)
  rm(methylDataRaw.filter10, envir=data.env)
  rm(methylDataRaw.filter10.rk, envir=data.env)
  rm(capture.region, envir=data.env)
}

func.env$runBetaRegression = function(){
  
  tempImagePath = paste0(dirname(meta.env$sample.info.file),"/", meta.env$temp.image.folder, "/")

  cat("Running beta regression across", meta.env$ncores, "instances ...\n")

  betaRegressionOnChromosome = function(chunk.name){

      inputFile   = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".Rdata")
      resultsFile = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".beta.Rdata")

      runRegression = function(obj){
        cat("Chunk", chunk.name, "with", nrow(obj), "rows\n")

        # To detect the CpG sites where the DNA methylation differs between case
        # and control samples we model the methylation within a beta regression with
        # the group as explanatory variable and use the Wald test to test if there is a
        # group effect
        # By setting type = "BR" the maximum likelihood with bias reduction is
        # called.  This is especially useful, when the sample size is small.

        ptm = proc.time()
        result = betaRegression(formula=~Group, 
                        link='probit', 
                        type='BR', 
                        object=obj,
                        mc.cores=meta.env$ncores)

        saveRDS(result, file = resultsFile)
        PrintTimeTaken(proc.time() - ptm)
        return()
      }

      loadChromosomeChunk = function(inputFile) {
        # Load the saved data chunk with the given name and run regression
        cat("Loading data chunk", chunk.name, "\n")
        return(readRDS(inputFile))
      }

      # Check if results exist
      if(!file.exists(resultsFile)){
        runRegression(loadChromosomeChunk(inputFile))
      }
  }

  saveChunks = function(chunk.data, chunk.name){
    # Save the given data to the temp dir
    tempFile = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".Rdata")
    saveRDS(chunk.data, file = tempFile)
  }
    
  loadResultsChunk = function(chunk.name) {
    # Load the saved data chunk with the given name and run regression
    tempFile = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".beta.Rdata")
    cat("Loading data chunk for chunk", chunk.name, "\n")
    return(readRDS(tempFile))
  }
   

  # Beta regression takes *forever* if run on a small number of cores. While the task
  # itself is embarassingly parallel, the limit was the memory in use
  # at the point mclapply gets invoked (~44Gb per instance with the methylDataRaw object loaded). Any object not
  # modified during a fork() should share memory between parent and child, so we need to 
  # minimise the size of the data being modified.
  #
  # This approach is to save separate data chunks, and process them serially. Also, it still 
  # takes a while, so we save the results of each chunk to file so we can skip them if 
  # the job gets cancelled and resumed.
  # Based on approach 3 in: https://lcolladotor.github.io/2013/11/14/Reducing-memory-overhead-when-using-mclapply/
  # Given results on nibbler, chunk size of 50000 rows should be ok for more than 40 cores.
  predictedMeth.split = split(data.env$predictedMeth, ceiling(seq_along(data.env$predictedMeth)/meta.env$chunk.size), drop=TRUE)

  chunk.names = c(1:length(predictedMeth.split))

  # Save out each chunk data
  cat("Saving", length(chunk.names) , "data chunks to", meta.env$tempImagePath, "\n")
  invisible(mapply(saveChunks, predictedMeth.split, chunk.names))

  rm(predictedMeth.split)
  gc()

  # Process each chunk in turn
  invisible(lapply(chunk.names, betaRegressionOnChromosome ))

  cat("Loading results\n")
  betaResults.split   = lapply(chunk.names, loadResultsChunk )
  
  # betaResults.split   = lapply(predictedMeth.split, betaRegressionOnChromosome )
  betaResults = do.call("rbind", betaResults.split)
  assign( "betaResults", betaResults, envir=data.env)
}

func.env$runNullBetaRegression = function(){
  cat("Running resampled beta model regression for null hypothesis across", meta.env$ncores, "instances ...\n")

  tempImagePath = paste0(dirname(meta.env$sample.info.file),"/", meta.env$temp.image.folder, "/")

  betaRegressionOnChunk = function(chunk.name){

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

      inputFile   = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".null.Rdata")
      resultsFile = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".null.beta.Rdata")

      runRegression = function(obj){
        cat("Chunk", chunk.name, "with", nrow(obj), "rows\n")

        # To detect the CpG sites where the DNA methylation differs between case
        # and control samples we model the methylation within a beta regression with
        # the group as explanatory variable and use the Wald test to test if there is a
        # group effect
        # By setting type = "BR" the maximum likelihood with bias reduction is
        # called.  This is especially useful, when the sample size is small.

        ptm = proc.time()
        result = betaRegression(formula=~group.null, 
                        link='probit', 
                        type='BR', 
                        object=obj,
                        mc.cores=meta.env$ncores)

        saveRDS(result, file = resultsFile)
        PrintTimeTaken(proc.time() - ptm)
        return()
      }

      loadChunk = function(inputFile) {
        # Load the saved data chunk with the given name and run regression
        cat("Loading data chunk", chunk.name, "\n")
        return(readRDS(inputFile))
      }

      # Check if results exist
      if(!file.exists(resultsFile)){
        runRegression(loadChunk(inputFile))
      }
  }

  saveChunks = function(chunk.data, chunk.name){
    # Save the given data to the temp dir
    tempFile = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".null.Rdata")
    saveRDS(chunk.data, file = tempFile)
  }
    
  loadResultsChunk = function(chunk.name) {
    # Load the saved data chunk with the given name and run regression
    tempFile = paste0(meta.env$tempImagePath, "Chunk", chunk.name, ".null.beta.Rdata")
    cat("Loading data chunk for chunk", chunk.name, "\n")
    return(readRDS(tempFile))
  }

  # Since there are only a small number of samples, use an even division.
  # The group null labels should cover the number of samples in your predictmeth object
  # Do not use real group assignments. I.e. have an even split of case and control between
  # samples.

  # Since there are 3 case and 3 control, we need to select 2 and 2.
  # TODO - make dynamic based on group column
  predictedMethNull = data.env$predictedMeth[,c(1, 3, 2, 5)]

  colData(predictedMethNull)$group.null=rep(c(1,2), nrow(colData(predictedMethNull))/2)

  cat("Selected samples for modelling:\n")
  cat(colData(predictedMethNull))

  predictedMethNull.split = split(predictedMethNull, ceiling(seq_along(predictedMethNull)/meta.env$chunk.size), drop=TRUE)
  chunk.names = c(1:length(predictedMethNull.split))

  # Save out each chunk data
  cat("Saving", length(chunk.names) , "data chunks to", meta.env$tempImagePath, "\n")
  invisible(mapply(saveChunks, predictedMethNull.split, chunk.names))

  rm(predictedMethNull.split)
  gc()

  # Process each chunk in turn
  invisible(lapply(chunk.names, betaRegressionOnChunk ))

  cat("Loading results\n")
  betaResultsNull.split   = lapply(chunk.names, loadResultsChunk )
  
  betaResultsNull = do.call("rbind", betaResultsNull.split)
  assign( "betaResultsNull", betaResultsNull, envir=data.env)

  # Estimate the variogram for the Z scores obtained for the resampled data
  vario = makeVariogram(betaResultsNull)

  assign( "vario", vario, envir=data.env)
  # save.image(paste0(tempImagePath, "Null_variogram.RData")) # save workspace

  # Plot the variogram
  png(file=paste0(meta.env$tempImagePath,"Null_variogram.png"),width = 480, height=480)
  plot(vario$variogram$v)
  dev.off()
}

func.env$smoothValues = function(){
  cat('Smoothing\n')
  ## Break, and call new script with desired sill value, or continue with default sill of 1?
  ## sill_value= readline("What is the value of sill to use (0 to 1)? please check graph ")
  sill_value= 1

  # It is necessary to smooth the variogram. Especially for greater h the variogram 
  # tends to oscillate strongly. This is the reason why the default bandwidth 
  # increases with increasing h. Nevertheless, the smoothed variogram may further 
  # increase or decrease after a horizontal part (sill). This is mostly due to the 
  # small number of observations for high distances. To wipe out this bias it is
  # useful to set the smoothed variogram to a fixed value above a certain h, 
  # usually the mean value of the horizontal part. If a smoothed value v.sm is
  # greater than sill for distance h_{range}, this v.sm and all other smoothed 
  # values with h > h_{range} are set to sill. Internally, the function lokerns 
  # from package lokerns is used for smoothing.

  vario.sm = smoothVariogram(data.env$vario, sill=as.numeric(sill_value))
  assign( "vario.sm", vario.sm, envir=data.env)

  # Plot the smoothed variogram
  png(file=paste0(meta.env$tempImagePath,"Smooth_variogram.png"),width = 480, height=480)
  plot(data.env$vario$variogram$v)
  lines(vario.sm$variogram[,c("h", "v.sm")], col = "red", lwd = 1.5)
  dev.off()

  # Replace the pValsList object (which consists of the test results of the
  # resampled data) by the test results of interest (for group effect):
  vario.aux = makeVariogram(betaResults, make.variogram=F)
  vario.sm$pValsList  = vario.aux$pValsList

  # vario.sm now contains the smoothed variogram under the Null hypothesis
  # together with the P values (and Z scores) from the Wald test, that the group
  # has no effect on methylation.  The correlation of the Z scores between two
  # locations in a cluster can now be estimated:
  locCor = estLocCor(vario.sm)

  # We test each CpG cluster for the presence of at least one differentially 
  # methylated location at <q>; what can be interpreted as the size-weighted FDR on
  # clusters:
  clusters.rej     = testClusters(locCor, FDR.cluster=0.1)

  # We then trim the rejected CpG clusters to remove the not differentially
  # methylated CpG sites at <q1>; what can be interpreted as the location-wise FDR:
  clusters.trimmed = trimClusters(clusters.rej, FDR.loc=0.05)

  # We can now define the boundaries of DMRs as rejected CpG sites within which rejected 
  # CpG sites solely are located.  Within the DMRs the distance between neighbored 
  # rejected CpG sites should not exceed max.dist base pairs (usually the same as for 
  # max.dist in clusterSites), otherwise, the DMR is splitted. DMRs are also splitted if 
  # the methylation difference switches from positive to negative, or vice versa, if 
  # diff.dir = TRUE. That way we ensure that within a DMR all CpG sites are hypermethylated, 
  # and hypomethylated respectively.
  DMRs    = findDMRs(clusters.trimmed, max.dist=100, diff.dir=T)
  assign( "DMRs", DMRs, envir=data.env)

  DMRs.df = data.frame(DMRs)
  assign( "DMRs.df", DMRs.df, envir=data.env)
}

func.env$findOverlapsWithGTF = function(){
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
  GTF     = import.gff(meta.env$gtf.file, format="gtf",feature.type="gene")
  bed.kit = import.bed(meta.env$target.region.file)

  assign( "GTF", GTF, envir=data.env)
  assign( "bed.kit", bed.kit, envir=data.env)

  tempoverlap1     = findOverlaps(data.env$DMRs, GTF)
  tempoverlap1.kit = findOverlaps(data.env$DMRs, bed.kit)
  tempoverlap2     = data.env$DMRs[queryHits(tempoverlap1)]
  tempoverlap2.kit = data.env$DMRs[queryHits((tempoverlap1.kit))]

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

  assign( "DMRs.annotated", DMRs.annotated, envir=data.env)
  assign( "DMRs.annotated.kit", DMRs.annotated.kit, envir=data.env)

  write.table(file=paste0(meta.env$tempImagePath, "DMRs_annotated_kit.tsv"),   sep="\t",row.names = F, DMRs.annotated.kit)
  write.table(file=paste0(meta.env$tempImagePath, "DMRs_annotated_genes.tsv"), sep="\t",row.names = F, DMRs.annotated)
}

func.env$plotDMRs = function(){

  outputDir = paste0(dirname(meta.env$sample.info.file),"/DMR_images/")
  if(!dir.exists(outputDir)){
    dir.create(outputDir)
  }

  for (i in 1:length(data.env$DMRs)){
    temp = subsetByOverlaps(predictedMeth, data.env$DMRs[i])
    t2   = rowRanges(temp) %>% as.data.frame
    t2   = cbind(index=as.numeric(rownames(t2)), t2)
    t2$seqnames2 = paste0(t2$seqnames, '_', t2$start)
    df = melt(methLevel(temp))

    annot = rowRanges(temp)
    annot$seqnames2 = paste0(seqnames(annot), '_', start(annot))
    annot = left_join(as.data.frame(annot), data.env$DMRs.annotated,by='seqnames2')
    annot = paste0(annot[complete.cases(annot),]$annotation, collapse= ' \n ')

    colnames(df) = c('index','Sample.Name','methLevel')
    df = left_join(df,as.data.frame(colData(temp)), by='Sample.Name')
    df = left_join(df,t2,by='index')

    positions = rowRanges(subsetByOverlaps(data.env$methylDataRaw.filter10.rk.clustered,data.env$DMRs[i])) %>% as.data.frame
    coverages = totalReads(subsetByOverlaps(data.env$methylDataRaw.filter10.rk.clustered,data.env$DMRs[i])) %>% as.data.frame
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

##################
#
# Run the analysis
#
##################

# The order in which functions should be run
func.env$func.order = c(func.env$readSamples, func.env$readTargetRegions, func.env$defineCpGClusters, func.env$cleanEnvironment, func.env$runBetaRegression,
  func.env$runNullBetaRegression, func.env$smoothValues, func.env$findOverlapsWithGTF, func.env$plotDMRs)
# The corresponding save points
func.env$func.names = c("Read_files.RData", "Filtered_methylation.RData", "Predicted_methylation.RData", "Cleaned.RData", "Beta_regression.RData",
  "Beta_regression_null.RData", "DMRs_found.RData", "Overlaps_with_GTF.RData", "DMRs.RData", "end.RData")

for( i in 1:length(func.env$func.order)){
  thisFile = func.env$func.names[i]
  nextFile = func.env$func.names[i+1]
  thisFunc = func.env$func.order[[i]]
  func.env$runOrSkip(nextFile, thisFile, thisFunc)
}