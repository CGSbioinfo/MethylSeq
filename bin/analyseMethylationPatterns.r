# Perform differential methylation analysis using BiSeq.
# Extended comments are mostly copied from the BiSeq guide:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/BiSeq/inst/doc/BiSeq.pdf

cat("Loading R packages ...\n")

suppressMessages(library(BiSeq))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(parallel))

meta.env = new.env() # parameters for the analysis
data.env = new.env() # data for the ongoing analysis - serialisable
func.env = new.env() # functions defined in this file

meta.env$sample.group.file  = commandArgs(TRUE)[1]
meta.env$cov.file.folder    = commandArgs(TRUE)[2]
meta.env$target.region.file = commandArgs(TRUE)[3] # Can be NA
meta.env$temp.folder.name   = commandArgs(TRUE)[4]
meta.env$ncores             = as.numeric(commandArgs(TRUE)[5])
meta.env$chunk.size         = as.numeric(commandArgs(TRUE)[6]) # number of rows per data chunk in beta regression
meta.env$gtf.file           = commandArgs(TRUE)[7]
meta.env$dmr.bandwidth      = as.numeric(commandArgs(TRUE)[8])
meta.env$sill               = as.numeric(commandArgs(TRUE)[9])
meta.env$tissue             = commandArgs(TRUE)[10] # choose the tissue to subset by

#' This function loads the external functions source file.
#' Detects the directory containing this script, and loads the function
#' file from the same directory. Allows this script to be invoked from a
#' different working directory without hardcoding paths.
#' @title Load external functions
#' @export
loadFunctionsFile = function(){
  cat("Loading functions ...\n")
  file.arg.name = "--file=" # implicit argument when RScript is invoked
  script.name   = sub(file.arg.name, "", commandArgs()[grep(file.arg.name, commandArgs())])
  script.folder = dirname(script.name)
  script.to.load = paste(sep="/", script.folder, "functions.r")
  source(script.to.load)
}

loadFunctionsFile()

##################
#
# Check validity of 
# file paths
#
##################

quit.if.not.exists(meta.env$sample.group.file, "Sample group file") 
quit.if.not.exists(meta.env$cov.file.folder, "Coverage file folder")
quit.if.not.exists(meta.env$gtf.file, "Gene annotation file")

if(meta.env$target.region.file=='NA'){
  meta.env$target.region.file=NA
}

if(!is.na(meta.env$target.region.file)){
  quit.if.not.exists(meta.env$target.region.file, "Target region file")
}

meta.env$temp.folder.name = slashTerminate(meta.env$temp.folder.name)

meta.env$temp.image.path = paste0(dirname(meta.env$sample.group.file),"/", meta.env$temp.folder.name)
ensure.dir.exists(meta.env$temp.image.path)

meta.env$server = system("hostname", intern = TRUE)

meta.env$log.file = paste0(meta.env$temp.image.path, format.Date(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".log.txt")
debug( meta.env$log.file, paste0("Running main on ", meta.env$server))
debug( meta.env$log.file, paste0("Platform: ", sessionInfo()$running ))
debug( meta.env$log.file, paste0("R version: ", getRversion()))
log.pkgs = function(pkg) debug( meta.env$log.file, paste0(pkg$Package, " - ", pkg$Version))
invisible(lapply(sessionInfo()$otherPkgs, log.pkgs))

##################
#
# Define the functions
#
##################

#' Check if the next analysis step output exists. If so, skip 
#' this analysis step. If not, check if this step's output exists
#' and if not, run the analysis.
#' @title Run or skip analysis step
#' @param nextStepTempFile the name of the temp file after this step
#' @param tempFile the name of the temp file for this step
#' @param runFunction the function to be run in this step
#' @export
func.env$runOrSkip = function(nextStepTempFile, tempFile, runFunction){
  if(!file.exists(paste0(meta.env$temp.image.path, nextStepTempFile))){
    tempFile = paste0(meta.env$temp.image.path, tempFile)
    if(file.exists(tempFile)){
      info( meta.env$log.file, paste0("Loading temporary analysis file... ", tempFile))
      runAndTime( f=function(){load(tempFile, envir = globalenv())} )
      info( meta.env$log.file, paste0("Loaded temporary analysis file ", tempFile))
    } else{
      runAndTime(runFunction)
      info( meta.env$log.file, paste0("Saving temporary analysis file... ", tempFile))

      # Only save the data values - ensures updated functions will not be overwritten 
      runAndTime(f=function(){save(data.env, file = tempFile)})
    }
  }
}

#' Read the methylation sample files
#' @title Read samples
#' @export
func.env$readSamples = function(){

  info( meta.env$log.file, paste0("Reading sample groups file '", meta.env$sample.group.file,"'..."))

  sampleGroups = read.csv(meta.env$sample.group.file, sep="\t", header=T, stringsAsFactors=F) %>%
    dplyr::filter(Tissue == meta.env$tissue) %>% dplyr::filter(Type == "Old")
  rownames(sampleGroups) = as.character(sampleGroups$Sample.Name)
  assign( "sampleGroups", sampleGroups, envir=data.env)
  info( meta.env$log.file, paste0("Found ", nrow(sampleGroups), " samples for ", meta.env$tissue))

  # Locate coverage files
  info( meta.env$log.file, "Locating cov files...")

  cov.file.list = list.files(path=meta.env$cov.file.folder,pattern=".cov",full.names = T)
  cov.file.list = grep(paste0(sampleGroups$Sample.Name,collapse='|'),
      cov.file.list, value=TRUE)

  write.samples = function(s) debug( meta.env$log.file, s)
  lapply(cov.file.list, write.samples)
  info( meta.env$log.file, "Reading cov files...")
  
  tryCatch({
    methylDataRaw = readBismark(cov.file.list, sampleGroups)
    info( meta.env$log.file, "Read cov files")
    assign( "methylDataRaw", methylDataRaw, envir=data.env)
  }, error = function(e){
    warn( meta.env$log.file, "Error reading cov files")
    warn( meta.env$log.file, e)
  })
  
}

#' Read the target region file and select these regions from methylation data
#' @title Read target regions
#' @export
func.env$readTargetRegions = function(){
  # Read target region annotations
  info( meta.env$log.file, paste0("Reading target region file ", meta.env$target.region.file))

  methylDataRaw.filter10 = filterByCov(data.env$methylDataRaw, minCov=10, global=F)

  # Do not filter if no target was provided
  if(is.na(meta.env$target.region.file)){
    methylDataRaw.filter10.rk = methylDataRaw.filter10
  } else {
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

    # subset the methylDataRaw object that overlaps the target regions
    methylDataRaw.filter10.rk = subsetByOverlaps(methylDataRaw.filter10, capture.region)
  }
 
  assign( "methylDataRaw.filter10.rk", methylDataRaw.filter10.rk, envir=data.env) # ignore capture region
  rm(methylDataRaw, envir=data.env) # no longer needed
  rm(methylDataRaw.filter10, envir=data.env) # no longer needed
}

#' Cluster CpGs and smooth the methylation
#' @title Define clusters
#' @export
func.env$defineCpGClusters = function(){
  info( meta.env$log.file, "Defining Cpg clusters")

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

  assign( "methylDataRaw.filter10.rk.clustered", data.clustered, envir=data.env)

  # In the smoothing step CpG sites with high coverages get high weights. To
  # reduce bias due to unusually high coverages we limit the coverage.
  # Default is the 90% quantile
  ind.cov = totalReads(data.clustered) > 0
  quant   = quantile(totalReads(data.clustered)[ind.cov],0.9)
  data.limited = limitCov(data.clustered, maxCov=quant)

  #   We then smooth the methylation values of CpG sites within the clusters with
  # the given bandwidth (h)
  info( meta.env$log.file, "Smoothing methylation...")

  predictedMeth = predictMeth(data.limited, h=meta.env$dmr.bandwidth, mc.cores=6)
  assign( "predictedMeth", predictedMeth, envir=data.env)

  rm(methylDataRaw.filter10.rk, envir=data.env)  #no longer needed
}

#' Run beta regression on a data chunk
#' @title Run regression on chunk
#' @param chunk.name the name of the data chunk to process
#' @param is.null T if this is the null model regression, F otherwise
#' @export
func.env$betaRegressionOnChunk = function(chunk.name, is.null){

    null.string = ifelse(is.null, ".null", "")

    inputFile   = paste0(meta.env$temp.image.path, "Chunk", chunk.name, null.string, ".Rdata")
    lockFile    = paste0(meta.env$temp.image.path, "Chunk", chunk.name, null.string, ".lck")
    resultsFile = paste0(meta.env$temp.image.path, "Chunk", chunk.name, null.string, ".beta.Rdata")

    runRegression = function(obj){
      info( meta.env$log.file, paste("Chunk", chunk.name, "with", nrow(obj), "rows", sep=" "))


      ptm = proc.time()

      if(is.null){
        result = betaRegression(formula=~group.null, 
                      link='probit', 
                      type='BR', 
                      object=obj,
                      mc.cores=meta.env$ncores)
      } else {
        result = betaRegression(formula=~Group, 
                        link='probit', 
                        type='BR', 
                        object=obj,
                        mc.cores=meta.env$ncores)
      }
      saveRDS(result, file = resultsFile)
      time = printTimeTaken(proc.time() - ptm)
      info(meta.env$log.file, paste(meta.env$server, meta.env$ncores, chunk.name, nrow(obj), time, sep = "\t"))
      return()
    }

    loadChromosomeChunk = function(inputFile) {
      # Load the saved data chunk with the given name and run regression
      info( meta.env$log.file, paste0("Loading data chunk ", chunk.name))
      return(readRDS(inputFile))
    }

    # Check if results exist
    if(!file.exists(resultsFile)){

      if(!file.exists(lockFile)){

        file.create(lockFile)
        runRegression(loadChromosomeChunk(inputFile))
        file.remove(lockFile)
      }
    }
}

#' Run beta regression on the samples
#' @title Run beta regression
#' @export
func.env$runBetaRegression = function(){
  
  info( meta.env$log.file, paste("Running beta regression across", meta.env$ncores, "instances...", sep=" "))

  # To detect the CpG sites where the DNA methylation differs between case
  # and control samples we model the methylation within a beta regression with
  # the group as explanatory variable and use the Wald test to test if there is a
  # group effect
  # By setting type = "BR" the maximum likelihood with bias reduction is
  # called.  This is especially useful, when the sample size is small.

  saveChunks = function(chunk.data, chunk.name){
    # Save the given data to the temp dir
    tempFile = paste0(meta.env$temp.image.path, "Chunk", chunk.name, ".Rdata")
    saveRDS(chunk.data, file = tempFile)
  }
    
  loadResultsChunk = function(chunk.name) {
    # Load the saved data chunk with the given name and run regression
    tempFile = paste0(meta.env$temp.image.path, "Chunk", chunk.name, ".beta.Rdata")
    info( meta.env$log.file, paste("Loading data chunk for chunk", chunk.name, sep=" "))
    return(readRDS(tempFile))
  }
   

  # Beta regression takes *forever* if run on a small number of cores. While the task
  # itself is embarassingly parallel, a limit is the memory in use at the point
  # mclapply gets invoked. Any object not modified during a UNIX fork() should share 
  # memory between parent and child, but in this case the methylDataRaw object is
  # copied if present (and it is the largest object in memory by far).
  #
  # 1) Minimise the size of the data being modified by removing unneeded variables from data.env.
  #    This allows more parallel instances within the node.
  # 
  # 2) Serialise the data into chunks. This allows them to be processed by multiple nodes.
  #    Based on approach 3 in: https://lcolladotor.github.io/2013/11/14/Reducing-memory-overhead-when-using-mclapply/
  #
  # Within a node, each chunk will take as long as the longest core time. When setting the chunk size, 
  # consider the tradeoff between loading a new chunk, and the cpu downtime waiting for the last core to complete.
  # In tests, a chunk size of 25000 was optimal for a dataset of ~250k rows of 6 samples at both 30 and 36 cores.
  # Before the chunking was added, the script ran for seven days on 6 cores before the job was killed 
  # when the cluster was shut down. Benchmarking 20 cores or lower has not been completed.
  nchunks     = ceiling(length(data.env$predictedMeth)/meta.env$chunk.size)
  chunk.names = c(1:nchunks)
  info( meta.env$log.file, paste("Regression running on", nchunks, "chunks", sep=" "))
  
  chunksSaved = function(){
    # This is a very basic test - it checks if the number of chunks present matches the expected
    # number of chunks. This will fail if chunk size is decreased between runs.
    length(list.files(path=meta.env$temp.folder.name, pattern="Chunk\\d+.Rdata",full.names = T))==nchunks
  }

  if(chunksSaved()){
    info( meta.env$log.file, paste("Reading existing", nchunks , "data chunks", sep=" "))
  } else {
    info( meta.env$log.file, paste("Saving", nchunks , "data chunks to", meta.env$temp.image.path, sep=" "))
    predictedMeth.split = split(data.env$predictedMeth, ceiling(seq_along(data.env$predictedMeth)/meta.env$chunk.size), drop=TRUE)
    invisible(mapply(saveChunks, predictedMeth.split, chunk.names))
  }

  rm(predictedMeth.split)
  info( meta.env$log.file, "Parallel script can now be invoked")

  # Process each chunk in turn
  invisible(lapply(chunk.names, func.env$betaRegressionOnChunk, F ))

  # At this point, you can invoke the parallelBetaRegression.r script on other node(s).
  # srun <options> /usr/bin/Rscript bin/parallelBetaRegression.r <temp.folder.name> <ncores> <is.null.regression>
  # Example:
  # srun -w calculon -c 36 /usr/bin/Rscript bin/parallelBetaRegression.r tmpImages_30core_25k_chunk/ 30 F

  # Wait for all lock files to be removed - halts until other nodes finish processing
  locksRemoved = function(){
    files = list.files(path=meta.env$temp.folder.name, pattern="Chunk\\d+.lck",full.names = T)
    length(files)==0
  }
  info( meta.env$log.file, "Waiting for other nodes to complete")
  while (!locksRemoved()) {
    Sys.sleep(1)
  }

  info( meta.env$log.file, "Loading results")
  betaResults.split   = lapply(chunk.names, loadResultsChunk )
  betaResults = do.call("rbind", betaResults.split)
  assign( "betaResults", betaResults, envir=data.env)

  # Remove the chunk files
  file.list = list.files(path=meta.env$temp.folder.name, pattern="Chunk\\d+.*.Rdata", full.names = T)
  file.remove(file.list)
}

#' Create a null model and run the null regression on the samples
#' @title Run null regression
#' @export
func.env$runNullBetaRegression = function(){
  info( meta.env$log.file, paste("Running resampled beta model regression for null hypothesis across", meta.env$ncores, "instances ...", sep=" "))

  #   # The aim is to detect CpG clusters containing at least one differentially methylated location.
  #   # To do so the P values p from the Wald tests are transformed to Z scores which are normally 
  #   # distributed under Null hypothesis (no group effect). As cluster test statistic a standardized
  #   # Z score average is used. To estimate the standard deviation of the Z scores we have to estimate 
  #   # the correlation and hence the variogram of methylation between two CpG sites within a cluster.
  #   # The estimation of the standard deviation requires that the distribution of the Z scores 
  #   # follows a standard normal distribution.  However, if methylation in both groups differs for 
  #   # many CpG sites the density distribution of P values shows a peak near 0.  To ensure that the 
  #   # P values are roughly uniformly distributed to get a variance of the Z scores that is Gaussian
  #   # with variance 1 we recommend to estimate the variogram (and hence the correlation of Z scores)
  #   # under the null hypothesis. To do so we model the beta regression again for resampled data.

  saveChunks = function(chunk.data, chunk.name){
    # Save the given data to the temp dir
    tempFile = paste0(meta.env$temp.image.path, "Chunk", chunk.name, ".null.Rdata")
    saveRDS(chunk.data, file = tempFile)
  }
    
  loadResultsChunk = function(chunk.name) {
    # Load the saved data chunk with the given name and run regression
    tempFile = paste0(meta.env$temp.image.path, "Chunk", chunk.name, ".null.beta.Rdata")
    info( meta.env$log.file, paste("Loading data chunk", chunk.name, sep=" "))
    return(readRDS(tempFile))
  }
  
  #' Choose the samples for null modelling. There must be an even distribution
  #' of samples into each null group. For example, if there are 3 case and 3 control, 
  #' select 2 of each.
  #' @title Choose null samples
  createNullSample = function(){

    sample.groups  = data.env$sampleGroups %>% distinct(Group) %>% unlist()
    sample.numbers = data.env$sampleGroups %>% group_by(Group) %>% summarise(n= n()) %>% select(n)
    sample.count   = min(ceiling(sample.numbers/2)) # Number of samples from each group to select 

    info( meta.env$log.file, paste("Choosing", sample.count, "samples from each group"))
    getNullSamplesForGroup = function(i){
      data.env$sampleGroups %>% filter(Group == sample.groups[i]) %>% head(sample.count) %>% select(Sample.Name) %>% unlist()
    }

    group1 = getNullSamplesForGroup(1)
    group2 = getNullSamplesForGroup(2)

    nullSamples = data.env$predictedMeth[, c(group1,group2)]

    colData(nullSamples)$group.null=rep(c(1,2), nrow(colData(nullSamples))/2)

    columns = as.data.frame(colData(nullSamples)[, c("Sample.Name", "Group", "group.null")])
    printNullGroup = function(i){
        null = columns %>% filter(group.null == i) %>% select(Sample.Name) %>% unlist()
        info( meta.env$log.file, paste0("Chosen samples for null group ", i, ": ", paste(null, collapse=", ")))
    }

    printNullGroup(1)
    printNullGroup(2)

    nullSamples
  }

  predictedMethNull = createNullSample()

  info( meta.env$log.file, "Selected samples for null modelling")

  nchunks     = ceiling(length(predictedMethNull)/meta.env$chunk.size)
  chunk.names = c(1:nchunks)
  
  chunksSaved = function(){
    # This is a very basic test - it checks if the number of chunks present matches the expected
    # number of chunks. This will fail if chunk size is decreased between runs.
    length(list.files(path=meta.env$temp.folder.name, pattern="Chunk\\d+.null.Rdata",full.names = T))==nchunks
  }

  if(chunksSaved()){
    info( meta.env$log.file, paste("Reading existing", nchunks , "data chunks", sep=" "))
  } else {
    info( meta.env$log.file, paste("Saving", nchunks , "data chunks to", meta.env$temp.image.path, sep=" "))
    predictedMethNull.split = split(predictedMethNull, ceiling(seq_along(predictedMethNull)/meta.env$chunk.size), drop=TRUE)
    invisible(mapply(saveChunks, predictedMethNull.split, chunk.names))
  }

  rm(predictedMethNull.split)
  info( meta.env$log.file, "Parallel script can now be invoked")

  # Process each chunk in turn
  invisible(lapply(chunk.names, func.env$betaRegressionOnChunk, T ))

  # At this point, you can invoke the parallelBetaRegression.r script on other node(s).
  # srun <options> /usr/bin/Rscript bin/parallelBetaRegression.r <temp.folder.name> <ncores> <is.null.regression>
  # Example:
  # srun -w calculon -c 20 /usr/bin/Rscript bin/parallelBetaRegression.r tmpImages/ 20 T

  # Wait for all lock files to be removed - halts until other nodes finish processing
  locksRemoved = function(){
    files = list.files(path=meta.env$temp.folder.name, pattern="Chunk\\d+.null.lck",full.names = T)
    length(files)==0
  }
  info( meta.env$log.file, "Waiting for other nodes to complete")
  while (!locksRemoved()) {
    Sys.sleep(5)
  }

  info( meta.env$log.file, "Loading results")  
  betaResultsNull.split   = lapply(chunk.names, loadResultsChunk )
  
  betaResultsNull = do.call("rbind", betaResultsNull.split)
  assign( "betaResultsNull", betaResultsNull, envir=data.env)

  info( meta.env$log.file, "Creating variogram...")  
  # Estimate the variogram for the Z scores obtained for the resampled data
  vario = makeVariogram(betaResultsNull)

  assign( "vario", vario, envir=data.env)

  # Plot the variogram
  png(file=paste0(meta.env$temp.image.path,"Null_variogram.png"),width = 480, height=480)
  plot(vario$variogram$v)
  dev.off()

  # Remove the chunk files
  file.list = list.files(path=meta.env$temp.folder.name, pattern="Chunk\\d+.*.Rdata", full.names = T)
  file.remove(file.list)
}

#' Quit the script so the variogram can be examined.
#' @title Save and quit 
#' @export
func.env$quitForVar = function(){
    save(data.env, file = paste0(meta.env$temp.image.path, "Save_6.Rdata"))
    info( meta.env$log.file, "Quitting so you can examine the variogram")
    info( meta.env$log.file, "Adjust the sill value in the analysis info and rerun step 6")  
    quit(save="no", status=0)
}

#' Smooth the variogram and identify DMRs
#' @title Smooth variogram
#' @export
func.env$smoothValues = function(){
  
  sill_value= 0.30
  info( meta.env$log.file, paste("Smoothing with sill value",sill_value,"..."))  
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
  png(file=paste0(meta.env$temp.image.path,"Smooth_variogram.png"),width = 480, height=480)
  plot(data.env$vario$variogram$v)
  lines(vario.sm$variogram[,c("h", "v.sm")], col = "red", lwd = 1.5)
  dev.off()

  # Replace the pValsList object (which consists of the test results of the
  # resampled data) by the test results of interest (for group effect):
  vario.aux = makeVariogram(data.env$betaResults, make.variogram=F)
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

  info( meta.env$log.file, 'Trimming clusters')  
  # Error in IRanges with negative width. Remove row 50345
  # idx <- GenomicRanges:::get_out_of_bound_index(ext_grn)
  # if (length(idx) != 0L)
  #     ext_grn <- ext_grn[-idx]

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
  info( meta.env$log.file, 'Finding DMRs')
  DMRs    = findDMRs(clusters.trimmed, max.dist=100, diff.dir=T)
  assign( "DMRs", DMRs, envir=data.env)

  DMRs.df = data.frame(DMRs)
  assign( "DMRs.df", DMRs.df, envir=data.env)
}

#' Find overlaps between the DMRs and genes, and export tables.
#' @title Find genes overlapping DMRs
#' @export
func.env$findOverlapsWithGTF = function(){
  info( meta.env$log.file, "Finding DMR overlaps with genes") 
  outputDir = paste0(dirname(meta.env$sample.group.file),"/Report/figure/DMR_results/")
  if(!dir.exists(outputDir)){
    info( meta.env$log.file, paste("Creating output directory", outputDir ))
    dir.create(outputDir)
  }

  merge_hits=function(x){
    data.frame(seqnames=x$seqnames[1], start=x$start[1], end=x$end[1], width=x$width[1],
                         median.p=x$median.p[1],
                         median.meth.group1=x$median.meth.group1[1],
                         median.meth.group2=x$median.meth.group2[1],
                         median.meth.diff=x$median.meth.diff[1],
                         annotation=paste(x$annotation, collapse=' | '))
  }



  GTF     = import.gff(meta.env$gtf.file, format="gtf",feature.type="gene")
  tempoverlap1     = findOverlaps(data.env$DMRs, GTF)
  tempoverlap2     = data.env$DMRs[queryHits(tempoverlap1)]
  
  tempoverlap2$annotation=paste(GTF[subjectHits(tempoverlap1)]$gene_id,
                                GTF[subjectHits(tempoverlap1)]$type,
                                GTF[subjectHits(tempoverlap1)]$gene_biotype,
                                GTF[subjectHits(tempoverlap1)]$gene_name, sep='; ')

  tempoverlap2     = as.data.frame(tempoverlap2)
  tempoverlap2$seqnames2     = paste0(tempoverlap2$seqnames,'_', tempoverlap2$start)
  
  DMRs.annotated     = tempoverlap2 %>% group_by(seqnames2) %>% do(merge_hits(.)) %>% as.data.frame
 
  write.table(file=paste0(outputDir, "DMRs_overlapping_annotated_genes.tsv"), sep="\t",row.names = F, DMRs.annotated)

  if(!is.na(meta.env$target.region.file)){
    bed.kit = import.bed(meta.env$target.region.file)
    tempoverlap1.kit = findOverlaps(data.env$DMRs, bed.kit)
    tempoverlap2.kit = data.env$DMRs[queryHits((tempoverlap1.kit))]
    tempoverlap2.kit$annotation = bed.kit[subjectHits(tempoverlap1.kit)]$name
    tempoverlap2.kit = data.frame(tempoverlap2.kit)
    tempoverlap2.kit$seqnames2 = paste0(tempoverlap2.kit$seqnames,'_', tempoverlap2.kit$start)
    DMRs.annotated.kit = tempoverlap2.kit %>% group_by(seqnames2) %>% do(merge_hits(.)) %>% as.data.frame
    write.table(file=paste0(outputDir, "DMRs_overlapping_target_regions.tsv"),   sep="\t",row.names = F, DMRs.annotated.kit)
  }

  info( meta.env$log.file, "Plotting DMRs") 
  
  if(!dir.exists(outputDir)){
    dir.create(outputDir)
  }

  plotDMR = function(dmr.data){

    temp = subsetByOverlaps(data.env$predictedMeth, dmr.data)
    t2   = rowRanges(temp) %>% as.data.frame
    t2   = cbind(index=as.numeric(rownames(t2)), t2)
    t2$seqnames2 = paste0(t2$seqnames, '_', t2$start)
    df = melt(methLevel(temp))

    annot = rowRanges(temp)
    annot$seqnames2 = paste0(seqnames(annot), '_', start(annot))
    annot = left_join(as.data.frame(annot), DMRs.annotated,by='seqnames2')
    annot = paste0(annot[complete.cases(annot),]$annotation, collapse= ' \n ')

    colnames(df) = c('index','Sample.Name','methLevel')
    df = left_join(df,as.data.frame(colData(temp)), by='Sample.Name')
    df = left_join(df,t2,by='index')

    positions = rowRanges(subsetByOverlaps(data.env$methylDataRaw.filter10.rk.clustered,  dmr.data)) %>% as.data.frame
    coverages = totalReads(subsetByOverlaps(data.env$methylDataRaw.filter10.rk.clustered, dmr.data)) %>% as.data.frame
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

  invisible(lapply(data.env$DMRs, plotDMR))
}

##################
#
# Run the analysis
#
##################

# The order in which functions should be run
func.env$func.order = c(func.env$readSamples, func.env$readTargetRegions, func.env$defineCpGClusters, func.env$runBetaRegression,
  func.env$runNullBetaRegression, func.env$quitForVar, func.env$smoothValues, func.env$findOverlapsWithGTF)

# The corresponding save points
func.env$func.names = c(paste0(meta.env$tissue, "_save_", 1:8, ".RData"), paste0(meta.env$tissue,"_save_9_end.RData"))

for( i in 1:length(func.env$func.order)){
  thisFile = func.env$func.names[i]
  nextFile = func.env$func.names[i+1]
  thisFunc = func.env$func.order[[i]]
  func.env$runOrSkip(nextFile, thisFile, thisFunc)
}

info( meta.env$log.file, "BiSeq analysis complete") 