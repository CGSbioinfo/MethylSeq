# Carry out beta regression on existing data chunk files. This can be used in conjunection
# with analyseMethylationPatterns.r to speed up the regression step
#
# Params:
#   1 - the temporary .Rdata directory
#   2 - the number of cores to use
#   3 - run standard regression (F) or null regression (T)
meta.env = new.env() # parameters for the analysis

meta.env$temp.image.folder  = commandArgs(TRUE)[1]
meta.env$ncores             = as.numeric(commandArgs(TRUE)[2])
meta.env$isnull             = commandArgs(TRUE)[3]
meta.env$log.file           = commandArgs(TRUE)[4]

#' This function loads the external functions file.
#' @title Load external functions
#' @export
loadFunctionsFile = function(){
  cat("Loading functions ...\n")
  file.arg.name = "--file="
  script.name   = sub(file.arg.name, "", commandArgs()[grep(file.arg.name, commandArgs())])
  script.folder = dirname(script.name)
  script.to.load = paste(sep="/", script.folder, "functions.r")
  source(script.to.load)
}

loadFunctionsFile()

install.missing(packages=c("dplyr", "parallel", biopackages=c("BiSeq")))

meta.env$server = system("hostname", intern = TRUE)

# meta.env$log.file = paste0(meta.env$temp.image.folder, "log.txt")
info( meta.env$log.file, paste("Running parallel on", meta.env$server))

meta.env$temp.image.folder = slashTerminate(meta.env$temp.image.folder)
quit.if.not.exists(meta.env$temp.image.folder, "Temp analysis folder")

##################
#
# Run
#
##################

betaRegressionOnChromosome = function(chunk.name){

    inputFile   = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".Rdata")
    lockFile    = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".lck") # Allow multiple nodes
    resultsFile = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".beta.Rdata")

    runRegression = function(obj){
      info( meta.env$log.file, paste("Chunk", chunk.name, "with", nrow(obj), "rows"))
      # cat("Chunk", chunk.name, "with", nrow(obj), "rows\n")
      ptm = proc.time()
      result = betaRegression(formula=~Group, 
                      link='probit', 
                      type='BR', 
                      object=obj,
                      mc.cores=meta.env$ncores)

      saveRDS(result, file = resultsFile)
      time = printTimeTaken(proc.time() - ptm)
      debug( meta.env$log.file, paste(meta.env$server, meta.env$ncores, chunk.name, nrow(obj), time, sep="\t"))
      return()
    }

    loadChromosomeChunk = function(inputFile) {
      # Load the saved data chunk with the given name and run regression
      info( meta.env$log.file, paste("Loading data chunk", chunk.name))
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

betaNullRegressionOnChromosome = function(chunk.name){

    inputFile   = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".null.Rdata")
    resultsFile = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".null.beta.Rdata")
    lockFile    = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".null.lck")

    runRegression = function(obj){
      info( meta.env$log.file, paste("Chunk", chunk.name, "with", nrow(obj), "rows"))
      ptm = proc.time()
      result = betaRegression(formula=~group.null, 
                      link='probit', 
                      type='BR', 
                      object=obj,
                      mc.cores=meta.env$ncores)

      saveRDS(result, file = resultsFile)
      time = printTimeTaken(proc.time() - ptm)
      debug( meta.env$log.file, paste(meta.env$server, meta.env$ncores, chunk.name, nrow(obj), time, sep="\t"))
      return()
    }

    loadChromosomeChunk = function(inputFile) {
      # Load the saved data chunk with the given name and run regression
      info( meta.env$log.file, paste("Loading data chunk", chunk.name))
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

pattern     = ifelse(meta.env$isnull, "Chunk\\d+.null.Rdata", "Chunk\\d+.Rdata")
file.list   = list.files(path=meta.env$temp.image.folder, pattern=pattern,full.names = T)
chunk.names = c(1:length(file.list))

if(meta.env$isnull){
  info( meta.env$log.file, paste("Running parallel null beta regression across", meta.env$ncores, "instances..."))
  invisible(lapply(chunk.names, betaNullRegressionOnChromosome ))

} else {
  info( meta.env$log.file, paste("Running parallel beta regression across", meta.env$ncores, "instances..."))
  invisible(lapply(chunk.names, betaRegressionOnChromosome ))
}


