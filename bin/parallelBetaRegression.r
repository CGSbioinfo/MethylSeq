# Carry out beta regression on existing data chunk files. This can be used in conjunection
# with analyseMethylationPatterns.r to speed up the regression step
suppressMessages(library(BiSeq))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))

meta.env = new.env() # parameters for the analysis

meta.env$temp.image.folder  = commandArgs(TRUE)[1]
meta.env$ncores             = as.numeric(commandArgs(TRUE)[2])
meta.env$chunk.size         = as.numeric(commandArgs(TRUE)[3])

cat("Running parallel beta regression across", meta.env$ncores, "instances ...\n")

betaRegressionOnChromosome = function(chunk.name){

    inputFile   = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".Rdata")
    lockFile    = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".lck") # Allow multiple nodes
    resultsFile = paste0(meta.env$temp.image.folder, "Chunk", chunk.name, ".beta.Rdata")

    runRegression = function(obj){
      cat("Chunk", chunk.name, "with", nrow(obj), "rows\n")

      ptm = proc.time()
      result = betaRegression(formula=~Group, 
                      link='probit', 
                      type='BR', 
                      object=obj,
                      mc.cores=meta.env$ncores)

      saveRDS(result, file = resultsFile)
      return()
    }

    loadChromosomeChunk = function(inputFile) {
      # Load the saved data chunk with the given name and run regression
      cat("Loading data chunk", chunk.name, "\n")
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

file.list = list.files(path=meta.env$temp.image.folder, pattern="Chunk\\d+.Rdata",full.names = T)
chunk.names = c(1:length(file.list))

# Process each chunk in turn
invisible(lapply(chunk.names, betaRegressionOnChromosome ))
