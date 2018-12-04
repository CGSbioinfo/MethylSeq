# Test pairwise loading of samples using BiSeq::readBismark
# If there are no errors, the samples are compatible
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(BiSeq))

source("src/functions.R")

log.file = "logs/pairwise.log"
done.file = "logs/pairwise.done"

sample.group.file = "sample_groups.txt"
cov.file.folder   = "methylExtraction/"

debug( log.file, paste0("Running main on ", system("hostname", intern = TRUE)))
debug( log.file, paste0("Platform: ", sessionInfo()$running ))
debug( log.file, paste0("R version: ", getRversion()))
log.pkgs = function(pkg) debug( log.file, paste0(pkg$Package, " - ", pkg$Version))
invisible(lapply(sessionInfo()$otherPkgs, log.pkgs))

sampleGroups = read.csv(sample.group.file, sep="\t", header=T, stringsAsFactors=F) %>%
  dplyr::filter(Tissue == "Cortex")
rownames(sampleGroups) = as.character(sampleGroups$Sample.Name)

testNames = expand.grid(sampleGroups$Sample.Name, sampleGroups$Sample.Name, sampleGroups$Sample.Name, stringsAsFactors=F) %>% 
  dplyr::filter(Var1!=Var2 & Var2!=Var3 & Var1!=Var3)
  
testNames = as.data.frame(t(apply(testNames, 1, sort))) 
names(testNames) = c("v1", "v2", "v3") 
testNames = testNames %>% 
  mutate(Comb = paste(v1, v2, v3, sep = '|')) %>%
  dplyr::select(Comb) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Comb) %>%
  tidyr::separate(Comb, c("Var1", "Var2", "Var3"), sep="\\|", remove=T)

info(log.file, paste0("Found ", nrow(testNames), " combinations to test"))

read.file = function(samples){
  # Locate coverage files
  cov.file.list = list.files(path=cov.file.folder, pattern=".cov", full.names = T)
  cov.file.list = grep(paste0(c(as.character(samples$Sample.Name)),collapse='|'),
                       cov.file.list, value=TRUE)
  methylDataRaw = readBismark(cov.file.list, samples)
}

subset.samples = function(s1, s2, s3){
  tryCatch({
    info(log.file, paste("Testing:", s1, s2, s3, sep="\t"))
    filt = sampleGroups %>% dplyr::filter(Sample.Name %in% c(s1, s2, s3))
    read.file(filt)
    cat(paste(s1, s2, s3, sep="\t"), "\n", file=done.file, append=TRUE)
  }, error = function(e){
    warn(log.file, "Error in", s1, "and", s2, " and ", s3,"\n")
    warn(log.file, e)
  })
}

invisible(mapply(subset.samples, testNames$Var1, testNames$Var2, testNames$Var3))


