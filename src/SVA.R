# Use SVA to account for cell composition in bisulfite sequencing data
# Based on tutorial at https://akhilesh362.wordpress.com/
# and the SVA manual https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

# Run as:
# $ RScript SVA.R <sample_group_file> <cov.folder> <temp.folder> <gtf.file> <tissue>

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

cat("Loading R packages ...\n")
install.missing(packages    = c("MASS", "ggplot2", "dplyr", "tidyr"),
                biopackages = c("sva", "limma", "BiSeq", "qvalue",  "rtracklayer"),
                repos       = c("https://cran.ma.imperial.ac.uk/", "https://www.stats.bris.ac.uk/R/"))

data.env = new.env() # data for the ongoing analysis
meta.env = new.env() # metadata for the ongoing analysis

meta.env$sample.group.file  = commandArgs(TRUE)[1]
meta.env$cov.file.folder    = commandArgs(TRUE)[2]
meta.env$temp.folder.name   = commandArgs(TRUE)[3]
meta.env$gtf.file           = commandArgs(TRUE)[4]
meta.env$tissue             = commandArgs(TRUE)[5] # choose the tissue to subset by

meta.env$temp.image.path = slash.terminate(paste0(dirname(meta.env$sample.group.file),"/", meta.env$temp.folder.name))
ensure.dir.exists(meta.env$temp.image.path)
meta.env$log.file = paste0(meta.env$temp.image.path, format.Date(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".SVA.log.txt")
log.session.info(meta.env$log.file)

debug(meta.env$log.file, paste("Command line parameters:", paste(commandArgs(TRUE), collapse=" ")))

default.min.beta.diff   = 0.05  # minimum required difference in methylation fraction
default.min.reads       = 5     # minimum number of reads in every sample
default.high.cov.filter = 0.99  # exclude loci with more than this proportion of the maximum read depth in any sample

meta.env$output.dir = paste0(dirname(meta.env$sample.group.file),"/Report/figure/SVA_results/")
ensure.dir.exists(meta.env$output.dir)

#' Load and filter data
#' 
#' Create the table for SVA input.  SVA requires a matrix with features
#' in the rows and samples in the columns. Null values are not permitted.
#' 
#' Filter the incoming data by readcount. Discard the CpGs with the top n% of total read coverage in each sample; this should
#' help remove artefacts due to PCR duplication
load.data = function(){
  min.reads  = default.min.reads
  cov.filter = default.high.cov.filter
  info( meta.env$log.file, paste("Excluding loci with more than", cov.filter,"read coverage in any sample"))
  info( meta.env$log.file, paste("Excluding loci with less than", min.reads, "reads in any sample"))
  
  APPLY_TO_ROWS = 1 # parameter for base::apply
  APPLY_TO_COLS = 2 # parameter for base::apply
  info( meta.env$log.file, paste0("Reading sample groups file '", meta.env$sample.group.file,"'..."))

  data.env$sampleGroups = read.csv(meta.env$sample.group.file, sep="\t", header=T, stringsAsFactors=F) %>%
    dplyr::filter(Tissue == meta.env$tissue)
  rownames(data.env$sampleGroups) = as.character(data.env$sampleGroups$Sample.Name)
  info( meta.env$log.file, paste0("Found ", nrow(data.env$sampleGroups), " samples for ", meta.env$tissue))

  # Locate coverage files
  info( meta.env$log.file, "Locating cov files...")

  cov.file.list = list.files(path=meta.env$cov.file.folder,pattern=".cov",full.names = T)
  cov.file.list = grep(paste0(data.env$sampleGroups$Sample.Name,collapse='|'),
      cov.file.list, value=TRUE)

  debug(meta.env$log.file, cov.file.list)

  info( meta.env$log.file, "Reading cov files...")
  
  tryCatch({
    methylDataRaw = readBismark(cov.file.list, data.env$sampleGroups)
    info( meta.env$log.file, "Read cov files")
    assign( "methylDataRaw", methylDataRaw, envir=data.env)
  }, error = function(e){
    warn( meta.env$log.file, "Error reading cov files")
    warn( meta.env$log.file, e)
    quit(save="no", status=1)
  })
  
  starting.reads = totalReads(methylDataRaw)
  names(starting.reads) = data.env$sampleGroups$Sample.Name
  
  # Plot the read counts for the given loci
  plotReadCounts = function(reads, file.name){
    tryCatch({
      info(meta.env$log.file, "Creating coverage plot")
      
      df = as.data.frame(reads)
      names(df) = data.env$sampleGroups$Sample.Name
      
      gathered = df %>% tidyr::gather(data.env$sampleGroups$Sample.Name, key = "sample_id", value = "total_reads")
      plot.file  = paste0(meta.env$temp.image.path, file.name ,".png")
      g = ggplot(gathered, aes(x=sample_id, y=total_reads))+
        geom_violin()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      debug(meta.env$log.file, paste("Saving coverage image to", plot.file))
      ggsave(file=plot.file, plot=g)
    }, error=function(e){
      warn(meta.env$log.file, "Error making coverage image")
      warn(meta.env$log.file, e)
    }, warning=function(e){
      warn(meta.env$log.file, "Issue making coverage image")
      warn(meta.env$log.file, e)
    })
  }
  
  plotReadCounts(starting.reads, paste0(meta.env$tissue, "_total_reads"))
  
  # In each col, find the nth quantile. Using method 8, as it is median-unbiased no matter the distribution -see ?quantile
  sample.quantiles = apply(starting.reads, APPLY_TO_COLS, function(col) quantile(col, cov.filter, type=8))

  debug(meta.env$log.file, paste0("Sample coverage filter: ", sample.quantiles))

  info(meta.env$log.file, "Creating predicates...")

  # Predicate for min reads per sample
  valid.min =  starting.reads[,1]>=min.reads
  for(i in 1:ncol(starting.reads)) { valid.min = valid.min & starting.reads[,i]>=min.reads }

  # Predicate for max reads per sample
  valid.max = starting.reads[,1]<sample.quantiles[1]
  for(i in 1:ncol(starting.reads)) { valid.max = valid.max & starting.reads[,i]<sample.quantiles[i] }

  # Predicate for loci which are constant in all samples
  valid.variances = apply(methReads(methylDataRaw), APPLY_TO_ROWS, function(row) var(row)!=0)


  # Combine predicates
  valid.all = valid.min & valid.max & valid.variances

  # Apply the predicate to the read data
  info(meta.env$log.file, "Applying predicates...")
  meth.reads  = methReads(methylDataRaw)[valid.all,TRUE]
  total.reads = totalReads(methylDataRaw)[valid.all,TRUE]
  info(meta.env$log.file,paste("Retained", nrow(total.reads),"of", nrow(starting.reads),"loci"))  
  plotReadCounts(total.reads, paste0(meta.env$tissue,"_filtered_reads.min-reads_",min.reads, ".max-cov_", cov.filter))
  
  # Ensure methylated fractions of zero and 1 are never possible
  data.env$b_values = (meth.reads+1)/ (total.reads+2)
  colnames(data.env$b_values) = data.env$sampleInfo$Sample.Name
  data.env$methylDataRaw = methylDataRaw[valid.all,TRUE] # keep the valid ranges in the data environent
}

prep.data = function(){
  # Mean centre the b-values to prevent fully methylated or unmethylated values
  # swamping the results
  centred_m = as.data.frame(data.env$b_values) %>% 
    mutate(meanrow = rowMeans(.)) %>% 
    mutate_all(funs(.-meanrow)) %>%
    dplyr::select(-meanrow)
  data.env$centred_m = as.matrix(centred_m)

  # Calculate the cg mean b_values per group
  data.env$groups       = unique(data.env$sampleGroups$Group)
  data.env$mean.group.1 = rowMeans(data.env$b_values[,data.env$sampleGroups$Group==data.env$groups[1]])
  data.env$mean.group.2 = rowMeans(data.env$b_values[,data.env$sampleGroups$Group==data.env$groups[2]])
  data.env$means        = cbind(data.env$mean.group.1, data.env$mean.group.2, data.env$mean.group.1-data.env$mean.group.2)

  # Create the real and null model matrix for estimation
  data.env$mod.real = model.matrix(~Group, data=data.env$sampleGroups) 
  data.env$mod.null = model.matrix(~1, data=data.env$sampleGroups)
}

run.sva = function(){
  info(meta.env$log.file, "Estimating surrogate variables")

  tryCatch({
    data.env$svobj = sva::sva(data.env$centred_m,data.env$mod.real,data.env$mod.null,n.sv=NULL)
    }, error = function(e){
      warn("Error calculating surrogates. Check if any locus b-values are perfectly correlated with group.")
      quit(save="no", status=1)
    })

  if(!exists("svobj", envir = data.env)){
    warn(meta.env$log.file, "SVA failed using BE method, estimating number of SVs using Leek method")
    n.sv  = sva::num.sv(data.env$centred_m,data.env$mod.real,method="leek")
    info(meta.env$log.file, paste0("Estimated number of surrogate variables is ", n.sv))
    if(n.sv>0){
      data.env$svobj = sva::sva(data.env$centred_m,data.env$mod.real,data.env$mod.null,n.sv=n.sv)
    } else {
      warn(meta.env$log.file, "No surrogate variables detected")
    }
  } else {
    info(meta.env$log.file, paste0("Estimated number of surrogate variables is ", data.env$svobj$n.sv))
  }
}

estimate.factors = function(){
  # Add the estimated factors to the model, and fit the new model
  info(meta.env$log.file, paste("Fitting",data.env$svobj$n.sv,"surrogate variables"))
  data.env$mod.real.sv = cbind(data.env$mod.real,data.env$svobj$sv)
  data.env$mod.null.sv = cbind(data.env$mod.null,data.env$svobj$sv)
  data.env$fit = limma::lmFit(data.env$centred_m, data.env$mod.real.sv, method="robust")

  # Use ebayes to calculate the test statistics 
  info(meta.env$log.file, "Calculating test statistics")
  fit.e1  = eBayes(data.env$fit)

  # Get FDR corrected differential probe list, with any probes below the pvalue threshold
   data.env$tab.all     = topTable(fit.e1, p.value=0.05, number=Inf, adjust = "fdr")
   rm(centred_m, envir=data.env) # clean up values no longer needed
}

#' Export a table of significant differentially methylated loci
#' 
#' Adds genomic ranges to the table of loci and exports to the given file.
#' Ignores any loci for which the difference in methylation between groups
#' is less than the given minimum.
#'
export.sig.loci = function(){
  info(meta.env$log.file, "Exporting significant loci")
  # Get the names of the significant rows and fetch the corresponding genomic ranges
  loci = data.env$tab.all
  min_beta_diff = default.min.beta.diff

  sig.rows   = as.numeric(rownames(loci))
  sig.ranges = data.env$methylDataRaw[sig.rows]
  
  # Build the ranges into a table
  data.env$result.table = data.frame( chr    = seqnames(sig.ranges),
                             start  = start(sig.ranges),
                             end    = end(sig.ranges),
                             strand = strand(sig.ranges))
  
  read.counts = totalReads(sig.ranges)
  total.read.column.names = paste0(data.env$sampleGroups$Sample.Name, "_total_reads")
  colnames(read.counts) = total.read.column.names
  
  # Get the proxy beta values for the significant samples
  sig.b_vals = data.env$b_values[sig.rows,]
  colnames(sig.b_vals) = data.env$sampleGroups$Sample.Name
  
  sig.means = data.env$means[sig.rows,]
  colnames(sig.means) = c( paste0("Mean_",data.env$groups[1]), paste0("Mean_",data.env$groups[2]), "Diff")
  
  # Mung it all together in one big table
  data.env$result.table = cbind(data.env$result.table, loci, sig.means, sig.b_vals, read.counts)
  
  # Filter by beta difference
  data.env$result.table = data.env$result.table %>% dplyr::filter( abs(Diff)>=min_beta_diff)

  filename = paste0( meta.env$output.dir, meta.env$tissue, "_SVA_loci" )
  table.file = paste0(filename, ".", default.min.reads, "_reads.", default.high.cov.filter, "_cov.", min_beta_diff, "_diff.csv")
  debug(meta.env$log.file, paste("Exporting table to", table.file))
  write.table(data.env$result.table, file=table.file, sep=",", row.names = F, col.names = T)
}


find.loci.in.genes = function(){
  info(meta.env$log.file, "Finding loci in genes")
  # Read gene annotation info
  gtf.granges = import.gff(meta.env$gtf.file, format="gtf",feature.type="gene")

  # TODO: check if either seqnames is prefixed with 'chr' - can cause mismatches

  # Create GRanges from SVA loci
  data.env$sva.granges = with(data.env$result.table, GRanges(seqnames = chr, 
                                       IRanges(start=start, 
                                               end=end),
                                       strand = strand,
                                       seqinfo=Seqinfo(seqlevels(gtf.granges), 
                                                       seqlengths(gtf.granges)),
                                       pval = P.Value))

  overlaps     = findOverlaps(data.env$sva.granges, gtf.granges)
  sva.overlaps = data.env$sva.granges[queryHits(overlaps)]
  gtf.overlaps = gtf.granges[subjectHits(overlaps)]
  
  sva.overlaps$annotation = with(gtf.overlaps, paste(gene_id, type, gene_biotype, gene_name, sep='; '))
  
  # Get each CpG annotated with the gene it is within
  cpgs.annotated = as.data.frame(sva.overlaps) %>% 
    group_by(seqnames, start, end, width, strand, pval) %>% 
    summarise(annot = paste(annotation, collapse="|")) 

  output.file  = paste0(meta.env$output.dir, meta.env$tissue, "_CpGs_in_genes.tsv")
  write.table(cpgs.annotated, file=output.file, sep="\t", row.names = F)
  
  # Get the unique gene list
  genes = as.data.frame(gtf.overlaps) %>% 
    select(seqnames, start, end, gene_id, gene_biotype, gene_name) %>%
    arrange(seqnames, start) %>% 
    distinct()
  output.file  = paste0(meta.env$output.dir, meta.env$tissue, "_unique_genes.tsv")
  write.table(genes, file=output.file, sep="\t", row.names = F)
}

compare.loci.to.biseq.dmrs = function(){
  info(meta.env$log.file, "Comparing loci to BiSeq DMR results")
  dmr.file = paste0(dirname(meta.env$sample.group.file),"/Report/figure/DMR_results/", meta.env$tissue, "_DMRs_overlapping_target_regions.tsv")
  
  if(!file.exists(dmr.file)){
    info(meta.env$log.file, paste("Expected DMR file was not found:", dmr.file))
    return()
  }

  dmr.data = read.csv(dmr.file, header = T, sep="\t", stringsAsFactors = F)
  dmr.granges = with(dmr.data, GRanges(seqnames = seqnames, 
                                     IRanges(start=start, 
                                             end=end),
                                     seqinfo=seqinfo(data.env$sva.granges)))

  dmrs.with.sva = subsetByOverlaps(data.env$sva.granges,dmr.granges)

  dmr.sva.table = as.data.frame(dmrs.with.sva) %>% arrange(seqnames, start)
  dmr.output.file  = paste0(meta.env$output.dir, meta.env$tissue, "_SVA overlapping_DMRs.tsv")
  write.table(dmr.sva.table, file=dmr.output.file, sep="\t", row.names = F)

  sva.with.dmrs = subsetByOverlaps(dmr.granges, data.env$sva.granges)

  sva.dmr.table = as.data.frame(sva.with.dmrs) %>% arrange(seqnames, start)
  sva.output.file  = paste0(meta.env$output.dir, meta.env$tissue, "_DMRs overlapping_SVA.tsv")
  write.table(sva.dmr.table, file=sva.output.file, sep="\t", row.names = F)
}

##################
#
# Run the analysis
#
##################

func.order = c(load.data, prep.data, run.sva, estimate.factors, export.sig.loci, find.loci.in.genes, compare.loci.to.biseq.dmrs)
func.names = c(paste0(meta.env$temp.image.path, meta.env$tissue, "_save_", 1:(length(func.order)+1), ".RData")) # The corresponding save points

# Go through the list of functions, testing each to see if it needs to be run based on whether a save file exists
run.functions = function(i)  run.or.skip(func.names[i+1], func.names[i], func.order[[i]], meta.env$log.file, data.env)

invisible(lapply(1:length(func.order), run.functions))

info(meta.env$log.file, "SVA analysis complete")