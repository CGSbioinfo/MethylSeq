# Common functions for the methylation sequencing 
# pipline.

#' Load the requested packages, installing them if needed.
#' @title Install missing packages
#' @param packages the R packages to install
#' @param biopackages the Bioconductor packages to install
#' @param repos the CRAN repositories to try. Passes through to install.packages
#' @examples
#' install.missing(packages=c("ggplot2", "dplyr"), biopackages=c("BiSeq", "rtracklayer"))
#' @export
install.missing = function(packages, biopackages=c(), repos = "https://cran.ma.imperial.ac.uk/" ) {
  
  source("https://bioconductor.org/biocLite.R")  
  installIfNeededFromBioconductor = function(pkg){
    if(pkg %in% rownames(installed.packages()) == FALSE) {
      biocLite(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
  }
  
  installIfNeeded = function(pkg){
    if(pkg %in% rownames(installed.packages()) == FALSE) {
      install.packages(pkg, dependencies = TRUE, repos = repos)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, warn.conflicts = FALSE))
  }

  tryCatch({
      lapply(biopackages, installIfNeededFromBioconductor)
      lapply(packages, installIfNeeded)
    }, error=function(e){
      cat("Error loading required libraries\n")
      message(e)
      quit()
    })
}

#' Log the details of the current session info to the given
#' file.
#' @title Log session info
#' @param log.file the file path to log to
#' @export
log.session.info = function(log.file){
  debug( log.file, paste0("Running on ", system("hostname", intern = TRUE)))
  debug( log.file, paste0("Platform: ", sessionInfo()$running ))
  debug( log.file, paste0("R version: ", getRversion()))
  log.pkgs = function(pkg) debug( log.file, paste0(pkg$Package, " - ", pkg$Version))
  invisible(lapply(sessionInfo()$otherPkgs, log.pkgs))
}

#' This function checks if a file path exists.
#' If the path does not exist, prints the given
#' message and quits with exit code 1
#' @title Check if the given path exists and quit if missing
#' @param path the file path to check
#' @param msg the error message to print if path is missing
#' @examples
#' dir = "/path/to/folder/"
#' quit.if.not.exists(dir)
#' setwd(dir)
#' @export
quit.if.not.exists = function(path, msg){
    if(!file.exists(path)){
        cat(msg, "'", path, "' not found\n" )
        quit(save="no", status=1)
    }
}

#' This function creates a directory if not present.
#' @title Check if the given path exists and create if needed
#' @param dir the directory path to check
#' @examples
#' dir = "/path/to/folder/"
#' ensure.dir.exists(dir)
#' @export
ensure.dir.exists = function(dir){
    if(!dir.exists(dir)){
        dir.create(dir)
    }
}

#' Formats the given relative \code{\link{proc.time()}} result, 
#' prints the time to screen, and returns the time in tab format.
#' @title Format a proc.time result
#' @param proc.time.object - a result of \code{\link{proc.time()}}
#' @return a tab separated string of days to seconds
#' @examples
#'
#' ptm = proc.time()
#' DoAFunction()
#' printTimeTaken(proc.time() - ptm)
#' @export
printTimeTaken = function(proc.time.object){
    days  = floor(proc.time.object / 86400)
    proc.time.object   = proc.time.object - (days*86400)
    hours = floor(proc.time.object / 3600)
    proc.time.object   = proc.time.object - (hours*3600)
    mins  = floor(proc.time.object / 60)
    secs  = proc.time.object - (mins*60)
    cat("Completed in", days[3], "days", hours[3], "hours", mins[3],"minutes and", secs[3], "seconds\n" )
    paste(days[3], hours[3], mins[3], secs[3], sep="\t")
}

#' Formats the given relative \code{\link{proc.time()}} result, 
#' prints the time to screen, and returns the time in tab format.
#' @title Format a proc.time result
#' @param proc.time.object - a result of \code{\link{proc.time()}}
#' @return a tab separated string of days to seconds
#' @examples
#'
#' ptm = proc.time()
#' DoAFunction()
#' print.time.taken(proc.time() - ptm, log.file)
#' @export
print.time.taken = function(proc.time.object, log.file){
    days  = floor(proc.time.object / 86400)
    proc.time.object   = proc.time.object - (days*86400)
    hours = floor(proc.time.object / 3600)
    proc.time.object   = proc.time.object - (hours*3600)
    mins  = floor(proc.time.object / 60)
    secs  = proc.time.object - (mins*60)
    debug(log.file, paste("Completed in", days[3], "days", hours[3], "hours", mins[3],"minutes and", secs[3], "seconds", sep=" ") )
    paste(days[3], hours[3], mins[3], secs[3], sep="\t")
}

#' Check if the next analysis step output exists. If so, skip 
#' this analysis step. If not, check if this step's output exists
#' and if not, run the analysis.
#' @title Run or skip analysis step
#' @param next.temp.file the name of the temp file after this step
#' @param temp.file the name of the temp file for this step
#' @param function.to.run the function to be run in this step
#' @param log.file the file to log output to
#' @param data.env the environment to save to file
#' @export
run.or.skip = function(next.temp.file, temp.file, function.to.run, log.file, data.env){
  if(!file.exists(next.temp.file)){
    # The next step has not been run. Has the current step been run?
    if(file.exists(temp.file)){
      # The current step has been run. Load the data.
      info(log.file, paste0("Loading temporary analysis file ", temp.file, "..."))
      run.and.time( f=function(){ load(temp.file, envir = globalenv()) }, log.file = log.file )
      info( log.file, paste0("Loaded temporary analysis file ", temp.file))
    } else{
      # The current step has not been run. Run it, and save the data.
      run.and.time(function.to.run, log.file)
      info( log.file, paste0("Saving temporary analysis file ", temp.file, "..."))
      # Only save the given data environment - ensures updated script functions will not be overwritten 
      run.and.time(f=function(){ save(data.env, file = temp.file) }, log.file = log.file)
    }
  }
}

#' Run the given function and print
#' the time it took to complete.
#' @title Run and time the given function
#'
#' @param f the function to run
#' @return the return value of f.
#' @examples
#' f = function(){
#'  sleep(10)
#' }
#' runAndTime(f)
#' @export
runAndTime = function(f){
    if(!is.function(f)){
        cat("A function was not supplied")
        return()
    }
    ptm = proc.time()
    r = f()
    printTimeTaken(proc.time() - ptm)
    r
}

#' Run the given function and print
#' the time it took to complete.
#' @title Run and time the given function
#'
#' @param f the function to run
#' @return the return value of f.
#' @examples
#' f = function(){
#'  sleep(10)
#' }
#' runAndTime(f)
#' @export
run.and.time = function(f, log.file){
    if(!is.function(f)){
        cat("A function was not supplied")
        return()
    }
    ptm = proc.time()
    r = f()
    print.time.taken(proc.time() - ptm, log.file)
    r
}


#' Print the size of the objects in the given environmnent
#' @param env the environment
#' @examples
#'  
#' printEnvironmentSize(globalenv())
#' @export
printEnvironmentSize = function(env){
    cat("Environment:\n")
    for ( o in ls(envir=env) ) {
      cat("\t", o, ": ", object.size(get(o, envir=env)), "\n")
    }
}

#' Log a message to console and file with the current time
#' @param file the target file
#' @param msg the message
#' @export
logToFile = function(file, msg){
    cat(msg, "\n")
    msg = gsub("^ ", "", paste(Sys.time(), msg, "\n", sep="\t"))
    cat(msg, file=file, append=TRUE)
}

#' Log a message (or list of messages) to file with log level INFO
#' @param file the target file
#' @param msg the message(s)
#' @export
info = function(file, msg){
    write.line = function(s) logToFile(file, paste("INFO", s, sep="\t"))
    invisible(lapply(msg, write.line))
}

#' Log a message (or list of messages) to file with log level DEBUG
#' @param file the target file
#' @param msg the message(s)
#' @export
debug = function(file, msg){
    write.line = function(s) logToFile(file, paste("DEBUG", s, sep="\t"))
    invisible(lapply(msg, write.line))
}

#' Log a message to file with log level WARN
#' @param file the target file
#' @param msg the message
#' @export
warn = function(file, msg){
    logToFile(file, paste("WARN", msg, sep="\t"))
}

#' Ensure a string ends with '/'
#' @param s the string
#' @return the string ending with /
#' @export
slash.terminate = function(s){
    ifelse(endsWith(s, "/"), s, paste0(s, "/"))
}
