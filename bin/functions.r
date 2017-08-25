# Common functions for the methylation sequencing 
# pipline.

#' This function checks if a file path exists.
#' If the path does not exist, prints the given
#' message and quits with exit code 1
#' @title Check if the given path exists and quit if missing
#' @param path the file path to check
#' @param msg the error message to print if path is missing
#' @examples
#' dir = "/path/to/folder/"
#' pathExistsOrQuit(dir)
#' setwd(dir)
#' @export
pathExistsOrQuit = function(path, msg){
    if(!file.exists(path)){
        cat(msg, "'", path, "' not found\n" )
        quit(save="no", status=1)
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
