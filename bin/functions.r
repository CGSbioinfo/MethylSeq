# Common functions for the methylation sequencing 
# pipline.

checkPath = function(path, msg){
    # Check if the given path exists.
    # If it does not, print the given
    # message and quit with exit code 1
    #
    # path - the file path to check
    # msg  - the error message to print
    if(!file.exists(path)){
        cat(msg, "'", path, "' not found\n" )
        quit(save="no", status=1)
    }
}

PrintTimeTaken = function(proc.time.object){
    # Format the given relative proc.time() object.
    #
    # proc.time.object - a result of proc.time()
    #
    # Example usage: 
    # ptm = proc.time()
    # DoAFunction()
    # PrintTimeTaken(proc.time() - ptm)
    days  = floor(proc.time.object / 86400)
    proc.time.object   = proc.time.object - (days*86400)
    hours = floor(proc.time.object / 3600)
    proc.time.object   = proc.time.object - (hours*3600)
    mins  = floor(proc.time.object / 60)
    secs  = proc.time.object - (mins*60)
    cat("Completed in", days[3], "days", hours[3], "hours", mins[3],"minutes and", secs[3], "seconds\n" )
}

RunAndTime = function(f){
    # Run the given function and print
    # the time it took to complete
    #
    # f - the function to run
    if(!is.function(f)){
        cat("A function was not supplied")
        return()
    }
    ptm = proc.time()
    f()
    PrintTimeTaken(proc.time() - ptm)
}

PrintEnvironmentSize = function(){
    # Print the size of the objects in the global
    # environmnent
    cat("Global environment:\n")
    for ( o in ls(envir=globalenv()) ) {
      cat("\t", o, ": ", object.size(get(o)), "\n")
    }
}
