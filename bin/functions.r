# Common functions for the methylation sequencing 
# pipline.


# Check if the given path exists.
# If it does not, print the given
# message and quit with exit code 1
#
# path - the file path to check
# msg  - the error message to print
checkPath = function(path, msg){
    if(!file.exists(path)){
        cat(msg, "'", path, "' not found\n" )
        quit(save="no", status=1)
    }
}

