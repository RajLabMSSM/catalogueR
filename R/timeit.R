# Time a process
timeit <- function(func, digits = 3) {
    start <- Sys.time()
    out <- func
    end <- Sys.time()
    # print(paste("Completed in",round(end-start,digits = digits),"seconds."))
    round(end - start, digits = digits)
    return(out)
}
