# Time a process
timeit <- function(func, digits = 3) {
    start <- Sys.time()
    out <- func
    end <- Sys.time()
    round(end - start, digits = digits)
    return(out)
}
