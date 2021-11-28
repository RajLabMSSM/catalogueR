# rbind a list of files
rbind.file.list <- function(file.list,
                            verbose = TRUE,
                            nCores = 4) {
    merged.dat <- parallel::mclapply(file.list, function(x) {
        messager(x, v = verbose)
        dat <- data.table::fread(x, nThread = 1)
        return(dat)
    }, mc.cores = nCores) %>% data.table::rbindlist(fill = TRUE)
    return(merged.dat)
}
