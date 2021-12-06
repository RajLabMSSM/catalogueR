progress_bar_check <- function(progress_bar,
                               verbose = TRUE) {
    if (progress_bar) {
        if (requireNamespace("pbmcapply")) {
            lapply_func <- pbmcapply::pbmclapply
        } else {
            messager(
                "++ R package `pbmcapply` not installed.",
                "Turning off `progress_bar`.", v = verbose)
            lapply_func <- parallel::mclapply
        }
    } else {
        lapply_func <- parallel::mclapply
    }
    return(lapply_func)
}
