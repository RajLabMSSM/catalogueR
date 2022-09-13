multithread_optimizer <- function(qtl_datasets,
                                  sumstats_paths,
                                  verbose = TRUE) {
    messager("+ Optimizing multi-threading...", v = verbose)
    if (length(qtl_datasets) > length(sumstats_paths)) {
        messager("++ Multi-threading across QTL datasets.", v = verbose)
        multithread_opts <- list(
            qtl = TRUE,
            loci = F,
            tabix = FALSE
        )
    } else {
        messager("++ Multi-threading across loci.", v = verbose)
        multithread_opts <- list(
            qtl = F,
            loci = TRUE,
            tabix = FALSE
        )
    }
    return(multithread_opts)
}
