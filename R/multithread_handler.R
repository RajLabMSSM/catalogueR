# Make sure there's not conflicting levels of parallelization
multithread_handler <- function(multithread_qtl,
                                multithread_loci,
                                multithread_tabix) {
    if (sum(c(multithread_qtl, multithread_loci, multithread_tabix)) > 1) {
        warning(
            "++ Only one multithreading option can be used at once. \n",
            "Setting: `multithread_qtl=T`, `multithread_loci=F`, `multithread_tabix=F`"
        )
        multithread_qtl <- T
        multithread_loci <- F
        multithread_tabix <- F
    }
    multithread_opts <- list(
        qtl = multithread_qtl,
        loci = multithread_loci,
        tabix = multithread_tabix
    )
    return(multithread_opts)
}
