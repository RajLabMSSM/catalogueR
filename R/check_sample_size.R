check_sample_size <- function(qtl.dat,
                              sample_size) {
    if (!"N" %in% colnames(qtl.dat)) {
        if (is.null(sample_size)) {
            stop_msg <- paste(
                "`N` column (sample size) not detected.",
                "Please provide `sample_size=` argument instead."
            )
            stop(stop_msg)
        } else {
            qtl.dat$N <- sample_size
        }
    }
    return(qtl.dat)
}
