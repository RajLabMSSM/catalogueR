# Make up a locus name based on coordinates
construct_locus_name <- function(query_dat,
                                 verbose = TRUE) {
    messager("++ Constructing locus name from coordinates", v = verbose)
    paste0("locus_chr", query_dat$CHR[1], "-",
           min(query_dat$POS), "-", max(query_dat$POS))
}
