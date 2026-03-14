# Make up a locus name based on coordinates
construct_locus_name <- function(query_dat,
                                 verbose = TRUE) {
    messager("++ Constructing locus name from coordinates", v = verbose)
    pos_col <- intersect(c("POS","BP"), colnames(query_dat))[1]
    if(is.na(pos_col)) stop("No position column (POS/BP) found in query_dat.")
    paste0("locus_chr", query_dat$CHR[1], "-",
           min(query_dat[[pos_col]]), "-", max(query_dat[[pos_col]]))
}
