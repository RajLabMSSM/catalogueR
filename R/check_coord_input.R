#' Check coordinates input
#' 
#' Make sure there's some way to specify genomic coordinates for tabix query.
#' @keywords internal
check_coord_input <- function(gwas_data,
                              chrom,
                              bp_lower,
                              bp_upper) {
    if (is.null(gwas_data) & 
        all(c(is.null(chrom), is.null(bp_lower), is.null(bp_upper)))) {
        stop_msg <- paste(
          "+ User must specify coordinates to fetch",
          "using either `gwas_data` OR `chrom`,`bp_lower`, and `bp_upper`")
        stop(stop_msg)
    }
}
