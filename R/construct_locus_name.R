# Make up a locus name based on coordinates
construct_locus_name <- function(gwas_data,
                                 verbose = TRUE) {
    messager("++ Constructing locus name from coordinates", v = verbose)
    paste0("locus_chr", gwas_data$CHR[1], "-", min(gwas_data$POS), "-", max(gwas_data$POS))
}
