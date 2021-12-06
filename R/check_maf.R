#' Check MAF
#' 
#' Check Minor Allele Frequency (MAF) column is 
#' prepared for running \pkg{coloc}. Must be numeric, <1, >0, 
#' and not contain NAs.
#' @keywords internal
check_maf <- function(gwas_shared,
                      eqtl_shared) {
    MAF <- maf.QTL <- NULL;
    if (!"MAF" %in% colnames(gwas_shared)) {
        messager(
            "WARNING: `MAF` column not provided in GWAS data.",
            "Borrowing MAF from QTL data instead."
        )
        gwas_shared$MAF <- eqtl_shared$maf.QTL
    }
    #### Ensure MAF is within expected range ####
    ## gwas_shared
    gwas_shared$MAF <- as.numeric(gwas_shared$MAF)
    gwas_shared <- subset(
        gwas_shared,
        MAF > 0 & MAF < 1 & (!is.na(MAF))
    )
    ## eqtl_shared
    eqtl_shared$maf.QTL <- as.numeric(eqtl_shared$maf.QTL)
    eqtl_shared <- subset(
        eqtl_shared,
        maf.QTL > 0 & maf.QTL < 1 & (!is.na(maf.QTL))
    )
    return(list(
        gwas_shared = gwas_shared,
        eqtl_shared = eqtl_shared
    ))
}
