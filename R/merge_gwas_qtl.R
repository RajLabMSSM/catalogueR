#' Merge GWAS with QTL
#' 
#' Merge GWAS data (query) and QTL data (results).
#' @family eQTL Catalogue
#' @keywords internal
#' @inheritParams eQTL_Catalogue.query
merge_gwas_qtl <- function(gwas_data,
                           qtl.subset,
                           verbose = TRUE) {
    A1 <- ref.QTL <- alt.QTL <- NULL;
    messager("++ Merging GWAS data and QTL query results.", v = verbose)
    gwas.qtl <- tryCatch(expr = {
        # Merging and allele flipping
        gwas.qtl <- data.table::merge.data.table(
            x = data.table::data.table(gwas_data),
            y = data.table::data.table(qtl.subset),
            # all.x  = TRUE,
            by.x = c("SNP"), # effect_allele
            by.y = c("rsid.QTL")
        ) %>%
            # subset(effect.is.ref|effect.is.alt) %>%
            data.table::data.table()
        if ("A1" %in% colnames(gwas.qtl) & "A2" %in% colnames(gwas.qtl)) {
            gwas.qtl <- gwas.qtl %>% dplyr::mutate(
                effect.is.ref = A1 == ref.QTL,
                effect.is.alt = A1 == alt.QTL
            )
        }
        return(gwas.qtl)
    }, error = function(e) {
        return(qtl.subset)
    })
    return(gwas.qtl)
}
