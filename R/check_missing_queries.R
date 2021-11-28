check_missing_queries <- function(qtl_datasets,
                                  loci,
                                  GWAS.QTL,
                                  verbose = TRUE) {
    messager("+ Checking results for missing queries...")
    GWASlocus__QTLid.queried <- unlist(lapply(qtl_datasets, function(x) {
        paste(loci, x, sep = "__")
    }))
    GWASlocus__QTLid.found <- unique(paste(GWAS.QTL$Locus.GWAS, GWAS.QTL$qtl_id, sep = "__"))

    missing_queries <- GWASlocus__QTLid.queried[!GWASlocus__QTLid.queried %in% GWASlocus__QTLid.found]
    if (length(missing_queries) > 0) {
        messager("+ QTL datasets with no hits/failed to pull data:", v = verbose)
        for (x in missing_queries) {
            messager("  +", x)
        }
    }
    return(missing_queries)
}
