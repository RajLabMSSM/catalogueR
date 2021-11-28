reconstruct_file_names <- function(coloc_QTLs) {
    gwas.qtl_files <- unique(paste0(coloc_QTLs$Locus.GWAS, "__", coloc_QTLs$qtl_id, ".tsv.gz"))
    return(gwas.qtl_files)
}
