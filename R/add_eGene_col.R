# Add eGene column
add_eGene_col <- function(qtl.dat) {
    # Remove any old col
    if ("eGene" %in% colnames(qtl.dat)) {
        qtl.dat <- subset(qtl.dat, select = -eGene)
    }
    qtl.dat <- qtl.dat %>% dplyr::mutate(eGene = ifelse(!is.na(gene.QTL) & gene.QTL != "", gene.QTL, molecular_trait_id.QTL))
    return(qtl.dat)
}
