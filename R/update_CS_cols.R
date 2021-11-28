update_CS_cols <- function(finemap_dat) {
    colnames(finemap_dat) <- gsub("*.Credible_Set$", ".CS", colnames(finemap_dat))
    return(finemap_dat)
}
