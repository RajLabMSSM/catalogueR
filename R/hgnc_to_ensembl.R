#' Convert HGNC gene symbols to ENSEMBL IDs
#'
#' @family utils
#' @examples
#' gene_symbols <- c("BDNF", "FOXP2", "BST1")
#' ensembl_ids <- hgnc_to_ensembl(gene_symbols)
hgnc_to_ensembl <- function(gene_symbols,
                            unique_only = TRUE,
                            verbose = TRUE) {
    messager("++ Converting: HGNC gene symbols ==> Ensembl IDs", v = verbose)
    if (unique_only) gene_symbols <- unique(gene_symbols)
    # columns(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
    gene_symbols[is.na(gene_symbols)] <- "NA"
    conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
        keys = gene_symbols,
        keytype = "SYMBOL",
        column = "GENEID"
    )
    return(conversion)
}
