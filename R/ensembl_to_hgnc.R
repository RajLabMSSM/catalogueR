#' Convert ENSEMBL IDs to HGNC gene symbols
#'
#' @family utils
#' @examples
#' ensembl_ids <- c("ENSG00000176697", "ENSG00000128573", "ENSG00000109743")
#' gene_symbols <- ensembl_to_hgnc(ensembl_ids = ensembl_ids)
ensembl_to_hgnc <- function(ensembl_ids,
                            unique_only = TRUE,
                            verbose = TRUE) {
    messager("++ Converting: Ensembl IDs ==> HGNC gene symbols", v = verbose)
    if (unique_only) ensembl_ids <- unique(ensembl_ids)
    ensembl_ids[is.na(ensembl_ids)] <- "NA"
    conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
        keys = ensembl_ids,
        keytype = "GENEID",
        column = "SYMBOL"
    )
    return(conversion)
}
