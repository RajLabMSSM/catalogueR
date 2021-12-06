#' Convert HGNC gene symbols to ENSEMBL IDs
#'
#' Convert HGNC gene symbols to ENSEMBL IDs using the
#' \link[EnsDb.Hsapiens.v75]{EnsDb.Hsapiens.v75} Bioconductor package.
#' 
#' @param gene_symbols Character vector of HGNC gene IDs.
#' @param unique_only Only query unique entries \code{gene_symbols}.
#' @param verbose Print messages.
#' 
#' @family utils
#' @export
#' @importFrom AnnotationDbi mapIds
#' @examples
#' gene_symbols <- c("BDNF", "FOXP2", "BST1")
#' ensembl_ids <- catalogueR::hgnc_to_ensembl(gene_symbols)
hgnc_to_ensembl <- function(gene_symbols,
                            unique_only = TRUE,
                            verbose = TRUE) {
    messager("++ Converting: HGNC gene symbols ==> Ensembl IDs", v = verbose)
    if (unique_only) gene_symbols <- unique(gene_symbols)
    # columns(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
    gene_symbols[is.na(gene_symbols)] <- "NA"
    conversion <- suppressWarnings(
        AnnotationDbi::mapIds(
            EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
            keys = gene_symbols,
            keytype = "SYMBOL",
            column = "GENEID")  
    )
    return(conversion)
}
