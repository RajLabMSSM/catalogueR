#' Convert ENSEMBL IDs to HGNC gene symbols
#'
#' Convert ENSEMBL IDs to HGNC gene symbols using the 
#' \link[EnsDb.Hsapiens.v75]{EnsDb.Hsapiens.v75} Bioconductor package.
#' 
#' @param ensembl_ids Character vector of Ensembl gene IDs.
#' @param unique_only Only query unique entries \code{ensembl_ids}.
#' @param verbose Print messages.
#' 
#' @family utils
#' @export
#' @importFrom AnnotationDbi mapIds
#' @examples
#' ensembl_ids <- c("ENSG00000176697", "ENSG00000128573", "ENSG00000109743")
#' gene_symbols <- catalogueR::ensembl_to_hgnc(ensembl_ids = ensembl_ids)
ensembl_to_hgnc <- function(ensembl_ids,
                            unique_only = TRUE,
                            verbose = TRUE) {
    messager("++ Converting: Ensembl IDs ==> HGNC gene symbols", v = verbose)
    if (unique_only) ensembl_ids <- unique(ensembl_ids)
    ensembl_ids[is.na(ensembl_ids)] <- "NA"
    conversion <- suppressWarnings(
        AnnotationDbi::mapIds(
            EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
            keys = ensembl_ids,
            keytype = "GENEID",
            column = "SYMBOL") 
    )
    return(conversion)
}
