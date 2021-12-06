#' Search eQTL Catalogue metadata
#'
#' Searches through multiple relevant metadata columns to find
#' eQTL Catalogue datasets
#' that match at least one of your substrings in a list.
#' All searches are case-insensitive.
#' If \code{qtl_search=NULL}, will return all available datasets.
#'
#' @inheritParams eQTL_Catalogue.query
#' @family eQTL Catalogue
#' @export
#' @examples
#' qtl_datasets <- eQTL_Catalogue.search_metadata(qtl_search = c(
#'     "Alasoo_2018",
#'     "monocyte"
#' ))
#' qtl_datasets.brain <- eQTL_Catalogue.search_metadata(qtl_search = "brain")
eQTL_Catalogue.search_metadata <- function(qtl_search = NULL,
                                           verbose = TRUE) {
    meta <- eQTL_Catalogue.list_datasets()
    if (is.null(qtl_search)) {
        messager("eQTL_Catalogue::",
            "Gathering data for all QTL Catalogue datasets...",
            v = verbose
        )
        qtl_datasets <- unique(meta$unique_id)
    } else {
        qtl_datasets <- grep(
            pattern = paste(qtl_search, collapse = "|"),
            x = meta$unique_id,
            value = TRUE,
            ignore.case = TRUE
        ) %>% unique()
    }
    return(qtl_datasets)
}
