#' Search eQTL Catalogue metadata
#'
#' Searches through multiple relevant metadata columns to find
#' eQTL Catalogue datasets
#' that match at least one of your substrings in a list.
#' All searches are case-insensitive.
#' If \code{qtl_search=NULL}, will return all available datasets.
#'
#' @inheritParams eQTLcatalogue_query
#' @family eQTL Catalogue
#' @export
#' @examples
#' qtl_datasets <- eQTLcatalogue_search_metadata(qtl_search = c(
#'     "Alasoo_2018",
#'     "monocyte"
#' ))
#' qtl_datasets.brain <- eQTLcatalogue_search_metadata(qtl_search = "brain")
eQTLcatalogue_search_metadata <- function(qtl_search = NULL,
                                           verbose = TRUE) {
    meta <- eQTLcatalogue_list_datasets()
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
        ) |> unique()
    }
    return(qtl_datasets)
}



#### Deprecation function #####
eQTL_Catalogue.search_metadata <- function(...){
  .Deprecated("eQTLcatalogue_search_metadata")
  eQTLcatalogue_search_metadata(...)
}
