#' Search eQTL Catalogue metadata
#'
#' Searches through multiple relevant metadata columns to find
#' eQTL Catalogue datasets
#' that match at least one of your substrings in a list.
#' All searches are case-insensitive.
#' If \code{qtl_search=NULL}, will return all available datasets.
#'
#' @param return_dataset_ids Return dataset IDs (default: \code{TRUE}) 
#' instead of human-readable dataset labels.
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
                                          return_dataset_ids=TRUE,
                                          verbose = TRUE) {
    meta <- eQTLcatalogue_list_datasets()
    if (is.null(qtl_search)) {
        messager("eQTL_Catalogue::",
            "Gathering data for all QTL Catalogue datasets...",
            v = verbose
        )
        if (return_dataset_ids) {
          qtl_datasets <- unique(meta$dataset_id)
        } else {
          qtl_datasets <- unique(meta$unique_label)
        }
        
        
    } else {
        if (return_dataset_ids){
          qtl_datasets <- meta[grep(
            pattern = paste(qtl_search, collapse = "|"),
            x = meta$unique_label, 
            ignore.case = TRUE
          )]$dataset_id |> unique()
          
        } else {
          qtl_datasets <- grep(
            pattern = paste(qtl_search, collapse = "|"),
            x = meta$unique_label,
            value = TRUE,
            ignore.case = TRUE
          ) |> unique()
        }
        
    }
    return(qtl_datasets)
}



#### Deprecation function #####
eQTL_Catalogue.search_metadata <- function(...){
  .Deprecated("eQTLcatalogue_search_metadata")
  eQTLcatalogue_search_metadata(...)
}
