#### ------------eQTL Catalogue------------####

#' eQTL Catalogue tabix header
#'
#'
#' The eQTL Catalogue tabix header (column names)
#'  is not always included in the file.
#' @source
#' \code{
#' eQTLcatalogue_header <- eQTLcatalogue_get_header(force_new_header = TRUE)
#' usethis::use_data(eQTLcatalogue_header, overwrite = TRUE)
#' }
#' @family eQTL Catalogue
#' @usage data("eQTLcatalogue_header")
"eQTLcatalogue_header"


#' eQTL Catalogue dataset metadata
#'
#' List of all queryable tabix-indexed eQTL Catalogue datasets
#' and their associated systems/tissues/cell types.
#' @source
#' \code{
#' meta <- eQTLcatalogue_list_datasets(force_new = TRUE) 
#' usethis::use_data(meta, overwrite = TRUE)
#' }
#' @family eQTL Catalogue
"meta"
