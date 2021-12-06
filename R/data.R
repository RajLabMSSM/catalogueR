#### ------------eQTL Catalogue------------####

#' eQTL Catalogue tabix header
#'
#'
#' The eQTL Catalogue tabix header (column names)
#'  is not always included in the file.
#' @source
#' \code{
#' eQTL_Catalogue.header <- tabix_header(force_new_header = TRUE)
#' usethis::use_data(eQTL_Catalogue.header, overwrite = TRUE)
#' }
#' @family eQTL Catalogue
#' @usage data("eQTL_Catalogue.header")
"eQTL_Catalogue.header"


#' eQTL Catalogue dataset metadata
#'
#' List of all queryable tabix-indexed eQTL Catalogue datasets
#' and their associated systems/tissues/cell types.
#' @source
#' \code{
#' meta <- eQTL_Catalogue.list_datasets(force_new = TRUE)
#' # Some paths in the metadata were originally wrong.
#' # Has since been corrected by authors.
#' # meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte",
#' #                                               "Fairfax_2014",ftp_path))
#' usethis::use_data(meta, overwrite = TRUE)
#' }
#' @family eQTL Catalogue
"meta"
