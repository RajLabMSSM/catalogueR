#' Fetch from eQTL Catalogue API
#' 
#' Sub-function of fetch_restAPI.
#' @keywords internal
#' @param link File link.
#' @param is_gwas Whether the file is a GWAS dataset.
#' @importFrom rlang is_empty
#' @importFrom tidyr unnest
fetch_from_eqtl_cat_API <- function(link,
                                    is_gwas = FALSE,
                                    timeout=5*60) {
  options(timeout = timeout)
  nullToNA <- function(x) {
      x[sapply(x, is.null)] <- NA
      return(x)
  }
  if (isTRUE(is_gwas)) {
      cols_to_nest <- c(
          "variant_id", "chromosome", "base_pair_location", "trait",
          "p_value", "ci_lower", "ci_upper", "beta",
          "effect_allele", "other_allele", "effect_allele_frequency",
          "odds_ratio", "study_accession", "code"
      )
  } else {
      cols_to_nest <- c(
          "study_id", "qtl_group", "rsid",
          "chromosome", "position", "pvalue", "condition_label",
          "tissue_label", "molecular_trait_id", "gene_id", "ac",
          "ref", "beta", "variant", "an", "median_tpm", "condition",
          "r2", "alt", "type", "maf", "tissue"
      )
  }
  is_paginated <- !stringr::str_detect(link, "paginate=False")
  # message("isPagined:", is_paginated)
  page <- 1
  merged_summaries <- data.frame()
  while (!is.null(link)) {
      # print(paste0("Fetching page #",page))
      api_raw_data <- jsonlite::fromJSON(link, 
                                         simplifyDataFrame = TRUE,
                                         flatten = TRUE)
      link <- api_raw_data$`_links`$`next`$href
      if (rlang::is_empty(api_raw_data$`_embedded`$associations)) {
          return(merged_summaries)
      }
      eqtl_raw_list_data <- do.call(
          rbind, lapply(api_raw_data$`_embedded`$associations, rbind))
      eqtl_data <- nullToNA(eqtl_raw_list_data) |>
          as.matrix() |>
          dplyr::as_tibble()
      if (isTRUE(is_paginated)) {
          eqtl_data <- dplyr::select(eqtl_data, -c("_links"))
      }
      eqtl_data <- tidyr::unnest(eqtl_data, 
                                 cols = dplyr::all_of(cols_to_nest))
      if (!is.null(link)) {
          page <- page + 1
      }
      merged_summaries <- merged_summaries |> rbind(eqtl_data)
  }
  return(data.table::data.table(merged_summaries[cols_to_nest]))
}
