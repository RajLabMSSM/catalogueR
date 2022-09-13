eQTLcatalogue_example_queries_local <- function(Rlib_path = NULL,
                                                locus_names = FALSE) {
  if (is.null(Rlib_path)) {
    cat_dir <- system.file("extdata/eQTL_Catalogue_queries",
                           package = "catalogueR"
    )
  } else {
    cat_dir <- file.path(
      Rlib_path,
      "catalogueR/extdata/eQTL_Catalogue_queries"
    )
  }
  sumstats_paths <- list.files(cat_dir,
                               pattern = "*.tsv.gz",
                               recursive = TRUE, full.names = TRUE
  )
  if (locus_names) {
    names(sumstats_paths) <- unlist(
      lapply(
        strsplit(basename(sumstats_paths), "_"),
        function(x) {
          x[1]
        }
      )
    )
  } else {
    names(sumstats_paths) <- gsub(".tsv|.gz", "", basename(sumstats_paths))
  }
  return(sumstats_paths)
}
