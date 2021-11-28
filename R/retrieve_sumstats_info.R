#' Merge coloc results with SNP-wise GWAS/QTL data
#'
#' @family coloc
#' @examples
#' \dontrun{
#' data("coloc_QTLs")
#' query_paths <- example_eQTL_Catalogue_query_paths()
#' output_dir <- unique(dirname(dirname(query_paths)))
#' }
retrieve_sumstats_info <- function(output_dir = "./catalogueR_queries",
                                   coloc_QTLs) {
    query_paths <- file.path(output_dir, qtl_id, )
    lapply(1:nrow(coloc_QTLs), function(snp,
                                        .output_dir = output_dir) {
        ROW <- coloc_QTLs[i, ]
        dat_path <- make_split_path(
            output_dir = .output_dir,
            qtl_id = ROW$qtl_id,
            loc = ROW$Locus.GWAS
        )
        dat <- data.table::fread(dat_path, nThread = 1)
    })
    file.path()
}
