#' Retrieve summary statistics info 
#' 
#' Merge coloc results with SNP-wise GWAS/QTL data
#' @family coloc
#' @keywords internal
#' @importFrom data.table fread
retrieve_sumstats_info <- function(output_dir = file.path(tempdir(),
                                                          "catalogueR_queries"),
                                   coloc_QTLs) {
    qtl_id <- 
    query_paths <- file.path(output_dir, qtl_id)
    lapply(seq_len(nrow(coloc_QTLs)), function(i,
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
