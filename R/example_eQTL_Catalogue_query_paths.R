#' Paths to example eQTL Catalogue query results
#'
#' Returns the paths to eQTL Catalogue query results stored within \emph{catalogueR}.
#' Each file is a merged data.table of the GWAS summary stats used to make the query,
#' and the eQTL Catalogue query results (which can contain data for multiple eGenes).
#'
#' GWAS data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#'
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @param Rlib_path This function will automatically find your Rlib path,
#' but you can override this by supplying it manually.
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @family Nalls23andMe_2019
#' These example files can be used
#' @examples
#' gwas.qtl_paths <- example_eQTL_Catalogue_query_paths()
example_eQTL_Catalogue_query_paths <- function(Rlib_path = NULL) {
    if (is.null(Rlib_path)) {
        cat_dir <- system.file("extdata/eQTL_Catalogue_queries", package = "catalogueR")
    } else {
        cat_dir <- file.path(Rlib_path, "catalogueR/extdata/eQTL_Catalogue_queries")
    }
    sumstats_paths <- list.files(cat_dir, pattern = "*.tsv.gz", recursive = TRUE, full.names = TRUE)
    locus_names <- unlist(lapply(strsplit(basename(sumstats_paths), "_"), function(x) {
        x[1]
    }))
    names(sumstats_paths) <- locus_names
    return(sumstats_paths)
}
