#' Paths to example summary stats
#'
#' Returns the paths to summary stats stored within \emph{catalogueR}.
#' Each file is the output of a locus that has been fine-mapping
#' using \emph{echolocatoR}.
#' Data originally comes from the Parkinson's disease GWAS
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
#'
#' @export
#' @examples
#' sumstats_paths <- example_sumstats_paths()
example_sumstats_paths <- function(Rlib_path = NULL) {
    if (is.null(Rlib_path)) {
        cat_dir <- system.file("extdata/Nalls23andMe_2019", package = "catalogueR")
    } else {
        cat_dir <- file.path(Rlib_path, "catalogueR/extdata/Nalls23andMe_2019")
    }
    sumstats_paths <- list.files(cat_dir,
        pattern = "*_subset.tsv.gz",
        recursive = TRUE,
        full.names = TRUE
    )
    locus_names <- unlist(lapply(
        strsplit(basename(sumstats_paths), "_"),
        function(x) {
            x[1]
        }
    ))
    names(sumstats_paths) <- locus_names
    return(sumstats_paths)
}
