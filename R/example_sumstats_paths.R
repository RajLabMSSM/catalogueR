#' Paths to example summary stats
#'
#' Returns the paths to summary stats stored within \emph{catalogueR}.
#' Each file is the output of a locus that has been fine-mapped
#' using \link[echolocatoR]{finemap_loci}.
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://doi.org/10.1016/S1474-4422(19)30320-5}{
#' Nalls et al., 2019 (The Lancet Neurology)}.
#'
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @param save_dir Folder to save example summary statistics in.
#' @param verbose Print messages. 
#' 
#' @source \url{https://doi.org/10.1016/S1474-4422(19)30320-5}
#' @family Nalls23andMe_2019
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom data.table fwrite
#' @examples
#' sumstats_paths <- catalogueR::example_sumstats_paths()
example_sumstats_paths <- function(save_dir = file.path(
                                       tempdir(),
                                       "GWAS/Nalls2019"
                                   ),
                                   verbose = FALSE) {
    ss <- list(
        "BST1" = echodata::BST1,
        "LRRK2" = echodata::LRRK2,
        "MEX3C" = echodata::MEX3C
    )
    messager("Writing example summary statistics.", v = verbose)
    sumstats_paths <- lapply(names(ss), function(locus,
                                                 .verbose = verbose) {
        out_path <- file.path(
            save_dir, locus,
            paste0(locus, "_Nalls23andMe_2019_subset.tsv.gz")
        )
        dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
        messager(locus, "==>", out_path, v = .verbose)
        data.table::fwrite(
            x = ss[[locus]],
            file = out_path
        )
        return(out_path)
    }) %>%
        `names<-`(names(ss)) %>%
        unlist()
    return(sumstats_paths)
}
