#' Example colocalization results
#'
#' Example colocalization results from running \link[catalogueR]{run_coloc}.
#' on GWAS summary stats from all loci in
#' \href{https://doi.org/10.1016/S1474-4422(19)30320-5}{Nalls23andMe_2019}.
#' These results were published in
#' \href{https://pubmed.ncbi.nlm.nih.gov/34617105/}{
#' Schilder & Raj (Human Molecular Genetics, 2021)}
#' 
#' @param full Download the non-filtered version of the colocalization results
#' (Default: \code{FALSE}).
#' @inheritParams get_data
#' 
#' @source
#' \code{
#' ##### Subset results ####
#' # Pre-processing
#' gwas.qtl_paths <- example_eQTL_Catalogue_query_paths()
#' coloc_QTLs <- run_coloc(gwas.qtl_paths = gwas.qtl_paths,
#'                         nThread = 4,
#'                         top_snp_only = TRUE,
#'                         save_path = "~/Desktop/coloc_results.tsv.gz")
#' # Import pre-processed results
#' URL <- file.path("https://github.com/RajLabMSSM/catalogueR/raw/master",
#'                  "data/coloc_QTLs.rda")
#' tmp <- file.path(tempdir(),basename(URL))
#' download.file(URL, tmp)
#' piggyback::pb_upload(file = tmp,
#'                      repo = "RajLabMSSM/catalogueR")
#'
#'
#' ##### Full results ####
#' URL <- file.path("https://github.com/RajLabMSSM/catalogueR/raw/master",
#'                  "data/coloc_QTLs_full.rda")
#' tmp <- file.path(tempdir(),basename(URL))
#' download.file(URL, tmp)
#' piggyback::pb_upload(file = tmp,
#'                      repo = "RajLabMSSM/catalogueR")
#' }
#' @family coloc
#' @export
#' @examples
#' coloc_QTLs <- catalogueR::get_coloc_QTLs()
get_coloc_QTLs <- function(storage_dir = tempdir(),
                           full = FALSE) {
    tmp <- get_data(
        fname = if (full) "coloc_QTLs_full.rda" else "coloc_QTLs.rda",
        storage_dir = storage_dir
    )
    obj <- load_rdata(tmp)
    return(obj)
}
