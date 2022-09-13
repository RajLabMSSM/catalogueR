#' Iterate queries to \emph{eQTL Catalogue}
#'
#' Determines which datasets to query using \code{qtl_search}.
#' Uses coordinates from stored summary stats files (e.g. GWAS)
#' to determine which regions to query from \emph{eQTL Catalogue}.
#' Each locus file can be stored separately,
#' or merged together to form one large file with all query results.
#'
#' @param sumstats_paths A list of paths to any number of summary stats files
#' whose coordinates you want to use to make queries to eQTL Catalogue.
#' If you wish to add custom names to the loci, simply add these as the
#' names of the path list
#'  (e.g. \code{c(BST1="<path>/<to>/<BST1_file>",
#'   LRRK2="<path>/<to>/<LRRK2_file>")}).
#'  Otherwise, loci will automatically named based on their min/max
#'  genomic coordinates.
#'
#' The minimum columns in these files required to make queries include:
#' \describe{
#' \item{SNP}{RSID of each SNP.}
#' \item{CHR}{Chromosome (can be in "chr12" or "12" format).}
#' \item{POS}{Genomic position of each SNP.}
#' \item{...}{Optional extra columns.}
#' }
#' @param output_dir The folder you want the merged gwas/qtl results to be
#' saved to (set to \code{NULL} to not save the results).
#' If \code{split_files=FALSE}, all query results will be merged into one and
#' saved as \emph{<output_dir>/eQTLcatalogue_tsv.gz}.
#' If \code{split_files=TRUE}, all query results will instead be split into
#'  smaller files and stored in \emph{<output_dir>/}.
#' @param qtl_search This function will automatically search for any datasets
#' that match your criterion.
#' For example, if you search "Alasoo_2018", it will query the datasets:
#' \itemize{
#' \item{Alasoo_2018.macrophage_naive}
#' \item{Alasoo_2018.macrophage_Salmonella}
#' \item{Alasoo_2018.macrophage_IFNg+Salmonella}
#' }
#' You can be more specific about which datasets you want to include,
#' for example by searching: "Alasoo_2018.macrophage_IFNg".
#' You can even search by tissue or condition type
#' (e.g. \code{c("blood","brain")})
#' and any QTL datasets containing those substrings (case-insensitive)
#' in their name or metadata will be queried too.
#' @param method Method for querying eQTL Catalogue:
#' \itemize{
#' \item{"REST" (default): }{Uses the REST API. Slow but can be used by anyone.}
#' \item{"tabix"}{Uses tabix \link[echotabix]{query}.
#' Fast, but requires the user to first get their IP address whitelisted 
#' by the EMBL-EBI server admin by putting in a request
#' \href{https://www.ebi.ac.uk/about/contact/support/}{here}.} 
#' }
#' \emph{Note}: "tabix" is about ~17x faster than the REST API,
#'  but is currently a far less reliable method than the REST API because 
#'  tabix tends to get blocked by eQTL Catalogue's firewall.
#' See \href{https://github.com/RajLabMSSM/catalogueR/issues/5}{here}
#' for more details.
#' @param multithread_tabix Multi-thread across within a single tabix file query
#'  (good when you have one-several large loci).
#' @param nThread The number of CPU cores you want to use to speed up your
#' queries through parallelization.
#' @param split_files Save the results as one file per QTL dataset
#' (with all loci within each file).
#' If this is set to \code{=TRUE}, then this function will return the list of
#' paths where these files were saved.
#' A helper function is provided to import and merge them back together in R.
#' If this is set to  \code{=FALSE}, then this function will instead return one
#' big merged \link[data.table]{data.table}
#' containing results from all QTL datasets and all loci.
#' \code{=FALSE} is not recommended when you have many large loci and/or many
#' QTL datasets,
#' because you can only fit so much data into memory.
#' @param quant_method eQTL Catalogue actually contains more than just
#'  eQTL data.
#' For each dataset, the following kinds of QTLs can be queried:
#' \describe{
#' \item{gene expression QTL}{\code{quant_method="ge"} (\emph{default})
#' or \code{quant_method="microarray"}, depending on the dataset.
#' \strong{catalogueR} will automatically select whichever option is available.}
#' \item{exon expression QTL}{\emph{*under construction*}
#' \code{quant_method="ex"}}
#' \item{transcript usage QTL}{\emph{*under construction*}
#' \code{quant_method="tx"}}
#' \item{promoter, splice junction and 3' end usage QTL}{
#' \emph{*under construction*}  \code{quant_method="txrev"}}
#' }
#' @param merge_with_gwas Whether you want to merge your QTL query results
#' with your GWAS data
#' (convenient, but takes up more storage).
#' @param force_new_subset By default, \strong{catalogueR} will use any
#' pre-existing files that match your query.
#' Set \code{force_new_subset=T} to override this and force a new query.
#' @param query_genome The genome build of your query coordinates
#' (e.g. \code{query_dat}).
#' If your coordinates are in \emph{hg19}, \strong{catalogueR} will
#' automatically lift them over
#' to \emph{hg38} (as this is the build that eQTL Catalogue uses).
#' @param conda_env Conda environment to search for tabix executable in. 
#' @param verbose Show more (\code{=TRUE}) or fewer (\code{=FALSE}) messages.
#' @family eQTL Catalogue
#'
#' @export
#' @importFrom data.table fwrite rbindlist
#' @examples
#' sumstats_paths <- echodata::get_Nalls2019_loci(limit_snps = 5)
#' GWAS.QTL <- catalogueR::eQTLcatalogue_query(
#'     sumstats_paths = sumstats_paths$BST1,
#'     qtl_search = "Alasoo_2018.macrophage_naive")
eQTLcatalogue_query <- function(sumstats_paths = NULL,
                                output_dir = file.path(
                                    tempdir(),
                                    "catalogueR_queries"
                                ),
                                qtl_search = NULL,  
                                # multithread_qtl = TRUE,
                                # multithread_loci = FALSE,
                                multithread_tabix = FALSE,
                                method = c("REST","tabix"),
                                quant_method = "ge",
                                split_files = TRUE,
                                merge_with_gwas = TRUE,
                                force_new_subset = FALSE,
                                query_genome = "hg19",
                                conda_env = "echoR_mini",
                                nThread = 1,
                                verbose = TRUE) {
    # Args hidden from user in favor of automatic optimizer
    # @param multithread_qtl Multi-thread across QTL datasets 
    # (good when you're querying lots of QTL datasets).
    # @param multithread_loci Multi-thread across loci
    # (good when you have lots of gwas loci).

    beta.QTL <- NULL;
    ## Be sure to close all connections first 
    ## (can clog up tabix queries?)
    base::closeAllConnections() 

    query.start <- Sys.time()
    #### Search metadata for matching datasets ####
    qtl_datasets <- eQTLcatalogue_search_metadata(
        qtl_search = qtl_search,
        verbose = verbose
    )
    #### Determine multi-threading level ####
    if (nThread==1) {
        multithread_qtl <-
            multithread_loci <-
            multithread_tabix <- FALSE
    } else {
        if (isTRUE(multithread_tabix)) {
            messager("++ Multi-threading within tabix.", v = verbose)
            multithread_qtl <- multithread_loci <- FALSE
        } else {
            multithread_opts <- multithread_optimizer(
                qtl_datasets = qtl_datasets,
                sumstats_paths = sumstats_paths
            )
            multithread_qtl <- multithread_opts$qtl
            multithread_loci <- multithread_opts$loci
            multithread_tabix <- multithread_opts$tabix
        }
    }
    #### Query eQTL Catalogue ####
    messager("eQTL_Catalogue:: Querying",
        formatC(length(qtl_datasets), big.mark = ","),
        "QTL datasets x", formatC(length(sumstats_paths), big.mark = ","),
        "GWAS loci",
        paste0("(", formatC(length(qtl_datasets) * length(sumstats_paths),
            big.mark = ","
        ), " total)"),
        v = verbose
    )
    #### Iterate over QTL datasets ####
    GWAS.QTL_all <- parallel::mclapply(qtl_datasets,
                                       function(qtl_id) {
        messager(qtl_id,v=verbose)
        GWAS.QTL <- eQTLcatalogue_iterate_fetch(
            sumstats_paths = sumstats_paths,
            output_dir = output_dir,
            qtl_id = qtl_id,
            quant_method = quant_method, 
            method = method,
            conda_env = conda_env,
            nThread = nThread,
            split_files = split_files,
            merge_with_gwas = merge_with_gwas,
            force_new_subset = force_new_subset,
            query_genome = query_genome,
            multithread_loci = multithread_loci,
            multithread_tabix = multithread_tabix,
            verbose = verbose
        )
        return(GWAS.QTL)
    }, mc.cores = if (multithread_qtl) nThread else 1)
    # END ITERATION OVER QTL_IDS

    #### Gather paths/results ####
    if (isTRUE(split_files)) {
        # Split protocol
        messager("++ Returning list of split files paths.", v = verbose)
        GWAS.QTL_all <- unlist(GWAS.QTL_all)
    } else {
        messager("++ Post-processing merged results.", v = verbose)
        # Merged protocol
        GWAS.QTL_all <- data.table::rbindlist(GWAS.QTL_all, fill = TRUE)
        # Remove NAs
        if ("beta.QTL" %in% colnames(GWAS.QTL_all)) {
            GWAS.QTL_all <- subset(GWAS.QTL_all, !is.na(beta.QTL))
        }
        #### Check for completeness ####
        if (!is.null(names(sumstats_paths))) {
            missing_queries <- check_missing_queries(
                qtl_datasets = qtl_datasets,
                loci = names(sumstats_paths),
                GWAS.QTL = GWAS.QTL_all,
                verbose = verbose
            )
        }
        #### Save ####
        if (!is.null(output_dir)) {
            output_file <- file.path(output_dir, "eQTLcatalogue_tsv.gz")
            messager("+ Saving merged query results ==>", output_file,
                v = verbose
            )
            dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
            data.table::fwrite(
                x = GWAS.QTL_all,
                file = output_file,
                nThread = nThread
            )
        }
    } # END MERGED PROTOCOL 
    check_dim(df = GWAS.QTL_all) 
    #### Clock it ####
    query.end <- Sys.time()
    messager(round(query.end - query.start, 1), v = verbose)
    return(GWAS.QTL_all)
}

#### Deprecation function #####
eQTL_Catalogue.query <- function(...){
  .Deprecated("eQTLcatalogue_query")
  eQTLcatalogue_query(...)
}

