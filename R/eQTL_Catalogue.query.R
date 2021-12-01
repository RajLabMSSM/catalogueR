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
#' saved to (set \code{output_dir=FALSE} if you don't want to save the results).
#' If \code{split_files=FALSE}, all query results will be merged into one and 
#' saved as \emph{<output_dir>/eQTL_Catalogue.tsv.gz}.
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
#' @param use_tabix Use Tabix to makes queries (Default: \emph{FALSE}). 
#' Tabix is about ~17x faster than the REST API, but is currently a far less 
#' reliable method than the REST API because it tends to get blocked by
#' eQTL Catalogue's firewall. 
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
#' @param genome_build The genome build of your query coordinates 
#' (e.g. \code{gwas_data}).
#' If your coordinates are in \emph{hg19}, \strong{catalogueR} will 
#' automatically lift them over
#' to \emph{hg38} (as this is the build that eQTL Catalogue uses).
#' @param progress_bar \code{progress_bar=T} allows progress to be monitored
#'  even when multithreading enabled.
#' Requires R package \code{\link{pbmcapply}}.
#' @param verbose Show more (\code{=TRUE}) or fewer (\code{=FALSE}) messages.
#' @family eQTL Catalogue
#'
#' @export
#' @importFrom data.table fwrite rbindlist
#' @examples
#' sumstats_paths <- catalogueR::example_sumstats_paths()
#' 
#' # Merged results
#' GWAS.QTL <- eQTL_Catalogue.query(sumstats_paths = sumstats_paths, 
#'                                  qtl_search = "Alasoo_2018",  
#'                                  merge_with_gwas = FALSE, 
#'                                  split_files = FALSE)
#' 
#' # Split results
#' gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths = sumstats_paths,
#'                                        qtl_search = "Alasoo_2018",
#'                                        merge_with_gwas = FALSE,
#'                                        progress_bar = TRUE)
#' GWAS.QTL <- gather_files(file_paths = gwas.qtl_paths)
eQTL_Catalogue.query <- function(sumstats_paths = NULL,
                                 output_dir = file.path(
                                     tempdir(),
                                     "catalogueR_queries"
                                 ),
                                 qtl_search = NULL,
                                 use_tabix = FALSE,
                                 conda_env = "echoR",
                                 nThread = 1,
                                 # multithread_qtl = TRUE,
                                 # multithread_loci = FALSE,
                                 multithread_tabix = FALSE,
                                 quant_method = "ge",
                                 infer_region = TRUE,
                                 split_files = TRUE,
                                 merge_with_gwas = TRUE,
                                 force_new_subset = FALSE,
                                 genome_build = "hg19",
                                 progress_bar = TRUE,
                                 verbose = TRUE) {
    # sumstats_paths <-  example_sumstats_paths()
    # merge_with_gwas=F; nThread=4; loci_names = basename(dirname(dirname(sumstats_paths))); qtl.id=qtl_datasets; multithread_qtl=T;  multithread_loci=F; quant_method="ge";   split_files=T; multithread_tabix=F; split_files=T;
    # qtl_search= c("Fairfax_2014") #c("ROSMAP","Alasoo_2018","Fairfax_2014", "Nedelec_2016","BLUEPRINT","HipSci.iPSC", "Lepik_2017","BrainSeq","TwinsUK","Schmiedel_2018", "blood","brain")
    # output_dir = "eQTL_Catalogue/Nalls23andMe_2019"; qtl.id="Alasoo_2018.macrophage_naive"; force_new_subset=F; progress_bar=F; verbose=T; genome_build="hg19";

    # Aargs hidden from user in favor of automatic optimizer
    # @param multithread_qtl Multi-thread across QTL datasets (good when you're querying lots of QTL datasets).
    # @param multithread_loci Multi-thread across loci (good when you have lots of gwas loci).

    ## Be doubly sure to close all connections first (can clog up tabix queries?)
    base::closeAllConnections()
    base::closeAllConnections()

    query.start <- Sys.time()
    # Setup
    cleanup_tbi(DIR = dirname(output_dir))
    lapply_func <- progress_bar_check(
        progress_bar = progress_bar,
        verbose = verbose
    )
    # Search metadata for matching datasets
    qtl_datasets <- eQTL_Catalogue.search_metadata(
        qtl_search = qtl_search,
        verbose = verbose
    )

    #### Determine multi-threading level ####
    if (multithread_tabix) {
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
    #### Query eQTL Catalogue ####
    messager("eQTL_Catalogue:: Querying",
             formatC(length(qtl_datasets), big.mark = ","), 
             "QTL datasets x", formatC(length(sumstats_paths), big.mark = ","),
             "GWAS loci",
        paste0("(", formatC(length(qtl_datasets) * length(sumstats_paths),
                            big.mark = ","), " total)"),
        v = verbose
    )
    #### Iterate over QTL datasets ####
    GWAS.QTL_all <- lapply_func(qtl_datasets, function(
        qtl_id,
       .sumstats_paths = sumstats_paths,
       .output_dir = output_dir,
       .quant_method = quant_method,
       .infer_region = infer_region,
       .use_tabix = use_tabix,
       .conda_env = conda_env,
       .nThread = nThread,
       .split_files = split_files,
       .merge_with_gwas = merge_with_gwas,
       .force_new_subset = force_new_subset,
       .genome_build = genome_build,
       .multithread_loci = multithread_loci,
       .multithread_tabix = multithread_tabix,
       .verbose = verbose) {
        message(qtl_id)
        GWAS.QTL <- data.table::data.table()
        # try({
        GWAS.QTL <- eQTL_Catalogue.iterate_fetch(
            sumstats_paths = .sumstats_paths,
            output_dir = .output_dir,
            qtl_id = qtl_id,
            quant_method = .quant_method,
            infer_region = .infer_region,
            use_tabix = .use_tabix,
            conda_env = .conda_env,
            nThread = .nThread,
            split_files = .split_files,
            merge_with_gwas = .merge_with_gwas,
            force_new_subset = .force_new_subset,
            genome_build = .genome_build,
            multithread_loci = .multithread_loci,
            multithread_tabix = .multithread_tabix,
            progress_bar = progress_bar,
            verbose = .verbose
        )
        # })

        return(GWAS.QTL)
        # END ITERATION OVER QTL_IDS
    }, mc.cores = if (multithread_qtl) nThread else 1)

    #### Gather paths/results ####
    if (split_files) {
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
        if (output_dir != FALSE) {
            output_file <- file.path(output_dir, "eQTL_Catalogue.tsv.gz")
            messager("+ Saving merged query results ==>", output_file, 
                     v = verbose)
            dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
            data.table::fwrite(GWAS.QTL_all,
                file = output_file,
                nThread = nThread
            )
        }
    } # END MERGED PROTOCOL

    check_dim(df = GWAS.QTL_all)
    # Clean
    cleanup_tbi(DIR = dirname(output_dir))
    # Clock it
    query.end <- Sys.time()
    messager(round(query.end - query.start, 1),v=verbose)
    return(GWAS.QTL_all)
}
