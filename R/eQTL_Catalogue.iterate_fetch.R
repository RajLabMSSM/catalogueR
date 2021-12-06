#' Iterate queries to \emph{eQTL Catalogue}
#'
#' Uses coordinates from stored summary stats files (e.g. GWAS)
#' to determine which regions to query from \emph{eQTL Catalogue}.
#' @source
#' \code{
#' sumstats_paths <-catalogueR::example_sumstats_paths()
#' qtl_id <- catalogueR::eQTL_Catalogue.list_datasets()$unique_id[1]
#' GWAS.QTL <- catalogueR:::eQTL_Catalogue.iterate_fetch(
#'     sumstats_paths = sumstats_paths,
#'      qtl_id = qtl_id,
#'      nThread = 1,
#'      split_files = FALSE,
#'      progress_bar = FALSE)
#' }
#' @inheritParams eQTL_Catalogue.query
#' @param progress_bar Show progress bar during parallelization across loci.
#'  \emph{WARNING!}: Progress bar (via \link[pbmcapply]{progressBar}) only works
#'   on Linux/Unix systems (e.g. mac) and NOT on Windows.
#' @family eQTL Catalogue
#' @keywords internal
#' @importFrom echotabix liftover
eQTL_Catalogue.iterate_fetch <- function(sumstats_paths,
                                         output_dir = file.path(
                                             tempdir(),
                                             "catalogueR_queries"
                                         ),
                                         qtl_id,
                                         method = c("REST","tabix","echotabix"),
                                         quant_method = "ge",
                                         infer_region = TRUE, 
                                         conda_env = "echoR",
                                         multithread_loci = TRUE,
                                         multithread_tabix = FALSE,
                                         nThread = 1,
                                         split_files = TRUE,
                                         merge_with_gwas = FALSE,
                                         force_new_subset = FALSE,
                                         progress_bar = FALSE,
                                         genome_build = "hg19",
                                         verbose = TRUE) {
    CHR <- BP <- NULL;
    # WARNING!: pbmclapply only worked on Linux/Unix systems (e.g. mac)
    # and NOT on Windows.
    lapply_func <- progress_bar_check(
        progress_bar = progress_bar,
        verbose = verbose
    )
    ########## -----# ITERATE ACROSS LOCI ---------------------
    GWAS.QTL <- lapply_func(seq_len(length(sumstats_paths)), 
                            function(i,
                                      .sumstats_paths = sumstats_paths,
                                      .qtl_id = qtl_id,
                                      .method = method,
                                      .quant_method = quant_method,
                                      .infer_region = infer_region,
                                      .nThread = nThread,
                                      .multithread_tabix = multithread_tabix, 
                                      .force_new_subset = force_new_subset,
                                      .genome_build = genome_build,
                                      .split_files = split_files,
                                      .merge_with_gwas = merge_with_gwas,
                                      .verbose = verbose) {

        #### Import GWAS data ####
        # Have to do it this way in order to get names of sumstats_paths
        loc_path <- .sumstats_paths[i]
        gwas_data <- data.table::fread(loc_path, nThread = 1)
        # get rid of "chr" just in case
        gwas_data$CHR <- gsub("chr", "", gwas_data$CHR)

        #### Name Locus ####
        if (!is.null(loc_path) && (!is.null(names(loc_path)))) {
            messager("++ Extracting locus name from `sumstats_paths` names.",
                v = .verbose
            )
            loc <- names(loc_path)
            if (!"Locus" %in% colnames(gwas_data)) {
                gwas_data <- cbind(
                    Locus = loc,
                    gwas_data
                )
            }
        } else {
            if ("Locus" %in% colnames(gwas_data)) {
                messager("++ Extracting locus name from GWAS file.",
                    v = .verbose
                )
                loc <- unique(gwas_data$Locus)[1]
            } else {
                loc <- construct_locus_name(gwas_data, verbose = .verbose)
                gwas_data <- cbind(Locus = loc, gwas_data)
            }
        }

        messager("_+_+_+_+_+_+_+_+_--- Locus: ", loc, " ---_+_+_+_+_+_+_+_+_",
            v = verbose
        )
        # Test if query file already exists
        split_path <- make_split_path(
            output_dir = output_dir,
            qtl_id = .qtl_id,
            loc = loc
        )
        if (file.exists(split_path) && .force_new_subset == FALSE) {
            messager("++ Using pre-existing file...", v = .verbose)
            if (split_files) {
                return(split_path)
            } else {
                gwas.qtl <- data.table::fread(split_path, nThread = 1)
                return(gwas.qtl)
            }
        } else {
            #### Convert from GRCh37 to GRCh38 ####
            if (.genome_build %in% c("hg19", "hg18")) { 
                gwas_data <- echotabix::liftover( 
                    sumstats_dt = gwas_data,
                    ref_genome = .genome_build,
                    convert_ref_genome = "HG38",
                    verbose = verbose
                ) %>%  
                    dplyr::rename(POS = BP) %>%
                    dplyr::mutate(
                        CHR = gsub("chr", "", CHR)
                    )
            }
            qtl.subset <- data.table::data.table()
            qtl.subset <- tryCatch(
                expr = {
                    eQTL_Catalogue.fetch(
                        unique_id = .qtl_id,
                        method = .method,
                        quant_method = .quant_method,
                        infer_region = .infer_region,
                        gwas_data = gwas_data, 
                        nThread = .nThread,
                        multithread_tabix = .multithread_tabix, 
                        conda_env = conda_env,
                        chrom = NULL,
                        bp_upper = NULL,
                        bp_lower = NULL,
                        add_qtl_id = TRUE,
                        convert_genes = TRUE,
                        verbose = .verbose
                    )
                },
                error = function(e) {
                    messager(e)
                    data.table::data.table()
                }
            )
            check_dim(df = qtl.subset)
            #### Merge results ####
            if (.merge_with_gwas) {
                gwas.qtl <- tryCatch(
                    expr = {
                        merge_gwas_qtl(
                            gwas_data = gwas_data,
                            qtl.subset = qtl.subset,
                            verbose = .verbose
                        )
                    },
                    error = function(e) {
                        messager(e)
                        qtl.subset
                    }
                )
            } else {
                gwas.qtl <- qtl.subset
            }
            check_dim(df = gwas.qtl)
            if (dim(gwas.qtl)[1] == 0) {
                messager("WARNING: Data dimensions are 0 x 0. Returning NULL")
                return(NULL)
            }
            gwas.qtl <- tryCatch(
                expr = {
                    #### Add locus name ####
                    messager("++ Adding `Locus.GWAS` column.", v = .verbose)
                    gwas.qtl <- cbind(
                        Locus.GWAS = loc,
                        gwas.qtl
                    )
                    check_dim(df = gwas.qtl)
                    #### Save split ####
                    if (.split_files) {
                        messager("++ Saving split file ==>", split_path,
                            v = .verbose
                        )
                        dir.create(dirname(split_path),
                            showWarnings = FALSE, recursive = TRUE
                        )
                        data.table::fwrite(gwas.qtl, split_path,
                            sep = "\t",
                            nThread = 1
                        )
                    }
                    gwas.qtl
                },
                error = function(e) {
                    message(e)
                    data.table::data.table()
                }
            )
            #### Return ####
            if (.split_files) {
                return(split_path)
            } else {
                return(gwas.qtl)
            }
        }
    }, mc.cores = if (multithread_loci) nThread else 1)
    ## END ITERATE ACROSS LOCI

    #### Return ####
    if (split_files) {
        messager("+ Returning list of split query results files.", v = verbose)
        return(unlist(GWAS.QTL))
    } else {
        check_dim(df = GWAS.QTL)
        messager("+ Returning merged data.table of query results.", v = verbose)
        messager(nrow(GWAS.QTL), "x", ncol(GWAS.QTL), v = verbose)
        GWAS.QTL <- data.table::rbindlist(GWAS.QTL, fill = TRUE)
        return(GWAS.QTL)
    }
    message(" ")
}
