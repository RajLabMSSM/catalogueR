#' Iterate queries to \emph{eQTL Catalogue}
#'
#' Uses coordinates from stored summary stats files (e.g. GWAS)
#' to determine which regions to query from \emph{eQTL Catalogue}.
#' @source
#' \code{
#' sumstats_paths <- echodata::get_Nalls2019_loci(limit_snps = 5)
#' qtl_id <- catalogueR::eQTLcatalogue_list_datasets()$unique_id[1]
#' GWAS.QTL <- catalogueR:::eQTLcatalogue_iterate_fetch(
#'     sumstats_paths = sumstats_paths,
#'      qtl_id = qtl_id,
#'      nThread = 1,
#'      split_files = FALSE)
#' }
#' @inheritParams eQTLcatalogue_query 
#' @inheritParams echotabix::query
#' @family eQTL Catalogue
#' 
#' @keywords internal
#' @importFrom echotabix liftover
eQTLcatalogue_iterate_fetch <- function(sumstats_paths,
                                        output_dir = file.path(
                                            tempdir(),
                                            "catalogueR_queries"
                                        ),
                                        qtl_id,
                                        method = c("REST","tabix"),
                                        quant_method = "ge",
                                        multithread_loci = TRUE,
                                        multithread_tabix = FALSE, 
                                        split_files = TRUE,
                                        merge_with_gwas = FALSE,
                                        force_new_subset = FALSE,
                                        query_genome = "hg19",
                                        conda_env = "echoR_mini",
                                        nThread = 1,
                                        verbose = TRUE) {
  
    #### ITERATE ACROSS LOCI ####
    GWAS.QTL <- parallel::mclapply(seq_len(length(sumstats_paths)), 
                                   function(i) {
        #### Import GWAS data ####
        # Have to do it this way in order to get names of sumstats_paths
        loc_path <- sumstats_paths[[i]]
        query_dat <- data.table::fread(loc_path)
        # get rid of "chr" just in case
        query_dat$CHR <- gsub("chr", "", query_dat$CHR)

        #### Name Locus ####
        if (!is.null(loc_path) && (!is.null(names(loc_path)))) {
            messager(
              "++ Extracting locus name from `sumstats_paths` names.",
              v = verbose
            )
            loc <- names(loc_path)
            if (!"Locus" %in% colnames(query_dat)) {
                query_dat <- cbind(
                    Locus = loc,
                    query_dat
                )
            }
        } else {
            if ("Locus" %in% colnames(query_dat)) {
                messager("++ Extracting locus name from GWAS file.",
                    v = verbose
                )
                loc <- unique(query_dat$Locus)[1]
            } else {
                loc <- construct_locus_name(query_dat,
                                            verbose = verbose)
                query_dat <- cbind(Locus = loc, query_dat)
            }
        }

        messager("_+_+_+_+_+_+_+_+_--- Locus: ", loc,
                 " ---_+_+_+_+_+_+_+_+_",
            v = verbose
        )
        # Test if query file already exists
        split_path <- make_split_path(
            output_dir = output_dir,
            qtl_id = qtl_id,
            loc = loc
        )
        if (file.exists(split_path) && force_new_subset == FALSE) {
            messager("++ Using pre-existing file...", v = verbose)
            if (split_files) {
                return(split_path)
            } else {
                gwas.qtl <- data.table::fread(split_path)
                return(gwas.qtl)
            }
        } else {
            #### Convert from GRCh37 to GRCh38 ####
            if (query_genome %in% c("hg19", "hg18")) { 
                query_dat <- echotabix::liftover( 
                    dat = query_dat,
                    query_genome = query_genome,
                    target_genome = "HG38",
                    style = "NCBI",
                    as_granges = FALSE,
                    verbose = verbose)
            } 
            qtl.subset <- tryCatch(
                expr = {
                    eQTLcatalogue_fetch(
                        unique_id = qtl_id,
                        query_granges = query_dat,
                        method = method,
                        quant_method = quant_method,
                        multithread_tabix = multithread_tabix, 
                        add_qtl_id = TRUE,
                        convert_genes = TRUE,
                        nThread = nThread,
                        conda_env = conda_env, 
                        verbose = verbose
                    )
                },
                error = function(e) {
                    messager(e);
                    data.table::data.table()
                }
            )
            check_dim(df = qtl.subset)
            #### Merge results ####
            if (merge_with_gwas) {
                gwas.qtl <- tryCatch(
                    expr = {
                        merge_gwas_qtl(
                            query_dat = query_dat,
                            qtl.subset = qtl.subset,
                            verbose = verbose
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
                messager("WARNING: Data dimensions are 0 x 0.",
                         "Returning NULL",v=verbose)
                return(NULL)
            }
            gwas.qtl <- tryCatch(
                expr = {
                    #### Add locus name ####
                    messager("++ Adding `Locus.GWAS` column.",
                             v = verbose)
                    gwas.qtl <- cbind(
                        Locus.GWAS = loc,
                        gwas.qtl
                    )
                    check_dim(df = gwas.qtl)
                    #### Save split ####
                    if (split_files) {
                        messager("++ Saving split file ==>", 
                                 split_path,v = verbose)
                        dir.create(dirname(split_path),
                            showWarnings = FALSE, recursive = TRUE
                        )
                        data.table::fwrite(gwas.qtl, split_path,
                            sep = "\t")
                    }
                    gwas.qtl
                },
                error = function(e) {
                    message(e)
                    data.table::data.table()
                }
            )
            #### Return ####
            if (split_files) {
                return(split_path)
            } else {
                return(gwas.qtl)
            }
        }
    }, mc.cores = if (multithread_loci) nThread else 1)
    ## END ITERATE ACROSS LOCI

    #### Return ####
    if (split_files) {
        messager("+ Returning list of split query results files.",
                 v = verbose)
        return(unlist(GWAS.QTL))
    } else {
        check_dim(df = GWAS.QTL)
        messager("+ Returning merged data.table of query results.",
                 v = verbose)
        messager(nrow(GWAS.QTL), "x", ncol(GWAS.QTL),
                 v = verbose)
        GWAS.QTL <- data.table::rbindlist(GWAS.QTL, fill = TRUE)
        return(GWAS.QTL)
    }
    message(" ")
}

#### Deprecation function #####
eQTL_Catalogue.iterate_fetch <- function(...){
  .Deprecated("eQTLcatalogue_iterate_fetch")
  eQTLcatalogue_iterate_fetch(...)
}
