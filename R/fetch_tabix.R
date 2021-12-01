#' 2. Query eQTL Catalogue datasets by region
#'
#' 2.1 Method 1: Tabix
#' Faster alternative to REST API.
#' @source
#' \code{
#' data("meta"); 
#' gwas_data <- echodata::BST1 
#' qtl.subset <- catalogueR:::fetch_tabix(unique_id=meta$unique_id[2],
#'                                        gwas_data=gwas_data,
#'                                        conda_env=NULL)
#' }
#' @inheritParams eQTL_Catalogue.query
#' @family eQTL Catalogue
#'
#' @keywords internal
#' @importFrom echotabix query_tabular 
fetch_tabix <- function(unique_id,
                        quant_method = "ge",
                        infer_region = TRUE,
                        gwas_data = NULL,
                        chrom = NULL,
                        bp_lower = NULL,
                        bp_upper = NULL,
                        is_gwas = FALSE,
                        nThread = 4,
                        conda_env = NULL,
                        verbose = TRUE) {
    # quant_method="ge"; infer_region=T;is_gwas=F; chrom=NULL; remove_tmp=F;  
  # add_chr=T; bp_lower=bp_upper=NULL; verbose=T; conda_env="echoR"; 
  # nThread=10; unique_id=meta$unique_id[2]; gwas_data=BST1;
    check_coord_input(
        gwas_data = gwas_data,
        chrom = chrom,
        bp_lower = bp_lower, 
        bp_upper = bp_upper
    )
    tabix.start <- Sys.time()
    # Get region
    if (infer_region & !is.null(gwas_data)) {
        messager("++ Inferring coordinates from gwas_data", v = verbose)
        chrom <- gsub("chr", "", unique(gwas_data$CHR))
        if (length(chrom) > 1) {
            stop("More than one chromosome detected.")
        }
        bp_lower <- min(gwas_data$POS, na.rm = TRUE)
        bp_upper <- max(gwas_data$POS, na.rm = TRUE)
    }
    region <- paste0(chrom, ":", bp_lower, "-", bp_upper)
    # messager("+ TABIX:: Querying region:", region)
    meta.sub <- choose_quant_method(
        ui = unique_id,
        qm = quant_method,
        verbose = verbose
    )
    # You have to get the tabix header separately
    ## since tabix has issues getting it from eQTL Catalogue.
    header <- tryCatch(expr = {
        tabix_header(
            tabix_path = meta.sub$ftp_path,
            force_new_header = FALSE
        )
    }, error = function(e) {
        tabix_header(
            tabix_path = meta.sub$ftp_path,
            force_new_header = FALSE
        )
    })

    # Run tabix
    #### Method 1 ####
    qtl.subset <- tabix_reader(
        tabix_path = meta.sub$ftp_path,
        region = region,
        conda_env = conda_env,
        nThread = nThread,
        verbose = verbose
    )
    #### Method 2 #### 
    # chrom <- 12; bp_lower=40189191; bp_upper=40334626
    # qtl.subset <- echotabix::query_tabular(fullSS_tabix = meta.sub$ftp_path,
    #                                        chrom = chrom,
    #                                        start_pos = bp_lower,
    #                                        end_pos = bp_upper,
    #                                        local = FALSE,
    #                                        verbose = TRUE)
    if (length(header) != ncol(qtl.subset)) {
        header <- tabix_header(
            tabix_path = meta.sub$ftp_path,
            force_new_header = TRUE
        )
    }

    colnames(qtl.subset) <- paste0(header, ".QTL")
    tabix.end <- Sys.time()
    messager("eQTL_Catalogue::", nrow(qtl.subset), "SNPs returned in",
        round(as.numeric(tabix.end - tabix.start), 1), "seconds.",
        v = verbose
    )
    return(qtl.subset)
}