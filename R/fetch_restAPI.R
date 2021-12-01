#' 2. Query eQTL Catalogue datasets by region
#'
#' 2.2 Method 2: RESTful API
#' Slower than tabix (unless you're only querying several specific SNPs).
#' @source 
#' \code{
#' data("meta")
#' gwas_data <- echodata::BST1 
#' qtl.subset <- fetch_restAPI(unique_id = meta$unique_id[1],
#'                             gwas_data = gwas_data)
#' }
#' @inheritParams eQTL_Catalogue.query
#' @family eQTL Catalogue
#' @keywords internal
#' @importFrom data.table data.table 
fetch_restAPI <- function(unique_id, # Alasoo_2018.macrophage_naive
                          quant_method = "ge",
                          infer_region = TRUE,
                          gwas_data = NULL,
                          chrom = NULL,
                          bp_lower = NULL,
                          bp_upper = NULL,
                          is_gwas = FALSE, # refers to the datasets being queried
                          size = NULL,
                          verbose = TRUE) {
    restAPI.start <- Sys.time()
    # Get region
    if (infer_region & !is.null(gwas_data)) {
        messager("+ Inferring coordinates from gwas_data", v = verbose)
        chrom <- gsub("chr", "", unique(gwas_data$CHR)[1])
        bp_lower <- min(gwas_data$POS)
        bp_upper <- max(gwas_data$POS)
    }
    meta.sub <- choose_quant_method(
        ui = unique_id,
        qm = quant_method,
        verbose = verbose
    )
    # gene_id <- hgnc_to_ensembl(gene_symbols = "LRRK2")
    link <- paste0(
        "http://www.ebi.ac.uk/eqtl/api/",
        "chromosomes/", chrom,
        "/associations?paginate=False",
        # Study name
        "&study=", meta.sub$study,
        # Condition
        "&qtl_group=", meta.sub$qtl_group,
        # ENSEMBL gene id
        # "&gene_id=",gene_id,
        # ENSEMBL molecular trait id
        # "&molecular_trait_id=",molecular_trait_id,
        # gene expression, transcript, etc.
        "&quant_method=", quant_method,
        # genomic position limits
        "&bp_lower=", bp_lower,
        "&bp_upper=", bp_upper,
        if (!is.null(size)) paste0("&size=", size) else ""
    )

    message("+ eQTL_Catalogue:: Fetching: ")
    message(link)
    qtl.subset <- fetch_from_eqtl_cat_API(link = link, is_gwas = is_gwas)
    colnames(qtl.subset) <- paste0(colnames(qtl.subset), ".QTL")
    restAPI.end <- Sys.time()
    messager("+ eQTL_Catalogue::", nrow(qtl.subset), 
             "SNPs returned in", 
             formatC(round(as.numeric(restAPI.end - restAPI.start), 1),
                     big.mark = ","), 
             "seconds.", v = verbose)
    return(data.table::data.table(qtl.subset))
}
