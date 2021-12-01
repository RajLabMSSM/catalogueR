#' 2. Query eQTL Catalogue datasets by region
#'
#' Choose between tabix (faster for larger queries)
#' or RESTful API (faster for small queries).
#' @inheritParams eQTL_Catalogue.query
#' @param bp_lower Minimum basepair position of the query window.
#' @param bp_upper Maxmimum basepair position of the query window.
#' @family eQTL Catalogue
#' @export
#' @importFrom data.table data.table
#' @examples
#' data("meta")
#' gwas_data <- echodata::BST1 
#' gwas.qtl <- eQTL_Catalogue.fetch(unique_id = meta$unique_id[1],
#'                                  gwas_data = gwas_data)
eQTL_Catalogue.fetch <- function(unique_id,
                                 quant_method = "ge",
                                 infer_region = TRUE,
                                 gwas_data = NULL,
                                 is_gwas = FALSE,
                                 nThread = 1,
                                 use_tabix = TRUE,
                                 chrom = NULL,
                                 bp_lower = NULL,
                                 bp_upper = NULL,
                                 multithread_tabix = FALSE,
                                 add_qtl_id = TRUE,
                                 convert_genes = TRUE,
                                 conda_env = "echoR",
                                 verbose = TRUE) {
    if (use_tabix) {
        # Tabix is about ~17x faster than the REST API.
        gwas.qtl <- fetch_tabix(
            unique_id = unique_id,
            quant_method = quant_method,
            infer_region = TRUE,
            gwas_data = gwas_data,
            chrom = chrom,
            bp_lower = bp_lower,
            bp_upper = bp_upper,
            is_gwas = FALSE,
            nThread = if (multithread_tabix) nThread else 1,
            conda_env = conda_env,
            verbose = verbose
        )
    } else {
        gwas.qtl <- fetch_restAPI(
            unique_id = unique_id,
            quant_method = quant_method,
            infer_region = TRUE,
            gwas_data = gwas_data,
            chrom = chrom,
            bp_lower = bp_lower,
            bp_upper = bp_upper,
            is_gwas = FALSE,
            verbose = verbose
        )
    }
    # Post=processing
    if (add_qtl_id) {
        gwas.qtl <- data.table::data.table(qtl_id = unique_id, gwas.qtl)
    }
    # Convert genes
    if (convert_genes) {
        try({
            gene_dict <- ensembl_to_hgnc(
                ensembl_ids = gwas.qtl$gene_id.QTL,
                verbose = verbose
            )
            gwas.qtl$gene.QTL <- gene_dict[gwas.qtl$molecular_trait_object_id.QTL]
        })
    }
    return(gwas.qtl)
}
