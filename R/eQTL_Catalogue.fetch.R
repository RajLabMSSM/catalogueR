#' Query eQTL Catalogue 
#'
#' Query eQTL Catalogue datasets with multiple methods options.
#' 
#' @param unique_id Unique eQTL Catalogue ID assigned in metadata
#'  ("unique_id" column in \code{data(meta)}). 
#' @param gwas_data \link[data.table]{data.table} of GWAS summary statistics.
#' @param add_qtl_id Add "qtl_id" (i.e. "unique_id") column to the query result.
#' @param convert_genes Convert Ensembl IDs to HGNC symbols.
#' @param chrom Chromosome of the query window.
#' @param bp_lower Minimum basepair position of the query window.
#' @param bp_upper Maxmimum basepair position of the query window.
#' @inheritParams eQTL_Catalogue.query
#' 
#' @family eQTL Catalogue
#' @export
#' @importFrom data.table data.table :=
#' @examples
#' data("meta")
#' GWAS.QTL_manual <- catalogueR:: eQTL_Catalogue.fetch(
#'     unique_id = meta$unique_id[1],  
#'     chrom = 8, 
#'     bp_lower = 21527069-500,
#'     bp_upper = 21527069+500)
eQTL_Catalogue.fetch <- function(unique_id,
                                 method = c("REST","tabix","echotabix"),
                                 quant_method = "ge",
                                 infer_region = TRUE,
                                 gwas_data = NULL, 
                                 nThread = 1, 
                                 chrom = NULL,
                                 bp_lower = NULL,
                                 bp_upper = NULL,
                                 multithread_tabix = FALSE,
                                 add_qtl_id = TRUE,
                                 convert_genes = TRUE,
                                 conda_env = "echoR",
                                 verbose = TRUE) {
    
    method <- tolower(method[1]) 
    use_tabix <- method %in% c("tabix","echotabix")
    if (use_tabix) {
        ## Tabix is about ~17x faster than the REST API,
        ## but gets blocked often by EBI server which thinks its an attack.
        gwas.qtl <- fetch_tabix(
            unique_id = unique_id,
            method = method,
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
    #### Post=processing ####
    if (add_qtl_id) {
        gwas.qtl <- data.table::data.table(
            qtl_id = unique_id,
            gwas.qtl
        )
    }
    #### Convert genes ####
    if (convert_genes) { 
        gene_dict <- ensembl_to_hgnc(
            ensembl_ids = gwas.qtl$gene_id.QTL,
            verbose = verbose
        ) 
        gwas.qtl$gene.QTL <- gene_dict[
            gwas.qtl$gene_id.QTL
        ]  
    }
    return(gwas.qtl)
}
