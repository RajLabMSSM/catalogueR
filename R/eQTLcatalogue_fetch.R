#' Query eQTL Catalogue 
#'
#' Query eQTL Catalogue datasets with multiple methods options.
#' 
#' @param dataset_id Unique eQTL Catalogue ID assigned in metadata
#'  ("dataset_id" column in \code{data(meta)}). 
#' @param query_dat \link[data.table]{data.table} of GWAS summary statistics.
#' @param add_dataset_id Add "dataset_id" (i.e. "dataset_id") column to the query result.
#' @param convert_genes Convert Ensembl IDs to HGNC symbols.
#' @param chrom Chromosome of the query window.
#' @param bp_lower Minimum basepair position of the query window.
#' @param bp_upper Maxmimum basepair position of the query window.
#' @inheritParams eQTLcatalogue_query
#' 
#' @family eQTL Catalogue
#' @export
#' @importFrom data.table data.table :=
#' @examples
#' data("meta")
#' query_granges <- echodata::BST1
#' GWAS.QTL_manual <- catalogueR::eQTLcatalogue_fetch(
#'   query_granges = query_granges,
#'   dataset_id = meta$dataset_id[1])
eQTLcatalogue_fetch <- function(dataset_id,
                                query_granges,
                                method = c("REST","tabix"),
                                quant_method = "ge", 
                                multithread_tabix = FALSE,
                                add_dataset_id = TRUE,
                                convert_genes = TRUE,
                                suffix = ".QTL",
                                timeout = 5*60,
                                conda_env = "echoR_mini",
                                nThread = 1, 
                                verbose = TRUE) { 
  dataset_id <- NULL;
  method <- tolower(method)[1]  
  if (method=="tabix") {
      ## Tabix is about ~17x faster than the REST API,
      ## but gets blocked often by EBI server which thinks its an attack.
    qtl.subset <- fetch_tabix(
          dataset_id = dataset_id,
          query_granges = query_granges, 
          quant_method = quant_method,  
          nThread = if (isTRUE(multithread_tabix)) nThread else 1,
          conda_env = conda_env,
          verbose = verbose)
  } else {
    qtl.subset <- fetch_restAPI(
          dataset_id = dataset_id,
          quant_method = quant_method, 
          query_granges = query_granges, 
          timeout = timeout,
          verbose = verbose)
  }
  #### Add suffix to all columns ####
  if(!is.null(suffix)){
    colnames(qtl.subset) <- paste0(colnames(qtl.subset), suffix)
  } 
  #### Post=processing ####
  if (isTRUE(add_dataset_id)) {
    qtl.subset[,dataset_id:=dataset_id]
  }
  #### Convert genes ####
  if (isTRUE(convert_genes)) { 
      gene_dict <- ensembl_to_hgnc(
          ensembl_ids = qtl.subset$gene_id.QTL,
          verbose = verbose) 
      qtl.subset$gene.QTL <- gene_dict[qtl.subset$gene_id.QTL]  
  }
  return(qtl.subset)
}

#### Deprecation function #####
eQTL_Catalogue.fetch <- function(...){
  .Deprecated("eQTLcatalogue_fetch")
  eQTLcatalogue_fetch(...)
}
