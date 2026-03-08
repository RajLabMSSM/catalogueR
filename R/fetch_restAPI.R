#' 2. Query eQTL Catalogue datasets by region
#'
#' 2.2 Method 2: RESTful API
#' Slower than tabix (unless you're only querying several specific SNPs).
#' @source
#' \code{
#' data("meta")
#' query_granges <- echotabix::construct_query(query_dat = echodata::BST1)
#' qtl.subset <- fetch_restAPI(dataset_id = meta$dataset_id[1],
#'                             query_granges = query_granges)
#' }
#' @inheritParams eQTLcatalogue_query
#' @inheritParams echotabix::construct_query
#' @family eQTL Catalogue
#' 
#' @keywords internal
#' @importFrom echotabix construct_query
#' @importFrom data.table data.table rbindlist
#' @importFrom GenomicRanges seqnames
fetch_restAPI <- function(dataset_id, # Alasoo_2018.macrophage_naive
                          query_granges,
                          quant_method = "ge", 
                          is_gwas = FALSE, # refers to the datasets being queried
                          size = NULL,
                          timeout = 5*60,
                          verbose = TRUE) {
  
    restAPI.start <- Sys.time()
    query_granges <- echotabix::construct_query(query_dat = query_granges,
                                                style = "NCBI",
                                                verbose = verbose)
    query_split <- GenomicRanges::split(
      x = query_granges,
      f = GenomicRanges::seqnames(query_granges),
      drop = TRUE
    )
    #### Chose quant method ####
    meta.sub <- choose_quant_method(
      di = dataset_id,
      qm = quant_method,
      verbose = verbose
    ) 
    #### Iterate over loci ####
    ### Via the main endpoint
    qtl <- lapply(query_split,
                  function(gr){
                    
      chrom <- as.character(unique(GenomicRanges::seqnames(gr)))[[1]]
      messager("Querying chromosome:",chrom,v=verbose)
      #### Via the chromosomes endpoint ####
      ## https://www.ebi.ac.uk/eqtl/api-docs/#get--eqtl-api-chromosomes-(string-chromosome) 
      link <- paste0(
        "https://www.ebi.ac.uk/eqtl/api/v3/",
        "chromosomes/",chrom,
        "/associations?paginate=False",
        # Study name
        "&study=", meta.sub$study_label,
        # Condition
        "&qtl_group=", meta.sub$sample_group,
        # ENSEMBL gene id
        # "&gene_id=",gene_id,
        # ENSEMBL molecular trait id
        # "&molecular_trait_id=",molecular_trait_id,
        # gene expression, transcript, etc.
        "&quant_method=", quant_method,
        # genomic position limits
        "&bp_lower=", min(GenomicRanges::start(gr)),
        "&bp_upper=", max(GenomicRanges::end(gr)),
        if (!is.null(size)) paste0("&size=", size) else ""
      ) 
      
      
      #### Via the top-level associations endpoint ####
      ## NOTE: Does not support bp_lower/upper params
      # link <- paste0(
      #   "https://www.ebi.ac.uk/eqtl/api/v3",
      #   "/associations?paginate=False",
      #   # Chromosome
      #   "&chromosome=", chrom,
      #   # Study name
      #   "&study=", meta.sub$study_label,
      #   # Condition
      #   "&qtl_group=", meta.sub$qtl_group,
      #   # ENSEMBL gene id
      #   # "&gene_id=",gene_id,
      #   # ENSEMBL molecular trait id
      #   # "&molecular_trait_id=",molecular_trait_id,
      #   # gene expression, transcript, etc.
      #   "&quant_method=", quant_method, 
      #   if (!is.null(size)) paste0("&size=", size) else ""
      # ) 
      # 
      # #### Via the datasets endpoint ####
      # ## NOTE: Does not support bp_lower/upper params
      # link <- paste0(
      #   "https://www.ebi.ac.uk/eqtl/api/v3",
      #   "/datasets/", meta.sub$dataset_id,
      #   "/associations?paginate=False",
      #   # Chromosome
      #   "&chromosome=", chrom,
      #   # Study name
      #   # "&study=", meta.sub$study_id,
      #   # Condition
      #   # "&qtl_group=", meta.sub$qtl_group,
      #   # ENSEMBL gene id
      #   # "&gene_id=",gene_id,
      #   # ENSEMBL molecular trait id
      #   # "&molecular_trait_id=",molecular_trait_id,
      #   # gene expression, transcript, etc.
      #   "&quant_method=", quant_method, 
      #   if (!is.null(size)) paste0("&size=", size) else ""
      # )
      # 
      messager("+ catalogueR:: Fetching:",link,v=verbose) 
      #### Post-process ####
      qtl.subset <- fetch_from_eqtl_cat_API(link = link, 
                                            is_gwas = is_gwas, 
                                            timeout = timeout) 
      #### Report ####
      restAPI.end <- Sys.time()
      messager("+ catalogueR::", 
               formatC(nrow(qtl.subset), big.mark = ","),
               "SNPs returned in",
               formatC(round(as.numeric(restAPI.end - restAPI.start), 1),
                       big.mark = ","),"seconds.",v = verbose
               )
      return(qtl.subset)
    }) |> data.table::rbindlist(fill = TRUE) 
    ### Return ####
    return(qtl)
}
