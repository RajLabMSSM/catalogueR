#' Query eQTL Catalogue:tabix/echotabix
#'
#' Query eQTL Catalogue datasets by region with tabix or \pkg{echotabix}.
#' Faster alternative to REST API.
#' @inheritParams eQTLcatalogue_query
#' @source \href{https://github.com/RajLabMSSM/catalogueR/issues/5}{
#' eQTL Catalogue blocking tabix requests}
#' @source
#' \code{
#' data("meta");
#' query_granges <- echotabix::construct_query(query_dat = echodata::BST1)
#' qtl.subset <- catalogueR:::fetch_tabix(unique_id=meta$unique_id[2],
#'                                        query_granges=query_granges)
#' }
#' @inheritParams eQTLcatalogue_query
#' @inheritParams echotabix::query
#' @family eQTL Catalogue 
#' @keywords internal
#' @importFrom echotabix query
fetch_tabix <- function(unique_id,
                        query_granges, 
                        query_method = c("rsamtools", 
                                         "variantannotation",
                                         "conda", 
                                         "seqminer"), 
                        quant_method = "ge",
                        nThread = 1,
                        conda_env = "echoR_mini",
                        verbose = TRUE) {
  # echoverseTemplate:::args2vars(catalogueR:::fetch_tabix)
  # echoverseTemplate:::source_all()
   
  #### Warning message ####
    msg <- paste(
      "WARNING: Querying eQTL Catalogue with tabix will only work",
      "if your IP address has been whitelisted by",
      "an EMBL-EBI server administrator.",
      "Otherwise, it will eventually be blocked.",
      "Please request access via this form:",
      "https://www.ebi.ac.uk/about/contact/support/")
    messager(msg,v=verbose)
    
    tabix.start <- Sys.time()  
    meta.sub <- choose_quant_method(
        ui = unique_id,
        qm = quant_method,
        verbose = verbose
    ) 
    #### Run tabix #### 
    qtl.subset <- echotabix::query(target_path = fix_ftp(meta.sub$ftp_path),
                                   query_granges = query_granges,
                                   query_method = query_method,
                                   verbose = verbose) 
    colnames(qtl.subset)[-1] <- eQTLcatalogue_get_header()
    tabix.end <- Sys.time()
    messager("catalogueR::", 
             formatC(nrow(qtl.subset),big.mark = ","),"SNPs returned in",
        round(as.numeric(tabix.end - tabix.start), 1), "seconds.",
        v = verbose
    )
    return(qtl.subset)
}
