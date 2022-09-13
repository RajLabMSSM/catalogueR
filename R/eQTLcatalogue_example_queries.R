#' Paths to example eQTL Catalogue query results
#'
#' Returns the paths to eQTL Catalogue query results stored within
#' \pkg{catalogueR}. Each file is a merged \link[data.table]{data.table}
#'  of the GWAS summary stats used to make
#' the query, and the eQTL Catalogue query results
#' (which can contain data for multiple eGenes).
#'
#' GWAS data originally comes from the Parkinson's disease GWAS
#' by \href{https://doi.org/10.1016/S1474-4422(19)30320-5}{
#' Nalls et al., 2019 (The Lancet Neurology)}.
#' 
#' @param fnames Character vector of file names to download.
#' @inheritParams echodata::get_data
#' 
#' @source \url{https://doi.org/10.1016/S1474-4422(19)30320-5} 
#' @family Nalls23andMe_2019
#' @returns Path to merged example GWAS-QTL summary statistics.
#' @source 
#' \code{
#' paths <- catalogueR:::eQTLcatalogue_example_queries_local()
#' fnames <- lapply(paths, function(x){
#'     piggyback::pb_upload(file = x, 
#'                          repo = "RajLabMSSM/catalogueR")
#'     return(basename(x))
#' })
#' }
#' @export
#' @importFrom piggyback pb_list
#' @importFrom echodata get_data
#' @examples
#' gwas.qtl_paths <- catalogueR::eQTLcatalogue_example_queries()
eQTLcatalogue_example_queries <- function(
    save_dir = tempdir(),
    fnames = c("BST1__Alasoo_2018.macrophage_IFNg+Salmonella.tsv.gz",
              "BST1__Alasoo_2018.macrophage_IFNg.tsv.gz",
              "BST1__Alasoo_2018.macrophage_naive.tsv.gz",          
              "BST1__Alasoo_2018.macrophage_Salmonella.tsv.gz",     
              "LRRK2__Alasoo_2018.macrophage_IFNg.tsv.gz",          
              "LRRK2__Alasoo_2018.macrophage_naive.tsv.gz",         
              "MEX3C__Alasoo_2018.macrophage_IFNg.tsv.gz",          
              "MEX3C__Alasoo_2018.macrophage_naive.tsv.gz")) {
    
    files <- piggyback::pb_list(repo = "RajLabMSSM/catalogueR")
    tsv.gz <- grep("\\.tsv.gz",files$file_name,value = TRUE)
    fnames_select <- fnames[fnames %in% tsv.gz]
    fnames_drop <- fnames[!fnames %in% tsv.gz]
    if(length(fnames_drop)>0){
        messager(length(fnames_drop),
                 "file(s) could not be found and will be ignored:\n",
                 paste("-",fnames_drop,collapse = "\n"))
    } 
    if(length(fnames_select)>0){
        messager(length(fnames_select),"file(s) will be downloaded.")
        out_paths <- lapply(fnames_select, function(x){
            group <- parse_gwas.qtl_path(gwas.qtl_path = x)
            final_dir <- file.path(save_dir,group) 
            echodata::get_data(fname = x, 
                               repo = "RajLabMSSM/catalogueR",
                               save_dir = final_dir)
        }) 
        out_paths <- unlist(out_paths)
        names(out_paths) <- gsub("\\.tsv\\.gz","",basename(out_paths))
        return(out_paths)
    } else {
        stop("0 files selected.")
    } 
}

#### Deprecation function #####
example_eQTL_Catalogue_query_paths <- function(...){
  .Deprecated("eQTLcatalogue_example_queries")
  eQTLcatalogue_example_queries(...)
}
