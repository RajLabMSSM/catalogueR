#### ------------eQTL Catalogue------------####

#' eQTL Catalogue tabix header
#' 
#' 
#' The eQTL Catalogue tabix header (column names)
#'  is not always included in the file. 
#' @source
#' \code{
#' eQTL_Catalogue.header <- tabix_header(force_new_header = TRUE)
#' usethis::use_data(eQTL_Catalogue.header, overwrite = TRUE)
#' }
#' @family eQTL Catalogue 
#' @usage data("eQTL_Catalogue.header")
"eQTL_Catalogue.header"


#' eQTL Catalogue dataset metadata
#'
#' List of all queryable tabix-indexed eQTL Catalogue datasets
#' and their associated systems/tissues/cell types.
#' @source 
#' \code{
#' meta <- eQTL_Catalogue.list_datasets(force_new = TRUE)
#' # Some paths in the metadata were originally wrong. 
#' # Has since been corrected by authors.
#' # meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte",
#' #                                               "Fairfax_2014",ftp_path))
#' usethis::use_data(meta, overwrite = TRUE)
#' }
#' @family eQTL Catalogue
"meta"


#' eQTL Catalogue query results
#'
#' Query: BST1 GWAS locus, Alasoo_2018.macrophage_IFNg QTL dataset.
#' Split output file from \link[catalogueR]{eQTL_Catalogue.query}.
#' @source 
#' \code{
#' root_dir <- "~/Desktop/catalogueR_queries/Alasoo_2018.macrophage_IFNg"
#' file_name <- "BST1__Alasoo_2018.macrophage_IFNg.tsv.gz"
#' BST1__Alasoo_2018.macrophage_IFNg <- data.table::fread(file.path(root_dir,
#'                                                                  file_name))
#' usethis::use_data(BST1__Alasoo_2018.macrophage_IFNg, overwrite = TRUE)
#' }
#' @family eQTL Catalogue 
"BST1__Alasoo_2018.macrophage_IFNg"


#' eQTL Catalogue query results
#'
#' Query: LRRK2 GWAS locus, Alasoo_2018.macrophage_IFNg QTL dataset.
#' Split output file from \link[catalogueR]{eQTL_Catalogue.query}.
#' @source 
#' \code{
#' root_dir <- "~/Desktop/catalogueR_queries/Alasoo_2018.macrophage_IFNg"
#' file_name <- "LRRK2__Alasoo_2018.macrophage_IFNg.tsv.gz"
#' LRRK2__Alasoo_2018.macrophage_IFNg <- data.table::fread(file.path(root_dir,
#'                                                                   file_name))
#' usethis::use_data(LRRK2__Alasoo_2018.macrophage_IFNg, overwrite = TRUE)
#' }
#' @family eQTL Catalogue 
"LRRK2__Alasoo_2018.macrophage_IFNg"


#' eQTL Catalogue query results
#'
#' Query: MEX3C GWAS locus, Alasoo_2018.macrophage_IFNg QTL dataset.
#' Split output file from \link[catalogueR]{eQTL_Catalogue.query}.
#' @source 
#' \code{
#' root_dir <- "~/Desktop/catalogueR_queries/Alasoo_2018.macrophage_IFNg"
#' file_name <- "MEX3C__Alasoo_2018.macrophage_IFNg.tsv.gz"
#' MEX3C__Alasoo_2018.macrophage_IFNg <- data.table::fread(file.path(root_dir,
#'                                                                   file_name))
#' usethis::use_data(MEX3C__Alasoo_2018.macrophage_IFNg, overwrite = TRUE)
#' }
#' @family eQTL Catalogue
"MEX3C__Alasoo_2018.macrophage_IFNg"


#### ----------coloc -------------####

#' Example colocalization results
#'
#' Example colocalization results from running \link[catalogueR]{run_coloc}.
#' on GWAS summary stats from all loci in
#'  \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls23andMe_2019}.
#' @source 
#' \code{
#' gwas.qtl_paths <- example_eQTL_Catalogue_query_paths()
#' coloc_QTLs <- run_coloc(gwas.qtl_paths = gwas.qtl_paths,
#'                         nThread = 4, 
#'                         top_snp_only = TRUE, 
#'                         save_path = "~/Desktop/coloc_results.tsv.gz")
#' usethis::use_data(coloc_QTLs, overwrite = TRUE)
#' }
#' @family coloc
"coloc_QTLs"


#' Example colocalization results
#'
#' Example colocalization results from running \link[catalogueR]{run_coloc}.
#' on GWAS summary stats from all loci in 
#' \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls23andMe_2019}.
#' @source 
#' \code{
#' library(dplyr)
#' root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' sub_dir <- "_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz"
#' coloc_QTLs_full <- data.table::fread(file.path(root_dir, sub_dir))
#' coloc_QTLs_full <- coloc_QTLs_full %>% 
#'   dplyr::rename(gene.QTL = eGene, qtl_id = qtl.id)
#' usethis::use_data(coloc_QTLs_full, overwrite = TRUE)
#' }
#' @family coloc 
"coloc_QTLs_full"

