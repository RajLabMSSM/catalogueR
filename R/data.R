


#### ------------Nalls23andMe_2019------------####


#' \emph{echolocatoR} output example (BST1 locus)
#'
#' An example results file after running
#' \code{\link{finemap_loci}} on the \emph{BST1} locus.
#'
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#'
#' @format data.table
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @family Nalls23andMe_2019
#' @examples
#' \dontrun{
#' root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' locus_dir <- "BST1/Multi-finemap/Multi-finemap_results.txt"
#' BST1 <- data.table::fread(file.path(root_dir, locus_dir))
#' BST1 <- update_CS_cols(finemap_dat = BST1)
#' BST1 <- find_consensus_SNPs(finemap_dat = BST1)
#' data.table::fwrite(BST1, "inst/extdata/Nalls23andMe_2019/BST1_Nalls23andMe_2019_subset.tsv.gz", sep = "\t")
#' usethis::use_data(BST1, overwrite = TRUE)
#' }
"BST1"




#' \emph{echolocatoR} output example (LRRK2 locus)
#'
#' An example results file after running
#' \code{\link{finemap_loci}} on the \emph{LRRK2} locus.
#'
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#'
#' @format data.table
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @family Nalls23andMe_2019
#' @examples
#' \dontrun{
#' root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' locus_dir <- "LRRK2/Multi-finemap/Multi-finemap_results.txt"
#' LRRK2 <- data.table::fread(file.path(root_dir, locus_dir))
#' LRRK2 <- update_CS_cols(finemap_dat = LRRK2)
#' LRRK2 <- find_consensus_SNPs(finemap_dat = LRRK2)
#' data.table::fwrite(LRRK2, "inst/extdata/Nalls23andMe_2019/LRRK2_Nalls23andMe_2019_subset.tsv.gz", sep = "\t")
#' usethis::use_data(LRRK2, overwrite = TRUE)
#' }
"LRRK2"




#' \emph{echolocatoR} output example (MEX3C locus)
#'
#' An example results file after running
#' \code{\link{finemap_loci}} on the \emph{MEX3C} locus.
#'
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#'
#' @format data.table
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @family Nalls23andMe_2019
#' @examples
#' \dontrun{
#' root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' locus_dir <- "MEX3C/Multi-finemap/Multi-finemap_results.txt"
#' MEX3C <- data.table::fread(file.path(root_dir, locus_dir))
#' MEX3C <- update_CS_cols(finemap_dat = MEX3C)
#' MEX3C <- find_consensus_SNPs(finemap_dat = MEX3C)
#' data.table::fwrite(MEX3C, "inst/extdata/Nalls23andMe_2019/MEX3C_Nalls23andMe_2019_subset.tsv.gz", sep = "\t")
#' usethis::use_data(MEX3C, overwrite = TRUE)
#' }
"MEX3C"




#### ------------eQTL Catalogue------------####


#' eQTL Catalogue tabix header
#'
#' @family eQTL Catalogue
#' @examples
#' \dontrun{
#' eQTL_Catalogue.header <- tabix_header(force_new_header = TRUE)
#' usethis::use_data(eQTL_Catalogue.header, overwrite = TRUE)
#' }
"eQTL_Catalogue.header"




#' eQTL Catalogue dataset metadata
#'
#' List of all queryable tabix-indexed eQTL Catalogue datasets
#' and their associated systems/tissues/cell types.
#' @family eQTL Catalogue
#' @examples
#' \dontrun{
#' meta <- eQTL_Catalogue.list_datasets(force_new = TRUE)
#' # Some paths in the metadata were originally wrong. Has since been corrected by authors.
#' ### meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte","Fairfax_2014",ftp_path))
#' usethis::use_data(meta, overwrite = TRUE)
#' }
"meta"




#' eQTL Catalogue query results
#'
#' Query: BST1 GWAS locus, Alasoo_2018.macrophage_IFNg QTL dataset.
#' Split output file from \code{eQTL_Catalogue.query()}.
#' @family eQTL Catalogue
#' @examples
#' \dontrun{
#' root_dir <- "~/Desktop/catalogueR_queries/Alasoo_2018.macrophage_IFNg"
#' file_name <- "BST1__Alasoo_2018.macrophage_IFNg.tsv.gz"
#' BST1__Alasoo_2018.macrophage_IFNg <- data.table::fread(file.path(root_dir, file_name))
#' usethis::use_data(BST1__Alasoo_2018.macrophage_IFNg, overwrite = TRUE)
#' }
"BST1__Alasoo_2018.macrophage_IFNg"





#' eQTL Catalogue query results
#'
#' Query: LRRK2 GWAS locus, Alasoo_2018.macrophage_IFNg QTL dataset.
#' Split output file from \code{eQTL_Catalogue.query()}.
#' @family eQTL Catalogue
#' @examples
#' \dontrun{
#' root_dir <- "~/Desktop/catalogueR_queries/Alasoo_2018.macrophage_IFNg"
#' file_name <- "LRRK2__Alasoo_2018.macrophage_IFNg.tsv.gz"
#' LRRK2__Alasoo_2018.macrophage_IFNg <- data.table::fread(file.path(root_dir, file_name))
#' usethis::use_data(LRRK2__Alasoo_2018.macrophage_IFNg, overwrite = TRUE)
#' }
"LRRK2__Alasoo_2018.macrophage_IFNg"





#' eQTL Catalogue query results
#'
#' Query: MEX3C GWAS locus, Alasoo_2018.macrophage_IFNg QTL dataset.
#' Split output file from \code{eQTL_Catalogue.query()}.
#' @family eQTL Catalogue
#' @examples
#' \dontrun{
#' root_dir <- "~/Desktop/catalogueR_queries/Alasoo_2018.macrophage_IFNg"
#' file_name <- "MEX3C__Alasoo_2018.macrophage_IFNg.tsv.gz"
#' MEX3C__Alasoo_2018.macrophage_IFNg <- data.table::fread(file.path(root_dir, file_name))
#' usethis::use_data(MEX3C__Alasoo_2018.macrophage_IFNg, overwrite = TRUE)
#' }
"MEX3C__Alasoo_2018.macrophage_IFNg"





#### ----------coloc -------------####


#' Example colocalization results
#'
#' Example colocalization results from running \code{\link{catalogueR::run_coloc}}
#' on GWAS summary stats from all loci in \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls23andMe_2019}.
#' @family coloc
#' @examples
#' \dontrun{
#' gwas.qtl_paths <- example_eQTL_Catalogue_query_paths()
#' coloc_QTLs <- run_coloc(gwas.qtl_paths = gwas.qtl_paths, nThread = 4, top_snp_only = TRUE, save_path = "~/Desktop/coloc_results.tsv.gz")
#' usethis::use_data(coloc_QTLs, overwrite = TRUE)
#' }
"coloc_QTLs"




#' Example colocalization results
#'
#' Example colocalization results from running \code{\link{catalogueR::run_coloc}}
#' on GWAS summary stats from all loci in \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls23andMe_2019}.
#' @family coloc
#' @examples
#' \dontrun{
#' library(dplyr)
#' root_dir <- "~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"
#' sub_dir <- "_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz"
#' coloc_QTLs_full <- data.table::fread(file.path(root_dir, sub_dir))
#' coloc_QTLs_full <- coloc_QTLs_full %>% dplyr::rename(gene.QTL = eGene, qtl_id = qtl.id)
#' usethis::use_data(coloc_QTLs_full, overwrite = TRUE)
#' }
"coloc_QTLs_full"
