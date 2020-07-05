
#------------Nalls23andMe_2019------------#



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
#' BST1 <- data.table::fread(file.path(root_dir,locus_dir))
#' BST1 <- update_CS_cols(finemap_dat=BST1) 
#' BST1 <- find_consensus_SNPs(finemap_dat=BST1)
#' data.table::fwrite(BST1,"inst/extdata/Nalls23andMe_2019/BST1_Nalls23andMe_2019_subset.tsv.gz", sep="\t")
#' usethis::use_data(BST1, overwrite = T)
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
#' LRRK2 <- data.table::fread(file.path(root_dir,locus_dir))
#' LRRK2 <- update_CS_cols(finemap_dat=LRRK2) 
#' LRRK2 <- find_consensus_SNPs(finemap_dat=LRRK2)
#' data.table::fwrite(LRRK2,"inst/extdata/Nalls23andMe_2019/LRRK2_Nalls23andMe_2019_subset.tsv.gz", sep="\t")
#' usethis::use_data(LRRK2, overwrite = T)
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
#' MEX3C <- data.table::fread(file.path(root_dir,locus_dir))
#' MEX3C <- update_CS_cols(finemap_dat=MEX3C) 
#' MEX3C <- find_consensus_SNPs(finemap_dat=MEX3C)
#' data.table::fwrite(MEX3C,"inst/extdata/Nalls23andMe_2019/MEX3C_Nalls23andMe_2019_subset.tsv.gz", sep="\t")
#' usethis::use_data(MEX3C, overwrite = T)
#' }
"MEX3C"




#------------eQTL Catalogue------------#



#' eQTL Catalogue dataset metadata
#' 
#' List of all queryable tabix-indexed eQTL Catalogue datasets
#' and their associated systems/tissues/cell types.
#' @family eQTL Catalogue
#' @examples 
#' \dontrun{
#' meta <- eQTL_Catalogue.list_datasets(force_new=T)
#' meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte","Fairfax_2014",ftp_path)) 
#' usethis::use_data(meta, overwrite=T)
#' }
"meta"