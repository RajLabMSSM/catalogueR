

#' List available eQTL datasets
#' 
#' Does some additional preprocessing of metadata to categorize tissue types.
#' @keywords metadata
#' @examples 
#' meta <- eQTL_Catalogue.list_datasets()
eQTL_Catalogue.list_datasets <- function(save_dir=F,
                                         force_new=F,
                                         verbose=F){ 
  printer <- function(..., v=T){if(v){print(paste(...))}}
  if(force_new==F){
    printer("+ Loading saved metadata.",v = verbose)
    data("meta")
  } else {
    printer("+ Downloading metadata from server.",v = verbose)
    # Main datasets
    URL <- "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"
    meta <- data.table::fread(URL, nThread = 4)
    # GTEx v8 datasets (standardized by not reprossed from raw data)
    URL_gtex <- "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv"
    meta_gtex <- data.table::fread(URL_gtex, nThread = 4)
    # Merged datasets
    meta <- rbind(meta, meta_gtex, fill=T)
    meta <- meta %>% dplyr::transmute(unique_id=paste0(study,".",qtl_group), !!!.) 
    # meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte","Fairfax_2014",ftp_path)) 
    if(save_dir!=F){
      meta_path <- file.path(save_dir,"eQTLcatalogue_tabix_ftp_paths.tsv")
      printer("Saving metadata ==>",meta_path)
      if(!dir.exists(dirname(meta_path)))dir.create(dirname(meta_path))
      data.table::fwrite(meta, meta_path, sep="\t")
    }
  }
  printer("++ eQTL Catalogue:: Currently contains",length(unique(meta$unique_id)),"QTL datasets from",
          length(unique(meta$study)),"studies across",length(unique(meta$tissue_label)),"tissues.",v=verbose)
  
  # Add group-level annotations 
  T_cells <- c("CD4+ T cell","CD8+ T cell","B cell","T cell","Tfh cell","Th17 cell","Th1 cell","Th2 cell","Treg naive","Treg memory")
  meta <- dplyr::mutate(meta, Tissue_group=ifelse(tissue_label %in% T_cells, "T-cell",tissue_label))
  
  # Add system-level annotations
  blood <- c("blood","macrophage","monocyte","CD16+ monocyte","neutrophil","NK cell","platelet",T_cells)
  CNS <- c("DLPFC") 
  meta$System <- ifelse(meta$tissue_label %in% blood,"Blood",ifelse(meta$tissue_label %in% CNS, "CNS", "Other"))
  return(meta)
}




#' Search eQTL Catalogue metadata
#' 
#' Searches through multiple relevant metadata columns to find eQTL Catalogue datasets
#' that match at least one of your substrings in a list. 
#' All searches are case-insensitive.
#' If \code{qtl_search=NULL}, will return all available datasets.  
#' 
#' @family eQTL Catalogue
#' @examples
#' qtl_datasets <- eQTL_Catalogue.search_metadata(qtl_search=c("Alasoo_2018","monocyte"))
#' qtl_datasets.brain <- eQTL_Catalogue.search_metadata(qtl_search="brain")
eQTL_Catalogue.search_metadata <- function(qtl_search=NULL,
                                           verbose=T){
  meta <- eQTL_Catalogue.list_datasets()
  if(is.null(qtl_search)){
    printer("eQTL_Catalogue:: Gathering data for all QTL Catalogue datasets...", v=verbose)
    qtl_datasets <- unique(meta$unique_id)
  }else {
    qtl_datasets <- grep(pattern = paste(qtl_search,collapse="|"),
                         x = meta$unique_id,
                         value = T,
                         ignore.case = T) %>% unique()
  } 
  return(qtl_datasets)
}



#' Annotate QTL datasets with metadata 
#' 
#' @keywords metadata 
#' @examples 
#' dat <- eQTL_Catalogue.annotate_tissues(dat)
eQTL_Catalogue.annotate_tissues <- function(dat,
                                            add_tissue_counts=F){ 
  meta <- eQTL_Catalogue.list_datasets() %>% 
    dplyr::select(unique_id, 
                  Study=study, 
                  Tissue=tissue_label, 
                  Tissue_group,
                  System) %>%
    dplyr::mutate(Tissue=gsub('CD16+ monocyte','monocyte', Tissue)) %>% 
    dplyr::mutate(Tissue=factor(Tissue, levels = unique(Tissue), ordered = T) ) %>%
    unique() 
  
  dat <- data.table::merge.data.table(dat %>% dplyr::rename(Study_id=Study), 
                                      meta, 
                                      by.x="qtl_id",
                                      by.y="unique_id" ) 
  # Add the number of datasets per Tissue as a col
  if(add_tissue_counts){
    dat_count <- dat %>% 
      dplyr::group_by(Tissue) %>% 
      dplyr::summarise(dataset_count=n_distinct(qtl_id, na.rm = T)) %>% 
      data.table::data.table() 
    countDict <-  setNames(dat_count$dataset_count, dat_count$Tissue)
    dat$dataset_count <- countDict[dat$Tissue]   #data.table::merge.data.table(x=dat, y=dat_count, by="Tissue")
    dat$Tissue_count <- paste0(dat$Tissue," (n=",dat$dataset_count,")")
    dat$Tissue_count <- factor(dat$Tissue_count, 
                               levels = unique(dat$Tissue_count), ordered = T) 
  } 
  return(dat)
}

