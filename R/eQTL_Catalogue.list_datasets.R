#' List available eQTL datasets
#'
#' Does some additional preprocessing of metadata to categorize tissue types.
#' 
#' @param save_dir Where to save the processed metadata.
#' @param force_new Re-import the metadata from GitHub instead of using 
#' metadata that comes included with \pkg{catalogueR}.
#' @param include_imported Include metadata from datasets that have not yet been
#' fully re-processed by eQTL Catalogue's standardized pipeline. 
#' @param verbose Print messages.
#' @keywords metadata
#' @export
#' @examples
#' meta <- catalogueR::eQTL_Catalogue.list_datasets()
eQTL_Catalogue.list_datasets <- function(save_dir = tempdir(),
                                         force_new = FALSE,
                                         include_imported = TRUE,
                                         verbose = FALSE) {
    if (force_new == FALSE) {
        messager("+ Loading saved metadata.", v = verbose)
        data("meta")
    } else {
        messager("+ Downloading metadata from server.", v = verbose)
        # Main datasets
        base_url <- file.path(
            "https://raw.githubusercontent.com/eQTL-Catalogue",
            "eQTL-Catalogue-resources/master"
        )
        URL <- file.path(base_url,"tabix/tabix_ftp_paths.tsv")
        meta <- data.table::fread(URL, nThread = 1)
        
        #### Import GTEX V8 metadata ####
        if(include_imported){
            messager("Including metadata from tabix_ftp_paths_imported.tsv",
                     v=verbose)
            meta <- tryCatch(expr = {
                # GTEx v8 metadata is in a separate file
                URL_gtex <- file.path(base_url,
                                      "tabix/tabix_ftp_paths_imported.tsv")
                meta_gtex <- data.table::fread(URL_gtex, nThread = 1)
                # Merged datasets
                rbind(meta, meta_gtex, fill = TRUE) %>%  
                    dplyr::mutate(unique_id = paste0(study, ".", qtl_group))
            }, error  = function(e) meta)
        } 
        # meta <- meta %>% dplyr::mutate(ftp_path= 
        #                                    gsub("Fairfax_2014_monocyte",
        #                                         "Fairfax_2014",ftp_path))
        if (save_dir != FALSE) {
            meta_path <- file.path(save_dir,
                                   "eQTLcatalogue_tabix_ftp_paths.tsv")
            messager("Saving metadata ==>", meta_path)
            if (!dir.exists(dirname(meta_path))) dir.create(dirname(meta_path))
            data.table::fwrite(meta, meta_path, sep = "\t")
        }
    }
    messager("++ eQTL Catalogue:: Currently contains",
             formatC(length(unique(meta$unique_id)), big.mark = ","), 
             "QTL datasets from",
        formatC(length(unique(meta$study)),big.mark = ","),
        "studies across",
        formatC(length(unique(meta$tissue_label)),big.mark = ","),
        "tissues.",
        v = verbose
    )
    #### Add group-level annotations ####
    T_cells <- c("CD4+ T cell", "CD8+ T cell", "B cell", "T cell", 
                 "Tfh cell", "Th17 cell", "Th1 cell", "Th2 cell",
                 "Treg naive", "Treg memory")
    meta <- dplyr::mutate(meta, 
                          Tissue_group = ifelse(tissue_label %in% T_cells, 
                                                "T-cell", tissue_label))
    #### Add system-level annotations ####
    blood <- c("blood", "macrophage", "monocyte", "CD16+ monocyte",
               "neutrophil", "NK cell", "platelet", T_cells)
    CNS <- c("DLPFC")
    meta$System <- ifelse(meta$tissue_label %in% blood, "Blood",
                          ifelse(meta$tissue_label %in% CNS, "CNS", "Other"))
    return(meta)
}
