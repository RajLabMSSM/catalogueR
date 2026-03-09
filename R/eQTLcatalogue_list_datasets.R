#' List available eQTL datasets
#'
#' Does some additional preprocessing of metadata to categorize tissue types.
#' @param save_dir Where to save the processed metadata.
#' @param force_new Re-import the metadata from GitHub instead of using
#' metadata that comes included with \pkg{catalogueR}.
#' @param include_imported Include metadata from datasets that have not yet been
#' fully re-processed by eQTL Catalogue's standardized pipeline.
#' @param verbose Print messages.
#' @keywords metadata
#' @source \href{https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources}{
#' eQTL-Catalogue GitHub repo}
#' 
#' @export
#' @importFrom data.table fread fwrite
#' @importFrom dplyr mutate
#' @examples
#' \dontrun{
#' meta <- catalogueR::eQTLcatalogue_list_datasets()
#' }
eQTLcatalogue_list_datasets <- function(save_dir = tempdir(),
                                        force_new = FALSE,
                                        include_imported = TRUE,
                                        verbose = FALSE) {
    study_label <- sample_group <- tissue_label <- NULL;
    
    if (force_new == FALSE) {
        messager("+ Loading saved metadata.", v = verbose)
        meta <- catalogueR::meta
    } else {
        messager("+ Downloading metadata from server.", v = verbose)
        # Main datasets
        base_url <- paste(
            "https://raw.githubusercontent.com/eQTL-Catalogue",
            "eQTL-Catalogue-resources/master",sep="/"
        )
        URL <- file.path(base_url, "tabix/tabix_ftp_paths.tsv")
        meta <- data.table::fread(URL) 
        if (save_dir != FALSE) {
            meta_path <- file.path(
                save_dir,
                "eQTLcatalogue_tabix_ftp_paths.tsv"
            )
            messager("Saving metadata ==>", meta_path,v=verbose)
            if (!dir.exists(dirname(meta_path))) dir.create(dirname(meta_path))
            data.table::fwrite(meta, meta_path, sep = "\t")
        }
    }
    messager("++ eQTL Catalogue:: Currently contains",
        formatC(length(unique(meta$unique_label)), big.mark = ","),
        "QTL datasets from",
        formatC(length(unique(meta$study_label)), big.mark = ","),
        "studies across",
        formatC(length(unique(meta$tissue_label)), big.mark = ","),
        "tissues.",
        v = verbose
    )
    #### Add group-level annotations ####
    T_cells <- c(
        "CD4+ T cell", "CD8+ T cell", "B cell", "T cell",
        "Tfh cell", "Th17 cell", "Th1 cell", "Th2 cell",
        "Treg naive", "Treg memory"
    )
    meta <- dplyr::mutate(meta,
        Tissue_group = ifelse(tissue_label %in% T_cells,
            "T-cell", tissue_label
        )
    )
    #### Add system-level annotations ####
    blood <- c(
        "blood", "macrophage", "monocyte", "CD16+ monocyte",
        "neutrophil", "NK cell", "platelet", T_cells
    )
    CNS <- c("DLPFC")
    meta$System <- ifelse(meta$tissue_label %in% blood, "Blood",
        ifelse(meta$tissue_label %in% CNS, "CNS", "Other")
    )
    
    # Create a human-readable unique label for each dataset
    meta <- meta |> dplyr::mutate(
      unique_label = paste0(study_label, ".", sample_group)
      )
    return(meta)
}

#### Deprecation function #####
eQTL_Catalogue.list_datasets <- function(...){
  .Deprecated("eQTLcatalogue_list_datasets")
  eQTLcatalogue_list_datasets(...)
}
