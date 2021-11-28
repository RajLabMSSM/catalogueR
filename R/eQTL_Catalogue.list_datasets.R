#' List available eQTL datasets
#'
#' Does some additional preprocessing of metadata to categorize tissue types.
#' @keywords metadata
#' @export
#' @examples
#' meta <- eQTL_Catalogue.list_datasets()
eQTL_Catalogue.list_datasets <- function(save_dir = F,
                                         force_new = F,
                                         verbose = FALSE) {
    messager <- function(..., v = TRUE) {
        if (v) {
            print(paste(...))
        }
    }
    if (force_new == FALSE) {
        messager("+ Loading saved metadata.", v = verbose)
        data("meta")
    } else {
        messager("+ Downloading metadata from server.", v = verbose)
        # Main datasets
        URL <- "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"
        meta <- data.table::fread(URL, nThread = 4)
        # # GTEx v8 datasets (metadata was previously in a separate file)
        # URL_gtex <- "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv"
        # meta_gtex <- data.table::fread(URL_gtex, nThread = 4)
        # Merged datasets
        # meta <- rbind(meta, meta_gtex, fill = TRUE)
        # meta <- meta %>% dplyr::transmute(unique_id=paste0(study,".",qtl_group), !!!.)
        meta <- meta %>% dplyr::mutate(unique_id = paste0(study, ".", qtl_group))
        # meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte","Fairfax_2014",ftp_path))
        if (save_dir != FALSE) {
            meta_path <- file.path(save_dir, "eQTLcatalogue_tabix_ftp_paths.tsv")
            messager("Saving metadata ==>", meta_path)
            if (!dir.exists(dirname(meta_path))) dir.create(dirname(meta_path))
            data.table::fwrite(meta, meta_path, sep = "\t")
        }
    }
    messager("++ eQTL Catalogue:: Currently contains", length(unique(meta$unique_id)), "QTL datasets from",
        length(unique(meta$study)), "studies across", length(unique(meta$tissue_label)), "tissues.",
        v = verbose
    )

    # Add group-level annotations
    T_cells <- c("CD4+ T cell", "CD8+ T cell", "B cell", "T cell", "Tfh cell", "Th17 cell", "Th1 cell", "Th2 cell", "Treg naive", "Treg memory")
    meta <- dplyr::mutate(meta, Tissue_group = ifelse(tissue_label %in% T_cells, "T-cell", tissue_label))

    # Add system-level annotations
    blood <- c("blood", "macrophage", "monocyte", "CD16+ monocyte", "neutrophil", "NK cell", "platelet", T_cells)
    CNS <- c("DLPFC")
    meta$System <- ifelse(meta$tissue_label %in% blood, "Blood", ifelse(meta$tissue_label %in% CNS, "CNS", "Other"))
    return(meta)
}
