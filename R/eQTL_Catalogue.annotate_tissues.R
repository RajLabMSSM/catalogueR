#' Annotate QTL datasets with metadata
#'
#' Annotate QTL datasets from \emph{eQTL Catalogue} with metadata.
#'
#' @keywords metadata
#' @importFrom dplyr %>% mutate select group_by summarise n_distinct rename
#' @importFrom data.table merge.data.table
#' @importFrom stats setNames
#' @examples
#' dat <- eQTL_Catalogue.annotate_tissues(dat)
eQTL_Catalogue.annotate_tissues <- function(dat,
                                            add_tissue_counts = FALSE) {
    meta <- eQTL_Catalogue.list_datasets() %>%
        dplyr::select(unique_id,
            Study = study,
            Tissue = tissue_label,
            Tissue_group,
            System
        ) %>%
        dplyr::mutate(Tissue = gsub("CD16+ monocyte", "monocyte", Tissue)) %>%
        dplyr::mutate(Tissue = factor(Tissue,
            levels = unique(Tissue),
            ordered = TRUE
        )) %>%
        unique()

    if (!"study" %in% colnames(dat)) {
        dat$Study <- "GWAS"
    }
    if (!"qtl_id" %in% colnames(dat)) {
        dat$qtl_id <- dat$qtl_ID
    }
    dat <- data.table::merge.data.table(dat %>% dplyr::rename(Study_id = Study),
        meta,
        by.x = "qtl_id",
        by.y = "unique_id"
    )
    # Add the number of datasets per Tissue as a col
    if (add_tissue_counts) {
        dat_count <- dat %>%
            dplyr::group_by(Tissue) %>%
            dplyr::summarise(dataset_count = dplyr::n_distinct(qtl_id,
                na.rm = TRUE
            )) %>%
            data.table::data.table()
        countDict <- stats::setNames(dat_count$dataset_count, dat_count$Tissue)
        dat$dataset_count <- countDict[dat$Tissue]
        # data.table::merge.data.table(x=dat, y=dat_count, by="Tissue")
        dat$Tissue_count <- paste0(dat$Tissue, " (n=", dat$dataset_count, ")")
        dat$Tissue_count <- factor(dat$Tissue_count,
            levels = unique(dat$Tissue_count),
            ordered = TRUE
        )
    }
    return(dat)
}
