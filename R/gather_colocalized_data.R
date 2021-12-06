#' Prepare data for coloc plot
#'
#' Use the coloc results (coloc_QTLs) to determine which full summary stats
#'  (gwas.qtl_paths) to plot.
#'  
#' @param gwas.qtl_paths Query results paths from 
#' \link[catalogueR]{eQTL_Catalogue.query}.
#' @param coloc_QTLs Colocalization results from \link[catalogueR]{run_coloc}. 
#' @param qtl_thresh QTL uncorrected p-value ("pvalues.QTL") threshold.
#' @param coloc_thresh Colocalization Posterior Probability threshold,
#'  using the formula:
#' \code{(PP.H3 + PP.H4 >= coloc_thresh) & (PP.H4 / PP.H3 >= 2}.
#' @param gwas_label Label for the GWAS subplot.
#' @param remove_extra_panes Remove SNPs from non-significant panes.
#' @param verbose Print messages. 
#' 
#' @family coloc
#' @export
merge_colocalized_data <- function(gwas.qtl_paths,
                                    coloc_QTLs = NULL,
                                    qtl_thresh = 1e-5,
                                    coloc_thresh = .8,
                                    gwas_label = "GWAS",
                                    remove_extra_panes = TRUE,
                                    verbose = TRUE) {
    PP.H4 <- PP.H3 <- Locus.GWAS <- qtl_id <- gene.QTL <- PP.Hyp4 <- 
        PP.H4.thresh <- pvalues.QTL <- POS <- Locus.QTL.eGene <- NULL;
    # Filter colocalization
    coloc_QTLs_sig <- coloc_QTLs %>%
        dplyr::mutate(
            PP.H4.thresh = ifelse(PP.H4 > coloc_thresh, PP.H4, NA),
            PP.Hyp4 = ifelse((PP.H3 + PP.H4 >= coloc_thresh) &
                                 (PP.H4 / PP.H3 >= 2), PP.H4, NA)
        ) %>%
        dplyr::mutate(Locus.QTL.eGene = paste(Locus.GWAS,
                                              qtl_id, gene.QTL, sep = "--"))
    if (!is.null(coloc_thresh)) {
        coloc_QTLs_sig <- subset(coloc_QTLs_sig,
                                 (!is.na(PP.Hyp4)) & (!is.na(PP.H4.thresh)))
    }
    if (nrow(coloc_QTLs_sig) == 0) {
        stop_msg <- paste(
            "No rows left after filtering coloc PP.",
            "Try using less stringent `coloc_thresh` (smaller, or =NULL).")
        stop(stop_msg)
    }

    # Filter by QTL p-value
    if (!is.null(qtl_thresh)) {
        coloc_QTLs_sig <- subset(coloc_QTLs_sig, pvalues.QTL < qtl_thresh)
    }
    if (nrow(coloc_QTLs_sig) == 0) {
        stop_msg2 <- paste(
            "No rows left after filtering QTL p-values.",
            "Try using less stringent `qtl_thresh` (larger, or =NULL).")
        stop(stop_msg2)
    }

    messager("++", dplyr::n_distinct(coloc_QTLs_sig$Locus.QTL.eGene), 
             "GWAS.QTL.eGene combinations will be plotted.", v = verbose)


    # Determine which GWAS.QTL files you actually need to import
    gwas.qtl_paths_sig <- reconstruct_file_names(coloc_QTLs = coloc_QTLs_sig)
    gwas.qtl_paths_select <- gwas.qtl_paths[
        basename(gwas.qtl_paths) %in% gwas.qtl_paths_sig]
    GWAS.QTL <- merge_files(file_paths = gwas.qtl_paths_select)

    # Filter out non-colocalized eGenes by unique ID
    plot_dat <- GWAS.QTL %>%
        dplyr::mutate(
            MB = POS / 1000000,
            Locus__eGene = paste(Locus.GWAS, gene.QTL, sep = "__"),
            Locus.QTL.eGene = paste(Locus.GWAS, qtl_id, gene.QTL, sep = "--"),
            GWAS.label = gwas_label
        )
    # Add coloc probs back into SNP-wise dataframe 
    # (only for tests that colocalized)
    h4_dict <- stats::setNames(
        coloc_QTLs_sig$PP.H4,
        coloc_QTLs_sig$Locus.QTL.eGene
    )
    plot_dat$PP.H4 <- h4_dict[plot_dat$Locus.QTL.eGene]

    # Which panes to remove if not colocalized 
    # but still appear in plot because it's a grid
    tmp_subset <- subset(plot_dat, 
                         Locus.QTL.eGene %in% unique(
                             coloc_QTLs_sig$Locus.QTL.eGene))
    if (remove_extra_panes) {
        # Remove SNPs from nonsig panes
        plot_dat <- tmp_subset
    } else {
        plot_dat <- subset(
            plot_dat,
            ((qtl_id %in% tmp_subset$qtl_id) & 
                 (Locus.GWAS %in% tmp_subset$Locus.GWAS)) |
                ((qtl_id %in% tmp_subset$qtl_id) &
                     (gene.QTL %in% tmp_subset$gene.QTL))
            # (gene.QTL %in% tmp_subset$gene.QTL)
        )
    }
    return(plot_dat)
}


gather_colocalized_data <- function(...){
    .Deprecated("merge_colocalized_data")
    merge_colocalized_data(...)
}
