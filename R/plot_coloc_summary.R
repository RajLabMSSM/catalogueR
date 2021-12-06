#' Plot summary of coloc tests
#'
#' Plots the output from \link[catalogueR]{run_coloc}.
#' 
#' @param label_snp_groups Whether to label SNP groups.
#' @param save_dir Directory to save plot in.
#' @param show_plot Whether to print the plot. 
#' @inheritParams merge_colocalized_data
#' 
#' @family coloc
#' @export
#' @importFrom dplyr %>% rename mutate n_distinct group_by summarise
#' @importFrom tidyr separate
#' @importFrom data.table data.table
#' @examples
#' coloc_QTLs <- catalogueR::get_coloc_QTLs()
#' gg_coloc <- catalogueR::plot_coloc_summary(
#'     coloc_QTLs = coloc_QTLs,
#'     coloc_thresh = .5
#' )
plot_coloc_summary <- function(coloc_QTLs,
                               coloc_thresh = .8,
                               label_snp_groups = FALSE,
                               save_dir = tempdir(),
                               show_plot = TRUE) {
    requireNamespace("ggplot2")
    
    Locus.GWAS <- qtl_id <- PP.H4 <- PP.H3 <- gene.QTL <- qtl_ID <- PP.Hyp4 <- 
        PP.H4.thresh <- eGene <- leadSNP <- Consensus_SNP <- SNP <- Support <- 
        Locus.eGene <- Consensus.sigQTL <- UCS.sigQTL <- leadGWAS.sigQTL <-
        NULL;
    coloc_dat <- subset(coloc_QTLs, !is.na(Locus.GWAS)) %>%
        dplyr::rename(qtl_ID = qtl_id) %>%
        dplyr::mutate(
            PP.H4.thresh = ifelse(PP.H4 >= coloc_thresh, PP.H4, NA),
            PP.Hyp4 = ifelse((PP.H3 + PP.H4 >= coloc_thresh) &
                (PP.H4 / PP.H3 >= 2), PP.H4, NA),
            Locus.eGene = paste0(Locus.GWAS, " (", gene.QTL, ")")
        ) %>%
        tidyr::separate(
            col = qtl_ID, into = c("group.QTL", "id"),
            sep = "\\.", remove = FALSE
        ) %>%
        # only plot colocalized loci-egene combinations
        # (otherwise, plot will be massive)
        subset((!is.na(PP.Hyp4) & !is.na(PP.H4.thresh)), .drop = FALSE) %>%
        data.table::data.table()
    if (nrow(coloc_dat) == 0) {
        stop_msg <- paste(
            "\n\n** No coloc tests reached the PPH4 threshold",
            "(", coloc_thresh, ").",
            "Try lowering `coloc_thresh`. **"
        )
        stop(stop_msg)
    }

    needed_cols <- c("leadSNP", "Consensus_SNP", "Support")
    if (label_snp_groups) {
        if (all(needed_cols %in% colnames(coloc_dat))) {
            coloc_dat <- coloc_dat %>%
                dplyr::group_by(Locus.GWAS, eGene, qtl_ID) %>%
                dplyr::summarise(
                    PP.Hyp4 = max(PP.Hyp4),
                    leadGWAS.sigQTL = sum(leadSNP, na.rm = TRUE),
                    Consensus.sigQTL = sum(Consensus_SNP, na.rm = TRUE),
                    UCS.sigQTL = dplyr::n_distinct(SNP[Support > 0],
                        na.rm = TRUE
                    )
                ) %>%
                subset(!is.na(PP.Hyp4)) %>%
                data.table::data.table()
        } else {
            missing_cols <- needed_cols[!needed_cols %in% colnames(coloc_dat)]
            messager(
                "++ Cannot label SNP groups. Missing columns: \n",
                paste(missing_cols, collapse = ", ")
            )
        }
    }

    coloc_dat <- eQTL_Catalogue.annotate_tissues(dat = coloc_dat)
    #### Heatmap ####
    gg_coloc <- ggplot2::ggplot(
        data = coloc_dat,
        ggplot2::aes(
            x = Locus.eGene,
            y = qtl_ID, fill = PP.Hyp4
        )
    ) +
        ggplot2::geom_tile(stat = "identity") +

        # scale_fill_continuous(limits=c(minPP4,1)) +
        ggplot2::scale_fill_gradient(
            na.value = "transparent",
            low = scales::alpha("blue", .7),
            high = scales::alpha("red", .7)
        ) +
        ggplot2::scale_y_discrete(drop = TRUE) +
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::facet_grid(
            facets = System ~ .,
            switch = "y", scales = "free", space = "free"
        ) +
        ggplot2::labs(
            title = paste0(
                "Colocalized GWAS loci x QTL loci (> ",
                coloc_thresh * 100, "% probability)\n"
            ),
            x = "GWAS Locus (QTL eGene)",
            y = "QTL Dataset\n",
            fill = "Colocalization\nProbability"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = .5),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 0),
            panel.border = ggplot2::element_rect(colour = "transparent"),
            strip.placement = "outside",
            strip.background = ggplot2::element_rect(fill = "transparent"),
            strip.text.y = ggplot2::element_text(angle = 90)
        )

    if (label_snp_groups) {
        gg_coloc <- gg_coloc +
            # Consensus SNP markers
            ggplot2::geom_point(
                data = subset(coloc_dat, Consensus.sigQTL > 0),
                ggplot2::aes(
                    x = Locus.eGene, y = qtl_ID,
                    color = "Consensus_SNP"
                ),
                color = "cyan2", shape = 16, size = 1.5,
                show.legend = FALSE
            ) +
            # UCS SNP markers
            ggplot2::geom_point(
                data = subset(coloc_dat, UCS.sigQTL > 0),
                ggplot2::aes(
                    x = Locus.eGene, y = qtl_ID,
                    color = "Union_Credible_Set"
                ),
                size = 3, shape = 5, color = "cyan2",
                show.legend = FALSE
            ) +
            # lead GWAS SNP markers
            ggplot2::geom_tile(
                data = subset(coloc_dat, leadGWAS.sigQTL > 0),
                ggplot2::aes(x = Locus.eGene, y = qtl_ID),
                fill = "transparent", color = "black", size = .7
            )
    }
    #### Print ####
    if (show_plot) print(gg_coloc)
    #### Save ####
    if (save_dir != FALSE) {
        save_path <- file.path(save_dir, paste0(
            "coloc_PP",
            coloc_thresh * 100, ".png"
        ))
        dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
        if (coloc_thresh == .99) {
            width <- 7
            height <- 9
        } else {
            width <- length(unique(coloc_dat$qtl_ID)) * .23 # .333333
            height <- length(unique(coloc_dat$Locus.eGene)) * .3 # 0.6428571
        }

        ggplot2::ggsave(
            filename = save_path,
            plot = gg_coloc,
            height = width, width = height, dpi = 400
        )
    }
    #### Return merged data ####
    return(gg_coloc)
}
