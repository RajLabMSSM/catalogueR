#' Plot summary of coloc tests
#'
#' Plots the output from \code{catalogueR::run_coloc}.
#' @family coloc
#' @examples
#' data("coloc_QTLs")
#' gg_coloc <- plot_coloc_summary(coloc_QTLs = coloc_QTLs, PP_thresh = .99, label_snp_groups = FALSE)
plot_coloc_summary <- function(coloc_QTLs,
                               PP_thresh = .8,
                               label_snp_groups = F,
                               save_dir = ".",
                               show_plot = TRUE) {
    # merged_DT <- openxlsx::read.xlsx("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/Nalls23andMe_2019_results.xlsx")
    # merged_DT <- echolocatoR::update_CS_cols(merged_DT)

    coloc_dat <- subset(coloc_QTLs, !is.na(Locus.GWAS)) %>%
        dplyr::rename(qtl_ID = qtl_id) %>%
        dplyr::mutate(
            PP.H4.thresh = ifelse(PP.H4 >= PP_thresh, PP.H4, NA),
            PP.Hyp4 = ifelse((PP.H3 + PP.H4 >= PP_thresh) & (PP.H4 / PP.H3 >= 2), PP.H4, NA),
            Locus.eGene = paste0(Locus.GWAS, " (", gene.QTL, ")")
        ) %>%
        tidyr::separate(col = qtl_ID, into = c("group.QTL", "id"), sep = "\\.", remove = FALSE) %>%
        # only plot colocalized loci-egene combinations (otherwise, plot will be massive)
        subset((!is.na(PP.Hyp4) & !is.na(PP.H4.thresh)), .drop = FALSE) %>%
        data.table::data.table()
    if (nrow(coloc_dat) == 0) {
        stop(
            "\n\n** No coloc tests reached the PPH4 threshold (", PP_thresh, ").",
            "Try lowering `PP_thresh`. **"
        )
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
                    UCS.sigQTL = n_distinct(SNP[Support > 0], na.rm = TRUE)
                ) %>%
                subset(!is.na(PP.Hyp4)) %>%
                data.table::data.table()
        } else {
            missing_cols <- needed_cols[!needed_cols %in% colnames(coloc_dat)]
            warning("++ Cannot label SNP groups. Missing columns: \n", paste(missing_cols, collapse = ", "))
        }
    }

    coloc_dat <- eQTL_Catalogue.annotate_tissues(dat = coloc_dat)
    # Heatmap
    gg_coloc <- ggplot(data = coloc_dat, aes(x = Locus.eGene, y = qtl_ID, fill = PP.Hyp4)) +
        # annotate("rect", xmin=-Inf, xmax=Inf, ymin=Tissue_group, ymax=Tissue_group) +
        # geom_rect(data = coloc_dat, aes(xmin=-Inf, xmax=Inf, ymin= qtl_ID, ymax=qtl_ID, color=Tissue_group), alpha=0.4, inherit.aes  = FALSE) +
        geom_tile(stat = "identity") +

        # scale_fill_continuous(limits=c(minPP4,1)) +
        scale_fill_gradient(na.value = "transparent", low = scales::alpha("blue", .7), high = scales::alpha("red", .7)) +
        scale_y_discrete(drop = TRUE) +
        scale_x_discrete(position = "top") +
        facet_grid(
            facets = System ~ .,
            switch = "y", scales = "free", space = "free"
        ) +
        labs(
            title = paste0("Colocalized GWAS loci x QTL loci (> ", PP_thresh * 100, "% probability)\n"),
            x = "GWAS Locus (QTL eGene)",
            y = "QTL Dataset\n",
            fill = "Colocalization\nProbability"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = .5),
            axis.text.x = element_text(angle = 45, hjust = 0),
            # axis.text.y = element_text(ifelse(unique(coloc_dat$Tissue_group) == 'T-cell', 'red', 'black')),
            panel.border = element_rect(colour = "transparent"),
            strip.placement = "outside",
            strip.background = element_rect(fill = "transparent"),
            strip.text.y = element_text(angle = 90)
        )

    if (label_snp_groups) {
        gg_coloc <- gg_coloc +
            # Consensus SNP markers
            geom_point(data = subset(coloc_dat, Consensus.sigQTL > 0), aes(x = Locus.eGene, y = qtl_ID, color = "Consensus_SNP"), color = "cyan2", shape = 16, size = 1.5, show.legend = FALSE) +
            # UCS SNP markers
            geom_point(data = subset(coloc_dat, UCS.sigQTL > 0), aes(x = Locus.eGene, y = qtl_ID, color = "Union_Credible_Set"), size = 3, shape = 5, color = "cyan2", show.legend = FALSE) +
            # lead GWAS SNP markers
            geom_tile(data = subset(coloc_dat, leadGWAS.sigQTL > 0), aes(x = Locus.eGene, y = qtl_ID), fill = "transparent", color = "black", size = .7)
        # geom_point(data = subset(coloc_dat, leadGWAS.sigQTL>0), aes(x=Locus.eGene, y=qtl_ID, color="Lead_GWAS_SNP"), size=5, shape=0, color="black", stroke=1, show.legend = FALSE)
    }
    if (show_plot) print(gg_coloc)


    if (save_dir != FALSE) {
        # save_dir = "./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC"
        save_path <- file.path(save_dir, paste0("coloc_PP", PP_thresh * 100, ".png"))
        dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
        if (PP_thresh == .99) {
            width <- 7
            height <- 9
        } else {
            width <- length(unique(coloc_dat$qtl_ID)) * .23 # .333333
            height <- length(unique(coloc_dat$Locus.eGene)) * .3 # 0.6428571
        }

        ggsave(save_path,
            plot = gg_coloc, height = width, width = height, dpi = 400
        )
    }
    # Return merged data
    return(gg_coloc)
}
