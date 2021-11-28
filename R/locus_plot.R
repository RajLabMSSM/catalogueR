#' Locus plot
#'
#' Plot colocalized GWAS and QTL results.
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr %>% group_by slice_min
#' @importFrom ggrepel geom_label_repel
#' @export
locus_plot <- function(coloc_QTLs,
                       locus_name,
                       LD_reference = "1KGphase3",
                       LD_results_dir = "data/COVID19_GWAS",
                       plot_results_dir = "plots",
                       force_new_LD = F,
                       plot_zoom = "1x",
                       plot_title = "eQTL Catalogue",
                       top_title = "GWAS",
                       top_x_lab = NULL,
                       show_plot = TRUE,
                       save_plot = TRUE,
                       dpi = 350,
                       height = 12,
                       width = 7) {
    requireNamespace("echolocatoR")
    #### Prepare plot data ####
    plot_dat <- prepare_plot_data(
        coloc_QTLs = coloc_QTLs,
        locus_name = locus_name,
        LD_reference = LD_reference,
        LD_results_dir = LD_results_dir,
        force_new_LD = force_new_LD
    )
    lead_qtls <- plot_dat %>%
        dplyr::group_by(Study, tissue_condition, gene.QTL) %>%
        dplyr::slice_min(order_by = pvalue.QTL, n = 1, with_ties = FALSE)

    #### Gene transcripts plot ####
    gg_genes <- echolocatoR::PLOT.transcript_model_track(finemap_dat = plot_dat) +
        theme(plot.margin = margin(rep(0, 4)))

    #### GWAS plot ####
    gg_gwas <- ggplot(plot_dat, aes(x = Mb, y = -log10(P), color = r2)) +
        geom_point(alpha = .7) +
        scale_color_gradient(low = "blue", high = "red", limits = c(0, 1), breaks = c(0, .5, 1)) +
        geom_point(
            data = subset(plot_dat, leadSNP), shape = 9, size = 3,
            aes(x = Mb, y = -log10(P))
        ) +
        ggrepel::geom_label_repel(
            data = subset(plot_dat, leadSNP)[1, ],
            fill = alpha("white", .5), seed = 2019,
            aes(x = Mb, y = -log10(P), label = SNP)
        ) +
        geom_hline(yintercept = -log10(5e-8), alpha = .5, linetype = "dashed") +
        labs(title = top_title, x = top_x_lab) +
        theme_bw() +
        theme(
            plot.margin = margin(rep(0, 4)),
            axis.text.x = element_blank()
        )

    #### QTL plots ####
    gg_qtl <- ggplot(subset(plot_dat, PP.H4 > .8), aes(x = Mb, y = -log10(pvalue.QTL), color = r2)) +
        geom_point(alpha = .7, show.legend = FALSE) +
        scale_color_gradient(low = "blue", high = "red", limits = c(0, 1), breaks = c(0, .5, 1)) +
        ### Highlight sig QTLs
        # geom_point(data = subset(plot_dat, sig_QTL),
        #            aes(x = Mb, y=-log10(pvalue.QTL)), size=4, color="cyan", shape=1) +
        ### Label lead QTLs
        geom_point(
            data = lead_qtls, shape = 9, size = 3,
            aes(x = Mb, y = -log10(pvalue.QTL))
        ) +
        ggrepel::geom_label_repel(
            data = lead_qtls,
            fill = alpha("white", .5), seed = 2019,
            aes(x = Mb, y = -log10(pvalue.QTL), label = SNP)
        ) +
        ### Facet
        facet_grid(facets = Study + tissue_condition + gene.QTL ~ .) +
        theme_bw() +
        theme(
            strip.background.y = element_rect(fill = "white"),
            strip.text.x = element_text(color = "black")
        ) +
        labs(title = plot_title) +
        theme(plot.margin = margin(rep(0.1, 4)))


    #### Merge plots ####
    plot_list <- list("genes" = gg_genes, "GWAS" = gg_gwas, "QTL" = gg_qtl)
    plot_list <- echolocatoR::PLOT.set_window_limits(
        TRKS = plot_list,
        finemap_dat = plot_dat,
        plot.zoom = plot_zoom
    )
    gg_merged <- patchwork::wrap_plots(plot_list) +
        patchwork::plot_layout(ncol = 1, heights = c(1 / 2, 1 / 4, 1)) +
        patchwork::plot_annotation(
            title = paste("Locus:", plot_dat$Locus[1]),
            tag_levels = letters
        )
    if (show_plot) print(gg_merged)

    if (save_plot) {
        ggsave(
            filename = file.path(
                plot_results_dir,
                paste("coloc.locus_plots", plot_dat$Locus[1], plot_zoom, "pdf", sep = ".")
            ),
            gg_merged, dpi = dpi, height = height, width = width
        )
    }
    return(gg_merged)
}
