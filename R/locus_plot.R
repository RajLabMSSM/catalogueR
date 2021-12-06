#' Locus plot
#'
#' Plot colocalized GWAS and QTL results.
#' @param coloc_QTLs Colocalization results from \link[catalogueR]{run_coloc}.
#' @param locus_name Locus name.
#' @param LD_results_dir Directory to store LD matrices in.
#' @param plot_results_dir Directory to store plots in.
#' @param plot_zoom How many times to zoom into plot.
#' @param QTL_title QTL subplot title.
#' @param GWAS_title GWAS subplot title.
#' @param GWAS_x_lab GWAS subplot x-axis label.
#' @param show_plot Whether to print the plot.
#' @param save_plot Whether to save the plot.
#' @inheritParams echoLD::load_or_create
#' @inheritParams ggplot2::ggsave
#' 
#' @importFrom patchwork wrap_plots plot_layout plot_annotation 
#' @importFrom dplyr %>% group_by slice_min
#' @importFrom ggrepel geom_label_repel
#' @export
locus_plot <- function(coloc_QTLs,
                       locus_name,
                       LD_reference = "1KGphase3",
                       LD_results_dir = file.path(tempdir(),"plots_LD"),
                       plot_results_dir = file.path(tempdir(),"plots"),
                       force_new_LD = FALSE,
                       plot_zoom = "1x",
                       QTL_title = "eQTL Catalogue",
                       GWAS_title = "GWAS",
                       GWAS_x_lab = NULL,
                       show_plot = TRUE,
                       save_plot = TRUE,
                       dpi = 350,
                       height = 12,
                       width = 7) {
    requireNamespace("echolocatoR")
    requireNamespace("ggplot2")
    
    Study <- tissue_condition <- gene.QTL <- pvalue.QTL <- Mb <- P <- r2 <- 
        leadSNP <- SNP <- PP.H4 <- NULL;
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
    gg_genes <- echolocatoR::PLOT.transcript_model_track(
        finemap_dat = plot_dat) +
        ggplot2::theme(plot.margin = ggplot2:: margin(rep(0, 4)))

    #### GWAS plot ####
    gg_gwas <- ggplot2::ggplot(plot_dat, 
                               ggplot2::aes(x = Mb, y = -log10(P),
                                            color = r2)) +
        ggplot2::geom_point(alpha = .7) +
        ggplot2::scale_color_gradient(low = "blue", high = "red",
                             limits = c(0, 1), breaks = c(0, .5, 1)) +
        ggplot2::geom_point(
            data = subset(plot_dat, leadSNP), shape = 9, size = 3,
            ggplot2::aes(x = Mb, y = -log10(P))
        ) +
        ggrepel::geom_label_repel(
            data = subset(plot_dat, leadSNP)[1, ],
            fill = ggplot2::alpha("white", .5), seed = 2019,
            ggplot2::aes(x = Mb, y = -log10(P), label = SNP)
        ) +
        ggplot2::geom_hline(yintercept = -log10(5e-8), alpha = .5,
                            linetype = "dashed") +
        ggplot2::labs(title = GWAS_title, x = GWAS_x_lab) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(rep(0, 4)),
            axis.text.x = ggplot2::element_blank()
        )
    #### QTL plots ####
    gg_qtl <- ggplot2::ggplot(subset(plot_dat, PP.H4 > .8),
                              ggplot2::aes(x = Mb, y = -log10(pvalue.QTL),
                                            color = r2)) +
        ggplot2::geom_point(alpha = .7, show.legend = FALSE) +
        ggplot2::scale_color_gradient(low = "blue", high = "red",
                             limits = c(0, 1), breaks = c(0, .5, 1)) +
        ### Highlight sig QTLs
        # geom_point(data = subset(plot_dat, sig_QTL),
                   # aes(x = Mb, y=-log10(pvalue.QTL)), 
                   #      size=4, color="cyan", shape=1) +
        ### Label lead QTLs
        ggplot2::geom_point(
            data = lead_qtls, shape = 9, size = 3,
            ggplot2::aes(x = Mb, y = -log10(pvalue.QTL))
        ) +
        ggrepel::geom_label_repel(
            data = lead_qtls,
            fill = ggplot2::alpha("white", .5), seed = 2019,
            ggplot2::aes(x = Mb, y = -log10(pvalue.QTL), label = SNP)
        ) +
        ### Facet
        ggplot2::facet_grid(facets = Study + tissue_condition + gene.QTL ~ .) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.background.y = ggplot2::element_rect(fill = "white"),
            strip.text.x = ggplot2::element_text(color = "black")
        ) +
        ggplot2::labs(title = QTL_title) +
        ggplot2::theme(plot.margin = ggplot2::margin(rep(0.1, 4)))


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
        ggplot2::ggsave(
            filename = file.path(
                plot_results_dir,
                paste("coloc.locus_plots", plot_dat$Locus[1],
                      plot_zoom, "pdf", sep = ".")
            ),
            gg_merged, dpi = dpi, height = height, width = width
        )
    }
    return(gg_merged)
}
